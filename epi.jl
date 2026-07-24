using KomaMRI
using KomaMRI.PulseDesigner: make_adc, make_delay, make_label, make_trapezoid
using Unitful


# EPI keeps its own readout timing because the ramp simultaneously constrains the
# bipolar Gx reversal and the Gy blip; it is not a Cartesian GRE readout ramp.
function _epi_readout_params(FOV, matrix, sys, BWpp; max_blip_steps=1)
    FOV_ro, FOV_pe = Float64.(FOV)
    N_ro = Int(matrix[1])

    dt_ideal = 1.0 / (N_ro * BWpp)
    N_dwell = max(round(Int, dt_ideal / sys.ADC_Δt), 1)
    dt = N_dwell * sys.ADC_Δt
    adc_center_shift = mod(dt / 2, sys.RF_Δt)

    Ga = 1.0 / (γ * dt * FOV_ro)
    Ga > sys.Gmax && error(
        "Readout Ga=$(round(Ga*1e3,digits=2)) mT/m exceeds " *
        "Gmax=$(round(sys.Gmax*1e3,digits=2)) mT/m - reduce BWpp")

    Ta_adc = dt * (N_ro - 1)
    Ta = ceil((Ta_adc + adc_center_shift) / sys.GR_Δt) * sys.GR_Δt

    ΔMy = 1.0 / (γ * FOV_pe)
    M_blip = max_blip_steps * ΔMy
    ζ_slew = sqrt(
        (Ga^2 + sqrt(Ga^4 + 4 * sys.Smax^2 * M_blip^2)) /
        (2 * sys.Smax^2)
    )
    ζ_adc = max(
        dt / 2 + sys.ADC_dead_time - adc_center_shift,
        dt / 2 + sys.ADC_dead_time + adc_center_shift - (Ta - Ta_adc),
    )
    ζ = ceil(
        max(ζ_slew, abs(M_blip) / sys.Gmax, ζ_adc) / sys.GR_Δt,
    ) * sys.GR_Δt

    return (; FOV_ro, FOV_pe, N_ro, dt, adc_center_shift, Ga, Ta_adc, Ta, ζ, ΔMy)
end


function _epi_line_groups(i_start, i_stop, n_shots, shot_partition)
    lines = collect(i_start:i_stop)
    if shot_partition == :interleaved
        return [collect(lines[shot:n_shots:end]) for shot in 1:n_shots]
    elseif shot_partition == :contiguous
        base_count, extra = divrem(length(lines), n_shots)
        groups = Vector{Vector{Int}}(undef, n_shots)
        offset = 1
        for shot in 1:n_shots
            count = base_count + (shot <= extra)
            groups[shot] = lines[offset:(offset + count - 1)]
            offset += count
        end
        return groups
    end
    error("shot_partition must be :interleaved or :contiguous (got $shot_partition)")
end


function _epi_common_lobe_timing(moments, sys)
    worst_moment = maximum(moment -> sqrt(sum(abs2, moment)), moments)
    flat, ramp = worst_moment > 0 ?
        _lobe_timing(worst_moment, sys) : (0.0, sys.GR_Δt)
    return (;
        flat,
        ramp,
        inv_area=worst_moment > 0 ? 1.0 / (flat + ramp) : 0.0,
        duration=flat + 2ramp,
    )
end


function _epi_gradient_lobe(moment, timing, sys)
    seq = Sequence(sys)
    addblock!(
        seq;
        x=make_trapezoid(; amplitude=moment[1] * timing.inv_area * u"T/m",
            flat_time=timing.flat * u"s", rise_time=timing.ramp * u"s",
            fall_time=timing.ramp * u"s", sys),
        y=make_trapezoid(; amplitude=moment[2] * timing.inv_area * u"T/m",
            flat_time=timing.flat * u"s", rise_time=timing.ramp * u"s",
            fall_time=timing.ramp * u"s", sys),
        z=make_trapezoid(; amplitude=moment[3] * timing.inv_area * u"T/m",
            flat_time=timing.flat * u"s", rise_time=timing.ramp * u"s",
            fall_time=timing.ramp * u"s", sys),
    )
    return seq
end


function _epi_adc_delay(params, polarity)
    return polarity > 0 ?
        params.ζ + params.adc_center_shift :
        params.ζ + params.Ta - params.Ta_adc - params.adc_center_shift
end


function _epi_readout_line(
    params,
    sys;
    polarity=1.0,
    line_index=0,
    navigator=false,
    average=0,
)
    reversed = Int(polarity < 0)
    RO = Sequence(sys)
    addblock!(
        RO;
        x=make_trapezoid(; amplitude=polarity * params.Ga * u"T/m",
            flat_time=params.Ta * u"s", rise_time=params.ζ * u"s",
            fall_time=params.ζ * u"s", sys),
        y=make_trapezoid(; amplitude=0.0u"T/m", flat_time=params.Ta * u"s",
            rise_time=params.ζ * u"s", fall_time=params.ζ * u"s", sys),
        z=make_trapezoid(; amplitude=0.0u"T/m", flat_time=params.Ta * u"s",
            rise_time=params.ζ * u"s", fall_time=params.ζ * u"s", sys)
    )
    adc_delay = _epi_adc_delay(params, polarity)
    adc = params.N_ro == 1 ?
        ADC(params.N_ro, params.Ta_adc, adc_delay) :
        make_adc(
            params.N_ro,
            params.dt * u"s";
            delay=(adc_delay - params.dt / 2) * u"s",
            sys,
        )
    RO.ADC = [adc]
    RO.EXT[1] = [
        make_label(:SET, :LIN, line_index),
        make_label(:SET, :REV, reversed),
        make_label(:SET, :SEG, reversed),
        make_label(:SET, :NAV, Int(navigator)),
        make_label(:SET, :AVG, average),
    ]
    return RO
end


function _epi_echo_group(
    lines,
    params,
    sys;
    full_precenter,
    prephaser_moment,
    prephaser_timing,
    rewinder_moment,
    rewinder_timing,
    rewind_m0,
)
    seq = _epi_gradient_lobe(prephaser_moment, prephaser_timing, sys)

    # Koma constructors and sequence concatenation reject negative delays.
    # Record the blips and mark them out-of-block only after assembly, when
    # no subsequent sequence copy is required.
    outboard_blip_blocks = Int[]
    for echo in eachindex(lines)
        polarity = isodd(echo) ? 1.0 : -1.0
        line_index = lines[echo] + full_precenter
        seq += _epi_readout_line(params, sys; polarity, line_index)

        if echo < length(lines)
            blip_moment = (lines[echo+1] - lines[echo]) * params.ΔMy
            BLIP = Sequence(sys)
            addblock!(
                BLIP;
                x=make_trapezoid(; amplitude=0.0u"T/m", flat_time=0.0u"s",
                    rise_time=params.ζ * u"s", fall_time=params.ζ * u"s", sys),
                y=make_trapezoid(; amplitude=blip_moment / params.ζ * u"T/m",
                    flat_time=0.0u"s", rise_time=params.ζ * u"s",
                    fall_time=params.ζ * u"s", sys),
                z=make_trapezoid(; amplitude=0.0u"T/m", flat_time=0.0u"s",
                    rise_time=params.ζ * u"s", fall_time=params.ζ * u"s", sys),
            )
            BLIP.DUR[1] = 0.0
            seq += BLIP
            push!(outboard_blip_blocks, length(seq))
        end
    end

    rewind_m0 &&
        (seq += _epi_gradient_lobe(rewinder_moment, rewinder_timing, sys))
    for block in outboard_blip_blocks
        seq.GR[2, block].delay = -params.ζ
    end
    return seq
end


function _epi_base(FOV, matrix, sys::Scanner, BWpp::Real;
    partial_fourier=1.0,
    n_shots=1,
    shot_partition=:interleaved,
    rewind_m0=true)

    length(FOV) == length(matrix) || error(
        "FOV and matrix must have the same length " *
        "(got $(length(FOV)) and $(length(matrix)))")
    length(FOV) == 2 || error("EPI readouts are only defined for 2D FOV/matrix inputs")
    all(N -> N > 0 && isinteger(N), matrix) ||
        error("matrix entries must be positive integers")
    FOV_ro, FOV_pe = Float64.(FOV)
    N_ro, N_pe = Int.(matrix)
    FOV_ro > 0 && FOV_pe > 0 || error("FOV entries must be positive")
    BWpp > 0 || error("BWpp must be positive")
    0.5 ≤ partial_fourier ≤ 1.0 || error(
        "partial_fourier must be between 0.5 and 1.0 (got $partial_fourier)")
    n_shots isa Integer || error("n_shots must be an integer")
    shot_partition in (:interleaved, :contiguous) || error(
        "shot_partition must be :interleaved or :contiguous (got $shot_partition)")

    # Phase partial Fourier uses integer line indices; this intentionally differs
    # from the half-sample centering used for Cartesian readout partial Fourier.
    full_precenter = N_pe ÷ 2
    N_precenter = round(Int, 2 * (partial_fourier - 0.5) * full_precenter)
    i_start = -N_precenter
    i_stop = N_pe - 1 - full_precenter
    N_acq = i_stop - i_start + 1
    1 ≤ n_shots ≤ N_acq || error("n_shots must be between 1 and $N_acq (got $n_shots)")

    line_groups = _epi_line_groups(i_start, i_stop, n_shots, shot_partition)
    max_blip_steps = maximum(
        lines -> length(lines) > 1 ? maximum(abs, diff(lines)) : 0,
        line_groups,
    )
    params = _epi_readout_params(FOV, matrix, sys, BWpp; max_blip_steps)
    (; dt, Ga, ζ, ΔMy) = params

    n_echo = N_ro ÷ 2
    M_ro = Ga * (ζ / 2 + params.adc_center_shift + n_echo * dt)
    Mx_pre = -M_ro
    readout_moment = Ga * (params.Ta + params.ζ)
    prephaser_moments = [
        (Mx_pre, first(lines) * ΔMy, 0.0)
        for lines in line_groups
    ]
    rewinder_moments = [
        let x_moment = Mx_pre + sum(
                isodd(echo) ? readout_moment : -readout_moment
                for echo in eachindex(lines)
            )
            (-x_moment, -last(lines) * ΔMy, 0.0)
        end
        for lines in line_groups
    ]
    prephaser_timing = _epi_common_lobe_timing(prephaser_moments, sys)
    rewinder_timing = _epi_common_lobe_timing(rewinder_moments, sys)

    function build_epi(shot)
        return _epi_echo_group(
            line_groups[shot],
            params,
            sys;
            full_precenter,
            prephaser_moment=prephaser_moments[shot],
            prephaser_timing,
            rewinder_moment=rewinder_moments[shot],
            rewinder_timing,
            rewind_m0,
        )
    end

    function build_navigator_readout(n_lines)
        n_lines isa Integer || error("n_lines must be an integer.")
        n_lines > 0 || error("n_lines must be positive.")

        prephaser_moment = (Mx_pre, 0.0, 0.0)
        line_moment = sum(
            isodd(line) ? readout_moment : -readout_moment
            for line in 1:n_lines
        )
        rewinder_moment = (-(Mx_pre + line_moment), 0.0, 0.0)
        navigator_prephaser_timing =
            _epi_common_lobe_timing([prephaser_moment], sys)
        navigator_rewinder_timing =
            _epi_common_lobe_timing([rewinder_moment], sys)

        seq = _epi_gradient_lobe(
            prephaser_moment,
            navigator_prephaser_timing,
            sys,
        )
        for line in 1:n_lines
            seq += _epi_readout_line(
                params,
                sys;
                polarity=isodd(line) ? 1.0 : -1.0,
                line_index=full_precenter,
                navigator=true,
                average=Int(line == n_lines),
            )
        end
        seq += _epi_gradient_lobe(
            rewinder_moment,
            navigator_rewinder_timing,
            sys,
        )
        append!(
            seq.EXT[end],
            [
                make_label(:SET, :REV, 0),
                make_label(:SET, :SEG, 0),
                make_label(:SET, :NAV, 0),
                make_label(:SET, :AVG, 0),
            ],
        )
        return seq
    end

    function center_time(shot)
        lines = line_groups[shot]
        center_line = argmin(abs.(lines))
        polarity = isodd(center_line) ? 1.0 : -1.0
        return prephaser_timing.duration +
            (center_line - 1) * (params.Ta + 2params.ζ) +
            _epi_adc_delay(params, polarity) +
            n_echo * params.dt
    end

    function echo_center_time(shot)
        n_lines = length(line_groups[shot])
        lower = (n_lines + 1) ÷ 2
        upper = (n_lines + 2) ÷ 2
        line_time(line) =
            prephaser_timing.duration +
            (line - 1) * (params.Ta + 2params.ζ) +
            _epi_adc_delay(params, isodd(line) ? 1.0 : -1.0) +
            n_echo * params.dt
        return (line_time(lower) + line_time(upper)) / 2
    end

    return (;
        readout=build_epi,
        center_time,
        echo_center_time,
        navigator_readout=build_navigator_readout,
        line_groups,
        prephaser_timing,
        rewinder_timing,
    )
end


function _epi_adc_center_time(seq, adc_block)
    adc = seq.ADC[adc_block]
    sample_spacing = adc.N == 1 ? 0.0 : adc.T / (adc.N - 1)
    return get_block_start_times(seq)[adc_block] +
        adc.delay +
        (adc.N ÷ 2) * sample_spacing
end


function _epi_shot_center_time(seq, lines)
    center_line = argmin(abs.(lines))
    adc_blocks = findall(block -> seq.ADC[block].N > 0, 1:length(seq))
    return _epi_adc_center_time(seq, adc_blocks[center_line])
end


function _epi_echo_center_time(seq)
    adc_blocks = findall(block -> seq.ADC[block].N > 0, 1:length(seq))
    lower = (length(adc_blocks) + 1) ÷ 2
    upper = (length(adc_blocks) + 2) ÷ 2
    return (
        _epi_adc_center_time(seq, adc_blocks[lower]) +
        _epi_adc_center_time(seq, adc_blocks[upper])
    ) / 2
end

# Koma does not currently expose this Pulseq-facing ADC check independently of
# its full timing validator, so the kernel keeps the small predicate local.
function _adc_events_fit_pulseq(seq, sys)
    for block in 1:length(seq)
        adc = seq.ADC[block]
        adc.N > 0 || continue
        dwell = adc.N == 1 ? adc.T : adc.T / (adc.N - 1)
        adc.delay >= dwell / 2 + sys.ADC_dead_time - 1e-12 || return false
        dur(adc, sys) <= seq.DUR[block] + 1e-12 || return false
    end
    return true
end


"""
    build_epi_navigator(excitation, epi_kernel; n_lines=3, post_delay)

Build a complete EPI navigator from an independent excitation and the readout
timing already designed by [`epi_readout_kernel`](@ref). The navigator acquires
`n_lines` alternating-polarity `ky=0` lines and explicitly rewinds the x
gradient moment, so either odd or even line counts are valid.

Navigator ADC blocks carry `NAV=1`, the zero-based center-line `LIN`, and
matching `REV`/`SEG` polarity state. `AVG` is zero except on the final
navigator line, where it is one. The terminal
rewinder resets `NAV`, `REV`, `SEG`, and `AVG` to zero so those states cannot
leak into a following readout. `post_delay` is a required nonnegative delay
between the navigator and subsequent imaging excitation, rounded up to the
block raster. The returned sequence can be prepended to an EPI acquisition or
passed as the `navigator` keyword to [`build_grase`](@ref) or
[`build_tse`](@ref).
"""
function build_epi_navigator(excitation, epi_kernel; n_lines=3, post_delay)
    hasproperty(epi_kernel, :navigator) ||
        error("The EPI kernel must provide `navigator`.")
    return epi_kernel.navigator(excitation; n_lines, post_delay)
end


"""
    epi_readout_kernel(FOV, matrix, sys; BWpp, kwargs...)

Design a two-dimensional Cartesian EPI readout kernel. Phase-encoding lines
are interleaved or divided into contiguous groups across shots. Each returned
shot is a self-contained EPI echo group with common worst-case prephaser,
blip, and rewinder timing. Readout polarity alternates while phase blips
overlap the readout ramps.

The returned named tuple contains:

- `readout(shot=1)`: build one shot with all gradients materialized inside
  Pulseq blocks;
- `navigator(excitation; n_lines=3, post_delay)`: build a complete navigator
  with an independent excitation and required post-navigator delay;
- `navigator_readout(n_lines=3)`: build only the alternating `ky=0` lines and
  x-moment rewinder;
- `lines`: centered phase-encoding indices acquired by each shot;
- `center_time(shot=1)`: analytically calculated time of the ADC sample nearest
  `(kx, ky) = (0, 0)`;
- `center_time(seq, shot=1)`: recover that time from a concrete readout;
- `echo_center_time(shot=1)`: temporal center of the EPI echo group, equal to
  the central line's `kx=0` sample for an odd-length group;
- `echo_center_time(seq, shot=1)`: recover that time from a concrete readout;
- `adc_timing_ok(seq)`: check ADC dead-time and block-fit constraints;
- `FOV` and `matrix`: validated two-dimensional design geometry;
- `BWpp`, `n_shots`, `shot_partition`, and `rewind_m0`: the requested design
  values.

ADC blocks carry zero-based `LIN` labels plus `REV`, `SEG`, `NAV`, and `AVG`
state. `FOV` and returned timing values use SI units.

# Keywords
- `BWpp`: Readout bandwidth per pixel. [`Hz/pixel`]
- `partial_fourier=1.0`: Acquired phase-encoding fraction. `ky=0` is always
  included.
- `n_shots=1`: Number of EPI echo groups.
- `shot_partition=:interleaved`: Use interleaved phase-encoding lines in each
  group. Pass `:contiguous` for adjacent line groups.
- `rewind_m0=true`: Append a simultaneous x/y rewinder so every group ends with
  zero gradient zeroth moment.
"""
function epi_readout_kernel(FOV, matrix, sys::Scanner;
    BWpp,
    partial_fourier=1.0,
    n_shots=1,
    shot_partition=:interleaved,
    rewind_m0=true,
)
    base = _epi_base(
        FOV,
        matrix,
        sys,
        BWpp;
        partial_fourier,
        n_shots,
        shot_partition,
        rewind_m0,
    )

    function validate_shot(shot)
        1 ≤ shot ≤ n_shots || error(
            "shot must be between 1 and $n_shots (got $shot)")
    end

    function readout(shot::Integer=1)
        validate_shot(shot)
        return materialize_gradients(base.readout(shot))
    end

    function center_time(seq::Sequence, shot::Integer=1)
        validate_shot(shot)
        return _epi_shot_center_time(seq, base.line_groups[shot])
    end

    function center_time(shot::Integer=1)
        validate_shot(shot)
        return base.center_time(shot)
    end

    function echo_center_time(seq::Sequence, shot::Integer=1)
        validate_shot(shot)
        return _epi_echo_center_time(seq)
    end

    function echo_center_time(shot::Integer=1)
        validate_shot(shot)
        return base.echo_center_time(shot)
    end

    function navigator(excitation; n_lines=3, post_delay)
        post_delay >= 0 || error("post_delay must be non-negative.")
        actual_post_delay = ceil_to_raster(post_delay, sys.DUR_Δt)
        seq = Sequence(sys)
        @addblock begin
            seq += excitation
            seq += base.navigator_readout(n_lines)
            actual_post_delay > 0 &&
                (seq += make_delay(actual_post_delay * u"s"))
        end
        return seq
    end

    return (;
        readout,
        navigator,
        navigator_readout=base.navigator_readout,
        center_time,
        echo_center_time,
        adc_timing_ok=seq -> _adc_events_fit_pulseq(seq, sys),
        lines=base.line_groups,
        FOV=Tuple(Float64.(FOV)),
        matrix=Tuple(Int.(matrix)),
        BWpp,
        n_shots,
        shot_partition,
        rewind_m0,
    )
end
