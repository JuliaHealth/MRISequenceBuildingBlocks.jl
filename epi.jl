using KomaMRI
using KomaMRI.PulseDesigner: make_adc, make_label, make_trapezoid
using Unitful


# EPI keeps its own readout timing because the ramp simultaneously constrains the
# bipolar Gx reversal and the Gy blip; it is not a Cartesian GRE readout ramp.
function _epi_readout_params(FOV, matrix, sys, BWpp; n_shots=1)
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
    M_blip = n_shots * ΔMy
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

function _epi_adc_delay(params, polarity)
    return polarity > 0 ?
        params.ζ + params.adc_center_shift :
        params.ζ + params.Ta - params.Ta_adc - params.adc_center_shift
end


function _epi_readout_line(params, sys; polarity=1.0, line_index=0)
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
        make_label(:SET, :NAV, 0),
        make_label(:SET, :AVG, 0),
    ]
    return RO
end


function _epi_base(FOV, matrix, sys::Scanner, BWpp::Real;
    partial_fourier=1.0,
    n_shots=1)

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

    # Phase partial Fourier uses integer line indices; this intentionally differs
    # from the half-sample centering used for Cartesian readout partial Fourier.
    full_precenter = N_pe ÷ 2
    N_precenter = round(Int, 2 * (partial_fourier - 0.5) * full_precenter)
    i_start = -N_precenter
    i_stop = N_pe - 1 - full_precenter
    N_acq = i_stop - i_start + 1
    1 ≤ n_shots ≤ N_acq || error("n_shots must be between 1 and $N_acq (got $n_shots)")

    params = _epi_readout_params(FOV, matrix, sys, BWpp; n_shots)
    (; dt, Ga, ζ, ΔMy) = params
    M_blip = n_shots * ΔMy
    G_blip = M_blip / ζ

    n_echo = N_ro ÷ 2
    M_ro = Ga * (ζ / 2 + params.adc_center_shift + n_echo * dt)

    Mx_pre = -M_ro
    Mz_pre = 0.0

    function shot_geometry(shot)
        i_first = i_start + shot - 1
        n_lines = length(i_first:n_shots:i_stop)
        return (; i_first, n_lines)
    end

    function design_shot_prephaser(i_first)
        My_pre_shot = i_first * ΔMy
        M_vec_pre_shot = sqrt(Mx_pre^2 + My_pre_shot^2 + Mz_pre^2)
        T_p_shot, ζ_p_shot = M_vec_pre_shot > 0 ?
            _lobe_timing(M_vec_pre_shot, sys) : (0.0, sys.GR_Δt)
        inv_area_p_shot = M_vec_pre_shot > 0 ? 1.0 / (T_p_shot + ζ_p_shot) : 0.0
        return (;
            y_moment=My_pre_shot,
            flat=T_p_shot,
            ramp=ζ_p_shot,
            inv_area=inv_area_p_shot,
            duration=T_p_shot + 2ζ_p_shot,
        )
    end

    function build_epi(shot)
        (; i_first, n_lines) = shot_geometry(shot)
        prephaser = design_shot_prephaser(i_first)

        PRE = Sequence(sys)
        addblock!(
            PRE;
            x=make_trapezoid(; amplitude=Mx_pre * prephaser.inv_area * u"T/m",
                flat_time=prephaser.flat * u"s", rise_time=prephaser.ramp * u"s",
                fall_time=prephaser.ramp * u"s", sys),
            y=make_trapezoid(; amplitude=prephaser.y_moment * prephaser.inv_area * u"T/m",
                flat_time=prephaser.flat * u"s", rise_time=prephaser.ramp * u"s",
                fall_time=prephaser.ramp * u"s", sys),
            z=make_trapezoid(; amplitude=Mz_pre * prephaser.inv_area * u"T/m",
                flat_time=prephaser.flat * u"s", rise_time=prephaser.ramp * u"s",
                fall_time=prephaser.ramp * u"s", sys)
        )

        seq = PRE
        # Koma constructors and sequence concatenation reject negative delays.
        # Record the blips and mark them out-of-block only after assembly, when
        # no subsequent sequence copy is required.
        outboard_blip_blocks = Int[]
        for line in 1:n_lines
            polarity = isodd(line) ? 1.0 : -1.0
            line_index = i_first + (line - 1) * n_shots + full_precenter
            RO = _epi_readout_line(params, sys; polarity, line_index)
            seq += RO

            if line < n_lines
                BLIP = Sequence(sys)
                addblock!(
                    BLIP;
                    x=make_trapezoid(; amplitude=0.0u"T/m", flat_time=0.0u"s",
                        rise_time=ζ * u"s", fall_time=ζ * u"s", sys),
                    y=make_trapezoid(; amplitude=G_blip * u"T/m", flat_time=0.0u"s",
                        rise_time=ζ * u"s", fall_time=ζ * u"s", sys),
                    z=make_trapezoid(; amplitude=0.0u"T/m", flat_time=0.0u"s",
                        rise_time=ζ * u"s", fall_time=ζ * u"s", sys)
                )
                BLIP.DUR[1] = 0.0
                seq += BLIP
                push!(outboard_blip_blocks, length(seq))
            end
        end

        for block in outboard_blip_blocks
            seq.GR[2, block].delay = -ζ
        end

        return seq
    end

    function center_time(shot)
        (; i_first) = shot_geometry(shot)
        prephaser = design_shot_prephaser(i_first)
        acquired_lines = collect(i_first:n_shots:i_stop)
        center_line = argmin(abs.(acquired_lines))
        polarity = isodd(center_line) ? 1.0 : -1.0
        return prephaser.duration +
            (center_line - 1) * (params.Ta + 2params.ζ) +
            _epi_adc_delay(params, polarity) +
            n_echo * params.dt
    end

    return (; readout=build_epi, center_time)
end


function _epi_shot_center_time(seq, matrix, partial_fourier, n_shots, shot)
    N_pe = Int(matrix[2])
    full_precenter = N_pe ÷ 2
    n_precenter = round(Int, 2 * (partial_fourier - 0.5) * full_precenter)
    i_start = -n_precenter
    i_stop = N_pe - 1 - full_precenter
    acquired_lines = collect(i_start + shot - 1:n_shots:i_stop)
    isempty(acquired_lines) && error("EPI shot $shot has no acquired lines.")
    center_line = argmin(abs.(acquired_lines))
    adc_blocks = findall(block -> seq.ADC[block].N > 0, 1:length(seq))
    adc = seq.ADC[adc_blocks[center_line]]
    sample_spacing = adc.N == 1 ? 0.0 : adc.T / (adc.N - 1)
    return get_block_start_times(seq)[adc_blocks[center_line]] +
        adc.delay +
        (adc.N ÷ 2) * sample_spacing
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
    epi_readout_kernel(FOV, matrix, sys; BWpp, kwargs...)

Design a two-dimensional Cartesian EPI readout kernel. Phase-encoding lines
are interleaved across shots, and each shot alternates readout polarity while
overlapping its phase blips with the readout ramps.

The returned named tuple contains:

- `readout(shot=1)`: build one shot with all gradients materialized inside
  Pulseq blocks;
- `center_time(shot=1)`: analytically calculated time of the ADC sample nearest
  `(kx, ky) = (0, 0)`;
- `center_time(seq, shot=1)`: recover that time from a concrete readout;
- `adc_timing_ok(seq)`: check ADC dead-time and block-fit constraints;
- `BWpp` and `n_shots`: the requested design values.

ADC blocks carry zero-based `LIN` labels plus `REV`, `SEG`, `NAV`, and `AVG`
state. `FOV` and returned timing values use SI units.

# Keywords
- `BWpp`: Readout bandwidth per pixel. [`Hz/pixel`]
- `partial_fourier=1.0`: Acquired phase-encoding fraction. `ky=0` is always
  included.
- `n_shots=1`: Number of interleaved EPI shots.
"""
function epi_readout_kernel(FOV, matrix, sys::Scanner;
    BWpp,
    partial_fourier=1.0,
    n_shots=1,
)
    base = _epi_base(
        FOV,
        matrix,
        sys,
        BWpp;
        partial_fourier,
        n_shots,
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
        return _epi_shot_center_time(seq, matrix, partial_fourier, n_shots, shot)
    end

    function center_time(shot::Integer=1)
        validate_shot(shot)
        return base.center_time(shot)
    end

    return (;
        readout,
        center_time,
        adc_timing_ok=seq -> _adc_events_fit_pulseq(seq, sys),
        BWpp,
        n_shots,
    )
end
