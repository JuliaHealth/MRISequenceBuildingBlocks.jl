using KomaMRI
using KomaMRI.PulseDesigner: make_delay, make_sinc_pulse, make_trapezoid
using Unitful


function _sequence_rf_center(seq)
    block_starts = get_block_start_times(seq)
    centers = Float64[]
    for block in 1:length(seq), coil in axes(seq.RF, 1)
        rf = seq.RF[coil, block]
        dur(rf) > 0 || continue
        push!(centers, block_starts[block] + rf.delay + rf_center(rf))
    end
    length(centers) == 1 || error(
        "Expected exactly one RF event, found $(length(centers)).")
    return only(centers)
end


function _gradient_amplitude_at(t, amplitude, x)
    x < first(t) && return 0.0
    x > last(t) && return 0.0

    i = searchsortedlast(t, x)
    i == 0 && return first(amplitude)
    i == length(t) && return last(amplitude)
    x == t[i] && return amplitude[i]

    fraction = (x - t[i]) / (t[i+1] - t[i])
    return (1 - fraction) * amplitude[i] + fraction * amplitude[i+1]
end


function _gradient_area_between(gradient, t0, t1)
    t1 <= t0 && return 0.0
    knot_times = times(gradient)
    knot_amplitudes = ampls(gradient)
    isempty(knot_times) && return 0.0

    start_time = max(t0, first(knot_times))
    stop_time = min(t1, last(knot_times))
    stop_time <= start_time && return 0.0

    integration_times = Float64[start_time]
    append!(integration_times, (
        t for t in knot_times if start_time < t < stop_time
    ))
    push!(integration_times, stop_time)

    return sum(
        (_gradient_amplitude_at(knot_times, knot_amplitudes, integration_times[i]) +
         _gradient_amplitude_at(knot_times, knot_amplitudes, integration_times[i+1])) *
        (integration_times[i+1] - integration_times[i]) / 2
        for i in 1:(length(integration_times) - 1)
    )
end


function _refocusing_bandwidth(slice_thickness, sys)
    gradient_limit = γ * sys.Gmax * slice_thickness
    approximate_b1_limit = 2γ * sys.B1
    return min(gradient_limit, approximate_b1_limit)
end


function _make_refocusing_events(
    slice_thickness,
    sys;
    bandwidth,
    time_bw_product,
    phase,
    apodization,
    rf_shift,
)
    requested_bandwidth = isnothing(bandwidth) ?
        _refocusing_bandwidth(slice_thickness, sys) : bandwidth
    requested_bandwidth > 0 || error("Refocusing bandwidth must be positive.")

    duration = ceil_to_raster(
        time_bw_product / requested_bandwidth,
        sys.GR_Δt,
    )
    duration > 0 || error("Refocusing duration must be positive.")

    while true
        rf, gradient, _ = make_sinc_pulse(
            π * u"rad";
            duration=duration * u"s",
            slice_thickness=slice_thickness * u"m",
            phase_offset=phase * u"rad",
            time_bw_product,
            apodization,
            use=Refocusing(),
            sys,
        )

        rf.delay += rf_shift
        gradient.T += ceil_to_raster(rf_shift, sys.GR_Δt)

        if !isnothing(bandwidth) || !isfinite(sys.B1)
            return rf, gradient, duration
        end

        peak_b1 = maximum(abs, ampls(rf))
        peak_b1 <= sys.B1 || peak_b1 ≈ sys.B1 || begin
            longer = ceil_to_raster(duration * peak_b1 / sys.B1, sys.GR_Δt)
            duration = max(longer, duration + sys.GR_Δt)
            continue
        end
        return rf, gradient, duration
    end
end


"""
    build_refocusing_block(slice_thickness, sys; kwargs...)

Build a single slice-selective 180° refocusing block without the usual slice
rephaser. The spin-echo builder accounts for the slice-gradient area in its
pre- and post-refocusing crusher lobes.

# Keywords
- `bandwidth=nothing`: Refocusing bandwidth. The default selects the shortest
  pulse that satisfies the scanner gradient and RF-amplitude limits. [`Hz`]
- `time_bw_product=4`: RF time-bandwidth product.
- `phase=π/2`: Refocusing RF phase. [`rad`]
- `apodization=0.46`: Cosine-window weight.
- `rf_shift=0.0`: RF-raster timing shift used to fit the requested echo time. [`s`]
"""
function build_refocusing_block(
    slice_thickness,
    sys;
    bandwidth=nothing,
    time_bw_product=4,
    phase=π / 2,
    apodization=0.46,
    rf_shift=0.0,
)
    slice_thickness > 0 || error("slice_thickness must be positive.")
    time_bw_product > 0 || error("time_bw_product must be positive.")
    0 <= apodization <= 1 || error("apodization must be between 0 and 1.")
    rf_shift >= 0 || error("rf_shift must be non-negative.")
    isapprox(rf_shift / sys.RF_Δt, round(rf_shift / sys.RF_Δt); rtol=0, atol=1e-9) ||
        error("rf_shift must be on the RF raster.")

    rf, gradient, _ = _make_refocusing_events(
        slice_thickness,
        sys;
        bandwidth,
        time_bw_product,
        phase,
        apodization,
        rf_shift,
    )

    seq = Sequence(sys)
    addblock!(seq, rf; z=gradient)
    seq.DUR[end] = ceil_to_raster(dur(seq[end], sys), sys.DUR_Δt)
    check_timing(seq, sys)
    check_hw_limits(seq, sys)
    return seq
end


function _z_lobe(target_area, sys)
    iszero(target_area) && return Sequence(sys)
    gradient = make_trapezoid(; area=target_area * u"T*s/m", sys)
    seq = Sequence(sys)
    addblock!(seq; z=gradient)
    return seq
end


function _crusher_pair(crusher_phase, slice_thickness, sys, refocus)
    desired_area = crusher_phase / (2π * γ * slice_thickness)
    refocusing_center = _sequence_rf_center(refocus)
    slice_gradient = refocus.GR[3, 1]
    pre_area = desired_area -
        _gradient_area_between(slice_gradient, 0.0, refocusing_center)
    post_area = desired_area -
        _gradient_area_between(slice_gradient, refocusing_center, dur(slice_gradient))
    return _z_lobe(pre_area, sys), _z_lobe(post_area, sys)
end


function _fit_refocused_train_timing(
    excitation_tail,
    packet_center,
    packet_duration,
    n_packets,
    slice_thickness,
    sys;
    crusher_phase,
    refocusing_phase,
    refocusing_bandwidth,
    refocusing_time_bw_product,
    refocusing_apodization,
)
    candidates = NamedTuple[]
    for rf_shift in 0:sys.RF_Δt:(sys.DUR_Δt - sys.RF_Δt)
        refocus = build_refocusing_block(
            slice_thickness,
            sys;
            bandwidth=refocusing_bandwidth,
            time_bw_product=refocusing_time_bw_product,
            phase=refocusing_phase,
            apodization=refocusing_apodization,
            rf_shift,
        )
        pre_crusher, post_crusher =
            _crusher_pair(crusher_phase, slice_thickness, sys, refocus)
        refocusing_center = _sequence_rf_center(refocus)

        refocusing_left = dur(pre_crusher) + refocusing_center
        refocusing_right =
            dur(refocus) - refocusing_center + dur(post_crusher)
        packet_left = packet_center
        packet_right = packet_duration - packet_center

        first_half = excitation_tail + refocusing_left
        before_packet = refocusing_right + packet_left
        after_packet = packet_right + refocusing_left
        target_half = n_packets == 1 ?
            max(first_half, before_packet) :
            max(first_half, before_packet, after_packet)

        first_delay =
            ceil_to_raster(target_half - first_half, sys.DUR_Δt)
        pre_packet_delay =
            ceil_to_raster(target_half - before_packet, sys.DUR_Δt)
        post_packet_delay = n_packets == 1 ? 0.0 :
            ceil_to_raster(target_half - after_packet, sys.DUR_Δt)

        actual_first_half = first_half + first_delay
        actual_before_packet = before_packet + pre_packet_delay
        actual_after_packet = n_packets == 1 ?
            actual_before_packet : after_packet + post_packet_delay
        timing_error = n_packets == 1 ?
            abs(actual_first_half - actual_before_packet) :
            max(
                abs(actual_first_half - actual_before_packet),
                abs(actual_after_packet - actual_before_packet),
            )
        maximum_spacing = n_packets == 1 ?
            actual_first_half + actual_before_packet :
            max(
                actual_first_half + actual_before_packet,
                actual_after_packet + actual_before_packet,
            )

        push!(candidates, (;
            refocus,
            pre_crusher,
            post_crusher,
            refocusing_center,
            first_delay,
            pre_packet_delay,
            post_packet_delay,
            timing_error,
            maximum_spacing,
        ))
    end

    best = sort(candidates; by=candidate -> (
        candidate.maximum_spacing,
        candidate.timing_error,
    ))[1]
    best.timing_error <= sys.DUR_Δt + 1e-12 || error(
        "No refocused-train timing solution centers every readout packet on a " *
        "spin echo within one block raster.")
    return best
end


function _fit_spin_echo_timing(
    excitation_tail,
    readout_center,
    slice_thickness,
    sys,
    TE;
    crusher_phase,
    refocusing_phase,
    refocusing_bandwidth,
    refocusing_time_bw_product,
    refocusing_apodization,
    min_post_delay,
)
    candidates = NamedTuple[]

    for rf_shift in 0:sys.RF_Δt:(sys.DUR_Δt - sys.RF_Δt)
        refocus = build_refocusing_block(
            slice_thickness,
            sys;
            bandwidth=refocusing_bandwidth,
            time_bw_product=refocusing_time_bw_product,
            phase=refocusing_phase,
            apodization=refocusing_apodization,
            rf_shift,
        )
        pre_crusher, post_crusher =
            _crusher_pair(crusher_phase, slice_thickness, sys, refocus)
        refocusing_center = _sequence_rf_center(refocus)

        pre_arm = excitation_tail + dur(pre_crusher) + refocusing_center
        post_arm = dur(refocus) - refocusing_center +
            dur(post_crusher) + readout_center

        if isnothing(TE)
            rastered_min_post_delay = ceil_to_raster(min_post_delay, sys.DUR_Δt)
            if pre_arm >= post_arm + rastered_min_post_delay
                target_pre_delay = 0.0
                target_post_delay = pre_arm - post_arm
            else
                target_pre_delay = post_arm + rastered_min_post_delay - pre_arm
                target_post_delay = rastered_min_post_delay
            end
        else
            target_pre_delay = TE / 2 - pre_arm
            target_post_delay = TE / 2 - post_arm
        end

        pre_delays = (
            floor_to_raster(max(target_pre_delay, 0.0), sys.DUR_Δt),
            ceil_to_raster(max(target_pre_delay, 0.0), sys.DUR_Δt),
        )
        post_delays = (
            floor_to_raster(max(target_post_delay, 0.0), sys.DUR_Δt),
            ceil_to_raster(max(target_post_delay, 0.0), sys.DUR_Δt),
        )

        for pre_delay in pre_delays, post_delay in post_delays
            target_pre_delay >= -1e-12 || continue
            target_post_delay >= -1e-12 || continue
            post_delay + 1e-12 >= min_post_delay || continue

            spin_echo_TE = 2 * (pre_arm + pre_delay)
            readout_TE = pre_arm + pre_delay + post_arm + post_delay
            error = isnothing(TE) ?
                abs(spin_echo_TE - readout_TE) :
                max(abs(spin_echo_TE - TE), abs(readout_TE - TE))
            error <= sys.DUR_Δt + 1e-12 || continue

            push!(candidates, (;
                refocus,
                pre_crusher,
                post_crusher,
                pre_delay,
                post_delay,
                actual_TE=readout_TE,
                spin_echo_TE,
                error,
            ))
        end
    end

    if isempty(candidates)
        requested = isnothing(TE) ? "the minimum TE" : "TE=$(TE) s"
        error(
            "No spin-echo timing solution places the readout center on the spin " *
            "echo within one block raster of $requested.")
    end

    best = isnothing(TE) ?
        sort(candidates; by=candidate -> (
            max(candidate.actual_TE, candidate.spin_echo_TE),
            candidate.error,
            candidate.post_delay,
        ))[1] :
        sort(candidates; by=candidate -> (
            candidate.error,
            candidate.post_delay,
        ))[1]
    return (;
        refocus=best.refocus,
        pre_crusher=best.pre_crusher,
        post_crusher=best.post_crusher,
        pre_delay=best.pre_delay,
        post_delay=best.post_delay,
        actual_TE=best.actual_TE,
    )
end


"""
    build_spin_echo(
        excitation,
        readout_kernel,
        sys,
        readout_args...;
        slice_thickness,
        TE=nothing,
        kwargs...,
    )

Build one spin echo from an existing excitation sequence and a readout kernel.
The kernel must provide `readout(args...)` and `center_time(args...)` callables.
The same positional `readout_args` are forwarded to both, so GRE encoding
indices and EPI shot numbers require no readout-specific logic here.

`epi_readout_kernel` materializes its out-of-block gradients before returning
the readout sequence, making it directly appendable here.

# Keywords
- `slice_thickness`: Refocusing slice thickness. [`m`]
- `TE=nothing`: Echo time measured from the excitation RF center. `nothing`
  selects the minimum feasible TE. [`s`]
- `crusher_phase=0.0`: Crusher phase accumulation across the slice. [`rad`]
- `refocusing_phase=π/2`: Refocusing RF phase. [`rad`]
- `refocusing_bandwidth=nothing`: Refocusing bandwidth. [`Hz`]
- `refocusing_time_bw_product=4`: Refocusing RF time-bandwidth product.
- `refocusing_apodization=0.46`: Refocusing cosine-window weight.
- `min_post_delay=0.0`: Minimum delay between the post-crusher and readout. [`s`]
"""
function build_spin_echo(
    excitation,
    readout_kernel,
    sys,
    readout_args...;
    slice_thickness,
    TE=nothing,
    crusher_phase=0.0,
    refocusing_phase=π / 2,
    refocusing_bandwidth=nothing,
    refocusing_time_bw_product=4,
    refocusing_apodization=0.46,
    min_post_delay=0.0,
)
    isnothing(TE) || TE > 0 || error("TE must be positive or nothing.")
    slice_thickness > 0 || error("slice_thickness must be positive.")
    min_post_delay >= 0 || error("min_post_delay must be non-negative.")

    readout = readout_kernel.readout(readout_args...)
    center = readout_kernel.center_time(readout_args...)
    0 <= center <= dur(readout) || error("Readout center must lie within the readout.")

    excitation_center = _sequence_rf_center(excitation)
    excitation_tail = dur(excitation) - excitation_center
    excitation_tail >= -1e-12 || error("Excitation RF center is outside its sequence.")
    excitation_tail = max(excitation_tail, 0.0)

    (; refocus, pre_crusher, post_crusher, pre_delay, post_delay, actual_TE) =
        _fit_spin_echo_timing(
            excitation_tail,
            center,
            slice_thickness,
            sys,
            TE;
            crusher_phase,
            refocusing_phase,
            refocusing_bandwidth,
            refocusing_time_bw_product,
            refocusing_apodization,
            min_post_delay,
        )

    seq = Sequence(sys)
    seq += excitation
    pre_delay > 0 && (seq += make_delay(pre_delay * u"s"))
    seq += pre_crusher
    seq += refocus
    seq += post_crusher
    post_delay > 0 && (seq += make_delay(post_delay * u"s"))
    seq += readout

    seq.DEF["TE"] = actual_TE
    return seq
end
