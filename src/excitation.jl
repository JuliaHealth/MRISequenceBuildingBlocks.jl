"""
    build_hard_pulse(flip_angle, sys; B1=sys.limits.B1, kwargs...)

Build a nonselective hard pulse with the requested flip angle [rad] and a
peak amplitude at most `min(B1, sys.limits.B1)` [T]. Choose the shortest duration
on the RF raster; block padding includes the scanner dead time and ring-down.
Other keywords, such as `use`, `phase_offset`, and `freq_offset`, are passed
to KomaMRI's `build_block_pulse`.
"""
function build_hard_pulse(flip_angle, sys; B1=sys.limits.B1, kwargs...)
    flip_angle, B1 = to_SI.((flip_angle, B1))
    isfinite(flip_angle) && !iszero(flip_angle) || error("flip_angle must be finite and nonzero [rad].")
    peak_B1 = min(B1, sys.limits.B1)
    isfinite(peak_B1) && peak_B1 > 0 || error("A positive finite B1 cap is required [T].")
    duration = ceil_to_raster(abs(flip_angle) / (2π * γ * peak_B1), sys.limits.RF_Δt)
    return build_block_pulse(flip_angle; duration=duration, sys, kwargs...)
end

"""
    load_rf_excitation(rf_file, sys; moment_compensation=false)

Load B1 [T], Gx/Gy and optional Gz [T/m], and scalar duration T [s] from
JLD2. Waveforms are uniform knots over [0,T]. Preserve the supplied waveforms
and their relative timing; the effective RF center is the final knot.
With `moment_compensation=true`, return two blocks: a gradient-moment
prewinder, then the supplied RF/gradients. The total physical gradient area
is zero on each axis. The default returns only the RF/gradient block.
"""
function load_rf_excitation(rf_file, sys; moment_compensation=false)
    B1, Gx, Gy, Gz, T = jldopen(rf_file, "r") do file
        gx = read(file, "Gx")
        (read(file, "B1"), gx, read(file, "Gy"),
            haskey(file, "Gz") ? read(file, "Gz") : zeros(length(gx)),
            read(file, "T"))
    end
    T isa Real && isfinite(T) && T > 0 ||
        error("T must be a positive finite scalar duration [s].")
    for (name, waveform) in (("B1", B1), ("Gx", Gx), ("Gy", Gy), ("Gz", Gz))
        waveform isa AbstractVector && length(waveform) ≥ 2 && all(isfinite, waveform) ||
            error("$name must be a finite waveform vector with at least two knots.")
    end

    pulse_delay = ceil_to_raster(sys.limits.RF_dead_time, sys.limits.GR_Δt)
    rf = RF(B1, T, 0.0, pulse_delay; center=T, use=Excitation())
    # Pulseq measures RF dead time to the sample edge, before the first knot.
    pulse_delay += ceil_to_raster(max(sys.limits.RF_dead_time - delay(rf, sys), 0.0), sys.limits.GR_Δt)
    rf.delay = pulse_delay
    excitation = Sequence(sys)
    addblock!(excitation, rf;
        x=Grad(Gx, T, 0.0, 0.0, pulse_delay),
        y=Grad(Gy, T, 0.0, 0.0, pulse_delay),
        z=Grad(Gz, T, 0.0, 0.0, pulse_delay))
    excitation.DUR[end] = ceil_to_raster(
        max(excitation.DUR[end], dur(rf, sys)), sys.limits.DUR_Δt)
    if moment_compensation
        excitation = build_moment_compensation(excitation, sys) + excitation
    end
    check_timing(excitation, sys)
    check_hw_limits(excitation, sys)
    return excitation
end

"""
    build_pencil_beam(flip_angle, diameter, FOV, matrix, sys; B1=sys.limits.B1, apodization=0.46)

Build an analytic XY pencil-beam excitation as two blocks: an M0-compensating
prewinder, then inward spiral gradients and RF. No ADC or post-rewinder is
included. The spiral starts from KomaMRI's single-interleaf Glover design.

Inputs use radians and metres; `matrix` sets nominal spiral resolution
`FOV / matrix`. The signed RF weighting is `2J₁(u)/u`, with
`u = π * diameter * |k|`, multiplied by the Archimedean-spiral density
weight `|G|*r/sqrt(r^2 + λ^2)`, where `r=|k|` and `λ=1/(2π*FOV)`.
`k` is computed from the remaining gradient area.
`diameter` is the nominal diameter of the cylindrical small-tip profile.
Finite k-space extent and larger flip angles modify that profile.
The radial window is `1 - a + a*cos(π*|k|/kmax)`, where `a=apodization`;
the default is Hamming apodization, and zero disables the window.

Normalize the on-axis RF integral to `flip_angle`, then stretch RF and gradients
together to satisfy gradient, slew, and `min(B1, sys.limits.B1)` limits.
RF resampling preserves its integral and peak amplitude bound without iteration.
`B1` is a peak amplitude cap in tesla. RF is sampled
on the RF raster; gradient ramps, RF dead time, and ring-down are included in
the second block. Its effective RF center is the final RF sample.

Reference: Coristine and Stuber (2014), doi:10.1186/1532-429X-16-S1-W34.
"""
function build_pencil_beam(flip_angle, diameter, FOV, matrix, sys; B1=sys.limits.B1, apodization=0.46)
    flip_angle, diameter, FOV = to_SI.((flip_angle, diameter, FOV))
    B1_limit = min(to_SI(B1), sys.limits.B1)
    isfinite(B1_limit) && B1_limit > 0 || error("A positive finite B1 cap is required [T].")
    isfinite(flip_angle) && !iszero(flip_angle) || error("flip_angle must be finite and nonzero [rad].")
    isfinite(diameter) && isfinite(FOV) && 0 < diameter <= FOV ||
        error("Require 0 < diameter <= FOV [m].")
    matrix isa Integer && matrix >= 2 || error("matrix must be an integer of at least 2.")
    0 <= apodization <= 1 || error("apodization must lie in [0,1].")
    all(x -> isfinite(x) && x > 0, (sys.limits.Gmax, sys.limits.Smax)) ||
        error("Finite positive gradient and slew limits are required.")

    spiral = spiral_base(FOV, matrix, sys; Nint=1)(0)
    reference = reverse(complex.(spiral.GR[1, 1].A, spiral.GR[2, 1].A))
    gradient_dt, rf_dt = sys.limits.GR_Δt, sys.limits.RF_Δt
    reference_duration = (length(reference) - 1) * gradient_dt
    remaining_k = [reverse(cumsum(reverse(
        (reference[1:end-1] + reference[2:end]) .* (γ * gradient_dt / 2)))); 0.0im]
    kmax = abs(first(remaining_k))
    samples = round(Int, reference_duration / rf_dt)
    rf_amplitudes = map(0:samples) do sample
        t = sample * rf_dt
        i = min(floor(Int, t / gradient_dt) + 1, length(reference) - 1)
        local_t = t - (i - 1) * gradient_dt
        slope = (reference[i + 1] - reference[i]) / gradient_dt
        k = remaining_k[i] - γ * (reference[i] * local_t + slope * local_t^2 / 2)
        g = reference[i] + slope * local_t
        u = π * diameter * abs(k)
        jinc = iszero(u) ? 1.0 : 2besselj1(u) / u
        window = 1 - apodization + apodization * cos(π * abs(k) / kmax)
        density = abs(g) * abs(k) / hypot(abs(k), 1 / (2π * FOV))
        jinc * window * density
    end
    cumulative_rf = [0.0; cumsum((rf_amplitudes[1:end-1] + rf_amplitudes[2:end]) .* (rf_dt / 2))]
    integral = last(cumulative_rf)
    isfinite(integral) && !iszero(integral) || error("Pencil-beam RF has zero or nonfinite integral.")
    rf_amplitudes .*= flip_angle / (2π * γ * integral)
    cumulative_rf .*= flip_angle / (2π * γ * integral)

    peak_gradient = maximum(g -> max(abs(real(g)), abs(imag(g))), reference)
    peak_slew = maximum(g -> max(abs(real(g)), abs(imag(g))), diff(reference)) / gradient_dt
    stretch = max(1.0, peak_gradient / sys.limits.Gmax, sqrt(peak_slew / sys.limits.Smax),
        maximum(abs, rf_amplitudes) / B1_limit)
    duration = ceil_to_raster(reference_duration * stretch, gradient_dt)
    stretch = duration / reference_duration

    # Resample gradients once on the hardware raster after time stretching.
    intervals = round(Int, duration / gradient_dt)
    gradient = map(0:intervals) do sample
        position = sample * (length(reference) - 1) / intervals
        i = min(floor(Int, position) + 1, length(reference) - 1)
        fraction = position - (i - 1)
        ((1 - fraction) * reference[i] + fraction * reference[i + 1]) / stretch
    end
    gradient[end] = 0.0im

    # Average the stretched RF over each sample's trapezoidal integration cell.
    # These averages preserve both the total area and the amplitude bound.
    samples = round(Int, duration / rf_dt)
    edges = [0.0; ((0:samples-1) .+ 0.5) .* rf_dt; duration]
    rf_areas = map(edges) do t
        position = t / stretch / rf_dt
        i = min(floor(Int, position) + 1, length(rf_amplitudes) - 1)
        local_t = (position - (i - 1)) * rf_dt
        slope = (rf_amplitudes[i + 1] - rf_amplitudes[i]) / rf_dt
        cumulative_rf[i] + rf_amplitudes[i] * local_t + slope * local_t^2 / 2
    end
    rf_amplitudes = diff(rf_areas) ./ diff(edges)

    rise = ceil_to_raster(max(abs(real(first(gradient))), abs(imag(first(gradient)))) /
        sys.limits.Smax, gradient_dt)
    # Explicit ramp knots preserve its area in Pulseq arbitrary gradients.
    ramp_samples = round(Int, rise / gradient_dt)
    ramp = [first(gradient) * sample / ramp_samples for sample in 0:ramp_samples-1]
    waveform = [ramp; gradient]
    gradient_delay = ceil_to_raster(sys.limits.RF_dead_time, gradient_dt)
    rf = RF(rf_amplitudes, duration, 0.0, gradient_delay + rise + rf_dt / 2;
        center=duration, use=Excitation())
    excitation = Sequence(sys)
    addblock!(excitation, rf;
        x=Grad(real.(waveform), rise + duration, 0.0, 0.0, gradient_delay),
        y=Grad(imag.(waveform), rise + duration, 0.0, 0.0, gradient_delay))
    excitation.DUR[end] = ceil_to_raster(max(dur(excitation), dur(rf, sys)), sys.limits.DUR_Δt)
    seq = build_moment_compensation(excitation, sys) + excitation
    check_timing(seq, sys)
    check_hw_limits(seq, sys)
    return seq
end
