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

Normalize the on-axis RF integral to `flip_angle`. The pulse is lengthened
as needed to satisfy gradient, slew, and `min(B1, sys.limits.B1)` limits.
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
    all(x -> isfinite(x) && x > 0, (sys.limits.Gmax, sys.limits.Smax, sys.limits.B1)) ||
        error("Finite positive gradient, slew, and B1 limits are required.")

    spiral = spiral_base(FOV, matrix, sys; Nint=1)(0)
    reference = reverse(complex.(spiral.GR[1, 1].A, spiral.GR[2, 1].A))
    gradient_dt, rf_dt = sys.limits.GR_Δt, sys.limits.RF_Δt
    reference_duration = (length(reference) - 1) * gradient_dt
    peak_gradient = maximum(g -> max(abs(real(g)), abs(imag(g))), reference)
    peak_slew = maximum(g -> max(abs(real(g)), abs(imag(g))), diff(reference)) / gradient_dt
    stretch = max(1.0, peak_gradient / sys.limits.Gmax, sqrt(peak_slew / sys.limits.Smax))
    duration = ceil_to_raster(reference_duration * stretch, gradient_dt)

    while true
        intervals = round(Int, duration / gradient_dt)
        stretch = duration / reference_duration
        # Resample the time-reversed gradients, preserving their k-space geometry.
        gradient = map(0:intervals) do sample
            position = sample * (length(reference) - 1) / intervals
            i = min(floor(Int, position) + 1, length(reference) - 1)
            fraction = position - (i - 1)
            ((1 - fraction) * reference[i] + fraction * reference[i + 1]) / stretch
        end
        gradient[end] = 0.0im
        remaining_k = [reverse(cumsum(reverse(
            (gradient[1:end-1] + gradient[2:end]) .* (γ * gradient_dt / 2)))); 0.0im]

        samples = round(Int, duration / rf_dt)
        B1 = map(1:samples) do sample
            t = (sample - 0.5) * rf_dt
            i = min(floor(Int, t / gradient_dt) + 1, length(gradient) - 1)
            local_t = t - (i - 1) * gradient_dt
            slope = (gradient[i + 1] - gradient[i]) / gradient_dt
            k = remaining_k[i] - γ * (gradient[i] * local_t + slope * local_t^2 / 2)
            g = gradient[i] + slope * local_t
            u = π * diameter * abs(k)
            jinc = iszero(u) ? 1.0 : 2besselj1(u) / u
            window = 1 - apodization + apodization * cos(π * abs(k) / abs(first(remaining_k)))
            density = abs(g) * abs(k) / hypot(abs(k), 1 / (2π * FOV))
            jinc * window * density
        end
        integral = rf_dt * (sum(B1) - (first(B1) + last(B1)) / 2)
        B1 .*= flip_angle / (2π * γ * integral)
        peak_B1 = maximum(abs, B1)
        if peak_B1 > B1_limit
            duration = max(duration + gradient_dt,
                ceil_to_raster(duration * peak_B1 / B1_limit, gradient_dt))
            continue
        end

        rise = ceil_to_raster(max(abs(real(first(gradient))), abs(imag(first(gradient)))) /
            sys.limits.Smax, gradient_dt)
        # Explicit ramp knots preserve its area in Pulseq arbitrary gradients.
        ramp_samples = round(Int, rise / gradient_dt)
        ramp = [first(gradient) * sample / ramp_samples for sample in 0:ramp_samples-1]
        waveform = [ramp; gradient]
        gradient_delay = ceil_to_raster(sys.limits.RF_dead_time, gradient_dt)
        rf_duration = (samples - 1) * rf_dt
        rf = RF(B1, rf_duration, 0.0, gradient_delay + rise + rf_dt / 2;
            center=rf_duration, use=Excitation())
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
end
