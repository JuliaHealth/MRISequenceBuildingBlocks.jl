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
