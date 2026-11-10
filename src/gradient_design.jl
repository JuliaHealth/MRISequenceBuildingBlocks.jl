"""
    build_moment_compensation(seq, sys; target_moment=(0.0, 0.0, 0.0))

Build simultaneous X/Y/Z trapezoids so their physical gradient areas plus
those of all blocks in `seq` equal `target_moment` [T*s/m]. Append as a
rewinder or prepend as a prewinder. The default cancels the zeroth moment.
RF-induced coherence reversals are not included in this area calculation.
"""
function build_moment_compensation(seq, sys; target_moment=(0.0, 0.0, 0.0))
    target_moment = to_SI.(target_moment)
    length(target_moment) == 3 || error("target_moment must contain X/Y/Z areas [T*s/m].")
    moments = ntuple(axis -> target_moment[axis] - sum(area, seq.GR[axis,:]; init=0.0), 3)
    compensation = Sequence(sys)
    addblock!(compensation;
        x=make_trapezoid(; area=moments[1] * γ, sys),
        y=make_trapezoid(; area=moments[2] * γ, sys),
        z=make_trapezoid(; area=moments[3] * γ, sys))
    return compensation
end

"""
    _lobe_timing(moment, sys)

Return the flat-top and ramp durations for the minimum-time trapezoid that
produces the nonnegative gradient `moment` while respecting the scanner
gradient, slew, and raster limits. Both returned values use seconds.
"""
function _lobe_timing(moment, sys)
    moment ≥ 0 || error("Gradient moment must be non-negative.")
    if moment ≤ sys.limits.Gmax^2 / sys.limits.Smax
        ramp_ideal = sqrt(moment / sys.limits.Smax)
        ramp_ceil = ceil(ramp_ideal / sys.limits.GR_Δt) * sys.limits.GR_Δt
        ramp_floor = floor(ramp_ideal / sys.limits.GR_Δt) * sys.limits.GR_Δt

        if ramp_floor > 0 && ramp_floor < ramp_ceil
            gradient = moment / (ramp_floor + sys.limits.GR_Δt)
            if gradient / ramp_floor ≤ sys.limits.Smax
                return sys.limits.GR_Δt, ramp_floor
            end
        end
        return 0.0, ramp_ceil
    else
        ramp = ceil(sys.limits.Gmax / sys.limits.Smax / sys.limits.GR_Δt) * sys.limits.GR_Δt
        flat = max(
            0.0,
            ceil((moment / sys.limits.Gmax - ramp) / sys.limits.GR_Δt) * sys.limits.GR_Δt,
        )
        return flat, ramp
    end
end
