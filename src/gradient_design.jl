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
