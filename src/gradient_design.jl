using KomaMRI

"""
    _lobe_timing(moment, sys)

Return the flat-top and ramp durations for the minimum-time trapezoid that
produces the nonnegative gradient `moment` while respecting the scanner
gradient, slew, and raster limits. Both returned values use seconds.
"""
function _lobe_timing(moment, sys)
    moment ≥ 0 || error("Gradient moment must be non-negative.")
    if moment ≤ sys.Gmax^2 / sys.Smax
        ramp_ideal = sqrt(moment / sys.Smax)
        ramp_ceil = ceil(ramp_ideal / sys.GR_Δt) * sys.GR_Δt
        ramp_floor = floor(ramp_ideal / sys.GR_Δt) * sys.GR_Δt

        if ramp_floor > 0 && ramp_floor < ramp_ceil
            gradient = moment / (ramp_floor + sys.GR_Δt)
            if gradient / ramp_floor ≤ sys.Smax
                return sys.GR_Δt, ramp_floor
            end
        end
        return 0.0, ramp_ceil
    else
        ramp = ceil(sys.Gmax / sys.Smax / sys.GR_Δt) * sys.GR_Δt
        flat = max(
            0.0,
            ceil((moment / sys.Gmax - ramp) / sys.GR_Δt) * sys.GR_Δt,
        )
        return flat, ramp
    end
end
