using KomaMRI
using KomaMRI.PulseDesigner: make_delay, make_sinc_pulse, make_trapezoid
using Unitful


"""
    build_centered_sinc_pulse(flip_angle, sys; kwargs...)

Build a slice-selective sinc excitation whose RF center lies exactly at the
middle of its first block. RF and slice-gradient events come from
`make_sinc_pulse`; this function only adds the symmetric block padding needed
around them. The slice rephaser returned by `make_sinc_pulse` is appended as a
second block. The RF use defaults to `Excitation()`, which Pulseq serializes as
`e`. Other arguments and keywords follow the Unitful `make_sinc_pulse`
interface.
"""
function build_centered_sinc_pulse(flip_angle, sys; use=Excitation(), kwargs...)
    rf, slice_gradient, slice_rephaser =
        make_sinc_pulse(flip_angle; sys, use, kwargs...)
    isnothing(slice_gradient) && error(
        "build_centered_sinc_pulse requires a slice_thickness.")

    block = Sequence(sys)
    addblock!(block, rf; z=slice_gradient)
    minimum_duration = ceil_to_raster(dur(block[end], sys), sys.DUR_Δt)
    center = rf.delay + rf_center(rf)
    half_duration = ceil_to_raster(
        max(center, minimum_duration - center),
        sys.DUR_Δt / 2,
    )

    shift = half_duration - center
    rf.delay += shift
    slice_gradient.delay += shift

    seq = Sequence(sys)
    addblock!(seq, rf; z=slice_gradient)
    seq.DUR[end] = 2 * half_duration
    addblock!(seq; z=slice_rephaser)
    return seq
end

"""
    rf_sinc(B1, T, sys, G, frequency_offset, apodization, time_bw_product)

Build a custom sinc RF pulse with simultaneous gradients on all three axes,
followed by a three-axis rephaser. Inputs use SI values: `B1` in tesla, `T` in
seconds, `G` in tesla per metre, and `frequency_offset` in hertz.

This is the low-level constructor used by [`slice_selective_sinc`](@ref). The
caller supplies a peak RF amplitude and gradient vector that already describe
the desired slice-selection geometry. The serialized Pulseq RF delay is at
least `sys.RF_dead_time` after converting from Koma's sample-center timing.
"""
function rf_sinc(B1, T, sys, G, Δf, a, TBP)
    B1 > 0 || error("B1 must be positive")
    T > 0 || error("T must be positive")
    TBP > 0 || error("TBP must be positive")
    length(G) == 3 || error("G must have exactly three components")
    0 ≤ a ≤ 1 || error("window coefficient a must be between 0 and 1")
    B1 ≤ sys.B1 || error("B1 exceeds the scanner RF-amplitude limit")
    maximum(abs, G) ≤ sys.Gmax || error("G exceeds the scanner gradient limit")
    t0 = T / TBP
    ζ = ceil(maximum(abs.(G)) / sys.Smax / sys.GR_Δt) * sys.GR_Δt
    t_rf = 0:sys.RF_Δt:T
    t = t_rf .- T / 2
    A_rf = B1 .* sinc.(t ./ t0) .* ((1 - a) .+ a .* cos.((2π .* t) ./ (TBP * t0)))
    T_rew = ceil(max((T - ζ) / 2, 0.0) / sys.GR_Δt) * sys.GR_Δt
    G_rew_amp = iszero(T_rew + ζ) ? zero.(G) : G .* (-(T + ζ) / (2 * (T_rew + ζ)))
    gradient_delay = max(0, sys.RF_dead_time - ζ)
    rf_event = RF(collect(A_rf), diff(t_rf), Δf, gradient_delay + ζ)
    rf_delay_shortfall = sys.RF_dead_time - delay(rf_event, sys)
    rf_delay_shortfall > 0 &&
        (rf_event.delay += ceil_to_raster(rf_delay_shortfall, sys.RF_Δt))
    G_ss = (
        x=make_trapezoid(; amplitude=G[1] * u"T/m", flat_time=T * u"s",
            rise_time=ζ * u"s", fall_time=ζ * u"s",
            delay=gradient_delay * u"s", sys),
        y=make_trapezoid(; amplitude=G[2] * u"T/m", flat_time=T * u"s",
            rise_time=ζ * u"s", fall_time=ζ * u"s",
            delay=gradient_delay * u"s", sys),
        z=make_trapezoid(; amplitude=G[3] * u"T/m", flat_time=T * u"s",
            rise_time=ζ * u"s", fall_time=ζ * u"s",
            delay=gradient_delay * u"s", sys),
    )
    G_rew = (
        x=make_trapezoid(; amplitude=G_rew_amp[1] * u"T/m",
            flat_time=T_rew * u"s", rise_time=ζ * u"s", fall_time=ζ * u"s", sys),
        y=make_trapezoid(; amplitude=G_rew_amp[2] * u"T/m",
            flat_time=T_rew * u"s", rise_time=ζ * u"s", fall_time=ζ * u"s", sys),
        z=make_trapezoid(; amplitude=G_rew_amp[3] * u"T/m",
            flat_time=T_rew * u"s", rise_time=ζ * u"s", fall_time=ζ * u"s", sys),
    )
    seq = Sequence(sys)
    ringdown = make_delay(
        ceil_to_raster(
            dur(rf_event) + sys.RF_ring_down_time,
            sys.DUR_Δt,
        ) * u"s",
    )
    @addblock seq += (rf_event, ringdown; G_ss...) + (; G_rew...)
    return seq
end


"""
    slice_selective_sinc(flip_angle, slice_thickness, sys; kwargs...)

Build a two-block slice-selective sinc excitation and slice rephaser. Inputs
use SI values: flip angle in radians, slice thickness in metres, bandwidth in
hertz, and frequency offset in hertz.

# Keywords
- `BW=nothing`: RF bandwidth. `nothing` selects the largest bandwidth allowed
  by both the gradient and RF-amplitude limits.
- `Δf=0.0`: RF frequency offset. [`Hz`]
- `a=0.46`: Raised-cosine window coefficient.
- `TBP=4`: RF time-bandwidth product.
"""
function slice_selective_sinc(ϕ, Δz, sys; BW=nothing, Δf=0.0, a=0.46, TBP=4)
    ϕ > 0 || error("flip angle must be positive")
    Δz > 0 || error("slice thickness must be positive")
    TBP > 0 || error("TBP must be positive")
    isnothing(BW) || BW > 0 || error("BW must be positive")
    0 ≤ a ≤ 1 || error("window coefficient a must be between 0 and 1")

    BW_Gz = γ * sys.Gmax * Δz
    BW_B1 = 2π * γ * sys.B1 / ϕ

    BW = if isnothing(BW)
        min(BW_Gz, BW_B1)
    else
        BW > BW_Gz && error(
            "BW=$(round(BW))Hz requires Gz=$(round(BW/(γ*Δz)*1e3,digits=2)) mT/m " *
            "exceeding Gmax=$(round(sys.Gmax*1e3,digits=2)) mT/m")
        BW > BW_B1 && error(
            "BW=$(round(BW))Hz requires B1=$(round(ϕ*BW/(2π*γ)*1e6,digits=2)) μT " *
            "exceeding sys.B1=$(round(sys.B1*1e6,digits=2)) μT")
        BW
    end

    T = ceil(TBP / BW / sys.GR_Δt) * sys.GR_Δt
    BW = TBP / T
    Gz = BW / (γ * Δz)
    B1 = ϕ * BW / (2π * γ)

    return rf_sinc(B1, T, sys, (0.0, 0.0, Gz), Δf, a, TBP)
end
