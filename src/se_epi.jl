"""
    shot = se_epi(excitation, epi_kernel, sys; slice_thickness, TE=nothing, kwargs...)

Plan SE-EPI shots as named, independently appendable sequence parts. For shot
`i`, `shot(i)` returns `delay_before_180`, `crushed_180`, `delay_before_epi`,
`epi_readout`, and `TE`. The excitation is supplied separately by the caller:

```julia
parts = shot(i)
@addblock seq += excitation + parts.delay_before_180 + parts.crushed_180 +
    parts.delay_before_epi + parts.epi_readout
```

`crushed_180` contains Z crusher → slice-selective 180° → Z crusher.
Readouts include any spoiler requested in `epi_readout_kernel`. TE targets
the EPI train midpoint, using `epi_kernel.echo_center_time`. With `TE=nothing`,
select a common minimum feasible TE across all shots. Other keywords follow
[`build_spin_echo`](@ref). Returned parts are copies, so disabling ADCs for a
dummy does not change subsequent shots. NAV, dummy ordering, and TR delays
are assembled by the caller.
"""
function se_epi(excitation, epi_kernel, sys; slice_thickness, TE=nothing, kwargs...)
    slice_thickness, TE = to_SI.((slice_thickness, TE))
    parts = _se_epi_shots(excitation, epi_kernel, sys; slice_thickness, TE, kwargs...)
    return i -> deepcopy(parts[i])
end

function _se_epi_shots(excitation, epi_kernel, sys;
    slice_thickness, TE=nothing, train_center=true, kwargs...)
    center_time = train_center ? epi_kernel.echo_center_time : epi_kernel.center_time
    train = (readout=epi_kernel.readout, center_time)
    make_parts(i, echo_time) = _spin_echo_parts(excitation, train, sys, i;
        slice_thickness, TE=echo_time, kwargs...)
    if isnothing(TE)
        minimum_parts = [make_parts(i, nothing) for i in eachindex(epi_kernel.lines)]
        if train_center
            minimum_TE = maximum(minimum_parts) do parts
                first_half = dur(excitation) - _sequence_rf_center(excitation) +
                    dur(parts.delay_before_180) + _sequence_rf_center(parts.crushed_180)
                2 * max(first_half, parts.TE - first_half)
            end
            TE = ceil_to_raster(minimum_TE, sys.limits.DUR_Δt)
        else
            # Preserve build_se_epi's original center and minimum-TE policy.
            length(minimum_parts) == 1 && return minimum_parts
            TE = ceil_to_raster(maximum(parts.TE for parts in minimum_parts) +
                sys.limits.DUR_Δt, sys.limits.DUR_Δt)
        end
    end
    return [merge(make_parts(i, TE), (;TE)) for i in eachindex(epi_kernel.lines)]
end


"""
    build_se_epi(excitation, epi_kernel, sys; slice_thickness,
        trigger_channel=nothing, post_shot_delay=nothing, kwargs...)

Build a complete single-shot or multishot spin-echo EPI acquisition. The EPI
kernel supplies all phase-encoding groups, and every group is acquired once in
the kernel's natural shot order. The kernel must come from
[`epi_readout_kernel`](@ref) with `rewind_m0=true`.

Every shot either begins with a physiological input trigger or ends with a
fixed delay block. `trigger_channel` and `post_shot_delay` are mutually
exclusive. With `TE=nothing`, a single-shot acquisition uses its minimum TE;
multishot acquisitions use a common TE one block raster above the largest
individually minimized shot TE. A requested `TE` is applied to every shot.
This preserves the original `epi_kernel.center_time` target (the line nearest
ky=0). Use [`se_epi`](@ref) for named parts timed to the EPI train midpoint.

# Keywords
- `slice_thickness`: Refocusing slice thickness and third FOV dimension. [`m`]
- `TE=nothing`: Common echo time. `nothing` derives a feasible common value. [`s`]
- `trigger_channel=nothing`: Place a `:physio1` or `:physio2` input trigger
  immediately before every shot.
- `post_shot_delay=nothing`: Append a fixed delay after every shot, including
  the final shot. Rounded up to the block raster. [`s`]
- `crusher_phase=0.0`: Crusher phase accumulation across the slice. [`rad`]
- `refocusing_phase=π/2`: Refocusing RF phase. [`rad`]
- `refocusing_bandwidth=nothing`: Refocusing bandwidth. [`Hz`]
- `refocusing_time_bw_product=4`: Refocusing RF time-bandwidth product.
- `refocusing_apodization=0.46`: Refocusing cosine-window weight.
"""
function build_se_epi(
    excitation,
    epi_kernel,
    sys;
    slice_thickness,
    TE=nothing,
    trigger_channel=nothing,
    post_shot_delay=nothing,
    crusher_phase=0.0,
    refocusing_phase=π / 2,
    refocusing_bandwidth=nothing,
    refocusing_time_bw_product=4,
    refocusing_apodization=0.46,
)
    slice_thickness, TE, post_shot_delay, crusher_phase, refocusing_phase, refocusing_bandwidth =
        to_SI.((slice_thickness, TE, post_shot_delay, crusher_phase, refocusing_phase, refocusing_bandwidth))
    for property in (:readout, :center_time, :lines, :FOV, :matrix, :rewind_m0)
        hasproperty(epi_kernel, property) ||
            error("The EPI kernel must provide `$property`.")
    end
    length(epi_kernel.FOV) == 2 && length(epi_kernel.matrix) == 2 ||
        error("SE-EPI requires a two-dimensional EPI kernel.")
    epi_kernel.rewind_m0 ||
        error("SE-EPI requires an EPI kernel with rewind_m0=true.")
    isempty(epi_kernel.lines) &&
        error("The EPI kernel contains no shots.")
    slice_thickness > 0 || error("slice_thickness must be positive.")
    isnothing(trigger_channel) == isnothing(post_shot_delay) &&
        error("Provide exactly one of trigger_channel or post_shot_delay.")
    isnothing(post_shot_delay) || post_shot_delay >= 0 ||
        error("post_shot_delay must be non-negative.")

    trigger = isnothing(trigger_channel) ? nothing :
        build_trigger(trigger_channel; sys)
    actual_post_shot_delay = isnothing(post_shot_delay) ? nothing :
        ceil_to_raster(post_shot_delay, sys.limits.DUR_Δt)
    shot_delay = isnothing(actual_post_shot_delay) || iszero(actual_post_shot_delay) ?
        nothing : make_delay(actual_post_shot_delay)

    shots = _se_epi_shots(
        excitation,
        epi_kernel,
        sys;
        slice_thickness,
        TE,
        train_center=false,
        crusher_phase,
        refocusing_phase,
        refocusing_bandwidth,
        refocusing_time_bw_product,
        refocusing_apodization,
    )

    seq = Sequence(sys)
    for parts in shots
        if !isnothing(trigger)
            @addblock seq += trigger
        end
        @addblock seq += excitation + parts.delay_before_180 + parts.crushed_180 +
            parts.delay_before_epi + parts.epi_readout
        if !isnothing(shot_delay)
            @addblock seq += shot_delay
        end
    end

    pop!(seq.DEF, "TE", nothing)
    seq.DEF["FOV"] = [epi_kernel.FOV..., Float64(slice_thickness)]
    seq.DEF["Nx"], seq.DEF["Ny"] = epi_kernel.matrix
    seq.DEF["Nz"] = 1
    return seq
end
