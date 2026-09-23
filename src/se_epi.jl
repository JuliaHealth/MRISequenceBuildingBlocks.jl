using KomaMRI
using KomaMRI.PulseDesigner: build_trigger, make_delay
using Unitful


function _se_epi_shots(
    excitation,
    epi_kernel,
    sys,
    slice_thickness,
    TE;
    crusher_phase,
    refocusing_phase,
    refocusing_bandwidth,
    refocusing_time_bw_product,
    refocusing_apodization,
)
    function build_shot(shot, echo_time)
        return build_spin_echo(
            excitation,
            epi_kernel,
            sys,
            shot;
            slice_thickness,
            TE=echo_time,
            crusher_phase,
            refocusing_phase,
            refocusing_bandwidth,
            refocusing_time_bw_product,
            refocusing_apodization,
        )
    end

    shots = eachindex(epi_kernel.lines)
    isnothing(TE) || return [build_shot(shot, TE) for shot in shots]

    minimum_shots = [build_shot(shot, nothing) for shot in shots]
    length(minimum_shots) == 1 && return minimum_shots

    common_TE = ceil_to_raster(
        maximum(shot.DEF["TE"] for shot in minimum_shots) + sys.DUR_Δt,
        sys.DUR_Δt,
    )
    return [build_shot(shot, common_TE) for shot in shots]
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
        ceil_to_raster(post_shot_delay, sys.DUR_Δt)
    shot_delay = isnothing(actual_post_shot_delay) || iszero(actual_post_shot_delay) ?
        nothing : make_delay(actual_post_shot_delay * u"s")

    shots = _se_epi_shots(
        excitation,
        epi_kernel,
        sys,
        slice_thickness,
        TE;
        crusher_phase,
        refocusing_phase,
        refocusing_bandwidth,
        refocusing_time_bw_product,
        refocusing_apodization,
    )

    seq = Sequence(sys)
    for shot in shots
        if !isnothing(trigger)
            @addblock seq += trigger
        end
        @addblock seq += shot
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
