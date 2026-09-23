using KomaMRI
using KomaMRI.PulseDesigner: build_trigger, make_label, make_trapezoid
using LinearAlgebra: cross, dot, normalize
using Unitful


"""
    bssfp_readout_kernel(FOV, matrix, sys, BWpp; fixed_area=(0, 0, 0))

Design a balanced Cartesian readout kernel for one, two, or three dimensions.
It starts from a temporally centered GRE readout and appends simultaneous
x/y/z rewinders so every phase-encoding line has the same timing.

The returned named tuple contains `readout(indices...)`, `center_time`, and a
linear vector of centered Cartesian `encodings`. Two- and three-dimensional ADC
blocks retain the `LIN` and `PAR` labels from [`gre_readout_kernel`](@ref).

`fixed_area` supplies the x/y/z moments that the symmetric pre- and post-readout
lobes must incorporate, for example to balance gradients in the excitation
block. All areas use `T*s/m`.
"""
function bssfp_readout_kernel(FOV, matrix, sys, BWpp;
    fixed_area=(0.0, 0.0, 0.0))

    gre_kernel = gre_readout_kernel(
        FOV,
        matrix,
        sys,
        BWpp;
        fixed_area,
        center_readout=true,
    )
    gre_line = gre_kernel.readout
    ndim = length(FOV)
    Ax, Ay, Az = fixed_area
    readout = _design_cartesian_readout(
        FOV[1], matrix[1], sys, BWpp; center_readout=true)

    M_x_pre = Ax - readout.prephaser_moment
    M_x_rew = Ax + readout.prephaser_moment - readout.total_moment
    Mx_worst = max(abs(M_x_pre), abs(M_x_rew))
    pe = ndim ≥ 2 ? _design_encoding_axis(FOV[2], matrix[2], Ay) : nothing
    par = ndim == 3 ? _design_encoding_axis(FOV[3], matrix[3], Az) : nothing
    My_worst = isnothing(pe) ? abs(Ay) : pe.worst_moment
    Mz_worst = isnothing(par) ? abs(Az) : par.worst_moment

    T_rew, ζ_rew = _lobe_timing(
        sqrt(Mx_worst^2 + My_worst^2 + Mz_worst^2),
        sys,
    )
    inv_rew_area = 1.0 / (T_rew + ζ_rew)
    G_x_rew = M_x_rew * inv_rew_area

    function bssfp_1D()
        seq = gre_line()
        @addblock seq += (
            x=make_trapezoid(; amplitude=G_x_rew * u"T/m", flat_time=T_rew * u"s",
                rise_time=ζ_rew * u"s", fall_time=ζ_rew * u"s", sys),
            y=make_trapezoid(; amplitude=Ay * inv_rew_area * u"T/m",
                flat_time=T_rew * u"s", rise_time=ζ_rew * u"s",
                fall_time=ζ_rew * u"s", sys),
            z=make_trapezoid(; amplitude=Az * inv_rew_area * u"T/m",
                flat_time=T_rew * u"s", rise_time=ζ_rew * u"s",
                fall_time=ζ_rew * u"s", sys),
        )
        return seq
    end

    function bssfp_2D(i)
        _check_encoding_index(i, pe, "PE")
        G_y_rew = (Ay - i * pe.ΔM) * inv_rew_area
        seq = gre_line(i)
        @addblock seq += (
            x=make_trapezoid(; amplitude=G_x_rew * u"T/m", flat_time=T_rew * u"s",
                rise_time=ζ_rew * u"s", fall_time=ζ_rew * u"s", sys),
            y=make_trapezoid(; amplitude=G_y_rew * u"T/m",
                flat_time=T_rew * u"s", rise_time=ζ_rew * u"s",
                fall_time=ζ_rew * u"s", sys),
            z=make_trapezoid(; amplitude=Az * inv_rew_area * u"T/m",
                flat_time=T_rew * u"s", rise_time=ζ_rew * u"s",
                fall_time=ζ_rew * u"s", sys),
        )
        return seq
    end

    function bssfp_3D(i, j)
        _check_encoding_index(i, pe, "PE")
        _check_encoding_index(j, par, "Partition")
        G_y_rew = (Ay - i * pe.ΔM) * inv_rew_area
        G_z_rew = (Az - j * par.ΔM) * inv_rew_area
        seq = gre_line(i, j)
        @addblock seq += (
            x=make_trapezoid(; amplitude=G_x_rew * u"T/m", flat_time=T_rew * u"s",
                rise_time=ζ_rew * u"s", fall_time=ζ_rew * u"s", sys),
            y=make_trapezoid(; amplitude=G_y_rew * u"T/m",
                flat_time=T_rew * u"s", rise_time=ζ_rew * u"s",
                fall_time=ζ_rew * u"s", sys),
            z=make_trapezoid(; amplitude=G_z_rew * u"T/m",
                flat_time=T_rew * u"s", rise_time=ζ_rew * u"s",
                fall_time=ζ_rew * u"s", sys),
        )
        return seq
    end

    build_readout = ndim == 1 ? bssfp_1D : ndim == 2 ? bssfp_2D : bssfp_3D
    encodings = if ndim == 1
        [()]
    elseif ndim == 2
        [(i,) for i in pe.lo:pe.hi]
    else
        [(i, j) for j in par.lo:par.hi for i in pe.lo:pe.hi]
    end

    return (;
        readout=build_readout,
        center_time=gre_kernel.center_time,
        encodings,
    )
end


"""
    build_cartesian_bssfp(excitation, readout_kernel, sys; kwargs...)

Build a Cartesian bSSFP sequence for one image. The readout train can be split
across heartbeats; every heartbeat consists of an optional trigger wait, a
linear flip-angle ramp, and at most `lines_per_trigger` acquired lines. RF and
ADC phase alternate by 180° on every TR and restart with the same phase at each
heartbeat. Ramp shots use that heartbeat's first encoding with ADC sampling and
acquisition labels disabled.

`excitation` must be a sequence whose first block and first RF coil contain the
pulse to repeat. Only that block is inserted; any excitation rephaser must be
represented through the readout kernel's fixed-area balance.

# Keywords
- `trigger=:physio1`: Pulseq trigger channel. Pass `nothing` to omit the trigger.
- `post_trigger_delay=0.0`: Post-trigger wait duration, rounded up to the block
  raster. [`s`]
- `n_ramp_shots=13`: Number of linear flip-angle preparation shots.
- `lines_per_trigger=nothing`: Maximum acquired lines per heartbeat. `nothing`
  puts the complete image after one trigger.
- `view_order=:linear`: For a two-dimensional kernel, acquire phase-encoding
  lines in `:linear` or `:center_out` order.
"""
function build_cartesian_bssfp(
    excitation,
    readout_kernel,
    sys;
    trigger=:physio1,
    post_trigger_delay=0.0,
    n_ramp_shots=13,
    lines_per_trigger=nothing,
    view_order=:linear,
)
    post_trigger_delay ≥ 0 || error("post_trigger_delay must be non-negative.")
    n_ramp_shots isa Integer && n_ramp_shots ≥ 1 ||
        error("n_ramp_shots must be a positive integer.")
    isempty(readout_kernel.encodings) && error("The readout train is empty.")
    isnothing(lines_per_trigger) || (
        lines_per_trigger isa Integer && lines_per_trigger ≥ 1
    ) || error("lines_per_trigger must be a positive integer or nothing.")
    view_order in (:linear, :center_out) ||
        error("view_order must be :linear or :center_out.")

    encodings = collect(readout_kernel.encodings)
    n_lines = length(encodings)
    lines_per_trigger = isnothing(lines_per_trigger) ?
        n_lines : lines_per_trigger
    if all(encoding -> length(encoding) == 1, encodings)
        line_groups = cartesian_line_order(
            only.(encodings),
            lines_per_trigger;
            view_order,
        )
        encoding_groups = [[(line,) for line in group] for group in line_groups]
    else
        view_order == :linear ||
            error("center-out ordering requires a two-dimensional kernel.")
        encoding_groups = [
            encodings[first:min(first + lines_per_trigger - 1, end)]
            for first in 1:lines_per_trigger:n_lines
        ]
    end
    length(encoding_groups) > 1 && isnothing(trigger) &&
        error("A trigger is required when the image spans multiple heartbeats.")

    trigger_block = if isnothing(trigger)
        nothing
    else
        trigger_duration = ceil_to_raster(post_trigger_delay, sys.DUR_Δt)
        build_trigger(
            trigger;
            duration=trigger_duration * u"s",
            sys,
        )
    end

    seq = Sequence(sys)
    @addblock begin
        for encoding_group in encoding_groups
            isnothing(trigger_block) || (seq += trigger_block)

            ramp_readout = readout_kernel.readout(first(encoding_group)...)
            for block in eachindex(ramp_readout.ADC)
                ramp_readout.ADC[block].N > 0 || continue
                ramp_readout.ADC[block] = ADC(0, 0.0)
                empty!(ramp_readout.EXT[block])
            end

            for shot in 1:n_ramp_shots
                phase = cispi(shot - 1)
                scale = shot / n_ramp_shots
                seq += rf_excitation(excitation, scale * phase) +
                    phase * ramp_readout
            end

            for (shot, encoding) in enumerate(encoding_group)
                phase = cispi(n_ramp_shots + shot - 1)
                seq += rf_excitation(excitation, phase) +
                    phase * readout_kernel.readout(encoding...)
            end
        end
    end

    return seq
end


function _bssfp_disable_adc!(seq)
    for block in eachindex(seq.ADC)
        seq.ADC[block].N > 0 || continue
        seq.ADC[block] = ADC(0, 0.0)
        empty!(seq.EXT[block])
    end
    return seq
end

function _bssfp_cine_shot(
    excitation,
    readout,
    cardiac_bin,
    rf_index;
    flip_scale=1.0,
)
    phase = cispi(rf_index - 1)
    shot = rf_excitation(excitation, flip_scale * phase) + phase * readout
    phase_label = make_label(:SET, :PHS, cardiac_bin - 1)
    for block in eachindex(shot.ADC)
        shot.ADC[block].N > 0 || continue
        push!(shot.EXT[block], phase_label)
    end
    return shot
end


"""
    build_cine_bssfp(excitation, readout_kernel, sys; cardiac_bins, RR=nothing, kwargs...)

Build a two- or three-dimensional Cartesian CINE bSSFP sequence. Every spatial
encoding is acquired in every cardiac bin. The requested `RR` is divided into
`cardiac_bins`, and each bin is rounded to the nearest whole number of complete
TRs. An incomplete final encoding group is padded with ADC-disabled dummy TRs.
`PHS` identifies the zero-based cardiac bin; the readout kernel supplies `LIN`
and, for a three-dimensional acquisition, `PAR`.

Preparation follows the CINE PC builders: a linear ADC-disabled flip-angle ramp
is followed by approximately `steady_state_duration` of full-flip ADC-disabled
shots. Triggered acquisitions repeat preparation after every trigger;
retrospective acquisitions prepare only once at sequence start. RF and receiver
phase alternate by 180 degrees continuously through preparation, acquisition,
dummy shots, and triggers.

`excitation` must be a sequence whose first block and first RF coil contain the
pulse to repeat. Only that block is inserted; any excitation rephaser must be
represented through the readout kernel's fixed-area balance.

# Required keyword
- `cardiac_bins`: Number of cardiac phases per R-R interval.

# Optional keywords
- `RR=nothing`: Approximate R-R interval. Required when `cardiac_bins > 1`;
  omit it for an unpadded single-phase acquisition. [`s`]
- `trigger_channel=:physio1`: Pulseq physiological trigger channel. Use
  `nothing` for continuous retrospective acquisition.
- `post_trigger_delay=0.0`: Delay included in the trigger block. [`s`]
- `n_ramp_shots=13`: ADC-disabled linear flip-angle ramp from `FA/13` through
  full flip angle.
- `steady_state_duration=0.3`: Approximate full-flip preparation following the
  ramp. [`s`]
- `view_order=:linear`: Phase-encoding order, `:linear` or `:center_out` in 2D.
  Three-dimensional filling is linear with `ky` varying fastest.
- `encoding_order=nothing`: Exact acquisition order as a complete permutation of
  the readout kernel's centered encoding tuples. Used in every cardiac bin;
  requires `view_order=:linear`.
"""
function build_cine_bssfp(
    excitation,
    readout_kernel,
    sys;
    cardiac_bins,
    RR=nothing,
    trigger_channel=:physio1,
    post_trigger_delay=0.0,
    n_ramp_shots=13,
    steady_state_duration=0.3,
    view_order=:linear,
    encoding_order=nothing,
)
    cardiac_bins isa Integer && cardiac_bins > 0 ||
        error("cardiac_bins must be a positive integer.")
    isnothing(RR) || RR > 0 || error("RR must be positive.")
    isnothing(RR) && cardiac_bins > 1 &&
        error("RR is required when cardiac_bins is greater than one.")
    post_trigger_delay >= 0 ||
        error("post_trigger_delay must be non-negative.")
    n_ramp_shots isa Integer && n_ramp_shots >= 0 ||
        error("n_ramp_shots must be a non-negative integer.")
    steady_state_duration >= 0 ||
        error("steady_state_duration must be non-negative.")
    view_order in (:linear, :center_out) ||
        error("view_order must be :linear or :center_out.")
    isnothing(encoding_order) || view_order == :linear ||
        error("view_order must be :linear when encoding_order is supplied.")

    kernel_encodings = collect(readout_kernel.encodings)
    spatial_encodings = isnothing(encoding_order) ?
        kernel_encodings : collect(encoding_order)
    isempty(spatial_encodings) && error("The readout train is empty.")
    if !isnothing(encoding_order)
        length(spatial_encodings) == length(kernel_encodings) &&
            Set(spatial_encodings) == Set(kernel_encodings) ||
            error("encoding_order must be a complete permutation of the readout encodings.")
    end
    encoding_dims = length(first(spatial_encodings))
    encoding_dims in (1, 2) ||
        error("CINE bSSFP requires a two- or three-dimensional readout kernel.")
    all(encoding -> length(encoding) == encoding_dims, spatial_encodings) ||
        error("Readout encodings must have a consistent dimensionality.")
    encoding_dims == 2 && view_order != :linear &&
        error("three-dimensional CINE bSSFP currently supports only linear view ordering.")

    readouts = Dict(
        encoding => readout_kernel.readout(encoding...)
        for encoding in spatial_encodings
    )
    representative = _bssfp_cine_shot(
        excitation,
        readouts[first(spatial_encodings)],
        1,
        1,
    )
    TR = dur(representative)
    steady_state_shots = steady_state_duration > 0 ?
        max(round(Int, steady_state_duration / TR), 1) : 0
    lines_per_bin = isnothing(RR) ? length(spatial_encodings) :
        max(round(Int, RR / (cardiac_bins * TR)), 1)
    actual_RR = isnothing(RR) ? nothing : cardiac_bins * lines_per_bin * TR
    encoding_groups = if encoding_dims == 1
        line_groups = cartesian_line_order(
            only.(spatial_encodings),
            lines_per_bin;
            view_order,
        )
        [[(line,) for line in group] for group in line_groups]
    else
        [
            spatial_encodings[first:min(first + lines_per_bin - 1, end)]
            for first in 1:lines_per_bin:length(spatial_encodings)
        ]
    end

    trigger_duration = ceil_to_raster(post_trigger_delay, sys.DUR_Δt)
    trigger = isnothing(trigger_channel) ? nothing : build_trigger(
        trigger_channel;
        duration=trigger_duration * u"s",
        sys,
    )

    sequence = Sequence(sys)
    rf_index = 1
    lead_time = nothing
    @addblock begin
        for encodings in encoding_groups
            isnothing(trigger) || (sequence += trigger)

            prepare = (n_ramp_shots > 0 || steady_state_shots > 0) &&
                (!isnothing(trigger) || rf_index == 1)
            if prepare
                for ramp_shot in 1:n_ramp_shots
                    shot = _bssfp_cine_shot(
                        excitation,
                        readouts[first(encodings)],
                        1,
                        rf_index;
                        flip_scale=ramp_shot / n_ramp_shots,
                    )
                    sequence += _bssfp_disable_adc!(shot)
                    rf_index = rf_index + 1
                end
                for _ in 1:steady_state_shots
                    shot = _bssfp_cine_shot(
                        excitation,
                        readouts[first(encodings)],
                        1,
                        rf_index,
                    )
                    sequence += _bssfp_disable_adc!(shot)
                    rf_index = rf_index + 1
                end
            end

            for cardiac_bin in 1:cardiac_bins
                for spatial_encoding in encodings
                    isnothing(lead_time) && (lead_time = dur(sequence))
                    sequence += _bssfp_cine_shot(
                        excitation,
                        readouts[spatial_encoding],
                        cardiac_bin,
                        rf_index,
                    )
                    rf_index = rf_index + 1
                end
                for _ in (length(encodings) + 1):lines_per_bin
                    shot = _bssfp_cine_shot(
                        excitation,
                        readouts[last(encodings)],
                        cardiac_bin,
                        rf_index,
                    )
                    sequence += _bssfp_disable_adc!(shot)
                    rf_index = rf_index + 1
                end
            end
        end
    end

    sequence.DEF["TR"] = TR
    isnothing(actual_RR) || (sequence.DEF["RR"] = actual_RR)
    sequence.DEF["CardiacBins"] = cardiac_bins
    sequence.DEF["LinesPerCardiacBin"] = lines_per_bin
    sequence.DEF["LeadTime"] = lead_time
    isnothing(trigger_channel) ||
        (sequence.DEF["TriggerChannel"] = String(trigger_channel))

    @info "CINE bSSFP timing" requested_RR=RR actual_RR TR cardiac_bins lines_per_bin
    return sequence
end


"""
    build_cine_bssfp(FOV, matrix, sys, BWpp; excitation, kwargs...)

Build a two- or three-dimensional Cartesian CINE bSSFP sequence from an
already-designed slice- or slab-selective excitation. `FOV` contains the
readout, phase, and encoded slice/slab dimensions in metres. The excitation's
second-block z gradient is included in the balanced readout design.

`slice_normal` is a scanner-coordinate normal; `slice_shift` is a signed
distance in metres along that normal. The sequence is rotated from a
right-handed logical readout/phase/slice basis into scanner coordinates. Other
keywords, including `encoding_order`, pass to the lower-level cine builder.
The supplied excitation is not modified.
"""
function build_cine_bssfp(
    FOV,
    matrix,
    sys,
    BWpp;
    excitation,
    slice_normal=(0.0, 0.0, 1.0),
    slice_shift=0.0,
    kwargs...,
)
    length(FOV) == 3 || error("FOV must have readout, phase, and slice dimensions.")
    length(matrix) in (2, 3) || error("matrix must be two- or three-dimensional.")
    length(slice_normal) == 3 || error("slice_normal must have three components.")
    length(excitation) >= 2 || error("excitation must include a slice rephaser block.")

    normal = collect(Float64.(slice_normal))
    all(isfinite, normal) && sum(abs2, normal) > 0 ||
        error("slice_normal must be a finite, nonzero vector.")
    normal = normalize(normal)
    reference_axis = zeros(3)
    reference_axis[argmin(abs.(normal))] = 1.0
    readout_direction = normalize(reference_axis - dot(reference_axis, normal) * normal)
    phase_direction = cross(normal, readout_direction)
    orientation = hcat(readout_direction, phase_direction, normal)

    shifted_excitation = deepcopy(excitation)
    shifted_excitation.RF[1].Δf += γ * shifted_excitation.GR[3, 1].A * slice_shift
    slice_rephaser_area = area(shifted_excitation.GR[3, 2])
    readout_FOV = length(matrix) == 2 ? FOV[1:2] : FOV
    readout = bssfp_readout_kernel(readout_FOV, matrix, sys, BWpp;
        fixed_area=(0.0, 0.0, slice_rephaser_area))
    sequence = orientation * build_cine_bssfp(shifted_excitation, readout, sys; kwargs...)

    sequence.DEF["FOV"] = collect(FOV)
    sequence.DEF["Nx"] = matrix[1]
    sequence.DEF["Ny"] = matrix[2]
    sequence.DEF["Nz"] = length(matrix) == 3 ? matrix[3] : 1
    sequence.DEF["BandwidthPerPixel"] = BWpp
    sequence.DEF["SliceNormal"] = normal
    sequence.DEF["SliceShift"] = slice_shift
    sequence.DEF["ReadoutDirection"] = readout_direction
    sequence.DEF["PhaseDirection"] = phase_direction
    return sequence
end
