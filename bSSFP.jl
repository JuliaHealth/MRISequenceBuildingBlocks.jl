using KomaMRI
using KomaMRI.PulseDesigner: build_trigger, make_trapezoid
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
