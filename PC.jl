using KomaMRI
using KomaMRI.PulseDesigner:
    build_trigger, make_delay, make_label, make_trapezoid
using Unitful


_pc_first_moment(venc) = 1 / (2γ * venc)

function _pc_bipolar_timing(first_moment, sys)
    ramp = ceil(cbrt(first_moment / (2sys.Smax)) / sys.GR_Δt) * sys.GR_Δt
    gradient = first_moment / (2ramp^2)
    gradient ≤ sys.Gmax && return (0.0, ramp)

    ramp = ceil(sys.Gmax / sys.Smax / sys.GR_Δt) * sys.GR_Δt
    flat = ceil(
        max((-3ramp + sqrt(ramp^2 + 4first_moment / sys.Gmax)) / 2, 0.0) /
        sys.GR_Δt,
    ) * sys.GR_Δt
    return (flat, ramp)
end

function _pc_velocity_module(axis, first_moment, flat, ramp, sys)
    duration = 2 * (flat + 2ramp)
    if iszero(first_moment)
        seq = Sequence(sys)
        @addblock seq += make_delay(duration * u"s")
        return seq
    end

    gradient = first_moment / ((flat + ramp) * (flat + 2ramp))
    gradient ≤ sys.Gmax || error(
        "Velocity-encoding gradient $(round(gradient * 1e3, digits=2)) mT/m " *
        "exceeds Gmax=$(round(sys.Gmax * 1e3, digits=2)) mT/m.")

    function gradient_event(amplitude)
        return make_trapezoid(;
            amplitude=amplitude * u"T/m",
            flat_time=flat * u"s",
            rise_time=ramp * u"s",
            fall_time=ramp * u"s",
            sys,
        )
    end

    seq = Sequence(sys)
    positive = gradient_event(gradient)
    negative = gradient_event(-gradient)
    if axis === :RO
        @addblock seq += (; x=positive) + (; x=negative)
    elseif axis === :PE
        @addblock seq += (; y=positive) + (; y=negative)
    elseif axis === :SS
        @addblock seq += (; z=positive) + (; z=negative)
    else
        error("velocity_axis must be :RO, :PE, or :SS.")
    end
    return seq
end

_pc_gradient_area(gradient) = gradient.A * (gradient.T + gradient.rise)
_pc_rf_phase(index, increment) = increment * (index - 1) * index / 2

function _pc_disable_adc!(seq)
    for block in eachindex(seq.ADC)
        seq.ADC[block].N > 0 || continue
        seq.ADC[block] = ADC(0, 0.0)
        empty!(seq.EXT[block])
    end
    return seq
end

function _pc_shot(
    excitation,
    readout,
    velocity_module,
    encoding,
    cardiac_bin,
    rf_index,
    rf_spoil_increment;
    flip_scale=1.0,
    write_set_label=true,
    write_phs_label=true,
)
    phase = cis(_pc_rf_phase(rf_index, rf_spoil_increment))
    shot = rf_excitation(excitation, flip_scale * phase) +
        velocity_module +
        phase * readout
    labels = Extension[]
    write_set_label &&
        pushfirst!(labels, make_label(:SET, :SET, encoding - 1))
    write_phs_label &&
        push!(labels, make_label(:SET, :PHS, cardiac_bin - 1))
    for block in eachindex(shot.ADC)
        shot.ADC[block].N > 0 || continue
        append!(shot.EXT[block], labels)
    end
    return shot
end


"""
    build_pc(FOV, matrix, sys; venc, cardiac_bins, RR, BWpp, kwargs...)

Build a beat-interleaved two-dimensional cine PC-GRE sequence.
`FOV=(readout, phase, slice_thickness)` uses metres and `matrix=(Nx, Ny)`.
Each Cartesian line is acquired in every cardiac bin with a reference
(`SET=0`) and one velocity-encoded (`SET=1`) acquisition on successive
acquisition windows. `PHS` identifies the zero-based cardiac bin and `LIN`
identifies the zero-based phase-encoding line.

The requested `RR` is divided into `cardiac_bins`, then each bin is rounded to
the nearest whole number of complete TRs. A final incomplete line group is
filled with ADC-disabled dummy shots so every triggered acquisition has the
same realized RR interval. RF spoiling remains continuous across triggers,
encoding beats, line groups, and dummy shots.

# Required keywords
- `venc`: Velocity aliasing limit. [`m/s`]
- `cardiac_bins`: Number of cardiac phases per trigger.
- `RR`: Approximate R-R interval. [`s`]
- `BWpp`: Requested ADC bandwidth per acquired readout pixel. [`Hz/pixel`]

# Optional keywords
- `velocity_axis=:SS`: Velocity-encoding direction: `:RO`, `:PE`, or `:SS`.
- `flip_angle=π/8`: Excitation flip angle. [`rad`]
- `rf_bandwidth=nothing`: Slice-selective RF bandwidth. [`Hz`]
- `spoil_phase=4π`: Slice-direction spoiler phase across one slice. [`rad`]
- `n_ramp_shots=10`: ADC-disabled linear flip-angle ramp after each trigger,
  or once at sequence start when `trigger_channel=nothing`.
- `trigger_channel=:physio1`: Pulseq physiological trigger channel. Use
  `nothing` for a continuous retrospective acquisition.
- `post_trigger_delay=0.0`: Delay included in the trigger block. [`s`]
- `view_order=:linear`: Phase-encoding order, `:linear` or `:center_out`.
- `rf_spoil_increment=deg2rad(117)`: RF-spoiling phase increment. [`rad`]
"""
function build_pc(
    FOV,
    matrix,
    sys;
    venc,
    cardiac_bins,
    RR,
    BWpp,
    velocity_axis=:SS,
    flip_angle=π / 8,
    rf_bandwidth=nothing,
    spoil_phase=4π,
    n_ramp_shots=10,
    trigger_channel=:physio1,
    post_trigger_delay=0.0,
    view_order=:linear,
    rf_spoil_increment=deg2rad(117),
    _encoding_axis=nothing,
)
    length(FOV) == 3 ||
        error("FOV must contain (readout, phase, slice_thickness).")
    length(matrix) == 2 || error("matrix must contain (Nx, Ny).")
    all(value -> value > 0, FOV) || error("FOV entries must be positive.")
    all(value -> value isa Integer && value > 0, matrix) ||
        error("matrix entries must be positive integers.")
    venc > 0 || error("venc must be positive.")
    cardiac_bins isa Integer && cardiac_bins > 0 ||
        error("cardiac_bins must be a positive integer.")
    RR > 0 || error("RR must be positive.")
    BWpp > 0 || error("BWpp must be positive.")
    velocity_axis in (:RO, :PE, :SS) ||
        error("velocity_axis must be :RO, :PE, or :SS.")
    isnothing(_encoding_axis) || _encoding_axis in (:REF, :RO, :PE, :SS) ||
        error("encoding axis must be :REF, :RO, :PE, or :SS.")
    flip_angle > 0 || error("flip_angle must be positive.")
    spoil_phase ≥ 0 || error("spoil_phase must be non-negative.")
    n_ramp_shots isa Integer && n_ramp_shots ≥ 0 ||
        error("n_ramp_shots must be a non-negative integer.")
    post_trigger_delay ≥ 0 ||
        error("post_trigger_delay must be non-negative.")
    view_order in (:linear, :center_out) ||
        error("view_order must be :linear or :center_out.")

    FOV = Tuple(Float64.(FOV))
    matrix = Tuple(Int.(matrix))
    excitation = slice_selective_sinc(
        flip_angle,
        FOV[3],
        sys;
        BW=rf_bandwidth,
    )
    fixed_area = ntuple(
        axis -> _pc_gradient_area(excitation.GR[axis, 2]),
        3,
    )
    spoiler_area = spoil_phase / (2π * γ * FOV[3])
    readout_kernel = sgre_base(
        FOV[1:2],
        matrix,
        sys,
        BWpp;
        fixed_area,
        spoil_area=(0.0, 0.0, spoiler_area),
    )

    first_moment = _pc_first_moment(venc)
    flat, ramp = _pc_bipolar_timing(first_moment, sys)
    reference_module = _pc_velocity_module(:SS, 0.0, flat, ramp, sys)
    write_set_label = isnothing(_encoding_axis)
    write_phs_label = cardiac_bins > 1
    velocity_modules = if write_set_label
        (
            reference_module,
            _pc_velocity_module(velocity_axis, first_moment, flat, ramp, sys),
        )
    elseif _encoding_axis === :REF
        (reference_module,)
    else
        (_pc_velocity_module(_encoding_axis, first_moment, flat, ramp, sys),)
    end

    ky = collect(-(matrix[2] ÷ 2):(matrix[2] - 1 - matrix[2] ÷ 2))
    readouts = Dict(line => readout_kernel(line) for line in ky)
    representative = _pc_shot(
        excitation,
        readouts[first(ky)],
        velocity_modules[1],
        1,
        1,
        1,
        rf_spoil_increment,
        write_set_label=write_set_label,
        write_phs_label=write_phs_label,
    )
    TR = dur(representative)
    lines_per_bin = max(round(Int, RR / (cardiac_bins * TR)), 1)
    actual_phase_interval = lines_per_bin * TR
    actual_RR = cardiac_bins * actual_phase_interval
    line_groups = cartesian_line_order(ky, lines_per_bin; view_order)

    trigger_duration = ceil_to_raster(post_trigger_delay, sys.DUR_Δt)
    trigger = isnothing(trigger_channel) ? nothing : build_trigger(
        trigger_channel;
        duration=trigger_duration * u"s",
        sys,
    )

    sequence = Sequence(sys)
    rf_index = 1
    @addblock begin
        for lines in line_groups, encoding in eachindex(velocity_modules)
            isnothing(trigger) || (sequence += trigger)

            if n_ramp_shots > 0 && (!isnothing(trigger) || rf_index == 1)
                ramp_start = min(1.0, (3π / 180) / flip_angle)
                for scale in range(ramp_start, 1.0; length=n_ramp_shots)
                    shot = _pc_shot(
                        excitation,
                        readouts[first(lines)],
                        reference_module,
                        1,
                        1,
                        rf_index,
                        rf_spoil_increment;
                        flip_scale=scale,
                        write_set_label,
                        write_phs_label,
                    )
                    sequence += _pc_disable_adc!(shot)
                    rf_index = rf_index + 1
                end
            end

            for cardiac_bin in 1:cardiac_bins
                for line in lines
                    sequence += _pc_shot(
                        excitation,
                        readouts[line],
                        velocity_modules[encoding],
                        encoding,
                        cardiac_bin,
                        rf_index,
                        rf_spoil_increment,
                        write_set_label=write_set_label,
                        write_phs_label=write_phs_label,
                    )
                    rf_index = rf_index + 1
                end
                for _ in (length(lines) + 1):lines_per_bin
                    shot = _pc_shot(
                        excitation,
                        readouts[last(lines)],
                        velocity_modules[encoding],
                        encoding,
                        cardiac_bin,
                        rf_index,
                        rf_spoil_increment,
                        write_set_label=write_set_label,
                        write_phs_label=write_phs_label,
                    )
                    sequence += _pc_disable_adc!(shot)
                    rf_index = rf_index + 1
                end
            end
        end
    end

    actual_BWpp = 1 / (matrix[1] * _design_cartesian_readout(
        FOV[1], matrix[1], sys, BWpp).dwell)
    sequence.DEF["FOV"] = collect(FOV)
    sequence.DEF["Nx"], sequence.DEF["Ny"] = matrix
    sequence.DEF["Nz"] = 1
    sequence.DEF["TR"] = TR
    sequence.DEF["RR"] = actual_RR
    sequence.DEF["CardiacBins"] = cardiac_bins
    sequence.DEF["LinesPerCardiacBin"] = lines_per_bin
    sequence.DEF["BandwidthPerPixel"] = actual_BWpp
    sequence.DEF["Venc"] = Float64(venc)
    sequence.DEF["VelocityEncodingAxis"] = String(
        isnothing(_encoding_axis) ? velocity_axis : _encoding_axis)
    sequence.DEF["RfSpoilIncrement"] = Float64(rf_spoil_increment)
    isnothing(trigger_channel) ||
        (sequence.DEF["TriggerChannel"] = String(trigger_channel))

    @info "PC timing" requested_RR=RR actual_RR TR cardiac_bins lines_per_bin requested_BWpp=BWpp actual_BWpp
    return sequence
end


"""
    build_pc_encoding_scans(FOV, matrix, sys; kwargs...)

Build four full-matrix PC-GRE scans with separate `REF`, `X`, `Y`, and `Z`
velocity encodings. The returned named tuple contains four sequence objects.
Because each file identifies its encoding, ADC events never carry `SET` labels.
They carry `PHS` only when more than one cardiac phase is requested.

All keywords are forwarded to [`build_pc`](@ref). Reference and velocity-
encoded scans use matched module timing. Each retrospective scan receives its
own initial ramp; triggered scans ramp after each trigger as usual.
"""
function build_pc_encoding_scans(
    FOV,
    matrix,
    sys;
    kwargs...,
)
    names = (:REF, :X, :Y, :Z)
    axes = (:REF, :RO, :PE, :SS)
    sequences = ntuple(length(axes)) do index
        build_pc(
            FOV,
            matrix,
            sys;
            kwargs...,
            _encoding_axis=axes[index],
        )
    end
    result = NamedTuple{names}(sequences)
    for (name, sequence) in pairs(result)
        sequence.DEF["Encoding"] = String(name)
    end
    return result
end
