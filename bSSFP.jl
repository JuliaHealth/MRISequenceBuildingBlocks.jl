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

Build a Cartesian bSSFP sequence consisting of an optional trigger wait, a
linear flip-angle ramp, and the kernel's complete linear readout train. RF and
ADC phase alternate by 180° on every shot. Ramp shots use the first Cartesian
encoding with ADC sampling and acquisition labels disabled.

`excitation` must be a sequence whose first block and first RF coil contain the
pulse to repeat. Only that block is inserted; any excitation rephaser must be
represented through the readout kernel's fixed-area balance.

# Keywords
- `trigger=:physio1`: Pulseq trigger channel. Pass `nothing` to omit the trigger.
- `post_trigger_delay=0.0`: Post-trigger wait duration, rounded up to the block
  raster. [`s`]
- `n_ramp_shots=13`: Number of linear flip-angle preparation shots.
"""
function build_cartesian_bssfp(
    excitation,
    readout_kernel,
    sys;
    trigger=:physio1,
    post_trigger_delay=0.0,
    n_ramp_shots=13,
)
    post_trigger_delay ≥ 0 || error("post_trigger_delay must be non-negative.")
    n_ramp_shots isa Integer && n_ramp_shots ≥ 1 ||
        error("n_ramp_shots must be a positive integer.")
    isempty(readout_kernel.encodings) && error("The readout train is empty.")

    seq = Sequence(sys)
    if !isnothing(trigger)
        trigger_duration = ceil_to_raster(post_trigger_delay, sys.DUR_Δt)
        seq += build_trigger(
            trigger;
            duration=trigger_duration * u"s",
            sys,
        )
    end

    first_encoding = first(readout_kernel.encodings)
    ramp_readout = readout_kernel.readout(first_encoding...)
    for block in eachindex(ramp_readout.ADC)
        ramp_readout.ADC[block].N > 0 || continue
        ramp_readout.ADC[block] = ADC(0, 0.0)
        empty!(ramp_readout.EXT[block])
    end

    for shot in 1:n_ramp_shots
        phase = cispi(shot - 1)
        scale = shot / n_ramp_shots
        seq += rf_excitation(excitation, scale * phase) + phase * ramp_readout
    end

    for (shot, encoding) in enumerate(readout_kernel.encodings)
        phase = cispi(n_ramp_shots + shot - 1)
        seq += rf_excitation(excitation, phase) +
            phase * readout_kernel.readout(encoding...)
    end

    return seq
end
