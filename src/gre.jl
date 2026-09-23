using KomaMRI
using KomaMRI.PulseDesigner: make_adc, make_label, make_trapezoid
using Unitful

function _readout_partial_fourier_geometry(N_ro, partial_fourier)
    0.5 ≤ partial_fourier ≤ 1.0 || error(
        "readout_partial_fourier must be between 0.5 and 1.0 (got $partial_fourier)")
    full_precenter = N_ro ÷ 2
    acquired_precenter = round(Int, 2 * (partial_fourier - 0.5) * full_precenter)
    kx_start = -acquired_precenter
    kx_stop = N_ro - 1 - full_precenter
    N_acq = kx_stop - kx_start + 1
    # Standard Cartesian indices always include kx=0. For an even full
    # acquisition this gives -N/2:N/2-1 rather than half-integer samples.
    n_precenter = -kx_start
    return (; n_precenter, N_acq)
end

function _design_cartesian_readout(FOV_ro, N_ro, sys, BWpp;
    partial_fourier=1.0,
    center_readout=false)

    (; n_precenter, N_acq) =
        _readout_partial_fourier_geometry(N_ro, partial_fourier)

    dwell_ideal = 1.0 / (N_acq * BWpp)
    dwell = max(round(Int, dwell_ideal / sys.ADC_Δt), 1) * sys.ADC_Δt
    gradient = 1.0 / (γ * dwell * FOV_ro)
    gradient > sys.Gmax && error(
        "Readout gradient=$(round(gradient*1e3,digits=2)) mT/m exceeds " *
        "Gmax=$(round(sys.Gmax*1e3,digits=2)) mT/m — reduce BWpp")

    adc_duration = dwell * (N_acq - 1)
    ramp = ceil(gradient / sys.Smax / sys.GR_Δt) * sys.GR_Δt
    adc_delay = max(ramp, dwell / 2 + sys.ADC_dead_time)
    gradient_duration = ceil((
        adc_duration +
        max(adc_delay + dwell / 2 + sys.ADC_dead_time - 2ramp, 0.0)
    ) / sys.GR_Δt) * sys.GR_Δt

    if center_readout
        block_duration = ceil_to_raster(
            max(gradient_duration + 2ramp,
                2 * (adc_delay + n_precenter * dwell)),
            sys.DUR_Δt,
        )
        adc_delay += block_duration / 2 -
            (adc_delay + n_precenter * dwell)
        gradient_duration = block_duration - 2ramp
    end

    block = Sequence(sys)
    adc = N_acq == 1 ?
        ADC(N_acq, adc_duration, adc_delay) :
        make_adc(
            N_acq,
            dwell * u"s";
            delay=(adc_delay - dwell / 2) * u"s",
            sys,
        )
    @addblock block += (
        adc,
        x=make_trapezoid(;
            amplitude=gradient * u"T/m",
            flat_time=gradient_duration * u"s",
            rise_time=ramp * u"s",
            fall_time=ramp * u"s",
            sys,
        ),
    )

    # Moment applied before the readout so the acquired samples land at their
    # requested kx positions. This is not generally half the total readout moment
    # because ADC and gradient durations are rasterized independently.
    prephaser_moment =
        gradient * (adc_delay - ramp / 2 + n_precenter * dwell)
    # Full trapezoid moment, used to design balanced rewinders and spoilers.
    total_moment = gradient * (gradient_duration + ramp)

    return (;
        block,
        N_acq,
        dwell,
        gradient,
        ramp,
        adc_duration,
        adc_delay,
        gradient_duration,
        prephaser_moment,
        total_moment,
        adc_center=adc.delay + n_precenter * dwell,
    )
end

function _design_encoding_axis(FOV, N, fixed_moment)
    ΔM = 1.0 / (γ * FOV)
    lo, hi = -(N ÷ 2), N - 1 - N ÷ 2
    # Symmetric bounds are deliberately conservative for even matrices.
    M_extent = (N ÷ 2) * ΔM
    worst_moment = max(
        abs(fixed_moment - M_extent),
        abs(fixed_moment + M_extent),
    )
    return (; FOV, N, ΔM, lo, hi, worst_moment)
end

function _check_encoding_index(index, axis, name)
    axis.lo ≤ index ≤ axis.hi || error(
        "$name index $index out of bounds [$(axis.lo), $(axis.hi)]")
end

function _set_cartesian_labels!(seq, lin; par=nothing)
    labels = [make_label(:SET, :LIN, lin)]
    isnothing(par) || push!(labels, make_label(:SET, :PAR, par))
    for block in eachindex(seq.ADC)
        seq.ADC[block].N > 0 || continue
        append!(seq.EXT[block], labels)
    end
    return seq
end


"""
    gre_readout_kernel(FOV, matrix, sys, BWpp; kwargs...)

Design a reusable one-, two-, or three-dimensional Cartesian GRE readout. The
returned named tuple contains:

- `readout(indices...)`: build one readout at the requested centered phase-
  encoding indices;
- `center_time(indices...)`: time of the `kx=0` ADC sample relative to the
  beginning of that readout;
- `FOV`, `matrix`, and `rewind_m0`: validated design geometry and rewinder
  state.

Two- and three-dimensional readouts carry zero-based `LIN` and `PAR` labels on
their ADC blocks. `FOV`, `fixed_area`, and all returned timing values use SI
units.

# Keywords
- `fixed_area=(0, 0, 0)`: Fixed x/y/z gradient moments incorporated into the
  prephaser design. [`T*s/m`]
- `readout_partial_fourier=1.0`: Fraction of readout samples acquired before
  and through the asymmetric positive side; `kx=0` is always included.
- `center_readout=false`: Symmetrically pad the readout block so `kx=0` occurs
  at its temporal midpoint.
- `rewind_m0=false`: Append a common-duration three-axis rewinder so every
  returned readout has zero gradient zeroth moment.
"""
function gre_readout_kernel(FOV, matrix, sys, BWpp;
    fixed_area=(0.0, 0.0, 0.0),
    readout_partial_fourier=1.0,
    center_readout=false,
    rewind_m0=false)

    # ── Input validation ───────────────────────────────────────────────────────
    length(FOV) == length(matrix) || error(
        "FOV and matrix must have the same length " *
        "(got $(length(FOV)) and $(length(matrix)))")
    ndim = length(FOV)
    1 ≤ ndim ≤ 3 || error("FOV and matrix must have 1, 2, or 3 elements (got $ndim)")
    all(value -> value > 0, FOV) || error("FOV entries must be positive")
    all(N -> N > 0 && isinteger(N), matrix) ||
        error("matrix entries must be positive integers")
    BWpp > 0 || error("BWpp must be positive")
    length(fixed_area) == 3 || error("fixed_area must have exactly 3 elements (x, y, z)")
    dimensions = Int.(matrix)

    Ax, Ay, Az = fixed_area
    readout = _design_cartesian_readout(
        FOV[1],
        dimensions[1],
        sys,
        BWpp;
        partial_fourier=readout_partial_fourier,
        center_readout,
    )
    RO = readout.block
    M_readout_prephaser = readout.prephaser_moment

    # ── Worst-case moment per axis ─────────────────────────────────────────────
    Mx_worst = abs(-M_readout_prephaser + Ax)
    pe = ndim ≥ 2 ? _design_encoding_axis(FOV[2], dimensions[2], Ay) : nothing
    par = ndim == 3 ? _design_encoding_axis(FOV[3], dimensions[3], Az) : nothing
    My_worst = isnothing(pe) ? abs(Ay) : pe.worst_moment
    Mz_worst = isnothing(par) ? abs(Az) : par.worst_moment

    M_vec = sqrt(Mx_worst^2 + My_worst^2 + Mz_worst^2)

    # ── Prephaser timing ───────────────────────────────────────────────────────
    T_p, ζ_p = _lobe_timing(M_vec, sys)
    inv_area = 1.0 / (T_p + ζ_p)
    G_x_pre = (-M_readout_prephaser + Ax) * inv_area
    G_z_pre = Az * inv_area
    readout_center = T_p + 2ζ_p + readout.adc_center

    Mx_rewind = M_readout_prephaser - Ax - readout.total_moment
    My_rewind_worst = isnothing(pe) ? abs(Ay) : max(
        abs(pe.lo * pe.ΔM + Ay),
        abs(pe.hi * pe.ΔM + Ay),
    )
    Mz_rewind_worst = isnothing(par) ? abs(Az) : max(
        abs(par.lo * par.ΔM + Az),
        abs(par.hi * par.ΔM + Az),
    )
    M_rewind = sqrt(Mx_rewind^2 + My_rewind_worst^2 + Mz_rewind_worst^2)
    T_r, ζ_r = rewind_m0 && M_rewind > 0 ?
        _lobe_timing(M_rewind, sys) : (0.0, sys.GR_Δt)
    inv_area_r = rewind_m0 && M_rewind > 0 ? 1.0 / (T_r + ζ_r) : 0.0

    function append_rewinder!(seq, y_moment, z_moment)
        rewind_m0 && M_rewind > 0 || return seq
        @addblock seq += (
            x=make_trapezoid(; amplitude=Mx_rewind * inv_area_r * u"T/m",
                flat_time=T_r * u"s", rise_time=ζ_r * u"s",
                fall_time=ζ_r * u"s", sys),
            y=make_trapezoid(; amplitude=y_moment * inv_area_r * u"T/m",
                flat_time=T_r * u"s", rise_time=ζ_r * u"s",
                fall_time=ζ_r * u"s", sys),
            z=make_trapezoid(; amplitude=z_moment * inv_area_r * u"T/m",
                flat_time=T_r * u"s", rise_time=ζ_r * u"s",
                fall_time=ζ_r * u"s", sys),
        )
        return seq
    end

    # ── Inner callables ────────────────────────────────────────────────────────
    function gre_1D()
        PRE = Sequence(sys)
        @addblock PRE += (
            x=make_trapezoid(; amplitude=G_x_pre * u"T/m", flat_time=T_p * u"s",
                rise_time=ζ_p * u"s", fall_time=ζ_p * u"s", sys),
            y=make_trapezoid(; amplitude=Ay * inv_area * u"T/m", flat_time=T_p * u"s",
                rise_time=ζ_p * u"s", fall_time=ζ_p * u"s", sys),
            z=make_trapezoid(; amplitude=G_z_pre * u"T/m", flat_time=T_p * u"s",
                rise_time=ζ_p * u"s", fall_time=ζ_p * u"s", sys),
        )
        seq = Sequence(sys)
        @addblock seq += PRE + RO
        return append_rewinder!(seq, -Ay, -Az)
    end

    function gre_2D(i)
        _check_encoding_index(i, pe, "PE")
        G_pe = (i * pe.ΔM + Ay) * inv_area
        PRE = Sequence(sys)
        @addblock PRE += (
            x=make_trapezoid(; amplitude=G_x_pre * u"T/m", flat_time=T_p * u"s",
                rise_time=ζ_p * u"s", fall_time=ζ_p * u"s", sys),
            y=make_trapezoid(; amplitude=G_pe * u"T/m", flat_time=T_p * u"s",
                rise_time=ζ_p * u"s", fall_time=ζ_p * u"s", sys),
            z=make_trapezoid(; amplitude=G_z_pre * u"T/m", flat_time=T_p * u"s",
                rise_time=ζ_p * u"s", fall_time=ζ_p * u"s", sys),
        )
        seq = Sequence(sys)
        @addblock seq += PRE + RO
        _set_cartesian_labels!(seq, i - pe.lo)
        return append_rewinder!(seq, -i * pe.ΔM - Ay, -Az)
    end

    function gre_3D(i, j)
        _check_encoding_index(i, pe, "PE")
        _check_encoding_index(j, par, "Partition")
        G_pe = (i * pe.ΔM + Ay) * inv_area
        G_par = (j * par.ΔM + Az) * inv_area
        PRE = Sequence(sys)
        @addblock PRE += (
            x=make_trapezoid(; amplitude=G_x_pre * u"T/m", flat_time=T_p * u"s",
                rise_time=ζ_p * u"s", fall_time=ζ_p * u"s", sys),
            y=make_trapezoid(; amplitude=G_pe * u"T/m", flat_time=T_p * u"s",
                rise_time=ζ_p * u"s", fall_time=ζ_p * u"s", sys),
            z=make_trapezoid(; amplitude=G_par * u"T/m", flat_time=T_p * u"s",
                rise_time=ζ_p * u"s", fall_time=ζ_p * u"s", sys),
        )
        seq = Sequence(sys)
        @addblock seq += PRE + RO
        _set_cartesian_labels!(seq, i - pe.lo; par=j - par.lo)
        return append_rewinder!(
            seq,
            -i * pe.ΔM - Ay,
            -j * par.ΔM - Az,
        )
    end

    build_readout = ndim == 1 ? gre_1D : ndim == 2 ? gre_2D : gre_3D
    center_time(_...) = readout_center
    return (;
        readout=build_readout,
        center_time,
        FOV=Tuple(Float64.(FOV)),
        matrix=Tuple(dimensions),
        rewind_m0,
    )
end


"""
    sgre_base(FOV, matrix, sys, BWpp; kwargs...)

Return a callable that builds one Cartesian spoiled-GRE readout, including its
phase encoders and a simultaneous three-axis rewind/spoiler lobe. The callable
accepts the same centered phase-encoding indices as [`gre_readout_kernel`](@ref).

# Keywords
- `spoil_area=(0, 0, 0)`: Requested x/y/z spoiler moments after accounting for
  the readout and phase encoders. [`T*s/m`]
- `fixed_area=(0, 0, 0)`: Fixed x/y/z moments forwarded to the GRE prephaser.
  [`T*s/m`]
- `readout_partial_fourier=1.0`: Readout partial-Fourier fraction.
"""
function sgre_base(FOV, matrix, sys, BWpp;
    spoil_area = (0.0, 0.0, 0.0),
    fixed_area=(0.0, 0.0, 0.0),
    readout_partial_fourier=1.0)

    length(spoil_area) == 3 || error("spoil_area must have exactly 3 elements (x, y, z)")

    gre_line = gre_readout_kernel(
        FOV,
        matrix,
        sys,
        BWpp;
        fixed_area,
        readout_partial_fourier,
    ).readout

    ndim = length(FOV)
    dimensions = Int.(matrix)
    readout = _design_cartesian_readout(
        FOV[1],
        dimensions[1],
        sys,
        BWpp;
        partial_fourier=readout_partial_fourier,
    )
    M_x_nom = -readout.prephaser_moment + readout.total_moment
    Sx, Sy, Sz = spoil_area
    Mx_spo = abs(Sx - M_x_nom)

    pe = ndim ≥ 2 ? _design_encoding_axis(FOV[2], dimensions[2], Sy) : nothing
    par = ndim == 3 ? _design_encoding_axis(FOV[3], dimensions[3], Sz) : nothing
    ΔMy = isnothing(pe) ? 0.0 : pe.ΔM
    ΔMz = isnothing(par) ? 0.0 : par.ΔM
    My_spo = isnothing(pe) ? abs(Sy) : pe.worst_moment
    Mz_spo = isnothing(par) ? abs(Sz) : par.worst_moment

    # ── Spoiler lobe timing ────────────────────────────────────────────────────
    M_vec_spo = sqrt(Mx_spo^2 + My_spo^2 + Mz_spo^2)
    T_s, ζ_s = M_vec_spo > 0 ? _lobe_timing(M_vec_spo, sys) : (0.0, sys.GR_Δt)
    inv_area_s = M_vec_spo > 0 ? 1.0 / (T_s + ζ_s) : 0.0
    G_x_spo = (Sx - M_x_nom) * inv_area_s      # constant across all lines

    # ── Inner callables ────────────────────────────────────────────────────────
    function sGRE_1D()
        G_y_spo = Sy * inv_area_s
        G_z_spo = Sz * inv_area_s
        seq = gre_line()
        @addblock seq += (
            x=make_trapezoid(; amplitude=G_x_spo * u"T/m", flat_time=T_s * u"s",
                rise_time=ζ_s * u"s", fall_time=ζ_s * u"s", sys),
            y=make_trapezoid(; amplitude=G_y_spo * u"T/m", flat_time=T_s * u"s",
                rise_time=ζ_s * u"s", fall_time=ζ_s * u"s", sys),
            z=make_trapezoid(; amplitude=G_z_spo * u"T/m", flat_time=T_s * u"s",
                rise_time=ζ_s * u"s", fall_time=ζ_s * u"s", sys),
        )
        return seq
    end

    function sGRE_2D(i)
        G_y_spo = (Sy - i * ΔMy) * inv_area_s
        G_z_spo = Sz * inv_area_s
        seq = gre_line(i)
        @addblock seq += (
            x=make_trapezoid(; amplitude=G_x_spo * u"T/m", flat_time=T_s * u"s",
                rise_time=ζ_s * u"s", fall_time=ζ_s * u"s", sys),
            y=make_trapezoid(; amplitude=G_y_spo * u"T/m", flat_time=T_s * u"s",
                rise_time=ζ_s * u"s", fall_time=ζ_s * u"s", sys),
            z=make_trapezoid(; amplitude=G_z_spo * u"T/m", flat_time=T_s * u"s",
                rise_time=ζ_s * u"s", fall_time=ζ_s * u"s", sys),
        )
        return seq
    end

    function sGRE_3D(i, j)
        G_y_spo = (Sy - i * ΔMy) * inv_area_s
        G_z_spo = (Sz - j * ΔMz) * inv_area_s
        seq = gre_line(i, j)
        @addblock seq += (
            x=make_trapezoid(; amplitude=G_x_spo * u"T/m", flat_time=T_s * u"s",
                rise_time=ζ_s * u"s", fall_time=ζ_s * u"s", sys),
            y=make_trapezoid(; amplitude=G_y_spo * u"T/m", flat_time=T_s * u"s",
                rise_time=ζ_s * u"s", fall_time=ζ_s * u"s", sys),
            z=make_trapezoid(; amplitude=G_z_spo * u"T/m", flat_time=T_s * u"s",
                rise_time=ζ_s * u"s", fall_time=ζ_s * u"s", sys),
        )
        return seq
    end

    return ndim == 1 ? sGRE_1D : ndim == 2 ? sGRE_2D : sGRE_3D
end
