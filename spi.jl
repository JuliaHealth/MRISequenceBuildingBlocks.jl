using KomaMRI
using KomaMRI.PulseDesigner:
    build_adc, build_block_pulse, build_delay, make_label, make_trapezoid
using Unitful


_spi_axis_range(N) = -(N ÷ 2):(N - 1 - N ÷ 2)


"""
    spi_linear_order(matrix)

Return every centered Cartesian `(kx, ky, kz)` index in linear order, with
`kx` changing fastest and `kz` slowest. An even dimension `N` uses the
conventional index range `-N/2:N/2-1`.
"""
function spi_linear_order(matrix)
    length(matrix) == 3 || error("SPI requires a three-dimensional matrix.")
    all(N -> N > 0 && isinteger(N), matrix) ||
        error("Matrix dimensions must be positive integers.")
    dimensions = Int.(matrix)
    return [
        (kx, ky, kz)
        for kz in _spi_axis_range(dimensions[3])
        for ky in _spi_axis_range(dimensions[2])
        for kx in _spi_axis_range(dimensions[1])
    ]
end


function _nearest_spiral_parameter(point, turns, shell_radius)
    point_radius = sqrt(sum(abs2, point))
    point_radius == 0 && return (parameter=0.0, distance2=shell_radius^2)
    unit_point = point ./ point_radius
    latitude = asin(clamp(unit_point[3], -1, 1))
    longitude = atan(unit_point[2], unit_point[1])

    function alignment(phi)
        theta = 2 * turns * (phi + π / 2)
        return cos(phi) * cos(theta) * unit_point[1] +
            cos(phi) * sin(theta) * unit_point[2] +
            sin(phi) * unit_point[3]
    end

    winding = round(Int,
        (2 * turns * (latitude + π / 2) - longitude) / (2π))
    best_phi = latitude
    best_alignment = alignment(latitude)
    golden_ratio = (sqrt(5.0) - 1) / 2
    for turn in (winding - 1):(winding + 1)
        same_longitude = (longitude + 2π * turn) / (2 * turns) - π / 2
        lower = max(-π / 2, same_longitude - π / (2 * turns))
        upper = min(π / 2, same_longitude + π / (2 * turns))
        lower <= upper || continue
        left = upper - golden_ratio * (upper - lower)
        right = lower + golden_ratio * (upper - lower)
        left_value = alignment(left)
        right_value = alignment(right)
        for _ in 1:24
            if left_value > right_value
                upper = right
                right = left
                right_value = left_value
                left = upper - golden_ratio * (upper - lower)
                left_value = alignment(left)
            else
                lower = left
                left = right
                left_value = right_value
                right = lower + golden_ratio * (upper - lower)
                right_value = alignment(right)
            end
        end
        phi = (lower + upper) / 2
        value = alignment(phi)
        if value > best_alignment
            best_phi = phi
            best_alignment = value
        end
    end

    distance2 = point_radius^2 + shell_radius^2 -
        2 * point_radius * shell_radius * best_alignment
    return (parameter=best_phi, distance2=max(distance2, 0.0))
end


"""
    spi_center_out_order(matrix; shell_thickness=1, spiral_spacing=1)

Return every centered Cartesian `(kx, ky, kz)` index in concentric spherical
shells. Shell thickness and spiral spacing are measured in units of `Δk`.

Shell `n` contains every point with radius in
`[n * shell_thickness, (n + 1) * shell_thickness)`. A spherical spiral is
placed at the middle of the shell, and the Cartesian points are ordered by
their nearest positions along that spiral. Consecutive shells traverse the
spiral in opposite directions. Completing a shell takes precedence over
minimizing step size, so shells intersecting the corners of an anisotropic or
cubic matrix can contain larger jumps.
"""
function spi_center_out_order(
    matrix;
    shell_thickness=1.0,
    spiral_spacing=1.0,
)
    shell_thickness > 0 || error("shell_thickness must be positive.")
    spiral_spacing > 0 || error("spiral_spacing must be positive.")

    shells = Dict{Int,Vector{NTuple{3,Int}}}()
    for point in spi_linear_order(matrix)
        radius = sqrt(sum(abs2, point))
        shell = floor(Int, radius / shell_thickness + 1e-12)
        push!(get!(shells, shell, NTuple{3,Int}[]), point)
    end

    order = NTuple{3,Int}[]
    for shell in sort!(collect(keys(shells)))
        shell_radius = (shell + 0.5) * shell_thickness
        turns = max(1, ceil(Int, π * shell_radius / spiral_spacing))
        assignments = [
            let projection =
                    _nearest_spiral_parameter(point, turns, shell_radius)
                (projection.parameter, projection.distance2, point)
            end
            for point in shells[shell]
        ]
        if isodd(shell)
            sort!(assignments; by=value -> (
                value[1], value[2], value[3][3], value[3][2], value[3][1]))
        else
            sort!(assignments; by=value -> (
                -value[1], value[2], value[3][3], value[3][2], value[3][1]))
        end
        append!(order, last.(assignments))
    end
    return order
end


"""
    build_spi(FOV, matrix, kspace_order, sys; kwargs...)

Build a non-selective, three-dimensional single-point imaging sequence. Each
entry in `kspace_order` is a centered integer `(kx, ky, kz)` index. An `N×3`
matrix is also accepted, with one point per row. Repeated points and incomplete
orders are allowed.

Every TR consists of a hard excitation, simultaneous 3D phase encoding, an FID
ADC with the gradients off, and one combined rewind/spoiler gradient block. All
pre-encoding and post-ADC gradients use timing designed for their respective
largest vector moments over the supplied order. Each ADC block stores the
zero-based Cartesian location as `LIN`, `PAR`, and `SLC`. The hard-pulse
duration is the shortest even number of RF raster samples that satisfies
`sys.B1`.

# Keywords
- `flip_angle=deg2rad(13)`: Hard-pulse flip angle. [`rad`]
- `adc_samples=48`: Number of FID samples at each k-space point; must be
  divisible by four for Pulseq export.
- `adc_duration=2e-3`: Requested full ADC window; dwell is rounded to the ADC
  raster. [`s`]
- `TR=nothing`: Repetition time. `nothing` uses the minimum TR. [`s`]
- `rf_spoil_increment=deg2rad(117)`: Quadratic RF-spoiling increment. [`rad`]
- `spoil_phase=4π`: z-spoiler phase across one partition voxel. [`rad`]
"""
function build_spi(
    FOV,
    matrix,
    kspace_order,
    sys;
    flip_angle=deg2rad(13),
    adc_samples=48,
    adc_duration=2e-3,
    TR=nothing,
    rf_spoil_increment=deg2rad(117),
    spoil_phase=4π,
)
    length(FOV) == 3 || error("SPI requires a three-dimensional FOV.")
    length(matrix) == 3 || error("SPI requires a three-dimensional matrix.")
    all(value -> value > 0, FOV) || error("FOV dimensions must be positive.")
    all(N -> N > 0 && isinteger(N), matrix) ||
        error("Matrix dimensions must be positive integers.")
    flip_angle > 0 || error("flip_angle must be positive.")
    sys.B1 > 0 || error("sys.B1 must be positive.")
    adc_samples > 0 && isinteger(adc_samples) ||
        error("adc_samples must be a positive integer.")
    adc_samples % 4 == 0 || error("adc_samples must be divisible by four.")
    adc_duration > 0 || error("adc_duration must be positive.")
    isnothing(TR) || TR > 0 || error("TR must be positive or nothing.")
    spoil_phase >= 0 || error("spoil_phase must be non-negative.")

    supplied_points = if kspace_order isa AbstractMatrix
        size(kspace_order, 2) == 3 ||
            error("A matrix kspace_order must have three columns.")
        eachrow(kspace_order)
    else
        kspace_order
    end
    points = NTuple{3,Int}[]
    for point in supplied_points
        length(point) == 3 || error("Each k-space point must have three indices.")
        all(isinteger, point) || error("K-space indices must be integers.")
        push!(points, (Int(point[1]), Int(point[2]), Int(point[3])))
    end
    isempty(points) && error("kspace_order must contain at least one point.")

    bounds = ntuple(axis -> _spi_axis_range(Int(matrix[axis])), 3)
    for point in points, axis in 1:3
        point[axis] in bounds[axis] || error(
            "K-space index $(point[axis]) on axis $axis is outside " *
            "$(first(bounds[axis])):$(last(bounds[axis])).")
    end

    moment(point) = ntuple(axis -> point[axis] / (γ * FOV[axis]), 3)
    largest_moment = maximum(point -> sqrt(sum(abs2, moment(point))), points)
    encode_flat, encode_rise = largest_moment > 0 ?
        _lobe_timing(largest_moment, sys) : (0.0, 0.0)
    adc_dwell = max(
        round(Int, adc_duration / adc_samples / sys.ADC_Δt),
        1,
    ) * sys.ADC_Δt

    voxel_z = FOV[3] / matrix[3]
    spoiler_area = spoil_phase / (2π * γ * voxel_z)
    rewind_moment(point) = let encoded = moment(point)
        (-encoded[1], -encoded[2], -encoded[3] + spoiler_area)
    end
    largest_rewind_moment = maximum(
        point -> sqrt(sum(abs2, rewind_moment(point))),
        points,
    )
    rewind_flat, rewind_rise = largest_rewind_moment > 0 ?
        _lobe_timing(largest_rewind_moment, sys) : (0.0, 0.0)
    minimum_rf_duration = isfinite(sys.B1) ?
        flip_angle / (2π * γ * sys.B1) : sys.RF_Δt
    # An even number of RF samples keeps the half-duration Pulseq center offset
    # on the RF raster.
    block_pulse_duration = max(
        ceil_to_raster(minimum_rf_duration, 2sys.RF_Δt),
        2sys.RF_Δt,
    )
    # Pulseq stores a block pulse as two center samples. Offset the Koma scalar
    # event so its serialized waveform still begins after the RF dead time.
    block_pulse_delay = sys.RF_dead_time + block_pulse_duration / 2

    function gradient_block(maximum_moment, flat_time, rise_time)
        amplitude = maximum_moment / (flat_time + rise_time)
        gradient = make_trapezoid(;
            amplitude=amplitude * u"T/m",
            flat_time=flat_time * u"s",
            rise_time=rise_time * u"s",
            fall_time=rise_time * u"s",
            sys,
        )
        block = Sequence(sys)
        @addblock block += (x=gradient, y=gradient, z=gradient)
        return block
    end

    excitation = build_block_pulse(
        flip_angle * u"rad";
        duration=block_pulse_duration * u"s",
        delay=block_pulse_delay * u"s",
        use=Excitation(),
        sys,
    )
    encoding = largest_moment > 0 ?
        gradient_block(largest_moment, encode_flat, encode_rise) : nothing
    adc = build_adc(adc_samples, adc_dwell * u"s"; sys)
    rewind = largest_rewind_moment > 0 ?
        gradient_block(largest_rewind_moment, rewind_flat, rewind_rise) : nothing

    minimum_TR = dur(excitation) + dur(adc)
    largest_moment > 0 && (minimum_TR += dur(encoding))
    largest_rewind_moment > 0 && (minimum_TR += dur(rewind))
    padding = if isnothing(TR)
        0.0
    else
        target_TR = ceil_to_raster(TR, sys.DUR_Δt)
        padding = target_TR - minimum_TR
        padding >= -1e-12 || error(
            "TR=$(TR) s is shorter than the minimum TR=$(minimum_TR) s.")
        padding
    end
    padding_block = padding > 1e-12 ? build_delay(padding * u"s"; sys) : nothing

    sequence = Sequence(sys)
    @addblock for (shot, point) in enumerate(points)
        phase = mod(rf_spoil_increment * (shot - 1) * shot / 2, 2π)
        excitation.RF[1].ϕ = phase
        sequence += excitation

        if largest_moment > 0
            amplitudes = moment(point) ./ (encode_flat + encode_rise)
            for axis in 1:3
                encoding.GR[axis, 1].A = amplitudes[axis]
            end
            sequence += encoding
        end

        adc.ADC[1].ϕ = phase
        adc.EXT[1] = [
            make_label(:SET, :LIN, point[1] - first(bounds[1])),
            make_label(:SET, :PAR, point[2] - first(bounds[2])),
            make_label(:SET, :SLC, point[3] - first(bounds[3])),
        ]
        sequence += adc

        if largest_rewind_moment > 0
            amplitudes = rewind_moment(point) ./ (rewind_flat + rewind_rise)
            for axis in 1:3
                rewind.GR[axis, 1].A = amplitudes[axis]
            end
            sequence += rewind
        end
        padding > 1e-12 && (sequence += padding_block)
    end
    return sequence
end
