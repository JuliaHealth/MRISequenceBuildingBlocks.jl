using KomaMRI
using KomaMRI.PulseDesigner: make_delay
using Unitful


function _grase_epi_packets(epi_kernel)
    for property in (:readout, :echo_center_time, :lines, :rewind_m0, :FOV, :matrix)
        hasproperty(epi_kernel, property) ||
            error("The EPI kernel must provide `$property`.")
    end
    epi_kernel.rewind_m0 ||
        error("GRASE requires EPI echo groups with rewind_m0=true.")

    line_groups = epi_kernel.lines
    isempty(line_groups) && error("The EPI kernel contains no echo groups.")
    group_length = length(first(line_groups))
    group_length > 0 || error("EPI echo groups must not be empty.")
    all(lines -> length(lines) == group_length, line_groups) ||
        error("GRASE requires equal-length EPI echo groups.")
    isodd(group_length) ||
        error("GRASE requires an odd number of lines in every EPI echo group.")

    zero_groups = findall(lines -> any(iszero, lines), line_groups)
    length(zero_groups) == 1 ||
        error("Exactly one EPI echo group must contain ky=0.")
    zero_group = only(zero_groups)
    findfirst(iszero, line_groups[zero_group]) == (group_length + 1) ÷ 2 ||
        error("ky=0 must be the middle line of its EPI echo group.")

    packets = [epi_kernel.readout(group) for group in eachindex(line_groups)]
    packet_duration = dur(first(packets))
    all(packet -> isapprox(dur(packet), packet_duration; rtol=0, atol=1e-12), packets) ||
        error("GRASE requires equal-duration EPI echo groups.")

    packet_center = epi_kernel.echo_center_time(1)
    all(eachindex(packets)) do group
        analytic = epi_kernel.echo_center_time(group)
        concrete = epi_kernel.echo_center_time(packets[group], group)
        isapprox(analytic, packet_center; rtol=0, atol=1e-12) &&
            isapprox(concrete, packet_center; rtol=0, atol=1e-12)
    end || error("GRASE requires a common EPI echo-group center time.")

    0 <= packet_center <= packet_duration ||
        error("The EPI echo-group center must lie within the group.")
    return (; packets, packet_center, packet_duration, zero_group)
end


"""
    build_grase(excitation, epi_kernel, sys; slice_thickness, kwargs...)

Build one two-dimensional GRASE acquisition from an excitation sequence and
the echo groups supplied by [`epi_readout_kernel`](@ref). One 180° refocusing
pulse is applied per EPI echo group, and each group's temporal center is placed
on the corresponding spin echo. Groups are acquired in the kernel's natural
order. The output definitions include the three-dimensional box
`FOV = [FOVx, FOVy, slice_thickness]` and `Nx`, `Ny`, and `Nz=1`.
A supplied navigator sequence is inserted once before the imaging excitation.

The EPI kernel must use `rewind_m0=true` and provide equal-duration groups with
the same odd number of lines. Exactly one group must contain `ky=0`, as its
middle line. The minimum feasible refocusing-train timing is derived from the
EPI bandwidth and packet timing, RF and crusher durations, scanner limits, and
event rasters. The resulting effective `TE` and `SpinEchoSpacing` are stored in
`seq.DEF`.

# Keywords
- `slice_thickness`: Refocusing slice thickness. [`m`]
- `navigator=nothing`: Complete navigator sequence, normally from
  [`build_epi_navigator`](@ref), inserted before the imaging excitation.
- `crusher_phase=0.0`: Crusher phase accumulation across the slice on each
  side of a refocusing pulse. [`rad`]
- `refocusing_phase=π/2`: Refocusing RF phase. [`rad`]
- `refocusing_bandwidth=nothing`: Refocusing bandwidth. `nothing` selects the
  shortest pulse satisfying scanner gradient and RF-amplitude limits. [`Hz`]
- `refocusing_time_bw_product=4`: Refocusing RF time-bandwidth product.
- `refocusing_apodization=0.46`: Refocusing cosine-window weight.
"""
function build_grase(
    excitation,
    epi_kernel,
    sys;
    slice_thickness,
    navigator=nothing,
    crusher_phase=0.0,
    refocusing_phase=π / 2,
    refocusing_bandwidth=nothing,
    refocusing_time_bw_product=4,
    refocusing_apodization=0.46,
)
    slice_thickness > 0 || error("slice_thickness must be positive.")
    crusher_phase >= 0 || error("crusher_phase must be non-negative.")

    (; packets, packet_center, packet_duration, zero_group) =
        _grase_epi_packets(epi_kernel)
    excitation_center = _sequence_rf_center(excitation)
    excitation_tail = dur(excitation) - excitation_center
    excitation_tail >= -1e-12 ||
        error("Excitation RF center is outside its sequence.")
    excitation_tail = max(excitation_tail, 0.0)

    timing = _fit_refocused_train_timing(
        excitation_tail,
        packet_center,
        packet_duration,
        length(packets),
        slice_thickness,
        sys;
        crusher_phase,
        refocusing_phase,
        refocusing_bandwidth,
        refocusing_time_bw_product,
        refocusing_apodization,
    )

    seq = Sequence(sys)
    refocusing_centers = Float64[]
    packet_centers = Float64[]
    @addblock begin
        isnothing(navigator) || (seq += navigator)
        excitation_start = dur(seq)
        seq += excitation
        timing.first_delay > 0 &&
            (seq += make_delay(timing.first_delay * u"s"))

        for group in eachindex(packets)
            seq += timing.pre_crusher
            refocusing_start = dur(seq)
            seq += timing.refocus
            push!(
                refocusing_centers,
                refocusing_start + timing.refocusing_center,
            )
            seq += timing.post_crusher
            timing.pre_packet_delay > 0 &&
                (seq += make_delay(timing.pre_packet_delay * u"s"))

            packet_start = dur(seq)
            seq += packets[group]
            push!(packet_centers, packet_start + packet_center)
            if group < length(packets) && timing.post_packet_delay > 0
                seq += make_delay(timing.post_packet_delay * u"s")
            end
        end
    end

    spin_echo_error = abs(
        2refocusing_centers[1] -
        (excitation_start + excitation_center) -
        packet_centers[1],
    )
    for group in 2:length(packets)
        spin_echo_error = max(
            spin_echo_error,
            abs(
                2refocusing_centers[group] -
                packet_centers[group-1] -
                packet_centers[group],
            ),
        )
    end
    spin_echo_error <= sys.DUR_Δt + 1e-12 || error(
        "Assembled GRASE packet centers miss their spin echoes by " *
        "$(spin_echo_error) s.")

    spin_echo_spacing = length(packets) == 1 ?
        packet_centers[1] - (excitation_start + excitation_center) :
        packet_centers[2] - packet_centers[1]
    all(
        spacing -> isapprox(spacing, spin_echo_spacing; rtol=0, atol=1e-12),
        diff(packet_centers),
    ) || error("Assembled GRASE spin-echo spacing is not uniform.")

    seq.DEF["TE"] =
        packet_centers[zero_group] - (excitation_start + excitation_center)
    seq.DEF["SpinEchoSpacing"] = spin_echo_spacing
    seq.DEF["FOV"] = [epi_kernel.FOV..., Float64(slice_thickness)]
    seq.DEF["Nx"], seq.DEF["Ny"] = epi_kernel.matrix
    seq.DEF["Nz"] = 1
    return seq
end
