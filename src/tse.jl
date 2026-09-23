using KomaMRI
using KomaMRI.PulseDesigner: build_trigger, make_delay
using Unitful


function _tse_view_order(Ny, echo_train_length, view_order, center_echo)
    n_shots = Ny ÷ echo_train_length
    lines = collect(-(Ny ÷ 2):(Ny - 1 - Ny ÷ 2))

    if view_order == :center_out
        ordered = sort(lines; by=line -> (abs(line), line > 0))
        return permutedims(reshape(ordered, n_shots, echo_train_length))
    end

    if iseven(echo_train_length)
        lines = circshift(lines, -cld(n_shots, 2))
    end
    order = permutedims(reshape(lines, n_shots, echo_train_length))
    zero_echo = findfirst(iszero, order)[1]
    return circshift(order, (center_echo - zero_echo, 0))
end


"""
    build_tse(excitation, line_kernel, sys; slice_thickness,
        echo_train_length, trigger_channel=nothing, TR=nothing, kwargs...)

Build a complete two-dimensional turbo spin-echo acquisition. Each shot begins
with `excitation` and acquires `echo_train_length` M0-refocused Cartesian lines,
one on each spin echo. The line kernel must come from `gre_readout_kernel` with
`rewind_m0=true`, and the echo-train length must divide `Ny`.

`view_order=:linear` distributes adjacent phase encodes across each echo number
and places `ky=0` on `center_echo`, which defaults to the middle echo.
`view_order=:center_out` assigns the central phase encodes to the first echo and
then proceeds outward; its center echo is fixed to one. Echo spacing and
effective TE are minimum feasible derived outputs stored in `seq.DEF`.
Every shot must either begin with a physiological input trigger or be padded to
a requested repetition time; `trigger_channel` and `TR` are mutually exclusive.
A supplied navigator sequence is inserted once at the start of the complete
acquisition, outside the per-shot trigger/TR loop.

# Keywords
- `slice_thickness`: Refocusing slice thickness and third FOV dimension. [`m`]
- `echo_train_length`: Number of spin echoes acquired after each excitation.
- `view_order=:linear`: `:linear` or `:center_out` phase-encoding order.
- `center_echo=nothing`: Echo receiving `ky=0`. Defaults to the middle echo for
  linear ordering and echo one for center-out ordering.
- `trigger_channel=nothing`: Place a `:physio1` or `:physio2` input trigger
  immediately before every excitation.
- `TR=nothing`: Excitation-center repetition time. Rounded up to the block
  raster and padded after every echo train. [`s`]
- `navigator=nothing`: Complete navigator sequence, normally from
  [`build_epi_navigator`](@ref), inserted once at sequence start.
- `crusher_phase=0.0`: Crusher phase accumulation across the slice on each side
  of a refocusing pulse. [`rad`]
- `refocusing_phase=π/2`: Refocusing RF phase. [`rad`]
- `refocusing_bandwidth=nothing`: Refocusing bandwidth. `nothing` selects the
  shortest pulse satisfying scanner gradient and RF-amplitude limits. [`Hz`]
- `refocusing_time_bw_product=4`: Refocusing RF time-bandwidth product.
- `refocusing_apodization=0.46`: Refocusing cosine-window weight.
"""
function build_tse(
    excitation,
    line_kernel,
    sys;
    slice_thickness,
    echo_train_length,
    view_order=:linear,
    center_echo=nothing,
    trigger_channel=nothing,
    TR=nothing,
    navigator=nothing,
    crusher_phase=0.0,
    refocusing_phase=π / 2,
    refocusing_bandwidth=nothing,
    refocusing_time_bw_product=4,
    refocusing_apodization=0.46,
)
    for property in (:readout, :center_time, :FOV, :matrix, :rewind_m0)
        hasproperty(line_kernel, property) ||
            error("The GRE line kernel must provide `$property`.")
    end
    length(line_kernel.FOV) == 2 && length(line_kernel.matrix) == 2 ||
        error("TSE requires a two-dimensional GRE line kernel.")
    line_kernel.rewind_m0 ||
        error("TSE requires a GRE line kernel with rewind_m0=true.")
    slice_thickness > 0 || error("slice_thickness must be positive.")
    echo_train_length isa Integer ||
        error("echo_train_length must be an integer.")
    echo_train_length > 0 || error("echo_train_length must be positive.")
    Ny = line_kernel.matrix[2]
    Ny % echo_train_length == 0 ||
        error("echo_train_length must divide Ny=$Ny.")
    view_order in (:linear, :center_out) ||
        error("view_order must be :linear or :center_out.")
    isnothing(trigger_channel) == isnothing(TR) &&
        error("Provide exactly one of trigger_channel or TR.")
    isnothing(TR) || TR > 0 || error("TR must be positive.")
    actual_TR = isnothing(TR) ? nothing : ceil_to_raster(TR, sys.DUR_Δt)
    trigger = isnothing(trigger_channel) ? nothing :
        build_trigger(trigger_channel; sys)

    center_echo = isnothing(center_echo) ?
        (view_order == :linear ? cld(echo_train_length, 2) : 1) :
        center_echo
    center_echo isa Integer || error("center_echo must be an integer.")
    1 <= center_echo <= echo_train_length ||
        error("center_echo must be between 1 and $echo_train_length.")
    view_order == :center_out && center_echo != 1 &&
        error("center_echo must be 1 for center-out ordering.")

    order = _tse_view_order(Ny, echo_train_length, view_order, center_echo)
    packets = map(line_kernel.readout, order)
    packet_duration = dur(first(packets))
    all(packet -> isapprox(dur(packet), packet_duration; rtol=0, atol=1e-12), packets) ||
        error("TSE requires equal-duration GRE line packets.")
    packet_center = line_kernel.center_time(first(order))
    all(line -> isapprox(
        line_kernel.center_time(line),
        packet_center;
        rtol=0,
        atol=1e-12,
    ), order) || error("TSE requires a common GRE line center time.")

    excitation_center = _sequence_rf_center(excitation)
    excitation_tail = dur(excitation) - excitation_center
    timing = _fit_refocused_train_timing(
        excitation_tail,
        packet_center,
        packet_duration,
        echo_train_length,
        slice_thickness,
        sys;
        crusher_phase,
        refocusing_phase,
        refocusing_bandwidth,
        refocusing_time_bw_product,
        refocusing_apodization,
    )

    seq = Sequence(sys)
    effective_TE = nothing
    echo_spacing = nothing
    @addblock begin
        isnothing(navigator) || (seq += navigator)
        for shot in axes(order, 2)
            shot_start = dur(seq)
            isnothing(trigger) || (seq += trigger)
            excitation_start = dur(seq)
            seq += excitation
            shot_excitation_center = excitation_start + excitation_center
            timing.first_delay > 0 &&
                (seq += make_delay(timing.first_delay * u"s"))

            refocusing_centers = Float64[]
            packet_centers = Float64[]
            for echo in axes(order, 1)
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
                seq += packets[echo, shot]
                push!(packet_centers, packet_start + packet_center)
                iszero(order[echo, shot]) &&
                    (effective_TE = packet_centers[end] - shot_excitation_center)
                if echo < echo_train_length && timing.post_packet_delay > 0
                    seq += make_delay(timing.post_packet_delay * u"s")
                end
            end

            spin_echo_error = abs(
                2refocusing_centers[1] -
                shot_excitation_center -
                packet_centers[1],
            )
            for echo in 2:echo_train_length
                spin_echo_error = max(
                    spin_echo_error,
                    abs(
                        2refocusing_centers[echo] -
                        packet_centers[echo-1] -
                        packet_centers[echo],
                    ),
                )
            end
            spin_echo_error <= sys.DUR_Δt + 1e-12 || error(
                "Assembled TSE line centers miss their spin echoes by " *
                "$(spin_echo_error) s.")

            shot_spacing = echo_train_length == 1 ?
                packet_centers[1] - shot_excitation_center :
                packet_centers[2] - packet_centers[1]
            all(
                spacing -> isapprox(spacing, shot_spacing; rtol=0, atol=1e-12),
                diff(packet_centers),
            ) || error("Assembled TSE echo spacing is not uniform.")
            if isnothing(echo_spacing)
                echo_spacing = shot_spacing
            else
                isapprox(shot_spacing, echo_spacing; rtol=0, atol=1e-12) ||
                    error("Assembled TSE echo spacing differs between shots.")
            end

            if !isnothing(actual_TR)
                shot_duration = dur(seq) - shot_start
                shot_duration <= actual_TR + 1e-12 || error(
                    "TR=$(actual_TR) s is shorter than the minimum TSE shot " *
                    "duration $(shot_duration) s.")
                shot_delay = actual_TR - shot_duration
                shot_delay > 0 && (seq += make_delay(shot_delay * u"s"))
            end
        end
    end

    isnothing(effective_TE) && error("TSE view order does not contain ky=0.")
    seq.DEF["TE"] = effective_TE
    seq.DEF["EchoSpacing"] = echo_spacing
    seq.DEF["EchoTrainLength"] = echo_train_length
    seq.DEF["ViewOrder"] = String(view_order)
    isnothing(actual_TR) || (seq.DEF["TR"] = actual_TR)
    isnothing(trigger_channel) ||
        (seq.DEF["TriggerChannel"] = String(trigger_channel))
    seq.DEF["FOV"] = [line_kernel.FOV..., Float64(slice_thickness)]
    seq.DEF["Nx"], seq.DEF["Ny"] = line_kernel.matrix
    seq.DEF["Nz"] = 1
    return seq
end
