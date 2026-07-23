using KomaMRI

const _TIMING_ATOL = 1e-12

struct _GradPiece
    t::Vector{Float64}
    A::Vector{Float64}
end

struct _QueuedPiece
    target::Int
    piece::_GradPiece
end


function _grad_knots(gr)
    if gr.A isa Number
        t = cumsum([gr.delay; gr.rise; gr.T; gr.fall])
        A = [gr.first; gr.A; gr.A; gr.last]
    elseif gr.T isa Number
        n_intervals = length(gr.A) - 1
        flat_times = n_intervals > 0 ? fill(gr.T / n_intervals, n_intervals) : Float64[]
        t = cumsum([gr.delay; gr.rise; flat_times; gr.fall])
        A = [gr.first; gr.A; gr.last]
    else
        t = cumsum([gr.delay; gr.rise; gr.T; gr.fall])
        A = [gr.first; gr.A; gr.last]
    end
    return _compact_piece(_GradPiece(Float64.(t), Float64.(A)))
end


function _compact_piece(piece)
    t = Float64[]
    A = Float64[]
    for i in eachindex(piece.t, piece.A)
        if !isempty(t) && isapprox(piece.t[i], t[end]; rtol=0, atol=_TIMING_ATOL)
            A[end] = piece.A[i]
        else
            push!(t, piece.t[i])
            push!(A, piece.A[i])
        end
    end
    return _GradPiece(t, A)
end


_is_active(piece) = any(!iszero, piece.A)


function _interp_piece(piece, x)
    t, A = piece.t, piece.A
    x < t[1] - _TIMING_ATOL && return 0.0
    x > t[end] + _TIMING_ATOL && return 0.0

    i = searchsortedlast(t, x)
    i == 0 && return A[1]
    i == length(t) && return A[end]
    isapprox(x, t[i]; rtol=0, atol=_TIMING_ATOL) && return A[i]

    lambda = (x - t[i]) / (t[i+1] - t[i])
    return (1 - lambda) * A[i] + lambda * A[i+1]
end


function _partial_grad(piece)
    piece = _compact_piece(piece)
    t, A = piece.t, piece.A
    amplitude_atol = 16eps(Float64) * max(maximum(abs, A), 1.0)
    for i in eachindex(A)
        abs(A[i]) <= amplitude_atol && (A[i] = 0.0)
    end
    length(t) >= 2 && _is_active(piece) || return Grad(0.0, 0.0)

    delay = t[1]
    T = diff(t)
    return Grad(A, T, 0.0, 0.0, delay, A[1], A[end])
end


function _clip_piece(piece, clip_start, clip_stop)
    t = Float64[0.0]
    A = Float64[_interp_piece(piece, clip_start)]

    for x in piece.t
        if clip_start < x < clip_stop
            push!(t, x - clip_start)
            push!(A, _interp_piece(piece, x))
        end
    end

    push!(t, clip_stop - clip_start)
    push!(A, _interp_piece(piece, clip_stop))
    return _compact_piece(_GradPiece(t, A))
end


function _combine_pieces(pieces)
    isempty(pieces) && return Grad(0.0, 0.0)
    length(pieces) == 1 && return _partial_grad(pieces[1])

    nknots = sum(length(piece.t) for piece in pieces)
    t = Vector{Float64}(undef, nknots)
    i = 1
    for piece in pieces, x in piece.t
        t[i] = x
        i += 1
    end
    unique!(sort!(t))

    A = Vector{Float64}(undef, length(t))
    for i in eachindex(t)
        A[i] = sum(_interp_piece(piece, t[i]) for piece in pieces)
    end
    return _partial_grad(_GradPiece(t, A))
end


function _split_gradient!(buffer, outboard, seq, block_edges, axis, block, naxes)
    gr = seq.GR[axis, block]
    piece = _grad_knots(gr)
    _is_active(piece) || return nothing

    event_start = block_edges[block] + piece.t[1]
    event_stop = block_edges[block] + piece.t[end]
    event_stop > block_edges[block] || error(
        "Gradient in block $block axis $axis ends before the block starts.")
    event_start < block_edges[block+1] || error(
        "Gradient in block $block axis $axis starts after the block ends.")
    event_start >= block_edges[1] - _TIMING_ATOL || error(
        "Gradient in block $block axis $axis starts before the sequence.")
    event_stop <= block_edges[end] + _TIMING_ATOL || error(
        "Gradient in block $block axis $axis ends after the sequence.")

    gr.delay >= 0 && piece.t[end] <= seq.DUR[block] + _TIMING_ATOL && return nothing

    outboard[axis, block] = true
    first_target = max(searchsortedlast(block_edges, event_start), 1)
    last_target = min(searchsortedfirst(block_edges, event_stop) - 1, length(block_edges) - 1)
    for target_block in first_target:last_target
        clip_start = max(event_start, block_edges[target_block])
        clip_stop = min(event_stop, block_edges[target_block+1])
        clip_stop > clip_start + _TIMING_ATOL || continue

        local_clip_start = clip_start - block_edges[block]
        local_clip_stop = clip_stop - block_edges[block]
        clipped = _clip_piece(piece, local_clip_start, local_clip_stop)
        shift = clip_start - block_edges[target_block]
        for i in eachindex(clipped.t)
            clipped.t[i] += shift
        end

        target = axis + (target_block - 1) * naxes
        push!(buffer, _QueuedPiece(target, clipped))
    end
    return nothing
end


function _collect_split_pieces(seq, block_edges, threaded)
    naxes, nblocks = size(seq.GR)
    nwork = naxes * nblocks
    outboard = fill(false, naxes, nblocks)

    if threaded && Base.Threads.nthreads() > 1 && nwork > 1
        buffers = Vector{Vector{_QueuedPiece}}(undef, nwork)
        Base.Threads.@threads :static for linear in 1:nwork
            axis = mod1(linear, naxes)
            block = div(linear - axis, naxes) + 1
            buffer = _QueuedPiece[]
            sizehint!(buffer, 2)
            _split_gradient!(buffer, outboard, seq, block_edges, axis, block, naxes)
            buffers[linear] = buffer
        end
    else
        buffers = [_QueuedPiece[]]
        sizehint!(buffers[1], nwork)
        for linear in 1:nwork
            axis = mod1(linear, naxes)
            block = div(linear - axis, naxes) + 1
            _split_gradient!(buffers[1], outboard, seq, block_edges, axis, block, naxes)
        end
    end

    queues = [_GradPiece[] for _ in 1:naxes, _ in 1:nblocks]
    for buffer in buffers, queued in buffer
        push!(queues[queued.target], queued.piece)
    end
    return queues, outboard
end


function _combine_block_axis(seq, queues, outboard, axis, block)
    pieces = queues[axis, block]
    if isempty(pieces) && !outboard[axis, block]
        piece = _grad_knots(seq.GR[axis, block])
        return _is_active(piece) ? copy(seq.GR[axis, block]) : Grad(0.0, 0.0)
    end

    if !outboard[axis, block]
        piece = _grad_knots(seq.GR[axis, block])
        _is_active(piece) && push!(pieces, piece)
    end
    return _combine_pieces(pieces)
end


function _combine_queues(seq, queues, outboard, threaded)
    naxes, nblocks = size(seq.GR)
    nwork = naxes * nblocks
    GR = Matrix{Grad}(undef, size(seq.GR))

    if threaded && Base.Threads.nthreads() > 1 && nwork > 1
        Base.Threads.@threads for linear in 1:nwork
            axis = mod1(linear, naxes)
            block = div(linear - axis, naxes) + 1
            GR[axis, block] = _combine_block_axis(seq, queues, outboard, axis, block)
        end
    else
        for linear in 1:nwork
            axis = mod1(linear, naxes)
            block = div(linear - axis, naxes) + 1
            GR[axis, block] = _combine_block_axis(seq, queues, outboard, axis, block)
        end
    end

    return GR
end


"""
    materialize_gradients(seq; threaded=Base.Threads.nthreads() > 1)

Return a new sequence with gradients clipped to block boundaries and overlapping
gradient pieces combined. RF and ADC events are copied elementwise; extension and
definition data are deep-copied. The input sequence is not modified.

Gradients may begin before or end after their nominal source block, but their
full waveform must remain within the sequence. Set `threaded=false` to force
serial processing.
"""
function materialize_gradients(seq; threaded=Base.Threads.nthreads() > 1)
    block_edges = get_block_start_times(seq)
    queues, outboard = _collect_split_pieces(seq, block_edges, threaded)
    GR = _combine_queues(seq, queues, outboard, threaded)

    return Sequence(
        GR,
        copy.(seq.RF),
        copy.(seq.ADC),
        copy(seq.DUR),
        deepcopy(seq.EXT),
        deepcopy(seq.DEF),
    )
end
