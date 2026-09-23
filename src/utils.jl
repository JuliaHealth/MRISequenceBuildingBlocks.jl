using KomaMRI

"""
    rf_excitation(excitation, scale; block=1, coil=1)

Return a copy of one excitation block with the selected RF event scaled and
marked as an excitation. A real `scale` changes RF amplitude; a complex scale
also changes RF phase. Gradients and other events in the block are preserved.

Only the selected block is returned. Rephasers or other blocks from a
multi-block input sequence must be handled separately by the caller.
"""
function rf_excitation(rf, scale; block=1, coil=1)
    1 ≤ block ≤ length(rf) || error("RF block $block is outside the sequence")
    1 ≤ coil ≤ size(rf.RF, 1) || error("RF coil $coil is outside the sequence")
    dur(rf.RF[coil, block]) > 0 || error("Selected block and coil contain no RF event")
    rf_event = copy(rf[block])
    rf_event.RF[coil, 1] = scale * rf_event.RF[coil, 1]
    rf_event.RF[coil, 1].use = Excitation()
    return rf_event
end
