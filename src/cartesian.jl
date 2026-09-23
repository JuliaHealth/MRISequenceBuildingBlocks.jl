"""
    cartesian_line_order(lines, lines_per_shot; view_order=:linear)

Partition centered Cartesian line indices into acquisition shots containing at
most `lines_per_shot` lines. `view_order=:linear` preserves the supplied order;
`:center_out` sorts lines by increasing distance from `k=0`, with the negative
line first at equal distance. The final shot may contain fewer lines.
"""
function cartesian_line_order(lines, lines_per_shot; view_order=:linear)
    lines_per_shot isa Integer ||
        error("lines_per_shot must be an integer.")
    lines_per_shot > 0 || error("lines_per_shot must be positive.")
    view_order in (:linear, :center_out) ||
        error("view_order must be :linear or :center_out.")

    ordered = collect(lines)
    isempty(ordered) && error("lines must contain at least one line.")
    all(line -> line isa Integer, ordered) ||
        error("lines must contain integer indices.")
    length(unique(ordered)) == length(ordered) ||
        error("lines must contain unique indices.")
    view_order == :center_out &&
        sort!(ordered; by=line -> (abs(line), line > 0))

    return [
        ordered[first:min(first + lines_per_shot - 1, end)]
        for first in 1:lines_per_shot:length(ordered)
    ]
end
