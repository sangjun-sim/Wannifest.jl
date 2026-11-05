module KpointMatch

using ..Validate: validate_atol, validate_points

export match_kpoints, require_same_kpoints

# Sorting x limits candidate comparisons; every coordinate still uses absolute tolerance.
function _reject_duplicates(available, order, xs, atol, label)
    close(i, j) = all(abs(available[d, i] - available[d, j]) <= atol for d in 1:3)
    for a in eachindex(order)
        for b in (a + 1):length(order)
            xs[b] - xs[a] > atol && break
            close(order[a], order[b]) && error("$label has duplicate or ambiguous k-points")
        end
    end
    return nothing
end

"""
    match_kpoints(requested, available; atol, label) -> Vector{Int}

Index in `available` (3×M) of every column of `requested` (3×N), both in the same
Cartesian units. Each requested point must match exactly once within `atol` per
coordinate; extra available points are allowed. No reciprocal-lattice
translation or symmetry operation is applied: k and k+G are different points
here, because plane waves are later paired by their Miller indices.
"""
function match_kpoints(requested, available; atol, label::AbstractString)::Vector{Int}
    validate_atol(atol)
    validate_points(requested, "requested k-points")
    validate_points(available, "$label k-points")
    order = sortperm(vec(available[1, :]))
    xs = available[1, order]
    _reject_duplicates(available, order, xs, atol, label)
    indices = Int[]
    for k in axes(requested, 2)
        first = searchsortedfirst(xs, requested[1, k] - atol)
        last = searchsortedlast(xs, requested[1, k] + atol)
        matches = [order[p] for p in first:last
                   if all(abs(available[d, order[p]] - requested[d, k]) <= atol for d in 1:3)]
        length(matches) == 1 ||
            error("$label must match requested k-point $k exactly once; found $(length(matches))")
        push!(indices, only(matches))
    end
    allunique(indices) || error("$label k-point matching is not one-to-one")
    return indices
end

"""Two k lists must agree point by point, in the same order."""
function require_same_kpoints(reference, other; atol, label::AbstractString)
    size(reference) == size(other) ||
        error("$label holds $(size(other, 2)) k-points, expected $(size(reference, 2))")
    for k in axes(reference, 2)
        all(abs(reference[d, k] - other[d, k]) <= atol for d in 1:3) ||
            error("$label k-point $k = $(other[:, k]) differs from $(reference[:, k]); the k order must match")
    end
    return nothing
end

end
