module OrbitalProjection

export ProjectionGroup, ProjectionSpec, weights_for_eigenvectors

struct ProjectionGroup
    label::String
    indices::Vector{Int}
    color::String
    atom_count::Int

    function ProjectionGroup(
        label::AbstractString,
        indices::Vector{Int},
        color::AbstractString,
        atom_count::Integer=1,
    )
        clean_indices = sort!(copy(indices))
        count = Int(atom_count)
        count > 0 || error("projection group atom_count must be positive")
        return new(String(label), clean_indices, String(color), count)
    end
end

struct ProjectionSpec
    groups::Vector{ProjectionGroup}
    num_wann::Int
    disjoint::Bool
    covers_all::Bool

    function ProjectionSpec(groups::Vector{ProjectionGroup}, num_wann::Integer)
        nw = Int(num_wann)
        nw > 0 || error("ProjectionSpec requires positive num_wann")
        isempty(groups) && error("ProjectionSpec requires at least one group")

        all_idx = Int[]
        labels = String[]
        for group in groups
            isempty(group.label) && error("projection group label cannot be empty")
            push!(labels, group.label)
            isempty(group.indices) && error("group '$(group.label)' has empty indices")
            any(i -> i < 1 || i > nw, group.indices) &&
                error("group '$(group.label)' has index out of [1, $nw]")
            length(unique(group.indices)) == length(group.indices) ||
                error("group '$(group.label)' has duplicate indices")
            append!(all_idx, group.indices)
        end

        length(unique(labels)) == length(labels) || error("duplicate group labels")
        unique_idx = unique(all_idx)
        disjoint = length(unique_idx) == length(all_idx)
        covers_all = length(unique_idx) == nw
        return new(groups, nw, disjoint, covers_all)
    end
end

function weights_for_eigenvectors(spec::ProjectionSpec, evecs::AbstractMatrix{<:Complex})
    size(evecs, 1) == spec.num_wann ||
        error("evecs has $(size(evecs, 1)) rows, expected num_wann=$(spec.num_wann)")

    nbands = size(evecs, 2)
    ngroups = length(spec.groups)
    weights = zeros(Float64, nbands, ngroups)
    @inbounds for (ig, group) in enumerate(spec.groups)
        for ib in 1:nbands
            weights[ib, ig] = sum(abs2, view(evecs, group.indices, ib))
        end
    end
    return weights
end

end
