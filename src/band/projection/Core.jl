module Projection

using ..LocalAxisRotation
using ..Model: RunConfig
using ..OrbitalProjection: ProjectionGroup, ProjectionSpec, weights_for_eigenvectors
using ..ProjectionModel: BandProjectionResult, ProjectedDosResult
using ..SpinLayout
using ..Win90Basis

export build_projection_spec, build_basis_rotation_transform
export projection_weights, band_projection_result, projected_dos_result

function _canonicalize_projection_spec(spec::ProjectionSpec, layout::Symbol)::ProjectionSpec
    index_map = SpinLayout.source_to_canonical_indices(spec.num_wann, layout)
    index_map == collect(1:spec.num_wann) && return spec
    groups = ProjectionGroup[
        ProjectionGroup(
            group.label,
            [index_map[i] for i in group.indices],
            group.color,
            group.atom_count,
        )
        for group in spec.groups
    ]
    return ProjectionSpec(groups, spec.num_wann)
end

function _canonicalize_basis_entries(entries::Vector{LocalAxisRotation.LocalBasisEntry}, num_wann::Integer, layout::Symbol)
    index_map = SpinLayout.source_to_canonical_indices(num_wann, layout)
    index_map == collect(1:Int(num_wann)) && return entries
    return LocalAxisRotation.LocalBasisEntry[
        LocalAxisRotation.LocalBasisEntry(index_map[entry.index], entry.site, entry.orbital, entry.spin)
        for entry in entries
    ]
end

function build_projection_spec(config::RunConfig, num_wann::Integer)::Union{Nothing, ProjectionSpec}
    projection = config.projection
    projection.enabled || return nothing

    spec = if projection.mode == :index_groups
        groups = ProjectionGroup[
            ProjectionGroup(group.label, group.indices, group.color)
            for group in projection.groups
        ]
        ProjectionSpec(groups, Int(num_wann))
    elseif projection.mode == :win_groups
        isnothing(projection.win_path) && error("projection.win is required for mode=\"win_groups\"")
        basis = Win90Basis.read_win_basis(projection.win_path; spin_layout=config.spin.layout)
        Win90Basis.build_projection_spec(basis, projection.groups)
    else
        error("Unsupported projection.mode: $(projection.mode)")
    end

    spec.num_wann == num_wann ||
        error("projection basis has num_wann=$(spec.num_wann), but hr has num_wann=$num_wann")
    return _canonicalize_projection_spec(spec, config.spin.layout)
end

function build_basis_rotation_transform(config::RunConfig, num_wann::Integer)::Union{Nothing, Matrix{ComplexF64}}
    rotation = config.projection.basis_rotation
    rotation.enabled || return nothing
    isempty(rotation.local_axes) && error("projection.basis_rotation.local_axes is required")

    entries = if config.projection.mode == :win_groups
        isnothing(config.projection.win_path) && error("projection.win is required for basis_rotation")
        basis = Win90Basis.read_win_basis(config.projection.win_path; spin_layout=config.spin.layout)
        LocalAxisRotation.basis_entries_from_win(basis)
    else
        error("projection.basis_rotation requires mode=\"win_groups\"")
    end
    length(entries) == num_wann ||
        error("basis rotation metadata has $(length(entries)) entries, but hr has num_wann=$num_wann")
    canonical_entries = _canonicalize_basis_entries(entries, num_wann, config.spin.layout)
    return LocalAxisRotation.build_rotation_transform(
        Int(num_wann),
        canonical_entries,
        rotation.local_axes;
        strict_t2g=rotation.strict_t2g,
        leakage_tol=rotation.leakage_tol,
    )
end

function projection_weights(
    spec::ProjectionSpec,
    eigenvectors;
    basis_transform::Union{Nothing, AbstractMatrix{<:Complex}}=nothing,
)
    if isnothing(basis_transform)
        return weights_for_eigenvectors(spec, eigenvectors)
    end
    size(basis_transform, 2) == size(eigenvectors, 1) ||
        error("basis rotation has $(size(basis_transform, 2)) columns, but evecs has $(size(eigenvectors, 1)) rows")
    rotated = basis_transform * eigenvectors
    return weights_for_eigenvectors(spec, rotated)
end

function band_projection_result(spec::ProjectionSpec, weights::Array{Float64, 3})
    return BandProjectionResult(
        [group.label for group in spec.groups],
        [group.color for group in spec.groups],
        weights,
        spec.disjoint,
        spec.covers_all,
    )
end

function projected_dos_result(
    spec::ProjectionSpec,
    pdos::Matrix{Float64};
    centers_of_mass::Union{Nothing, Vector{Float64}}=nothing,
)
    size(pdos, 2) == length(spec.groups) ||
        error("projected DOS has $(size(pdos, 2)) groups, expected $(length(spec.groups))")
    if !isnothing(centers_of_mass)
        length(centers_of_mass) == length(spec.groups) ||
            error("projected DOS centers have $(length(centers_of_mass)) groups, expected $(length(spec.groups))")
    end
    atom_counts = [group.atom_count for group in spec.groups]
    normalized = copy(pdos)
    for (ig, atom_count) in enumerate(atom_counts)
        normalized[:, ig] ./= atom_count
    end
    centers = isnothing(centers_of_mass) ? fill(NaN, length(spec.groups)) : centers_of_mass
    return ProjectedDosResult(
        [group.label for group in spec.groups],
        [group.color for group in spec.groups],
        atom_counts,
        normalized,
        centers,
        spec.disjoint,
        spec.covers_all,
    )
end

end
