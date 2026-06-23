module ProjectionModel

using ..LocalAxisRotation

export ProjectionGroupConfig, ProjectionBasisRotationConfig, ProjectionConfig
export BandProjectionResult, ProjectedDosResult
export disabled_basis_rotation_config, disabled_projection_config

struct ProjectionGroupConfig
    label::String
    color::String
    indices::Vector{Int}
    species::Vector{String}
    sites::Vector{Int}
    site_labels::Vector{String}
    orbitals::Vector{String}
    orbital_shells::Vector{String}
    spin::String
end

struct ProjectionBasisRotationConfig
    enabled::Bool
    local_axes::Vector{LocalAxisRotation.AxisSpec}
    strict_t2g::Bool
    leakage_tol::Float64
end

disabled_basis_rotation_config() =
    ProjectionBasisRotationConfig(false, LocalAxisRotation.AxisSpec[], false, 1.0e-8)

struct ProjectionConfig
    enabled::Bool
    mode::Symbol
    color_group::Vector{String}
    plot_style::Symbol
    colorbar_colormap::String
    colorbar_colors::Vector{String}
    circle_max_size::Float64
    circle_stroke_width::Float64
    weights_data::String
    projected_bands_plot::String
    pdos_data::String
    pdos_plot::String
    projected_combined_plot::String
    win_path::Union{Nothing, String}
    groups::Vector{ProjectionGroupConfig}
    basis_rotation::ProjectionBasisRotationConfig
end

ProjectionConfig(
    enabled::Bool,
    mode::Symbol,
    color_group::Vector{String},
    plot_style::Symbol,
    colorbar_colormap::String,
    colorbar_colors::Vector{String},
    circle_max_size::Float64,
    circle_stroke_width::Float64,
    weights_data::String,
    projected_bands_plot::String,
    pdos_data::String,
    pdos_plot::String,
    projected_combined_plot::String,
    win_path::Union{Nothing, String},
    groups::Vector{ProjectionGroupConfig},
) = ProjectionConfig(
    enabled,
    mode,
    color_group,
    plot_style,
    colorbar_colormap,
    colorbar_colors,
    circle_max_size,
    circle_stroke_width,
    weights_data,
    projected_bands_plot,
    pdos_data,
    pdos_plot,
    projected_combined_plot,
    win_path,
    groups,
    disabled_basis_rotation_config(),
)

function disabled_projection_config()
    return ProjectionConfig(
        false,
        :index_groups,
        String[],
        :colorbar,
        "viridis",
        String[],
        9.0,
        1.0,
        "",
        "",
        "",
        "",
        "",
        nothing,
        ProjectionGroupConfig[],
        disabled_basis_rotation_config(),
    )
end

struct BandProjectionResult
    labels::Vector{String}
    colors::Vector{String}
    weights::Array{Float64, 3}
    disjoint::Bool
    covers_all::Bool
end

struct ProjectedDosResult
    labels::Vector{String}
    colors::Vector{String}
    atom_counts::Vector{Int}
    pdos::Matrix{Float64}
    centers_of_mass::Vector{Float64}
    disjoint::Bool
    covers_all::Bool
end

ProjectedDosResult(
    labels::Vector{String},
    colors::Vector{String},
    atom_counts::Vector{Int},
    pdos::Matrix{Float64},
    disjoint::Bool,
    covers_all::Bool,
) = ProjectedDosResult(
    labels,
    colors,
    atom_counts,
    pdos,
    fill(NaN, length(labels)),
    disjoint,
    covers_all,
)

ProjectedDosResult(
    labels::Vector{String},
    colors::Vector{String},
    pdos::Matrix{Float64},
    disjoint::Bool,
    covers_all::Bool,
) = ProjectedDosResult(labels, colors, fill(1, length(labels)), pdos, disjoint, covers_all)

end
