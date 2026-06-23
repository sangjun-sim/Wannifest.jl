module Model

using ..ProjectionModel:
    BandProjectionResult,
    ProjectedDosResult,
    ProjectionBasisRotationConfig,
    ProjectionConfig,
    ProjectionGroupConfig,
    disabled_basis_rotation_config,
    disabled_projection_config
using ..ObservablesModel: OamConfig, OamOrbitalSelection, SamConfig
using ..ObservablesModel: disabled_oam_config, disabled_sam_config
using ..WannierTypes

export HrData, WsvecTable, RunConfig
export InputFiles, OutputFiles, DosConfig, EnergyConfig, PlotConfig
export BandPlotConfig, DosPlotConfig, CombinedPlotConfig, SpinConfig
export ProjectionGroupConfig, ProjectionBasisRotationConfig, ProjectionConfig
export BandProjectionResult, ProjectedDosResult
export OamConfig, OamOrbitalSelection, SamConfig
export BandResult, DosResult, BandRunResult

const HrData = WannierTypes.HrBlocks
const WsvecTable = WannierTypes.WsvecTable

struct InputFiles
    hr_path::String
    wsvec_path::Union{Nothing, String}
    kpoints_path::String
    structure_path::String
end

struct OutputFiles
    bands_data::String
    dos_data::String
    bands_plot::String
    dos_plot::String
    combined_plot::String
end

struct DosConfig
    mesh::NTuple{3, Int}
    shift::NTuple{3, Float64}
    sigma::Float64
    npts::Int
    emin::Union{Nothing, Float64}
    emax::Union{Nothing, Float64}
end

struct EnergyConfig
    shift::Float64
end

struct PlotConfig
    interactive::Bool
    size::Tuple{Int, Int}
    energy_range::Tuple{Float64, Float64}
    font_size::Int
    targets::Vector{Symbol}
end

PlotConfig(
    interactive::Bool,
    size::Tuple{Int, Int},
    energy_range::Tuple{Float64, Float64},
    font_size::Int,
) = PlotConfig(interactive, size, energy_range, font_size, Symbol[])

struct BandPlotConfig
    linewidth::Float64
    colors::Vector{String}
    ylabel::String
end

struct DosPlotConfig
    linewidth::Float64
    color::String
    xlabel::String
    ylabel::String
end

struct CombinedPlotConfig
    dos_width_ratio::Float64
    dos_xlabel::String
    dos_ylabel::String
end

struct SpinConfig
    enabled::Bool
    layout::Symbol
    colors::Tuple{String, String}
end

SpinConfig(enabled::Bool, colors::Tuple{String, String}) = SpinConfig(enabled, :qe, colors)

struct RunConfig
    mode::Symbol
    files::InputFiles
    output::OutputFiles
    dos::DosConfig
    energy::EnergyConfig
    plot::PlotConfig
    band_plot::BandPlotConfig
    dos_plot::DosPlotConfig
    combined_plot::CombinedPlotConfig
    spin::SpinConfig
    projection::ProjectionConfig
    oam::OamConfig
    sam::SamConfig
    hermiticity_tol::Float64
    verbose::Bool
end

RunConfig(
    mode::Symbol,
    files::InputFiles, output::OutputFiles, dos::DosConfig, energy::EnergyConfig,
    plot::PlotConfig, band_plot::BandPlotConfig, dos_plot::DosPlotConfig,
    combined_plot::CombinedPlotConfig, spin::SpinConfig,
    projection::ProjectionConfig,
    hermiticity_tol::Float64,
    verbose::Bool,
) = RunConfig(mode, files, output, dos, energy, plot, band_plot, dos_plot,
    combined_plot, spin, projection, disabled_oam_config(), disabled_sam_config(), hermiticity_tol, verbose)

RunConfig(
    mode::Symbol,
    files::InputFiles, output::OutputFiles, dos::DosConfig, energy::EnergyConfig,
    plot::PlotConfig, band_plot::BandPlotConfig, dos_plot::DosPlotConfig,
    combined_plot::CombinedPlotConfig, spin::SpinConfig,
    hermiticity_tol::Float64,
    verbose::Bool,
) = RunConfig(mode, files, output, dos, energy, plot, band_plot, dos_plot,
    combined_plot, spin, disabled_projection_config(), disabled_oam_config(), disabled_sam_config(), hermiticity_tol, verbose)

struct BandResult
    kpoints_frac::Vector{Vector{Float64}}
    distances::Vector{Float64}
    eigenvalues::Matrix{Float64}
    tick_positions::Vector{Float64}
    tick_labels::Vector{String}
    segment_ranges::Vector{UnitRange{Int}}
    is_physical_distance::Bool
    projection::Union{Nothing, BandProjectionResult}
    oam::Union{Nothing, Array{Float64, 3}}
    sam::Union{Nothing, Array{Float64, 3}}
end

BandResult(
    kpoints_frac::Vector{Vector{Float64}},
    distances::Vector{Float64}, eigenvalues::Matrix{Float64},
    tick_positions::Vector{Float64}, tick_labels::Vector{String},
    segment_ranges::Vector{UnitRange{Int}},
    is_physical_distance::Bool,
) = BandResult(kpoints_frac, distances, eigenvalues, tick_positions, tick_labels,
    segment_ranges, is_physical_distance, nothing, nothing, nothing)

BandResult(
    kpoints_frac::Vector{Vector{Float64}},
    distances::Vector{Float64}, eigenvalues::Matrix{Float64},
    tick_positions::Vector{Float64}, tick_labels::Vector{String},
    segment_ranges::Vector{UnitRange{Int}},
    is_physical_distance::Bool,
    projection::Union{Nothing, BandProjectionResult},
) = BandResult(kpoints_frac, distances, eigenvalues, tick_positions, tick_labels,
    segment_ranges, is_physical_distance, projection, nothing, nothing)

BandResult(
    kpoints_frac::Vector{Vector{Float64}},
    distances::Vector{Float64}, eigenvalues::Matrix{Float64},
    tick_positions::Vector{Float64}, tick_labels::Vector{String},
    segment_ranges::Vector{UnitRange{Int}},
    is_physical_distance::Bool,
    projection::Union{Nothing, BandProjectionResult},
    oam::Union{Nothing, Array{Float64, 3}},
) = BandResult(kpoints_frac, distances, eigenvalues, tick_positions, tick_labels,
    segment_ranges, is_physical_distance, projection, oam, nothing)

function _trapezoidal_integral(values::AbstractVector{<:Real}, coordinates::AbstractVector{<:Real})
    length(values) == length(coordinates) ||
        error("trapezoidal integration values and coordinates must have the same length")
    length(values) >= 2 || error("trapezoidal integration requires at least two points")

    integral = 0.0
    @inbounds for i in 1:(length(values) - 1)
        width = Float64(coordinates[i + 1]) - Float64(coordinates[i])
        integral += 0.5 * width * (Float64(values[i]) + Float64(values[i + 1]))
    end
    return integral
end

function _density_center_of_mass(energies::AbstractVector{<:Real}, density::AbstractVector{<:Real})
    denominator = _trapezoidal_integral(density, energies)
    denominator == 0.0 && return NaN

    numerator = 0.0
    @inbounds for i in 1:(length(energies) - 1)
        width = Float64(energies[i + 1]) - Float64(energies[i])
        left = Float64(energies[i]) * Float64(density[i])
        right = Float64(energies[i + 1]) * Float64(density[i + 1])
        numerator += 0.5 * width * (left + right)
    end
    return numerator / denominator
end

function _combined_dos_center_of_mass(
    energies::Vector{Float64},
    dos::Vector{Float64},
    dos_down::Union{Nothing, Vector{Float64}},
)
    isnothing(dos_down) && return _density_center_of_mass(energies, dos)
    length(dos_down) == length(dos) ||
        error("spin-down DOS has $(length(dos_down)) points, expected $(length(dos))")
    return _density_center_of_mass(energies, dos .+ dos_down)
end

struct DosResult
    energies::Vector{Float64}
    dos::Vector{Float64}
    dos_down::Union{Nothing, Vector{Float64}}
    integral::Float64
    center_of_mass::Float64
    num_bands::Int
    window_is_auto::Bool
    projected::Union{Nothing, ProjectedDosResult}
end

DosResult(
    energies::Vector{Float64},
    dos::Vector{Float64},
    dos_down::Union{Nothing, Vector{Float64}},
    integral::Float64,
    num_bands::Int,
    window_is_auto::Bool,
) = DosResult(
    energies,
    dos,
    dos_down,
    integral,
    _combined_dos_center_of_mass(energies, dos, dos_down),
    num_bands,
    window_is_auto,
    nothing,
)

DosResult(
    energies::Vector{Float64},
    dos::Vector{Float64},
    dos_down::Union{Nothing, Vector{Float64}},
    integral::Float64,
    num_bands::Int,
    window_is_auto::Bool,
    projected::Union{Nothing, ProjectedDosResult},
) = DosResult(
    energies,
    dos,
    dos_down,
    integral,
    _combined_dos_center_of_mass(energies, dos, dos_down),
    num_bands,
    window_is_auto,
    projected,
)

struct BandRunResult
    config::RunConfig
    band_result::Union{Nothing, BandResult}
    dos_result::Union{Nothing, DosResult}
    hermiticity_ok::Bool
end

end
