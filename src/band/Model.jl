module Model

using ..ProjectionModel:
    BandProjectionResult,
    ProjectedDosResult,
    ProjectionBasisRotationConfig,
    ProjectionConfig,
    ProjectionGroupConfig
using ..ObservablesModel: OamConfig, OamOrbitalSelection, SamConfig
using ..WannierTypes

export HrData, WsvecTable, RunConfig
export InputFiles, OutputFiles, DosConfig, EnergyConfig, PlotConfig
export BandPlotConfig, DosPlotConfig, CombinedPlotConfig, SpinConfig, ReferenceConfig
export ProjectionGroupConfig, ProjectionBasisRotationConfig, ProjectionConfig
export BandProjectionResult, ProjectedDosResult
export OamConfig, OamOrbitalSelection, SamConfig
export BandResult, DosResult, BandRunResult, BandComparison
export disabled_reference_config

const HrData = WannierTypes.HrBlocks
const WsvecTable = WannierTypes.WsvecTable

struct InputFiles
    hr_path::String
    hr_down_path::Union{Nothing, String}
    wsvec_path::Union{Nothing, String}
    wsvec_down_path::Union{Nothing, String}
    kpoints_path::String
    structure_path::String
end

struct OutputFiles
    bands_data::String
    dos_data::String
    bands_plot::String
    dos_plot::String
    combined_plot::String
    compare_data::String
    compare_plot::String
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
    margins::NTuple{4, Float64}
    targets::Vector{Symbol}
    fermi_line::Bool
    # `nothing` resolves per plot: 0.0 for shifted band/DOS plots and the
    # reference-derived level for QE comparison plots.
    fermi_level::Union{Nothing, Float64}
end

struct BandPlotConfig
    linewidth::Float64
    colors::Vector{String}
    ylabel::String
    legend::Bool
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

struct ReferenceConfig
    enabled::Bool
    source::Symbol
    bands_xml_path::Union{Nothing, String}
    scf_xml_path::Union{Nothing, String}
    k_labels::Vector{String}
    align::Symbol
    window::Union{Nothing, Tuple{Float64, Float64}}
    color::String
    color_down::String
    linewidth::Float64
end

disabled_reference_config() = ReferenceConfig(
    false,
    :qe_xml,
    nothing,
    nothing,
    String[],
    :none,
    nothing,
    "gray60",
    "gray60",
    0.8,
)

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
    reference::ReferenceConfig
    hermiticity_tol::Float64
    verbose::Bool
end

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

struct BandComparison
    qe_eigenvalues::Matrix{Float64}
    n_in_window::Vector{Int}
    rms_eV::Vector{Float64}
    max_eV::Vector{Float64}
    max_at_k::Vector{Vector{Float64}}
    overall_rms_eV::Float64
end

BandResult(
    kpoints_frac::Vector{Vector{Float64}},
    distances::Vector{Float64}, eigenvalues::Matrix{Float64},
    tick_positions::Vector{Float64}, tick_labels::Vector{String},
    segment_ranges::Vector{UnitRange{Int}},
    is_physical_distance::Bool,
) = BandResult(kpoints_frac, distances, eigenvalues, tick_positions, tick_labels,
    segment_ranges, is_physical_distance, nothing, nothing, nothing)

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

struct BandRunResult
    config::RunConfig
    band_result::Union{Nothing, BandResult}
    band_result_down::Union{Nothing, BandResult}
    dos_result::Union{Nothing, DosResult}
    hermiticity_ok::Bool
    comparison::Union{Nothing, BandComparison}
    comparison_down::Union{Nothing, BandComparison}
end

end
