module Model

using ..WannierTypes

export HrData, WsvecTable
export RunFiles, PlaneConfig, EnergyConfig, PlotConfig, OutputConfig, RunConfig
export EnergySurfaceResult, ContourRunResult

const HrData = WannierTypes.HrBlocks
const WsvecTable = WannierTypes.WsvecTable

struct RunFiles
    hr_path::String
    structure_path::String
    wsvec_path::Union{Nothing, String}
end

struct PlaneConfig
    x_axis::Int
    y_axis::Int
    fixed_axis::Int
    fixed_value::Float64
    range_x::Tuple{Float64, Float64}
    range_y::Tuple{Float64, Float64}
    mesh::Tuple{Int, Int}
end

struct EnergyConfig
    shift::Float64
    bands::Vector{Int}
end

struct PlotConfig
    mode::Symbol
    interactive::Bool
    size::Tuple{Int, Int}
    energy_range::Tuple{Float64, Float64}
    colormap::String
    contour_levels::Int
end

struct OutputConfig
    output_dir::String
end

struct RunConfig
    files::RunFiles
    plane::PlaneConfig
    energy::EnergyConfig
    plot::PlotConfig
    output::OutputConfig
    spin_layout::Symbol
    hermiticity_tol::Float64
    verbose::Bool
end

struct EnergySurfaceResult
    x_axis::Vector{Float64}
    y_axis::Vector{Float64}
    kpoints_frac::Vector{Vector{Float64}}
    bands::Vector{Int}
    energies::Array{Float64, 3}
end

struct ContourRunResult
    config::RunConfig
    surface::EnergySurfaceResult
    data_path::String
    plot_handles::Vector{Any}
    hermiticity_ok::Bool
end

end
