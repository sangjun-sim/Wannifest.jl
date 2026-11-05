module Model

using LinearAlgebra

using ..WannierTypes

export RKey, WsvecEntry, WsvecTable, CenterTable, OrbitalSpec
export RunConfig, HrModel, SupercellIndex
export SupercellBuildReport, SupercellBuildResult, SuperhamRunResult

const RKey = WannierTypes.RKey
const WsvecEntry = WannierTypes.WsvecEntry
const WsvecTable = WannierTypes.WsvecTable

struct CenterTable
    centers_frac::Matrix{Float64}   # 3 x num_wann
    centers_cart::Matrix{Float64}   # 3 x num_wann
    labels::Vector{String}
    source::String
    mode::Symbol                    # :wannier | :manual_centers | :atomic_assumption
end

struct OrbitalSpec
    label::String
    center_frac::NTuple{3, Float64}
    atom_index::Union{Nothing, Int}
end

struct RunConfig
    hr_path::String
    structure_path::String
    win_path::Union{Nothing, String}
    wsvec_path::Union{Nothing, String}
    centres_path::Union{Nothing, String}
    output_hr::Union{Nothing, String}
    spin_layout::Symbol
    strict_geometry::Bool
    geometry_mode::Symbol           # :none | :wsvec | :manual_centers | :atomic_assumption
    supercell_matrix::Matrix{Int}
    manual_num_wann::Union{Nothing, Int}
    orbital_specs::Vector{OrbitalSpec}
end

struct HrModel
    blocks::WannierTypes.HrBlocks
    lattice::Matrix{Float64}        # 3x3, columns
    reciprocal::Matrix{Float64}     # 3x3, includes 2π
    wsvec::Union{Nothing, WsvecTable}
    centers::Union{Nothing, CenterTable}
end

function HrModel(
    header::AbstractString,
    lattice::Matrix{Float64},
    reciprocal::Matrix{Float64},
    num_wann::Integer,
    hoppings::Dict{RKey, Matrix{ComplexF64}},
    ndegen::Dict{RKey, Int},
    wsvec::Union{Nothing, WsvecTable},
    centers::Union{Nothing, CenterTable};
    normalization::Symbol=:raw,
)
    blocks = WannierTypes.HrBlocks(
        header,
        num_wann,
        length(hoppings),
        hoppings,
        ndegen,
        normalization,
    )
    return HrModel(blocks, Matrix{Float64}(lattice), Matrix{Float64}(reciprocal), wsvec, centers)
end

function Base.getproperty(model::HrModel, name::Symbol)
    if name === :header
        return getfield(getfield(model, :blocks), :header)
    elseif name === :num_wann
        return getfield(getfield(model, :blocks), :num_wann)
    elseif name === :nrpts
        return getfield(getfield(model, :blocks), :nrpts)
    elseif name === :hoppings
        return getfield(getfield(model, :blocks), :hoppings)
    elseif name === :ndegen
        return getfield(getfield(model, :blocks), :ndegen)
    end
    return getfield(model, name)
end

function Base.propertynames(::HrModel, private::Bool=false)
    public = (:header, :lattice, :reciprocal, :num_wann, :nrpts, :hoppings, :ndegen, :wsvec, :centers)
    return private ? (:blocks, public...) : public
end

struct SupercellIndex
    S::Matrix{Int}
    reps::Matrix{Int}               # 3 x Nsc
    multiplicity::Int
end

struct SupercellBuildReport
    wsvec_input::Bool
    wsvec_output_policy::Symbol     # :none | :dropped | :rebuilt
    center_output_policy::Symbol    # :none | :propagated
end

struct SupercellBuildResult
    model::HrModel
    report::SupercellBuildReport
end

struct SuperhamRunResult
    config::RunConfig
    model::HrModel
    super_model::HrModel
    output_hr::Union{Nothing, String}
    geometry_source::String
    build_report::SupercellBuildReport
    primitive_hermiticity_error::Float64
    multiplicity::Int
    folded_spectrum_diff::Float64
    primitive_eigenvalues::Vector{Float64}
    primitive_wsvec_eigenvalues::Union{Nothing, Vector{Float64}}
    supercell_eigenvalues::Vector{Float64}
end

end
