module CrystalCells

using LinearAlgebra

export CrystalCell
export natoms, species_counts, species_labels
export wrap_fractional, reciprocal_lattice

struct CrystalCell
    comment::String
    lattice::Matrix{Float64}          # 3x3, columns are a, b, c
    frac_positions::Matrix{Float64}   # 3xN
    species_names::Vector{String}
    species_ids::Vector{Int}
    source::String
end

natoms(cell::CrystalCell) = size(cell.frac_positions, 2)

function species_counts(cell::CrystalCell)::Vector{Int}
    counts = zeros(Int, length(cell.species_names))
    for id in cell.species_ids
        counts[id] += 1
    end
    return counts
end

species_labels(cell::CrystalCell) = [cell.species_names[id] for id in cell.species_ids]

function wrap_fractional(x::AbstractVector{<:Real}; tol::Float64=1e-10)
    y = Vector{Float64}(undef, length(x))
    y .= mod.(Float64.(x), 1.0)
    y[isapprox.(y, 1.0; atol=tol)] .= 0.0
    y[isapprox.(y, 0.0; atol=tol)] .= 0.0
    return y
end

function wrap_fractional(xs::AbstractMatrix{<:Real}; tol::Float64=1e-10)
    ys = Matrix{Float64}(undef, size(xs)...)
    ys .= mod.(Float64.(xs), 1.0)
    ys[isapprox.(ys, 1.0; atol=tol)] .= 0.0
    ys[isapprox.(ys, 0.0; atol=tol)] .= 0.0
    return ys
end

reciprocal_lattice(lattice::AbstractMatrix{<:Real}) = 2π .* inv(Matrix{Float64}(lattice))'

end
