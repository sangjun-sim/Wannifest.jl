module Transform

using LinearAlgebra
using ..PoscarIO: StructureCell, natoms, wrap_fractional
using ..SupercellGeometry

export parse_direction_matrix
export supercell_multiplicity
export build_supercell

function parse_direction_matrix(a::NTuple{3, Int}, b::NTuple{3, Int}, c::NTuple{3, Int})
    M = Matrix{Int}(hcat(collect(a), collect(b), collect(c)))
    rank(Matrix{Float64}(M)) == 3 || error("Direction matrix must be full-rank")
    return M
end

supercell_multiplicity(M::AbstractMatrix{<:Integer}) =
    SupercellGeometry.supercell_multiplicity(M)

function _fractional_distance(A::Matrix{Float64}, x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    d = Float64.(x) .- Float64.(y)
    d .-= round.(d)
    return norm(A * d)
end

function _deduplicate_positions(A::Matrix{Float64}, xs::Vector{Vector{Float64}}, ids::Vector{Int}; tol::Float64=1e-8)
    kept_xs = Vector{Vector{Float64}}()
    kept_ids = Int[]

    for (x, id) in zip(xs, ids)
        duplicate = false
        for (y, yid) in zip(kept_xs, kept_ids)
            if id == yid && _fractional_distance(A, x, y) <= tol
                duplicate = true
                break
            end
        end
        if !duplicate
            push!(kept_xs, x)
            push!(kept_ids, id)
        end
    end

    return kept_xs, kept_ids
end

function build_supercell(base::StructureCell, M::Matrix{Int}; tol::Float64=1e-8)
    geom = SupercellGeometry.from_user_matrix(M; tol=tol)
    A_new = SupercellGeometry.lattice_from_primitive(base.lattice, geom)
    M_inv = inv(Matrix{Float64}(geom.lattice_matrix))

    positions = Vector{Vector{Float64}}()
    ids = Int[]
    for atom_index in 1:natoms(base)
        x_old = base.frac_positions[:, atom_index]
        for rep in eachcol(geom.reps)
            x_new = wrap_fractional(M_inv * (x_old + rep))
            push!(positions, x_new)
            push!(ids, base.species_ids[atom_index])
        end
    end

    unique_positions, unique_ids = _deduplicate_positions(A_new, positions, ids; tol=tol)
    expected_atoms = natoms(base) * geom.multiplicity
    length(unique_positions) == expected_atoms || error("Expected $expected_atoms atoms, found $(length(unique_positions))")

    frac_positions = hcat(unique_positions...)
    return StructureCell(
        "$(base.comment) [supercell]",
        A_new,
        wrap_fractional(frac_positions),
        copy(base.species_names),
        unique_ids,
        base.source,
    )
end

end
