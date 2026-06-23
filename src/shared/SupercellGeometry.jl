module SupercellGeometry

using LinearAlgebra

export SupercellGeometryData
export from_user_matrix, lattice_from_primitive, primitive_k_base
export supercell_multiplicity

struct SupercellGeometryData
    input_matrix::Matrix{Int}
    lattice_matrix::Matrix{Int}
    reps::Matrix{Int}
    k_shifts::Matrix{Float64}
    multiplicity::Int
end

supercell_multiplicity(M::AbstractMatrix{<:Integer}) =
    round(Int, abs(det(Matrix{Float64}(M))))

function _half_open_contains(u::AbstractVector{<:Real}; tol::Float64=1e-10)
    return all((-tol .<= u) .& (u .< 1.0 - tol))
end

function _representative_bounds(M::Matrix{Int})
    corners = Vector{Vector{Float64}}()
    for a in (0.0, 1.0), b in (0.0, 1.0), c in (0.0, 1.0)
        push!(corners, Matrix{Float64}(M) * Float64[a, b, c])
    end
    stacked = hcat(corners...)
    mins = floor.(Int, minimum(stacked; dims=2)[:]) .- 1
    maxs = ceil.(Int, maximum(stacked; dims=2)[:]) .+ 1
    return mins, maxs
end

function _coset_representatives(M::Matrix{Int}; tol::Float64=1e-10)
    M_float = Matrix{Float64}(M)
    M_inv = inv(M_float)
    mins, maxs = _representative_bounds(M)
    reps_list = NTuple{3, Int}[]

    for n1 in mins[1]:maxs[1], n2 in mins[2]:maxs[2], n3 in mins[3]:maxs[3]
        n = Float64[n1, n2, n3]
        u = M_inv * n
        if _half_open_contains(u; tol=tol)
            push!(reps_list, (n1, n2, n3))
        end
    end

    expected = supercell_multiplicity(M)
    length(reps_list) == expected ||
        error("Failed to construct $expected coset representatives; found $(length(reps_list))")

    reps = Matrix{Int}(undef, 3, length(reps_list))
    for (j, rep) in enumerate(reps_list)
        reps[:, j] = collect(rep)
    end
    return reps
end

function _reciprocal_shifts(lattice_matrix::Matrix{Int}; tol::Float64=1e-10)
    dual_reps = _coset_representatives(Matrix{Int}(transpose(lattice_matrix)); tol=tol)
    primitive_map = transpose(inv(Matrix{Float64}(lattice_matrix)))
    shifts = Matrix{Float64}(undef, 3, size(dual_reps, 2))
    for j in axes(dual_reps, 2)
        shifts[:, j] = primitive_map * Float64.(dual_reps[:, j])
    end
    return shifts
end

function from_user_matrix(S_rows::AbstractMatrix{<:Integer}; tol::Float64=1e-10)
    size(S_rows) == (3, 3) || error("supercell matrix must be 3x3")
    input_matrix = Matrix{Int}(S_rows)
    rank(Matrix{Float64}(input_matrix)) == 3 || error("supercell matrix must be full-rank")

    lattice_matrix = Matrix{Int}(transpose(input_matrix))
    reps = _coset_representatives(lattice_matrix; tol=tol)
    k_shifts = _reciprocal_shifts(lattice_matrix; tol=tol)
    return SupercellGeometryData(
        input_matrix,
        lattice_matrix,
        reps,
        k_shifts,
        supercell_multiplicity(lattice_matrix),
    )
end

lattice_from_primitive(A::Matrix{Float64}, geom::SupercellGeometryData) =
    A * Matrix{Float64}(geom.lattice_matrix)

function primitive_k_base(geom::SupercellGeometryData, K::AbstractVector{<:Real})
    length(K) == 3 || error("supercell k-point must have length 3")
    return vec(transpose(inv(Matrix{Float64}(geom.lattice_matrix))) * Float64.(K))
end

end
