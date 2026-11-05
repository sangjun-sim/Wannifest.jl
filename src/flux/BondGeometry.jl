module BondGeometry

using LinearAlgebra

using ..Model: FluxBasisEntry, PairRecord
using ..WannierTypes: RKey

export enumerate_pairs, pair_distance, displacement_frac, displacement_cart
export cartesian_point

function displacement_frac(bi::FluxBasisEntry, bj::FluxBasisEntry, R::RKey)::NTuple{3, Float64}
    return (
        Float64(R[1]) + bj.center_frac[1] - bi.center_frac[1],
        Float64(R[2]) + bj.center_frac[2] - bi.center_frac[2],
        Float64(R[3]) + bj.center_frac[3] - bi.center_frac[3],
    )
end

function _matvec3(A::AbstractMatrix{<:Real}, v::NTuple{3, Float64})::NTuple{3, Float64}
    x = A * [v[1], v[2], v[3]]
    return (Float64(x[1]), Float64(x[2]), Float64(x[3]))
end

displacement_cart(A::AbstractMatrix{<:Real}, bi::FluxBasisEntry, bj::FluxBasisEntry, R::RKey) =
    _matvec3(A, displacement_frac(bi, bj, R))

cartesian_point(A::AbstractMatrix{<:Real}, frac::NTuple{3, Float64}) = _matvec3(A, frac)

function pair_distance(A::AbstractMatrix{<:Real}, bi::FluxBasisEntry, bj::FluxBasisEntry, R::RKey)::Float64
    return norm(collect(displacement_cart(A, bi, bj, R)))
end

function _cluster_distances(values::Vector{Float64}, tol::Float64)::Vector{Float64}
    isempty(values) && return Float64[]
    sorted = sort(values)
    clusters = Float64[sorted[1]]
    for value in sorted[2:end]
        abs(value - clusters[end]) <= tol && continue
        push!(clusters, value)
    end
    return clusters
end

function _shell_index(distance::Float64, clusters::Vector{Float64}, tol::Float64)::Int
    distance <= tol && return 0
    idx = findfirst(cluster -> abs(distance - cluster) <= tol, clusters)
    isnothing(idx) && error("Failed to assign shell for distance $distance")
    return idx
end

function enumerate_pairs(
    basis::Vector{FluxBasisEntry},
    lattice::AbstractMatrix{<:Real};
    search_bounds::RKey,
    distance_tol::Float64,
)::Vector{PairRecord}
    distances = Float64[]
    for rx in -search_bounds[1]:search_bounds[1]
        for ry in -search_bounds[2]:search_bounds[2]
            for rz in -search_bounds[3]:search_bounds[3]
                R = (rx, ry, rz)
                for bi in basis, bj in basis
                    dist = pair_distance(lattice, bi, bj, R)
                    dist > distance_tol && push!(distances, dist)
                end
            end
        end
    end

    clusters = _cluster_distances(distances, distance_tol)
    pairs = PairRecord[]
    for rx in -search_bounds[1]:search_bounds[1]
        for ry in -search_bounds[2]:search_bounds[2]
            for rz in -search_bounds[3]:search_bounds[3]
                R = (rx, ry, rz)
                for bi in basis, bj in basis
                    dist = pair_distance(lattice, bi, bj, R)
                    push!(pairs, PairRecord(bi.index, bj.index, R, dist, _shell_index(dist, clusters, distance_tol)))
                end
            end
        end
    end
    return pairs
end

end
