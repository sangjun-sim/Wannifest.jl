module FluxCompiler

using ..BasisSource: entry_label
using ..BondGeometry
using ..HrHermiticity: complete_hermiticity!, reverse_key
using ..Model: FluxBasisEntry, FluxConfig, FluxEdge, FluxEndpoint, FluxRow, PairRecord
using ..WannierTypes: RKey

export apply_flux_terms!, complete_hermiticity!

function _get_hopping!(hops::Dict{RKey, Matrix{ComplexF64}}, R::RKey, nw::Int)
    return get!(hops, R) do
        zeros(ComplexF64, nw, nw)
    end
end

function _accumulate_flux_pair!(
    hops::Dict{RKey, Matrix{ComplexF64}},
    pair::PairRecord,
    value::ComplexF64,
    num_wann::Int,
)
    H = _get_hopping!(hops, pair.R, num_wann)
    H[pair.i, pair.j] += value

    Rm = reverse_key(pair.R)
    Hm = _get_hopping!(hops, Rm, num_wann)
    Hm[pair.j, pair.i] += conj(value)
    return nothing
end

function _entries_by_index(basis::Vector{FluxBasisEntry})
    return Dict(entry.index => entry for entry in basis)
end

function _selector_indices(basis::Vector{FluxBasisEntry}, label::AbstractString)::Vector{Int}
    site_indices = [entry.index for entry in basis if entry.site_label == label]
    !isempty(site_indices) && return site_indices
    species_indices = [entry.index for entry in basis if entry.species == label]
    !isempty(species_indices) && return species_indices
    error("No Wannier orbitals found for site label or species '$label'")
end

function _endpoint_indices(endpoint::FluxEndpoint, basis::Vector{FluxBasisEntry}, num_wann::Int)::Vector{Int}
    if endpoint.value isa Int
        idx = endpoint.value::Int
        1 <= idx <= num_wann || error("Wannier index $idx is outside 1:$num_wann")
        return [idx]
    end
    return _selector_indices(basis, endpoint.value::String)
end

function _row_index_pairs(row::FluxRow, basis::Vector{FluxBasisEntry}, num_wann::Int)
    from_ids = _endpoint_indices(row.from, basis, num_wann)
    to_ids = _endpoint_indices(row.to, basis, num_wann)
    return Set((i, j) for i in from_ids for j in to_ids)
end

_pair_key(i::Int, j::Int, R::RKey) = (i, j, R[1], R[2], R[3])

function _physical_pair_key(pair::PairRecord)
    direct = _pair_key(pair.i, pair.j, pair.R)
    Rm = reverse_key(pair.R)
    reverse = _pair_key(pair.j, pair.i, Rm)
    return direct <= reverse ? direct : reverse
end

function _start_finish_frac(from_entry::FluxBasisEntry, to_entry::FluxBasisEntry, R::RKey)
    start = from_entry.center_frac
    finish = (
        Float64(R[1]) + to_entry.center_frac[1],
        Float64(R[2]) + to_entry.center_frac[2],
        Float64(R[3]) + to_entry.center_frac[3],
    )
    return start, finish
end

function _edge_from_pair(
    row::FluxRow,
    pair::PairRecord,
    basis_by_index,
    lattice::AbstractMatrix{<:Real},
)::FluxEdge
    from_entry = basis_by_index[pair.i]
    to_entry = basis_by_index[pair.j]
    start_frac, finish_frac = _start_finish_frac(from_entry, to_entry, pair.R)
    start_cart = BondGeometry.cartesian_point(lattice, start_frac)
    finish_cart = BondGeometry.cartesian_point(lattice, finish_frac)
    return FluxEdge(
        row.nn,
        pair.i,
        pair.j,
        pair.R,
        row.value,
        entry_label(from_entry),
        entry_label(to_entry),
        start_frac,
        finish_frac,
        start_cart,
        finish_cart,
    )
end

function apply_flux_terms!(
    hops::Dict{RKey, Matrix{ComplexF64}},
    basis::Vector{FluxBasisEntry},
    pairs::Vector{PairRecord},
    lattice::AbstractMatrix{<:Real},
    config::FluxConfig,
    num_wann::Int,
)::Vector{FluxEdge}
    basis_by_index = _entries_by_index(basis)
    edges = FluxEdge[]
    seen_physical_pairs = Set{Tuple{Int, Int, Int, Int, Int}}()

    for term in config.terms
        for row in term.rows
            wanted = _row_index_pairs(row, basis, num_wann)
            seen_row_pairs = Set{Tuple{Int, Int, Int, Int, Int}}()
            matches = 0
            for pair in pairs
                pair.shell == row.nn || continue
                pair.R == row.R || continue
                (pair.i, pair.j) in wanted || continue
                physical_key = _physical_pair_key(pair)
                physical_key in seen_row_pairs && continue
                physical_key in seen_physical_pairs && error(
                    "Duplicate flux term for Hermitian-equivalent pair " *
                    "i=$(pair.i), j=$(pair.j), R=$(pair.R)",
                )
                push!(seen_row_pairs, physical_key)
                push!(seen_physical_pairs, physical_key)
                _accumulate_flux_pair!(hops, pair, row.value, num_wann)
                push!(edges, _edge_from_pair(row, pair, basis_by_index, lattice))
                matches += 1
            end
            matches > 0 || error(
                "No flux pair matched nn=$(row.nn), from=$(row.from.value), " *
                "to=$(row.to.value), R=$(row.R)",
            )
        end
    end
    return edges
end

end
