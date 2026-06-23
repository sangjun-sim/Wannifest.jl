module AtomicOam

using LinearAlgebra

using ..AtomicAngularMomentum: D_ORDER, EG_ORDER, P_ORDER, T2G_ORDER
using ..AtomicAngularMomentum: d_operators, p_operators
using ..BasisLabelNormalize: canonical_orbital
using ..Win90OrbitalTokens: expand_orbital_token

export OamOperators, build_l_operators, oam_expectations

struct OamOperators
    Lx::Matrix{ComplexF64}
    Ly::Matrix{ComplexF64}
    Lz::Matrix{ComplexF64}
    L2::Matrix{ComplexF64}
    l2_mode::Symbol
end
function _subblock(ops, labels::Vector{String})
    idx = [findfirst(==(label), D_ORDER) for label in labels]
    any(isnothing, idx) && error("internal OAM subblock requested non-d orbital labels")
    rows = Int.(idx)
    return (ops[1][rows, rows], ops[2][rows, rows], ops[3][rows, rows])
end

function _selected_orbitals(token::AbstractString)
    clean = canonical_orbital(token)
    if clean == "eg"
        return copy(EG_ORDER)
    elseif clean == "t2g"
        return copy(T2G_ORDER)
    end
    return [canonical_orbital(label) for label in expand_orbital_token(clean)]
end

function _entries_by_site(entries)
    by_site = Dict{String, Vector{Any}}()
    for entry in entries
        push!(get!(by_site, entry.site, Any[]), entry)
    end
    return by_site
end

function _entries_by_spin(entries)
    by_spin = Dict{String, Vector{Any}}()
    for entry in entries
        push!(get!(by_spin, entry.spin, Any[]), entry)
    end
    return by_spin
end

function _entries_by_orbital(entries, site::AbstractString, spin::AbstractString)
    by_orbital = Dict{String, Any}()
    for entry in entries
        orbital = canonical_orbital(entry.orbital)
        haskey(by_orbital, orbital) && error("site $site spin $spin has duplicate orbital $orbital")
        by_orbital[orbital] = entry
    end
    return by_orbital
end

function _operator_block(labels::Vector{String})
    label_set = Set(labels)
    if label_set == Set(["s"])
        z = zeros(ComplexF64, 1, 1)
        return labels, (z, z, z)
    elseif label_set == Set(P_ORDER)
        return P_ORDER, p_operators()
    elseif label_set == Set(D_ORDER)
        return D_ORDER, d_operators()
    elseif label_set == Set(T2G_ORDER)
        return T2G_ORDER, _subblock(d_operators(), T2G_ORDER)
    elseif label_set == Set(EG_ORDER)
        z = zeros(ComplexF64, 2, 2)
        @warn "eg-only OAM block is zero because P_eg L P_eg vanishes in the retained basis."
        return EG_ORDER, (z, z, z)
    elseif label_set == Set(["pz"])
        z = zeros(ComplexF64, 1, 1)
        @warn "pz-only OAM block is set to zero by retained-basis convention; physical full-shell L2 would be 2."
        return ["pz"], (z, z, z)
    end

    p_overlap = intersect(label_set, Set(P_ORDER))
    d_overlap = intersect(label_set, Set(D_ORDER))
    if !isempty(p_overlap)
        missing = setdiff(Set(P_ORDER), label_set)
        error("incomplete p shell for OAM selection; missing $(sort!(collect(missing)))")
    elseif !isempty(d_overlap)
        missing_t2g = setdiff(Set(T2G_ORDER), label_set)
        missing_d = setdiff(Set(D_ORDER), label_set)
        error(
            "unsupported partial d shell for OAM selection; " *
            "missing t2g=$(sort!(collect(missing_t2g))), missing full d=$(sort!(collect(missing_d)))",
        )
    end
    error("unsupported OAM orbital selection $(sort(labels))")
end

function _place_block!(
    Lx::Matrix{ComplexF64},
    Ly::Matrix{ComplexF64},
    Lz::Matrix{ComplexF64},
    labels::Vector{String},
    ops,
    by_orbital,
    selected::Set{Int},
    site::AbstractString,
    spin::AbstractString,
)
    missing = [label for label in labels if !haskey(by_orbital, label)]
    isempty(missing) || error("site $site spin $spin lacks OAM orbital(s): $(join(missing, ", "))")
    indices = [by_orbital[label].index for label in labels]
    overlap = [idx for idx in indices if idx in selected]
    isempty(overlap) || error("OAM selections overlap at Wannier index/indices $(join(overlap, ", "))")
    union!(selected, indices)
    Lx[indices, indices] .= ops[1]
    Ly[indices, indices] .= ops[2]
    Lz[indices, indices] .= ops[3]
    return nothing
end

function build_l_operators(num_wann::Integer, entries, selections)
    nw = Int(num_wann)
    all(entry -> 1 <= entry.index <= nw, entries) ||
        error("OAM basis index outside hr range 1:$nw")

    Lx = zeros(ComplexF64, nw, nw)
    Ly = zeros(ComplexF64, nw, nw)
    Lz = zeros(ComplexF64, nw, nw)
    by_site = _entries_by_site(entries)
    selected = Set{Int}()

    for selection in selections
        site_entries = get(by_site, selection.site, nothing)
        isnothing(site_entries) && error("OAM selection references unknown site label $(selection.site)")
        labels = _selected_orbitals(selection.orbital_shell)
        canonical_labels, ops = _operator_block(labels)
        for (spin, spin_entries) in sort!(collect(_entries_by_spin(site_entries)); by=first)
            by_orbital = _entries_by_orbital(spin_entries, selection.site, spin)
            _place_block!(Lx, Ly, Lz, canonical_labels, ops, by_orbital, selected, selection.site, spin)
        end
    end

    L2 = Lx * Lx + Ly * Ly + Lz * Lz
    return OamOperators(Lx, Ly, Lz, L2, :operator_components)
end

function _expectation(psi, op::AbstractMatrix{<:Complex}; imaginary_tol::Float64)
    value = dot(psi, op * psi)
    abs(imag(value)) <= imaginary_tol ||
        @warn "OAM expectation has imaginary residual $(imag(value)); returning real part."
    return real(value)
end

function oam_expectations(ops::OamOperators, evecs; imaginary_tol::Float64=1.0e-10)
    nbands = size(evecs, 2)
    values = zeros(Float64, nbands, 5)
    for ib in 1:nbands
        psi = @view evecs[:, ib]
        lx = _expectation(psi, ops.Lx; imaginary_tol=imaginary_tol)
        ly = _expectation(psi, ops.Ly; imaginary_tol=imaginary_tol)
        lz = _expectation(psi, ops.Lz; imaginary_tol=imaginary_tol)
        l2 = _expectation(psi, ops.L2; imaginary_tol=imaginary_tol)
        values[ib, :] .= (lx, ly, lz, sqrt(lx^2 + ly^2 + lz^2), l2)
    end
    return values
end

end
