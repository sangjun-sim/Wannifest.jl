module AtomicSpin

using LinearAlgebra

using ..BasisLabelNormalize: canonical_orbital, canonical_spin

export SpinOperators, build_spin_operators, spin_expectations

struct SpinOperators
    Sx::Matrix{ComplexF64}
    Sy::Matrix{ComplexF64}
    Sz::Matrix{ComplexF64}
    S2::Matrix{ComplexF64}
end

const SX2 = ComplexF64[0 0.5; 0.5 0]
const SY2 = ComplexF64[0 -0.5im; 0.5im 0]
const SZ2 = ComplexF64[0.5 0; 0 -0.5]

function _expectation(psi, op::AbstractMatrix{<:Complex}; imaginary_tol::Float64)
    value = dot(psi, op * psi)
    abs(imag(value)) <= imaginary_tol ||
        @warn "Spin expectation has imaginary residual $(imag(value)); returning real part."
    return real(value)
end

function _basis_groups(num_wann::Int, entries)
    indices = [entry.index for entry in entries]
    sort(indices) == collect(1:num_wann) ||
        error("SAM basis metadata must cover all Wannier indices 1:$num_wann")

    groups = Dict{Tuple{String, String}, Dict{String, Int}}()
    for entry in entries
        1 <= entry.index <= num_wann || error("SAM basis index outside hr range 1:$num_wann")
        spin = canonical_spin(entry.spin)
        spin in ("unpolarized", "any") &&
            error("band.sam requires spinful up/dn basis metadata; found spin=$spin")
        spin in ("up", "dn") ||
            error("band.sam requires spin to be up or dn, got '$spin'")
        key = (String(entry.site), canonical_orbital(entry.orbital))
        by_spin = get!(groups, key, Dict{String, Int}())
        haskey(by_spin, spin) &&
            error("duplicate SAM basis entry for site=$(key[1]) orbital=$(key[2]) spin=$spin")
        by_spin[spin] = entry.index
    end
    isempty(groups) && error("band.sam requested but basis has no spinful up/dn pairs")
    return groups
end

function build_spin_operators(num_wann::Integer, entries)
    nw = Int(num_wann)
    nw > 0 || error("SAM basis requires positive num_wann")

    Sx = zeros(ComplexF64, nw, nw)
    Sy = zeros(ComplexF64, nw, nw)
    Sz = zeros(ComplexF64, nw, nw)

    for ((site, orbital), by_spin) in _basis_groups(nw, entries)
        haskey(by_spin, "up") && haskey(by_spin, "dn") || error(
            "band.sam requires complete up/dn pair for site=$site orbital=$orbital",
        )
        length(by_spin) == 2 || error("internal SAM basis group has unexpected spin entries")
        indices = [by_spin["up"], by_spin["dn"]]
        Sx[indices, indices] .= SX2
        Sy[indices, indices] .= SY2
        Sz[indices, indices] .= SZ2
    end

    S2 = Sx * Sx + Sy * Sy + Sz * Sz
    return SpinOperators(Sx, Sy, Sz, S2)
end

function spin_expectations(ops::SpinOperators, evecs; imaginary_tol::Float64=1.0e-10)
    nbands = size(evecs, 2)
    values = zeros(Float64, nbands, 5)
    for ib in 1:nbands
        psi = @view evecs[:, ib]
        sx = _expectation(psi, ops.Sx; imaginary_tol=imaginary_tol)
        sy = _expectation(psi, ops.Sy; imaginary_tol=imaginary_tol)
        sz = _expectation(psi, ops.Sz; imaginary_tol=imaginary_tol)
        s2 = _expectation(psi, ops.S2; imaginary_tol=imaginary_tol)
        values[ib, :] .= (sx, sy, sz, sqrt(sx^2 + sy^2 + sz^2), s2)
    end
    return values
end

end
