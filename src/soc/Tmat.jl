module Tmat

using LinearAlgebra

using ..Model: WannierChannel, KpointIndex, OnsiteDiagnostic
using ..Validate: validate_tmat_inputs, validate_points
using ..WannierTypes: RKey

export tmat_k, tmat_r, tmat_unitarity_error, expected_permutation
export tmat_onsite_diagnostic, warn_onsite_diagnostic

const WEIGHT0_WARN = 0.9       # below this, R = 0 no longer dominates T(R)
const PERM0_WARN = 0.5         # |T(0)| this far from P while still unitary: wrong spin layout
const UNITARITY0_SMALL = 0.1

"""
    tmat_k(up, down, soc, V, index) -> Array{ComplexF64, 3}

`T(k) = blockdiag(Ũ_up, Ũ_down)† V Ũ_soc` as `(nw_up+nw_down) × nw_soc × nk`, in the
k order of the up channel; rows are the up block, then the down block. The zero
blocks of blockdiag are never formed.
"""
function tmat_k(up::WannierChannel, down::WannierChannel, soc::WannierChannel, V::Array{ComplexF64, 3}, index::KpointIndex)
    validate_tmat_inputs(up, down, soc, V, index)
    nb_up, nw_up = up.num_bands, up.num_wann
    nk = size(V, 3)
    T = Array{ComplexF64}(undef, nw_up + down.num_wann, soc.num_wann, nk)
    for k in 1:nk
        U_soc = view(soc.U, :, :, index.soc[k])
        T[1:nw_up, :, k] = view(up.U, :, :, k)' * view(V, 1:nb_up, :, k) * U_soc
        T[(nw_up + 1):end, :, k] = view(down.U, :, :, index.down[k])' * view(V, (nb_up + 1):size(V, 1), :, k) * U_soc
    end
    return T
end

"""
    tmat_r(Tk, kfrac, Rs) -> Array{ComplexF64, 3}

`T(R) = (1/N) Σ_k e^{+2πi k·R} T(k) = ⟨w^SP_R|w^SOC_0⟩`, one block per entry of `Rs`
(integer lattice coordinates; `kfrac` is 3×N fractional). The sign is opposite to
Wannier90 `H(R) = (1/N) Σ_k e^{−2πi k·R} H(k) = ⟨w_0|H|w_R⟩`, so the hr-convention
matrix `⟨w^SP_0|w^SOC_R⟩` is `T(−R)`.
"""
function tmat_r(Tk::Array{ComplexF64, 3}, kfrac::AbstractMatrix{<:Real}, Rs::Vector{RKey})
    validate_points(kfrac, "Wannier k-points")
    nk = size(Tk, 3)
    size(kfrac, 2) == nk || error("T(k) holds $nk k-points but the k list holds $(size(kfrac, 2))")
    TR = zeros(ComplexF64, size(Tk, 1), size(Tk, 2), length(Rs))
    for (r, R) in enumerate(Rs), k in 1:nk
        phase = cis(2pi * (kfrac[1, k] * R[1] + kfrac[2, k] * R[2] + kfrac[3, k] * R[3])) / nk
        @views TR[:, :, r] .+= phase .* Tk[:, :, k]
    end
    return TR
end

_unitarity(T::AbstractMatrix) = opnorm(T * T' - I)

"""`max_k ‖T(k)T(k)† − I‖₂`, or `nothing` when T is not square (`nw_sp ≠ nw_soc`)."""
function tmat_unitarity_error(Tk::Array{ComplexF64, 3})
    size(Tk, 1) == size(Tk, 2) || return nothing
    return maximum(_unitarity(view(Tk, :, :, k)) for k in axes(Tk, 3))
end

"""
Expected structure of T(0) when the three `.win` files share their projections:
row i (up) ↔ SOC column 2i−1 and row nw_up+i (down) ↔ SOC column 2i, because a
Wannier90 v3 spinor run interleaves ↑↓ per orbital. `nothing` unless
`nw_up == nw_down` and `nw_soc == 2 nw_up`.
"""
function expected_permutation(nw_up::Integer, nw_down::Integer, nw_soc::Integer)
    nw_up == nw_down && nw_soc == 2nw_up || return nothing
    P = zeros(Float64, nw_up + nw_down, nw_soc)
    for i in 1:nw_up
        P[i, 2i - 1] = 1.0
        P[nw_up + i, 2i] = 1.0
    end
    return P
end

"""
    tmat_onsite_diagnostic(Tk, T0, nw_up, nw_down) -> OnsiteDiagnostic

How good the R = R′ approximation is. By Parseval `Σ_R ‖T(R)‖² = (1/N) Σ_k ‖T(k)‖²`,
so `weight0 = ‖T(0)‖²_F / Σ_R ‖T(R)‖²_F` needs no other R. `unitarity0 = ‖T(0)T(0)† − I‖₂`
(square T only) and `perm0 = ‖|T(0)| − P‖_max` (when P is defined).
"""
function tmat_onsite_diagnostic(Tk::Array{ComplexF64, 3}, T0::AbstractMatrix, nw_up::Integer, nw_down::Integer)
    size(T0) == size(Tk)[1:2] || error("T(0) must have the shape of T(k)")
    total = sum(abs2, Tk) / size(Tk, 3)
    weight0 = total > 0 ? sum(abs2, T0) / total : 0.0
    unitarity0 = size(T0, 1) == size(T0, 2) ? _unitarity(T0) : nothing
    P = expected_permutation(nw_up, nw_down, size(T0, 2))
    perm0 = isnothing(P) ? nothing : maximum(abs.(abs.(T0) .- P))
    return OnsiteDiagnostic(weight0, unitarity0, perm0)
end

"""Warn about a weak R = 0 term or a suspected spin-layout mismatch; return the number of warnings."""
function warn_onsite_diagnostic(diagnostic::OnsiteDiagnostic)::Int
    count = 0
    if diagnostic.weight0 < WEIGHT0_WARN
        @warn "weight0 = $(diagnostic.weight0) < $WEIGHT0_WARN: T(R = 0) does not dominate. SP and SOC Wannier centres may sit in different unit cells (Wannier90 does not guarantee the home cell); request more R with tmat_R."
        count += 1
    end
    perm0, unitarity0 = diagnostic.perm0, diagnostic.unitarity0
    if !isnothing(perm0) && !isnothing(unitarity0) && perm0 > PERM0_WARN && unitarity0 < UNITARITY0_SMALL
        @warn "perm0 = $perm0 although T(0) is nearly unitary (unitarity0 = $unitarity0): the SOC columns may follow the old block layout (all ↑, then all ↓) instead of the Wannier90 v3 interleaved order, or the projections differ."
        count += 1
    end
    return count
end

end
