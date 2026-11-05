module WannierEigensystem

using LinearAlgebra

using ..WannierKspace
using ..WannierTypes: HrBlocks, WsvecTable

export solve_kpoint, solve_kpoint_values

function _hermiticity_residual(Hk::AbstractMatrix{<:Complex})
    denom = opnorm(Hk)
    residual = norm(Hk - Hk')
    return denom == 0.0 ? residual : residual / denom
end

function _real_sorted_eigenpairs(Hk::Matrix{ComplexF64}; hermiticity_tol::Float64=1e-8)
    residual = _hermiticity_residual(Hk)
    residual <= hermiticity_tol || error(
        "H(k) is not Hermitian within tolerance " *
        "(norm(H-H')/opnorm(H) = $residual, tolerance = $hermiticity_tol).",
    )
    F = eigen(Hermitian(Hk))
    return F.values, F.vectors
end

function _real_sorted_eigenvalues(Hk::Matrix{ComplexF64}; hermiticity_tol::Float64=1e-8)
    residual = _hermiticity_residual(Hk)
    residual <= hermiticity_tol || error(
        "H(k) is not Hermitian within tolerance " *
        "(norm(H-H')/opnorm(H) = $residual, tolerance = $hermiticity_tol).",
    )
    return eigvals(Hermitian(Hk))
end

function solve_kpoint(
    hr::HrBlocks,
    k_frac::AbstractVector{<:Real};
    wsvec::Union{Nothing, WsvecTable}=nothing,
    hermiticity_tol::Float64=1e-8,
)
    Hk = isnothing(wsvec) ?
        WannierKspace.hamiltonian_k_plain(hr, k_frac) :
        WannierKspace.hamiltonian_k_wsvec(hr, wsvec, k_frac)
    return _real_sorted_eigenpairs(Hk; hermiticity_tol=hermiticity_tol)
end

function solve_kpoint_values(
    hr::HrBlocks,
    k_frac::AbstractVector{<:Real};
    wsvec::Union{Nothing, WsvecTable}=nothing,
    hermiticity_tol::Float64=1e-8,
)
    Hk = isnothing(wsvec) ?
        WannierKspace.hamiltonian_k_plain(hr, k_frac) :
        WannierKspace.hamiltonian_k_wsvec(hr, wsvec, k_frac)
    return _real_sorted_eigenvalues(Hk; hermiticity_tol=hermiticity_tol)
end

end
