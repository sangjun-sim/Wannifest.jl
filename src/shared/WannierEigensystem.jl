module WannierEigensystem

using LinearAlgebra

using ..WannierKspace
using ..WannierTypes: HrBlocks, WsvecTable

export EigenBundle, solve_kpoint, solve_kpoint_values, solve_kpoints

struct EigenBundle
    kpoints_frac::Vector{Vector{Float64}}
    eigenvalues::Matrix{Float64}
    eigenvectors::Vector{Matrix{ComplexF64}}
end

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
    values = real.(F.values)
    perm = sortperm(values)
    return Float64.(values[perm]), Matrix{ComplexF64}(F.vectors[:, perm])
end

function _real_sorted_eigenvalues(Hk::Matrix{ComplexF64}; hermiticity_tol::Float64=1e-8)
    residual = _hermiticity_residual(Hk)
    residual <= hermiticity_tol || error(
        "H(k) is not Hermitian within tolerance " *
        "(norm(H-H')/opnorm(H) = $residual, tolerance = $hermiticity_tol).",
    )
    return Float64.(sort(real.(eigvals(Hermitian(Hk)))))
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

function solve_kpoints(
    hr::HrBlocks,
    kpoints::Vector{<:AbstractVector{<:Real}};
    wsvec::Union{Nothing, WsvecTable}=nothing,
    hermiticity_tol::Float64=1e-8,
)
    nk = length(kpoints)
    eigenvalues = Matrix{Float64}(undef, nk, hr.num_wann)
    eigenvectors = Vector{Matrix{ComplexF64}}(undef, nk)

    for ik in 1:nk
        evals, evecs = solve_kpoint(
            hr,
            kpoints[ik];
            wsvec=wsvec,
            hermiticity_tol=hermiticity_tol,
        )
        eigenvalues[ik, :] .= evals
        eigenvectors[ik] = evecs
    end

    return EigenBundle([Float64.(k) for k in kpoints], eigenvalues, eigenvectors)
end

end
