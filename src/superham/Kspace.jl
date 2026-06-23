module Kspace

using LinearAlgebra

using ..Model: HrModel
using ..WannierEigensystem
using ..WannierKspace

export hamiltonian_k, eigenvalues_k, hamiltonian_k_ws, eigenvalues_k_ws

hamiltonian_k(model::HrModel, k_frac::AbstractVector{<:Real}) =
    WannierKspace.hamiltonian_k_plain(model.blocks, k_frac)

function eigenvalues_k(
    model::HrModel,
    k_frac::AbstractVector{<:Real};
    hermiticity_tol::Float64=1e-8,
)
    evals, _ = WannierEigensystem.solve_kpoint(
        model.blocks,
        k_frac;
        hermiticity_tol=hermiticity_tol,
    )
    return evals
end

function hamiltonian_k_ws(model::HrModel, k_frac::AbstractVector{<:Real})
    isnothing(model.wsvec) && throw(ArgumentError("hamiltonian_k_ws requires model.wsvec"))
    return WannierKspace.hamiltonian_k_wsvec(model.blocks, model.wsvec, k_frac)
end

function eigenvalues_k_ws(
    model::HrModel,
    k_frac::AbstractVector{<:Real};
    hermiticity_tol::Float64=1e-8,
)
    isnothing(model.wsvec) && throw(ArgumentError("eigenvalues_k_ws requires model.wsvec"))
    evals, _ = WannierEigensystem.solve_kpoint(
        model.blocks,
        k_frac;
        wsvec=model.wsvec,
        hermiticity_tol=hermiticity_tol,
    )
    return evals
end

end
