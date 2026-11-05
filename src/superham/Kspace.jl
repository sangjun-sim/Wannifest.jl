module Kspace

using ..Model: HrModel
using ..WannierEigensystem

export eigenvalues_k, eigenvalues_k_ws

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
