module Validate

using LinearAlgebra

using ..PairChecks: check_hr_pair_symmetry
using ..Kspace: eigenvalues_k, eigenvalues_k_ws
using ..Model: HrModel
using ..SupercellHam: build_supercell_index, supercell_multiplicity
using ..SupercellGeometry
using ..WannierHrIO: normalized_hoppings

export validate_hermiticity, validate_supercell_size, validate_folded_spectrum

function validate_hermiticity(model::HrModel; atol::Float64=1e-10)
    return check_hr_pair_symmetry(normalized_hoppings(model.blocks); atol=atol)
end

function validate_supercell_size(model::HrModel, S::Matrix{Int}, super_model::HrModel)
    mult = SupercellGeometry.from_user_matrix(S).multiplicity
    expected = mult * model.num_wann
    actual = super_model.num_wann
    actual == expected || error("Supercell orbital count mismatch: expected $expected, got $actual")
    return (multiplicity=mult, expected_orbitals=expected, actual_orbitals=actual)
end

function _primitive_k_shifts(S::Matrix{Int})
    geom = SupercellGeometry.from_user_matrix(S)
    shifts = Vector{Vector{Float64}}()
    for j in axes(geom.k_shifts, 2)
        push!(shifts, vec(geom.k_shifts[:, j]))
    end
    return shifts
end

function validate_folded_spectrum(model::HrModel, super_model::HrModel, S::Matrix{Int};
                                  kpoints::Vector{Vector{Float64}}=[Float64[0.13, 0.21, 0.07]],
                                  atol::Float64=1e-8)
    geom = SupercellGeometry.from_user_matrix(S)
    shifts = _primitive_k_shifts(S)
    max_diff = 0.0

    for K in kpoints
        length(K) == 3 || error("Each k-point must have length 3")
        base = SupercellGeometry.primitive_k_base(geom, K)
        eval_sc = eigenvalues_k(super_model, K)
        folded = Float64[]
        for shift in shifts
            primitive_evals = isnothing(model.wsvec) ?
                eigenvalues_k(model, base .+ shift) :
                eigenvalues_k_ws(model, base .+ shift)
            append!(folded, primitive_evals)
        end
        eval_fold = sort(folded)
        length(eval_sc) == length(eval_fold) || error("Spectrum size mismatch at K=$K")
        diff = maximum(abs.(eval_sc .- eval_fold))
        max_diff = max(max_diff, diff)
    end

    max_diff <= atol || error("Folded-spectrum validation failed with max diff $max_diff > $atol")
    return max_diff
end

end
