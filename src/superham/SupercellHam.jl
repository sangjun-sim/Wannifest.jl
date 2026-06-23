module SupercellHam

using LinearAlgebra

using ..CrystalCells: wrap_fractional
using ..Model: CenterTable, HrModel, RKey, SupercellBuildReport, SupercellBuildResult, SupercellIndex
using ..SupercellGeometry

export supercell_multiplicity, build_supercell_index, fold_hoppings
export propagate_centers, build_supercell, build_supercell_model

supercell_multiplicity(M::AbstractMatrix{<:Integer}) =
    SupercellGeometry.supercell_multiplicity(M)

function propagate_centers(model::HrModel, reps::Matrix{Int}, A_sc::Matrix{Float64})
    isnothing(model.centers) && return nothing

    primitive = model.centers::CenterTable
    nw = model.num_wann
    nsc = size(reps, 2)
    M_float = A_sc \ model.lattice
    centers_frac = Matrix{Float64}(undef, 3, nw * nsc)
    centers_cart = Matrix{Float64}(undef, 3, nw * nsc)
    labels = Vector{String}(undef, nw * nsc)

    for alpha in 1:nsc
        rep = Float64.(reps[:, alpha])
        for n in 1:nw
            idx = (alpha - 1) * nw + n
            x_prim = primitive.centers_frac[:, n] .+ rep
            x_sc = wrap_fractional(M_float * x_prim)
            centers_frac[:, idx] = x_sc
            centers_cart[:, idx] = A_sc * x_sc
            labels[idx] = primitive.labels[n]
        end
    end

    return CenterTable(centers_frac, centers_cart, labels, primitive.source, primitive.mode)
end

function build_supercell_index(geom::SupercellGeometry.SupercellGeometryData)::SupercellIndex
    return SupercellIndex(copy(geom.lattice_matrix), copy(geom.reps), geom.multiplicity)
end

function build_supercell_index(S::Matrix{Int}; tol::Float64=1e-10)::SupercellIndex
    return build_supercell_index(SupercellGeometry.from_user_matrix(S; tol=tol))
end

function _integer_supercell_shift(S::Matrix{Int}, d::RKey; tol::Float64=1e-10)
    x = inv(Matrix{Float64}(S)) * Float64[d[1], d[2], d[3]]
    xr = round.(Int, x)
    return isapprox(x, xr; atol=tol) ? (xr[1], xr[2], xr[3]) : nothing
end

_hr_degeneracy_scale(model::HrModel, R::RKey) =
    model.blocks.normalization == :raw ? inv(Float64(model.ndegen[R])) : 1.0

function _add_plain_block!(
    sc_hoppings::Dict{RKey, Matrix{ComplexF64}},
    model::HrModel,
    idx::SupercellIndex,
    S::Matrix{Int},
    R::RKey,
    H::Matrix{ComplexF64},
    nw_sc::Int;
    tol::Float64,
)
    nw = model.num_wann
    scale = _hr_degeneracy_scale(model, R)
    for alpha in 1:idx.multiplicity, beta in 1:idx.multiplicity
        d = (
            R[1] - idx.reps[1, beta] + idx.reps[1, alpha],
            R[2] - idx.reps[2, beta] + idx.reps[2, alpha],
            R[3] - idx.reps[3, beta] + idx.reps[3, alpha],
        )
        L = _integer_supercell_shift(S, d; tol=tol)
        isnothing(L) && continue

        Hsc = get!(sc_hoppings, L) do
            zeros(ComplexF64, nw_sc, nw_sc)
        end

        I = ((alpha - 1) * nw + 1):(alpha * nw)
        J = ((beta - 1) * nw + 1):(beta * nw)
        Hsc[I, J] .+= scale .* H
    end
    return nothing
end

function _add_wsvec_block!(
    sc_hoppings::Dict{RKey, Matrix{ComplexF64}},
    model::HrModel,
    idx::SupercellIndex,
    S::Matrix{Int},
    R::RKey,
    H::Matrix{ComplexF64},
    nw_sc::Int;
    tol::Float64,
)
    nw = model.num_wann
    wsvec = model.wsvec
    isnothing(wsvec) && error("internal error: _add_wsvec_block! requires model.wsvec")
    for alpha in 1:idx.multiplicity, beta in 1:idx.multiplicity
        for j in axes(H, 2), i in axes(H, 1)
            entry = get(wsvec.table, (R, i, j), nothing)
            isnothing(entry) && error("Missing wsvec entry for (R, i, j)=($R, $i, $j)")
            entry.n_shift > 0 || error("Invalid wsvec entry for (R, i, j)=($R, $i, $j): n_shift must be positive")
            weight = inv(Float64(entry.n_shift))
            for is in 1:entry.n_shift
                shift = @view entry.shifts[:, is]
                d = (
                    R[1] + shift[1] - idx.reps[1, beta] + idx.reps[1, alpha],
                    R[2] + shift[2] - idx.reps[2, beta] + idx.reps[2, alpha],
                    R[3] + shift[3] - idx.reps[3, beta] + idx.reps[3, alpha],
                )
                L = _integer_supercell_shift(S, d; tol=tol)
                isnothing(L) && continue

                Hsc = get!(sc_hoppings, L) do
                    zeros(ComplexF64, nw_sc, nw_sc)
                end
                row = (alpha - 1) * nw + i
                col = (beta - 1) * nw + j
                Hsc[row, col] += weight * H[i, j]
            end
        end
    end
    return nothing
end

function fold_hoppings(model::HrModel, idx::SupercellIndex; tol::Float64=1e-10)
    nw_sc = model.num_wann * idx.multiplicity
    sc_hoppings = Dict{RKey, Matrix{ComplexF64}}()

    @inbounds for (R, H) in model.hoppings
        if isnothing(model.wsvec)
            _add_plain_block!(sc_hoppings, model, idx, idx.S, R, H, nw_sc; tol=tol)
        else
            _add_wsvec_block!(sc_hoppings, model, idx, idx.S, R, H, nw_sc; tol=tol)
        end
    end
    return sc_hoppings
end

function build_supercell(
    model::HrModel,
    geom::SupercellGeometry.SupercellGeometryData;
    strict_geometry::Bool=false,
    tol::Float64=1e-10,
)::SupercellBuildResult
    strict_geometry && (isnothing(model.wsvec) || isnothing(model.centers)) &&
        error("strict_geometry=true requires both wsvec and center information")

    idx = build_supercell_index(geom)
    nw_sc = model.num_wann * idx.multiplicity
    sc_hoppings = fold_hoppings(model, idx; tol=tol)
    A_sc = SupercellGeometry.lattice_from_primitive(model.lattice, geom)
    B_sc = 2π .* inv(A_sc)'
    ndegen_sc = Dict{RKey, Int}(R => 1 for R in keys(sc_hoppings))
    centers_sc = propagate_centers(model, idx.reps, A_sc)

    super_model = HrModel(
        "supercell(" * model.header * ")",
        A_sc,
        B_sc,
        nw_sc,
        sc_hoppings,
        ndegen_sc,
        nothing,
        centers_sc,
    )
    report = SupercellBuildReport(
        !isnothing(model.wsvec),
        isnothing(model.wsvec) ? :none : :dropped,
        isnothing(centers_sc) ? :none : :propagated,
    )
    return SupercellBuildResult(super_model, report)
end

function build_supercell_model(model::HrModel, S::Matrix{Int}; strict_geometry::Bool=false, tol::Float64=1e-10)::HrModel
    geom = SupercellGeometry.from_user_matrix(S; tol=tol)
    return build_supercell(model, geom; strict_geometry=strict_geometry, tol=tol).model
end

end
