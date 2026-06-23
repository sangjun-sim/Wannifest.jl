module WsvecIO

using ..Model: CenterTable, HrModel, WsvecTable
using ..WannierWsvecGenerate
using ..WannierWsvecIO

export read_wsvec, attach_wsvec, generate_wsvec

read_wsvec(
    path::AbstractString;
    num_wann::Union{Nothing, Int}=nothing,
    spin_layout=:qe,
)::WsvecTable = WannierWsvecIO.read_wsvec(path; num_wann=num_wann, spin_layout=spin_layout)

function attach_wsvec(model::HrModel, ws::WsvecTable; require_unit_ndegen::Bool=false)
    WannierWsvecGenerate.assert_wsvec_usable(model.blocks, ws; require_unit_ndegen=require_unit_ndegen)
    return HrModel(model.blocks, model.lattice, model.reciprocal, ws, model.centers)
end

function generate_wsvec(
    model::HrModel,
    mp_grid::NTuple{3, Int};
    centers::Union{Nothing, CenterTable}=model.centers,
    kwargs...,
)::WsvecTable
    isnothing(centers) && error("generate_wsvec requires center information")
    return WannierWsvecGenerate.generate_wsvec(
        model.blocks,
        model.lattice,
        mp_grid,
        centers.centers_frac;
        centers_cart=centers.centers_cart,
        kwargs...,
    )
end

end
