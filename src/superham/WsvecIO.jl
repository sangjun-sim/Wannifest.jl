module WsvecIO

using ..Model: HrModel, WsvecTable
using ..WannierWsvecIO

export read_wsvec, attach_wsvec

read_wsvec(
    path::AbstractString;
    num_wann::Union{Nothing, Int}=nothing,
    spin_layout=:qe,
)::WsvecTable = WannierWsvecIO.read_wsvec(path; num_wann=num_wann, spin_layout=spin_layout)

function attach_wsvec(model::HrModel, ws::WsvecTable)
    return HrModel(model.blocks, model.lattice, model.reciprocal, ws, model.centers)
end

end
