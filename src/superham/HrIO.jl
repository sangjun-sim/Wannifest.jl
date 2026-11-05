module HrIO

using LinearAlgebra

using ..Model: HrModel
using ..WannierHrIO

export read_hr

function read_hr(path::AbstractString; lattice::Matrix{Float64}, spin_layout=:qe)::HrModel
    size(lattice) == (3, 3) || error("lattice must be 3x3")
    reciprocal = 2π .* inv(lattice)'
    blocks = WannierHrIO.read_hr(path; spin_layout=spin_layout)
    return HrModel(blocks, Matrix{Float64}(lattice), reciprocal, nothing, nothing)
end

end
