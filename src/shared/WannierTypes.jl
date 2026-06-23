module WannierTypes

export RKey, HrBlocks, WsvecEntry, WsvecTable, is_normalized

const RKey = NTuple{3, Int}

struct HrBlocks
    header::String
    num_wann::Int
    nrpts::Int
    hoppings::Dict{RKey, Matrix{ComplexF64}}
    ndegen::Dict{RKey, Int}
    normalization::Symbol

    function HrBlocks(
        header::AbstractString,
        num_wann::Integer,
        nrpts::Integer,
        hoppings::Dict{RKey, Matrix{ComplexF64}},
        ndegen::Dict{RKey, Int},
        normalization::Symbol,
    )
        normalization in (:raw, :normalized) ||
            error("Unsupported hr normalization: $normalization")
        return new(String(header), Int(num_wann), Int(nrpts), hoppings, ndegen, normalization)
    end
end

struct WsvecEntry
    n_shift::Int
    shifts::Matrix{Int}
end

struct WsvecTable
    table::Dict{Tuple{RKey, Int, Int}, WsvecEntry}
end

is_normalized(hr::HrBlocks) = hr.normalization == :normalized

end
