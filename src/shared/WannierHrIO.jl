module WannierHrIO

using ..PairChecks: pair_dict_max_error
using ..SpinLayout
using ..WannierTypes: HrBlocks, RKey

export read_hr, normalized_hoppings, pair_hermiticity_error

function _get_hopping!(dict::Dict{RKey, Matrix{ComplexF64}}, R::RKey, nw::Int)
    return get!(dict, R) do
        zeros(ComplexF64, nw, nw)
    end
end

function _read_ndegen(io::IO, nrpts::Int)
    values = Int[]
    while length(values) < nrpts
        eof(io) && error("Unexpected EOF while reading ndegen")
        line = strip(readline(io))
        isempty(line) && continue
        append!(values, parse.(Int, split(line)))
    end
    resize!(values, nrpts)
    return values
end

function _parse_hr_line(line::AbstractString)
    fields = split(strip(line))
    length(fields) == 7 || error("Malformed hr line: $line")
    R = (parse(Int, fields[1]), parse(Int, fields[2]), parse(Int, fields[3]))
    m = parse(Int, fields[4])
    n = parse(Int, fields[5])
    value = ComplexF64(parse(Float64, fields[6]), parse(Float64, fields[7]))
    return R, m, n, value
end

function read_hr(path::AbstractString; spin_layout=SpinLayout.DEFAULT_LAYOUT)::HrBlocks
    open(path, "r") do io
        eof(io) && error("Empty hr file: $(abspath(path))")
        header = strip(readline(io))
        num_wann = parse(Int, strip(readline(io)))
        nrpts = parse(Int, strip(readline(io)))
        ndegen_vec = _read_ndegen(io, nrpts)
        index_map = SpinLayout.source_to_canonical_indices(num_wann, spin_layout)

        hoppings = Dict{RKey, Matrix{ComplexF64}}()
        ndegen = Dict{RKey, Int}()

        for ir in 1:nrpts
            block_R = nothing
            block_deg = ndegen_vec[ir]
            for _ in 1:(num_wann * num_wann)
                eof(io) && error("Unexpected EOF while reading hopping block $ir from $(abspath(path))")
                R, m, n, value = _parse_hr_line(readline(io))
                if isnothing(block_R)
                    block_R = R
                    haskey(hoppings, R) && error("Duplicate hopping block for R=$R in $(abspath(path))")
                    ndegen[R] = block_deg
                elseif R != block_R
                    error("Encountered mixed R block in $(abspath(path)); expected $block_R, got $R")
                end
                H = _get_hopping!(hoppings, R, num_wann)
                1 <= m <= num_wann || error("hr row index $m is outside 1:$num_wann in $(abspath(path))")
                1 <= n <= num_wann || error("hr column index $n is outside 1:$num_wann in $(abspath(path))")
                H[index_map[m], index_map[n]] = value
            end
        end

        return HrBlocks(header, num_wann, nrpts, hoppings, ndegen, :raw)
    end
end

function normalized_hoppings(hr::HrBlocks)
    if hr.normalization == :normalized
        return hr.hoppings
    end
    return Dict{RKey, Matrix{ComplexF64}}(
        R => hr.hoppings[R] ./ hr.ndegen[R] for R in keys(hr.hoppings)
    )
end

pair_hermiticity_error(hr::HrBlocks) = pair_dict_max_error(normalized_hoppings(hr))

end
