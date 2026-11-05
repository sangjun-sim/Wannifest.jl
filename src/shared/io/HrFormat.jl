module HrFormat

using Printf

export write_hr_blocks_raw, write_hr_blocks_normalized

const RKey = NTuple{3, Int}

function _write_ndegen(io::IO, ndegen::Vector{Int})
    for start in 1:15:length(ndegen)
        stop = min(start + 14, length(ndegen))
        for d in ndegen[start:stop]
            @printf(io, "%5d", d)
        end
        write(io, "\n")
    end
end

function _write_block(io::IO, R::RKey, H::AbstractMatrix{ComplexF64})
    nw = size(H, 1)
    for n in 1:nw
        for m in 1:nw
            v = H[m, n]
            @printf(io, "%5d%5d%5d%5d%5d % .12f % .12f\n",
                    R[1], R[2], R[3], m, n, real(v), imag(v))
        end
    end
end

"""
    write_hr_blocks_raw(path, header, num_wann, hoppings, ndegen)

Write raw Wannier90-compatible blocks. `hoppings[R]` is interpreted as raw block data, and
`ndegen[R]` is emitted verbatim. The file format matches canonical Wannier90-style output:

- `_write_ndegen`: 15 integers per line with `%5d`
- `_write_block`: `%5d%5d%5d%5d%5d % .12f % .12f\\n`
- outer loop `n`, inner loop `m`, emitted index order `(m, n)`
"""
function write_hr_blocks_raw(
    path::AbstractString,
    header::AbstractString,
    num_wann::Int,
    hoppings::Dict{RKey, Matrix{ComplexF64}},
    ndegen::Dict{RKey, Int},
)
    Rlist = sort!(collect(keys(hoppings)))
    mkpath(dirname(abspath(path)))

    open(path, "w") do io
        println(io, header)
        println(io, num_wann)
        println(io, length(Rlist))
        _write_ndegen(io, [ndegen[R] for R in Rlist])
        for R in Rlist
            H = hoppings[R]
            size(H) == (num_wann, num_wann) || error("Block for R=$R does not have size ($num_wann, $num_wann)")
            _write_block(io, R, H)
        end
    end
    return path
end

"""
    write_hr_blocks_normalized(path, header, num_wann, hoppings)

Write blocks that are already normalized by primitive `ndegen`, i.e. `raw / ndegen_primitive`.
This is the v1 superham-export path. The emitted supercell file uses `ndegen_sc = 1` for every
block as a semantic marker. Do not use this helper for primitive raw Wannier90 export.
"""
function write_hr_blocks_normalized(
    path::AbstractString,
    header::AbstractString,
    num_wann::Int,
    hoppings::Dict{RKey, Matrix{ComplexF64}},
)
    ndegen = Dict{RKey, Int}(R => 1 for R in keys(hoppings))
    return write_hr_blocks_raw(path, header, num_wann, hoppings, ndegen)
end

end
