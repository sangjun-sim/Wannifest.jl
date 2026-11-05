module BinaryIO

export SocTensor, write_tensor, read_tensor, HEADER_BYTES

# Little-endian layout, stdlib only. Every offset is a multiple of 8:
#   0                UInt8[8]    magic "WNFSOC01"
#   8                UInt8[16]   name, NUL padded ("vmat" | "tmat_k" | "tmat_r")
#   24               Int64 × 3   n1 (rows), n2 (columns), nblock (nk or nR)
#   48               Int64       label kind (0 = Float64 fractional k, 1 = Int64 R)
#   56               Int64       reserved (0)
#   64               Float64/Int64[3, nblock]    labels
#   64 + 24·nblock   ComplexF64[n1, n2, nblock]  data, column-major
const MAGIC = Vector{UInt8}(codeunits("WNFSOC01"))
const HEADER_BYTES = 64
const NAME_BYTES = 16
const LABEL_FRACTIONAL_K = 0
const LABEL_LATTICE_R = 1
const HOST_IS_LITTLE_ENDIAN = ENDIAN_BOM == 0x04030201

struct SocTensor
    name::String
    data::Array{ComplexF64, 3}
    labels::Union{Matrix{Float64}, Matrix{Int64}}
end

_words(A::Array) = reinterpret(UInt64, vec(A))

# htol is the identity on little-endian hosts, so the array is written as it is there.
_write_le(io::IO, A::Array) = HOST_IS_LITTLE_ENDIAN ? write(io, A) : write(io, htol.(_words(A)))

function _read_le(io::IO, ::Type{T}, dims) where {T}
    A = Array{T}(undef, dims)
    read!(io, A)
    HOST_IS_LITTLE_ENDIAN || (_words(A) .= ltoh.(_words(A)))
    return A
end

"""
    write_tensor(path, name, data, labels) -> path

Write `data` (`n1 × n2 × nblock`) with one 3-vector label per block: fractional k
for a floating-point `labels` matrix, integer R for an integer one.
"""
function write_tensor(path::AbstractString, name::AbstractString, data::Array{ComplexF64, 3}, labels::AbstractMatrix)
    isascii(name) && 0 < ncodeunits(name) <= NAME_BYTES || error("Tensor name must be 1–$NAME_BYTES ASCII characters: $name")
    all(>(0), size(data)) || error("Tensor $name must not have an empty dimension")
    size(labels) == (3, size(data, 3)) || error("Tensor $name needs a 3×$(size(data, 3)) label matrix")
    kind, stored = eltype(labels) <: Integer ? (LABEL_LATTICE_R, Matrix{Int64}(labels)) :
                                               (LABEL_FRACTIONAL_K, Matrix{Float64}(labels))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        write(io, MAGIC)
        write(io, codeunits(name), zeros(UInt8, NAME_BYTES - ncodeunits(name)))
        for value in (size(data)..., kind, 0)
            write(io, htol(Int64(value)))
        end
        _write_le(io, stored)
        _write_le(io, data)
    end
    return String(path)
end

"""
    read_tensor(path) -> SocTensor

Inverse of `write_tensor`. A wrong magic, header value or file length is an error.
"""
function read_tensor(path::AbstractString)::SocTensor
    open(path, "r") do io
        filesize(io) >= HEADER_BYTES || error("$path: shorter than the $HEADER_BYTES-byte tensor header")
        read(io, length(MAGIC)) == MAGIC || error("$path: not a Wannifest soc tensor (bad magic)")
        name = String(rstrip(String(read(io, NAME_BYTES)), '\0'))
        n1, n2, nblock, kind, reserved = [ltoh(read(io, Int64)) for _ in 1:5]
        n1 > 0 && n2 > 0 && nblock > 0 || error("$path: invalid tensor dimensions $n1×$n2×$nblock")
        kind in (LABEL_FRACTIONAL_K, LABEL_LATTICE_R) && reserved == 0 || error("$path: invalid tensor header")
        expected = big(HEADER_BYTES) + 24big(nblock) + 16big(n1) * n2 * nblock
        filesize(io) == expected || error("$path: file size $(filesize(io)) does not match the header ($expected bytes)")
        labels = kind == LABEL_LATTICE_R ? _read_le(io, Int64, (3, nblock)) : _read_le(io, Float64, (3, nblock))
        return SocTensor(name, _read_le(io, ComplexF64, (n1, n2, nblock)), labels)
    end
end

end
