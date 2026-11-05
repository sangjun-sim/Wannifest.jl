module QeWfcIO

using ..Model: QeWfcHeader, QeWavefunction
using ..QeXmlSchema: QeWfcMetadata
using ..Validate

export qe_save_dir, wavefunction_files, sp_wfc_files, soc_wfc_files
export read_wfc_header, read_wfc, expand_gamma

const HEADER_BYTES = 44

qe_save_dir(outdir::AbstractString, prefix::AbstractString) = joinpath(outdir, prefix * ".save")

"""Collected `<stem><k>.dat` files of one save directory, for k = 1:nk."""
function wavefunction_files(savedir::AbstractString, stem::AbstractString, nk::Integer)
    isdir(savedir) || error("QE save directory not found: $(abspath(savedir))")
    names = readdir(savedir)
    pattern = Regex("^" * stem * "([1-9][0-9]*)\\.dat\$")
    files = Dict{Int, String}()
    for name in names
        matched = match(pattern, name)
        isnothing(matched) && continue
        files[parse(Int, matched.captures[1])] = joinpath(savedir, name)
    end
    if isempty(files) && any(name -> endswith(name, ".hdf5") || endswith(name, ".h5"), names)
        error("HDF5 wavefunctions are unsupported; supply QE non-HDF5 wfc*.dat files")
    end
    sort!(collect(keys(files))) == collect(1:nk) ||
        error("Expected $stem indices 1:$nk in $savedir; wavefunctions may be missing or stale")
    return [files[k] for k in 1:nk]
end

sp_wfc_files(savedir, nk) = (wavefunction_files(savedir, "wfcup", nk), wavefunction_files(savedir, "wfcdw", nk))
soc_wfc_files(savedir, nk) = wavefunction_files(savedir, "wfc", nk)

_read_i32(io, swap) = Int(swap ? bswap(read(io, Int32)) : read(io, Int32))

function _read_f64(io, swap)
    bits = read(io, UInt64)
    return reinterpret(Float64, swap ? bswap(bits) : bits)
end

"""One Fortran sequential record: 4-byte length markers around `nbytes` of payload."""
function _record(io, nbytes, swap)
    _read_i32(io, swap) == nbytes || error("Unexpected QE Fortran record size; unsupported ABI or corrupt file")
    nbytes <= filesize(io) - position(io) - 4 || error("Truncated QE Fortran record")
    bytes = read(io, nbytes)
    _read_i32(io, swap) == nbytes || error("QE Fortran record marker mismatch")
    return bytes
end

function _detect_swap(io, path)
    filesize(io) >= 4 || error("Truncated QE wavefunction file: $path")
    marker = read(io, Int32)
    swap = marker != HEADER_BYTES
    (swap ? bswap(marker) : marker) == HEADER_BYTES ||
        error("Unsupported QE binary ABI: expected 4-byte record markers and 44-byte header in $path")
    seekstart(io)
    return swap
end

function _read_header(io, swap)
    header = IOBuffer(_record(io, HEADER_BYTES, swap))
    ik = _read_i32(header, swap)
    kpoint = [_read_f64(header, swap) for _ in 1:3]
    ispin = _read_i32(header, swap)
    logical = _read_i32(header, swap)
    logical in (0, 1, -1) || error("Unsupported Fortran LOGICAL encoding")
    return QeWfcHeader(ik, kpoint, ispin, logical != 0, _read_f64(header, swap))
end

"""Read only the first record (`ik`, Cartesian k in bohr⁻¹, spin, gamma flag, scale factor)."""
function read_wfc_header(path::AbstractString)::QeWfcHeader
    open(path, "r") do io
        return _read_header(io, _detect_swap(io, path))
    end
end

function _miller_matrix(bytes::Vector{UInt8}, swap::Bool, igwx::Int)
    raw = reinterpret(Int32, bytes)
    miller = Matrix{Int}(undef, 3, igwx)
    @inbounds for i in eachindex(miller)
        miller[i] = Int(swap ? bswap(raw[i]) : raw[i])
    end
    return miller
end

function _fill_complex!(dest::AbstractArray{ComplexF64}, bytes::Vector{UInt8}, swap::Bool)
    raw = reinterpret(UInt64, bytes)
    @inbounds for i in eachindex(dest)
        re, im = raw[2i - 1], raw[2i]
        swap && ((re, im) = (bswap(re), bswap(im)))
        dest[i] = complex(reinterpret(Float64, re), reinterpret(Float64, im))
    end
    return dest
end

"""Rebuild the full sphere of a scalar Gamma-only file from C(−G) = conj(C(G))."""
function expand_gamma(miller, coefficients)
    validate_gamma(miller, coefficients)
    keys = [Tuple(g) for g in eachcol(miller)]
    nonzero = findall(!=((0, 0, 0)), keys)
    return hcat(miller, -miller[:, nonzero]), cat(coefficients, conj.(coefficients[nonzero, :, :]); dims=1)
end

"""
    read_wfc(path, meta; npol, expected_ik, ispin=nothing) -> QeWavefunction

Read one collected QE `wfc*.dat` file (non-HDF5, either byte order). Records:
`ik, xk, ispin, gamma_only, scalef` / `ngw, igwx, npol, nbnd` / `b1, b2, b3` /
`mill(3, igwx)` / one record of `npol*igwx` complex coefficients per band, with
the first `igwx` values being spinor component 1.
"""
function read_wfc(path::AbstractString, meta::QeWfcMetadata; npol::Integer, expected_ik::Integer, ispin=nothing)
    open(path, "r") do io
        swap = _detect_swap(io, path)
        header = _read_header(io, swap)
        validate_wfc_header(header, meta.gamma_only; expected_ik, ispin, path)
        dimensions = IOBuffer(_record(io, 16, swap))
        ngw, igwx, polarizations, nbands = [_read_i32(dimensions, swap) for _ in 1:4]
        validate_wfc_dimensions(ngw, igwx, polarizations, nbands; npol, path)
        # Bound dimensions against the actual file before allocating coefficient arrays.
        remaining = big(80) + 8 + 12big(igwx) + big(nbands) * (8 + 16big(igwx) * npol)
        remaining == filesize(io) - position(io) || error("Wavefunction dimensions/file length mismatch in $path")
        basis = IOBuffer(_record(io, 72, swap))
        reciprocal = reshape([_read_f64(basis, swap) for _ in 1:9], 3, 3)
        validate_wfc_reciprocal(reciprocal, meta.cell)
        miller = _miller_matrix(_record(io, 12igwx, swap), swap, igwx)
        validate_miller(miller)
        coefficients = Array{ComplexF64}(undef, igwx, npol, nbands)
        for band in 1:nbands
            _fill_complex!(view(coefficients, :, :, band), _record(io, 16igwx * npol, swap), swap)
        end
        validate_coefficients(coefficients)
        eof(io) || error("Unexpected trailing data in $path")
        if header.gamma_only
            validate_gamma_kpoint(header.kpoint)
            miller, coefficients = expand_gamma(miller, coefficients)
        end
        return QeWavefunction(header, miller, coefficients, reciprocal, ngw, igwx, abspath(path))
    end
end

end
