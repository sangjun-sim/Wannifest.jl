module Validate

using LinearAlgebra

using ..Model: QeWfcHeader, QeWavefunction, WannierChannel, num_kpoints

export validate_atol, validate_points
export validate_wfc_header, validate_wfc_dimensions, validate_wfc_reciprocal
export validate_miller, validate_coefficients, validate_gamma, validate_gamma_kpoint
export validate_spin_pair, validate_same_lattice, miller_index, validate_tmat_inputs

# Semantic checks shared by the QE readers and the V/T calculators. Successful checks return nothing.

function validate_atol(atol)
    atol isa Real && !(atol isa Bool) && isfinite(atol) && atol >= 0 ||
        error("atol must be a finite nonnegative real number")
    return nothing
end

function validate_points(points, label::AbstractString)
    points isa AbstractMatrix && size(points, 1) == 3 && size(points, 2) > 0 ||
        error("$label must be a nonempty 3×N matrix")
    all(x -> x isa Real && !(x isa Bool) && isfinite(x), points) ||
        error("$label must contain finite real coordinates")
    return nothing
end

function validate_wfc_header(header::QeWfcHeader, gamma_only::Bool; expected_ik, ispin=nothing, path="wavefunction")
    header.ik > 0 && all(isfinite, header.kpoint) || error("Invalid k-point header in $path")
    header.ik == expected_ik ||
        error("Wavefunction header k index $(header.ik) does not match its filename ($expected_ik): $path")
    header.gamma_only == gamma_only || error("XML/binary gamma_only mismatch in $path")
    header.scale_factor == 1.0 || error("Only pw.x wavefunctions with scale_factor=1 are supported: $path")
    isnothing(ispin) || header.ispin == ispin || error("Wrong spin channel in $path")
    return nothing
end

function validate_wfc_dimensions(ngw, igwx, polarizations, nbands; npol, path="wavefunction")
    ngw > 0 && igwx > 0 && nbands > 0 || error("Invalid wavefunction dimensions in $path")
    polarizations == npol || error("Expected npol=$npol in $path")
    return nothing
end

function validate_wfc_reciprocal(reciprocal, cell)
    all(isfinite, reciprocal) && isapprox(cell' * reciprocal, 2pi * I; atol=1e-7, rtol=1e-7) ||
        error("Binary reciprocal lattice is inconsistent with the XML cell")
    return nothing
end

function validate_miller(miller)
    length(Set(Tuple(g) for g in eachcol(miller))) == size(miller, 2) || error("Duplicate G vectors")
    return nothing
end

function validate_coefficients(coefficients)
    all(isfinite, coefficients) || error("Nonfinite wavefunction coefficients")
    return nothing
end

function validate_gamma(miller, coefficients)
    size(coefficients, 2) == 1 || error("Gamma-only spinors are unsupported")
    keys = [Tuple(g) for g in eachcol(miller)]
    stored = Set(keys)
    for i in findall(!=((0, 0, 0)), keys)
        Tuple(-miller[:, i]) in stored && error("Gamma file contains both G and -G")
    end
    zero_index = findfirst(==((0, 0, 0)), keys)
    if !isnothing(zero_index)
        all(z -> abs(imag(z)) <= 1e-10, coefficients[zero_index, :, :]) ||
            error("Gamma G=0 coefficients must be real")
    end
    return nothing
end

function validate_gamma_kpoint(kpoint)
    all(x -> abs(x) <= 1e-12, kpoint) || error("Gamma-only file has nonzero k-point")
    return nothing
end

function validate_spin_pair(up::QeWavefunction, down::QeWavefunction)
    isapprox(up.header.kpoint, down.header.kpoint; atol=1e-10, rtol=0) ||
        error("SP up/down k-points differ: $(up.source) vs $(down.source)")
    isapprox(up.reciprocal, down.reciprocal; atol=1e-10, rtol=0) ||
        error("SP up/down reciprocal lattices differ: $(up.source) vs $(down.source)")
    return nothing
end

"""SP and SOC must share one cell: plane waves are matched by their Miller indices."""
function validate_same_lattice(sp_cell, soc_cell; atol)
    all(abs.(sp_cell .- soc_cell) .<= atol) ||
        error("SP and SOC calculations use different cells; V is only defined for a common lattice")
    return nothing
end

"""Check one wavefunction and return the Miller triplet => row lookup used for the overlaps."""
function miller_index(wave::QeWavefunction, npol::Integer, label::AbstractString)
    miller, coefficients = wave.miller, wave.coefficients
    size(miller, 1) == 3 || error("$label.miller must be a 3×nG integer matrix")
    size(miller, 2) > 0 || error("$label has an empty plane-wave basis")
    size(coefficients, 1) == size(miller, 2) && size(coefficients, 2) == npol && size(coefficients, 3) > 0 ||
        error("$label.coefficients must be nG×$npol×nband complex data with nonempty bands")
    all(isfinite, coefficients) || error("$label.coefficients must be finite")
    indices = Dict{NTuple{3, Int}, Int}()
    for g in axes(miller, 2)
        key = (miller[1, g], miller[2, g], miller[3, g])
        haskey(indices, key) && error("$label has duplicate Miller indices $key")
        indices[key] = g
    end
    return indices
end

function _validate_channel(channel::WannierChannel)
    all(isfinite, channel.U) || error("$(channel.label) U matrices must contain finite values")
    size(channel.U) == (channel.num_bands, channel.num_wann, num_kpoints(channel)) ||
        error("$(channel.label) U matrices do not match num_bands×num_wann×nk")
    return nothing
end

"""Dimension and sanity checks for T(k) = blockdiag(Ũ_up, Ũ_down)† V Ũ_soc."""
function validate_tmat_inputs(up::WannierChannel, down::WannierChannel, soc::WannierChannel, V, index)
    foreach(_validate_channel, (up, down, soc))
    nk = num_kpoints(up)
    num_kpoints(down) == nk && num_kpoints(soc) == nk ||
        error("up, down and SOC Wannier k-point counts must match")
    ndims(V) == 3 && size(V) == (up.num_bands + down.num_bands, soc.num_bands, nk) ||
        error("V must be $(up.num_bands + down.num_bands) SP band rows × $(soc.num_bands) SOC band columns × $nk k-points")
    all(isfinite, V) || error("V must contain finite values")
    for (name, values) in (("down", index.down), ("soc", index.soc))
        length(values) == nk && sort(values) == 1:nk ||
            error("$name k-point correspondence must be a permutation of 1:$nk")
    end
    return nothing
end

end
