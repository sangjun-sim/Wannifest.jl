module Vmat

using ..Model: QeWavefunction, WfcFiles, KpointIndex
using ..QeWfcIO: read_wfc, read_wfc_header
using ..QeXmlSchema: QeWfcMetadata
using ..Validate: miller_index, validate_spin_pair

export UP_WEIGHTS, DOWN_WEIGHTS, vmat_at_k, wfc_kpoints, calculate_vmat

# ⟨↑| and ⟨↓| along z select spinor components 1 and 2. A tilted quantisation axis
# only changes these weights: (cos(θ/2), e^{−iφ} sin(θ/2)) and (−e^{+iφ} sin(θ/2), cos(θ/2)).
const UP_WEIGHTS = (1.0 + 0.0im, 0.0 + 0.0im)
const DOWN_WEIGHTS = (0.0 + 0.0im, 1.0 + 0.0im)

function _check_kept(kept, wave::QeWavefunction, label)
    nbands = size(wave.coefficients, 3)
    all(band -> 1 <= band <= nbands, kept) ||
        error("$label keeps a band outside the $nbands bands stored in $(wave.source)")
    return nothing
end

"""
Rows `row_offset+1 : row_offset+length(kept_sp)` of V:
`V[m, n] = Σ_G conj(C_SP[G, m]) (w₁ C_SOC[G, 1, n] + w₂ C_SOC[G, 2, n])`.
Plane waves are paired by Miller triplet; one missing from either cutoff contributes zero.
"""
function _overlap!(V, row_offset, scalar::QeWavefunction, spinor::QeWavefunction, spinor_basis, weights, kept_sp, kept_soc)
    scalar_rows, spinor_rows = Int[], Int[]
    for g in axes(scalar.miller, 2)
        row = get(spinor_basis, (scalar.miller[1, g], scalar.miller[2, g], scalar.miller[3, g]), 0)
        row == 0 && continue
        push!(scalar_rows, g)
        push!(spinor_rows, row)
    end
    bra = adjoint(scalar.coefficients[scalar_rows, 1, kept_sp])
    rows = (row_offset + 1):(row_offset + length(kept_sp))
    for (component, weight) in enumerate(weights)
        iszero(weight) && continue
        V[rows, :] .+= weight .* (bra * spinor.coefficients[spinor_rows, component, kept_soc])
    end
    return V
end

"""
    vmat_at_k(up, down, spinor, kept_up, kept_down, kept_soc) -> Matrix{ComplexF64}

`V[m, n] = ⟨ψ^SP_m(k)|ψ^SOC_n(k)⟩` at one k-point: kept up bands, then kept down
bands as rows, kept SOC bands as columns. No normalization or orthogonalization.
"""
function vmat_at_k(up::QeWavefunction, down::QeWavefunction, spinor::QeWavefunction,
                   kept_up::Vector{Int}, kept_down::Vector{Int}, kept_soc::Vector{Int})
    miller_index(up, 1, "SP up")
    miller_index(down, 1, "SP down")
    spinor_basis = miller_index(spinor, 2, "SOC spinor")
    _check_kept(kept_up, up, "up")
    _check_kept(kept_down, down, "down")
    _check_kept(kept_soc, spinor, "soc")
    V = zeros(ComplexF64, length(kept_up) + length(kept_down), length(kept_soc))
    _overlap!(V, 0, up, spinor, spinor_basis, UP_WEIGHTS, kept_up, kept_soc)
    _overlap!(V, length(kept_up), down, spinor, spinor_basis, DOWN_WEIGHTS, kept_down, kept_soc)
    return V
end

"""Header-only pass: Cartesian k (bohr⁻¹) of every file, as 3×nfiles."""
function wfc_kpoints(files::Vector{String})
    kpoints = Matrix{Float64}(undef, 3, length(files))
    for (i, path) in enumerate(files)
        header = read_wfc_header(path)
        header.ik == i || error("Wavefunction header k index $(header.ik) does not match its filename ($i): $path")
        kpoints[:, i] = header.kpoint
    end
    return kpoints
end

"""
    calculate_vmat(files, sp_meta, soc_meta, index, kept_up, kept_down, kept_soc; progress=false)

V(k) in the k order of the up Wannier channel, as `(nb_up+nb_down) × nb_soc × nk`.
Only the wavefunctions of one k-point are in memory at a time, and QE k-points
that the Wannier list does not use are never read.
"""
function calculate_vmat(files::WfcFiles, sp_meta::QeWfcMetadata, soc_meta::QeWfcMetadata, index::KpointIndex,
                        kept_up::Vector{Int}, kept_down::Vector{Int}, kept_soc::Vector{Int}; progress::Bool=false)
    nk = length(index.sp_file)
    V = Array{ComplexF64}(undef, length(kept_up) + length(kept_down), length(kept_soc), nk)
    for k in 1:nk
        i, j = index.sp_file[k], index.soc_file[k]
        up = read_wfc(files.up[i], sp_meta; npol=1, expected_ik=i, ispin=1)
        down = read_wfc(files.down[i], sp_meta; npol=1, expected_ik=i, ispin=2)
        validate_spin_pair(up, down)
        spinor = read_wfc(files.soc[j], soc_meta; npol=2, expected_ik=j)
        V[:, :, k] = vmat_at_k(up, down, spinor, kept_up, kept_down, kept_soc)
        progress && (k == nk || k % max(1, nk ÷ 10) == 0) && println(stderr, "  V(k): $k / $nk k-points")
    end
    return V
end

end
