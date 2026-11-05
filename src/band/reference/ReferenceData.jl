module ReferenceData

export ReferenceBandData, assert_lattice_match

"""
    ReferenceBandData

DFT band data used as the reference in `[band.reference]` comparisons,
independent of the code that produced it. `kpoints_frac` are reduced
coordinates in the reciprocal basis of `real_lattice` (columns are lattice
vectors, in angstrom). `eigenvalues` is `nk × nbands` in eV; `eigenvalues_down`
holds the second collinear channel when `lsda` is true.
"""
struct ReferenceBandData
    kpoints_frac::Vector{Vector{Float64}}
    eigenvalues::Matrix{Float64}
    eigenvalues_down::Union{Nothing, Matrix{Float64}}
    fermi_energy::Float64
    real_lattice::Matrix{Float64}
    lsda::Bool
    noncolin::Bool
end

ReferenceBandData(
    kpoints_frac::Vector{Vector{Float64}},
    eigenvalues::Matrix{Float64},
    fermi_energy::Float64,
    real_lattice::Matrix{Float64},
    lsda::Bool,
    noncolin::Bool,
) = ReferenceBandData(kpoints_frac, eigenvalues, nothing, fermi_energy, real_lattice, lsda, noncolin)

function assert_lattice_match(reference::ReferenceBandData, lattice; rtol::Float64=1e-6, atol::Float64=1e-8)
    isapprox(reference.real_lattice, lattice.real_lattice; rtol=rtol, atol=atol) && return nothing
    error("reference lattice and run.structure lattice differ; do not compare bands from different cells")
end

end
