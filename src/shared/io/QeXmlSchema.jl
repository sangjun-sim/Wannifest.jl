module QeXmlSchema

using EzXML
using LinearAlgebra: det

export HARTREE_TO_EV, BOHR_TO_ANGSTROM
export node_name, node_text, node_floats, descendants, only_descendant, xml_bool
export qe_cell_bohr, qe_real_lattice, qe_reciprocal_lattice, lsda_band_counts, split_lsda
export QeWfcMetadata, QeEigenvalues, qe_wfc_metadata, qe_eigenvalues

const HARTREE_TO_EV = 27.211386245988
const BOHR_TO_ANGSTROM = 0.529177210903
const SCHEMA_FILE = "data-file-schema.xml"

node_name(node) = last(split(EzXML.nodename(node), ':'))
node_text(node) = strip(EzXML.nodecontent(node))
node_floats(node) = parse.(Float64, split(replace(node_text(node), r"[dD]" => "e")))

function descendants(node, name::AbstractString, result=EzXML.Node[])
    for child in EzXML.eachelement(node)
        node_name(child) == name && push!(result, child)
        descendants(child, name, result)
    end
    return result
end

function only_descendant(node, name::AbstractString; path::AbstractString=name)
    nodes = descendants(node, name)
    length(nodes) == 1 || error("QE XML must contain exactly one $path, found $(length(nodes))")
    return only(nodes)
end

function xml_bool(node)
    value = lowercase(node_text(node))
    value in ("true", ".true.", "t", "1") && return true
    value in ("false", ".false.", "f", "0") && return false
    error("Malformed QE boolean value: $(node_text(node))")
end

"""Cell vectors as columns, in bohr."""
function qe_cell_bohr(root)
    structure = only_descendant(root, "atomic_structure")
    cell = only_descendant(structure, "cell"; path="atomic_structure/cell")
    vectors = [node_floats(only_descendant(cell, "a$i"; path="atomic_structure/cell/a$i")) for i in 1:3]
    all(length(v) == 3 for v in vectors) || error("QE XML cell vectors must have three entries")
    # QE data-file-schema.xml stores cell vectors directly in bohr; `alat` is
    # metadata for other quantities and must not scale these values again.
    return hcat(vectors...)
end

qe_real_lattice(root) = BOHR_TO_ANGSTROM .* qe_cell_bohr(root)

function qe_reciprocal_lattice(root)
    basis = only_descendant(root, "basis_set")
    reciprocal = only_descendant(basis, "reciprocal_lattice"; path="basis_set/reciprocal_lattice")
    vectors = [node_floats(only_descendant(reciprocal, "b$i"; path="reciprocal_lattice/b$i")) for i in 1:3]
    all(length(v) == 3 for v in vectors) || error("QE XML reciprocal vectors must have three entries")
    return hcat(vectors...)
end

function _optional_band_count(band_structure, name::AbstractString)
    nodes = descendants(band_structure, name)
    isempty(nodes) && return nothing
    length(nodes) == 1 || error("QE lsda: expected at most one $name, found $(length(nodes))")
    count = parse(Int, node_text(only(nodes)))
    count >= 0 || error("QE lsda: $name must be non-negative")
    return count
end

function lsda_band_counts(band_structure, total::Integer)
    nbnd_up = _optional_band_count(band_structure, "nbnd_up")
    nbnd_down = _optional_band_count(band_structure, "nbnd_dw")
    if isnothing(nbnd_up) && isnothing(nbnd_down)
        iseven(total) || error("QE lsda: cannot split $total eigenvalues equally without nbnd_up/nbnd_dw")
        return total ÷ 2, total ÷ 2
    elseif isnothing(nbnd_up)
        return total - nbnd_down, nbnd_down
    elseif isnothing(nbnd_down)
        return nbnd_up, total - nbnd_up
    end
    return nbnd_up, nbnd_down
end

"""Split `rows` (nk × total eigenvalues, up block first) into the two LSDA channels."""
function split_lsda(rows::Matrix{Float64}, band_structure)
    total = size(rows, 2)
    nbnd_up, nbnd_down = lsda_band_counts(band_structure, total)
    nbnd_up >= 0 && nbnd_down >= 0 && nbnd_up + nbnd_down == total ||
        error("QE lsda: nbnd_up($nbnd_up) + nbnd_dw($nbnd_down) != eigenvalues 길이($total)")
    return rows[:, 1:nbnd_up], rows[:, (nbnd_up + 1):total]
end

struct QeWfcMetadata
    cell::Matrix{Float64}             # bohr, columns a1 a2 a3
    alat::Float64                     # bohr
    gamma_only::Bool
    nk::Int
    nbnd::Int                         # bands per wfc file (LSDA: spin-up channel)
    nbnd_down::Union{Nothing, Int}    # LSDA spin-down channel
end

struct QeEigenvalues
    kcart::Matrix{Float64}                  # 3×nk, Cartesian bohr⁻¹
    up::Matrix{Float64}                     # nbnd×nk, eV
    down::Union{Nothing, Matrix{Float64}}   # LSDA only
end

function _output_node(savedir::AbstractString)
    path = joinpath(savedir, SCHEMA_FILE)
    isfile(path) || error("QE save directory has no $SCHEMA_FILE: $(abspath(savedir))")
    return only_descendant(EzXML.root(EzXML.readxml(path)), "output"; path="output")
end

function _check_spin_settings(output, soc::Bool)
    magnetization = only_descendant(output, "magnetization"; path="output/magnetization")
    flag(name) = xml_bool(only_descendant(magnetization, name; path="magnetization/$name"))
    lsda, noncolin, spinorbit = flag("lsda"), flag("noncolin"), flag("spinorbit")
    expected = soc ? (!lsda && noncolin && spinorbit) : (lsda && !noncolin && !spinorbit)
    expected || error(soc ? "Expected noncolin=true, spinorbit=true, lsda=false in the SOC QE XML" :
                            "Expected lsda=true, noncolin=false, spinorbit=false in the SP QE XML")
    return nothing
end

function _alat(output)
    structure = only_descendant(output, "atomic_structure")
    haskey(structure, "alat") || error("QE XML atomic_structure has no alat attribute")
    alat = parse(Float64, replace(structure["alat"], r"[dD]" => "e"))
    isfinite(alat) && alat > 0 || error("Invalid QE alat: $alat")
    return alat
end

function _ks_blocks(band_structure)
    nk = parse(Int, node_text(only_descendant(band_structure, "nks"; path="band_structure/nks")))
    nk > 0 || error("Invalid number of k-points in the QE XML: $nk")
    blocks = descendants(band_structure, "ks_energies")
    length(blocks) == nk || error("QE XML declares nks=$nk but holds $(length(blocks)) ks_energies blocks")
    return blocks
end

_block_eigenvalues(block) = node_floats(only_descendant(block, "eigenvalues"; path="ks_energies/eigenvalues"))

"""
    qe_wfc_metadata(savedir; soc) -> QeWfcMetadata

Read what the binary wavefunction reader needs from `data-file-schema.xml` and
check the spin settings: LSDA for `soc=false`, noncollinear spin-orbit for `soc=true`.
"""
function qe_wfc_metadata(savedir::AbstractString; soc::Bool)::QeWfcMetadata
    output = _output_node(savedir)
    _check_spin_settings(output, soc)
    cell = qe_cell_bohr(output)
    all(isfinite, cell) && abs(det(cell)) > 0 || error("Invalid or singular QE cell in $(abspath(savedir))")
    gamma_only = xml_bool(only_descendant(only_descendant(output, "basis_set"), "gamma_only"))
    soc && gamma_only && error("Gamma-only spinor storage is unsupported; use explicit K_POINTS")
    band_structure = only_descendant(output, "band_structure"; path="output/band_structure")
    blocks = _ks_blocks(band_structure)
    total = length(_block_eigenvalues(first(blocks)))
    nbnd, nbnd_down = soc ? (total, nothing) : lsda_band_counts(band_structure, total)
    return QeWfcMetadata(cell, _alat(output), gamma_only, length(blocks), nbnd, nbnd_down)
end

"""
    qe_eigenvalues(savedir; soc) -> QeEigenvalues

Eigenvalues in eV (`nbnd × nk`) with the k-points in Cartesian bohr⁻¹, in XML order,
which is the numbering of the `wfc*.dat` files.
"""
function qe_eigenvalues(savedir::AbstractString; soc::Bool)::QeEigenvalues
    output = _output_node(savedir)
    band_structure = only_descendant(output, "band_structure"; path="output/band_structure")
    blocks = _ks_blocks(band_structure)
    scale = 2pi / _alat(output)
    kcart = Matrix{Float64}(undef, 3, length(blocks))
    eigenvalues = Vector{Vector{Float64}}()
    for (ik, block) in enumerate(blocks)
        k = node_floats(only_descendant(block, "k_point"; path="ks_energies/k_point"))
        length(k) == 3 || error("QE k_point must have three entries")
        kcart[:, ik] = scale .* k
        push!(eigenvalues, HARTREE_TO_EV .* _block_eigenvalues(block))
    end
    all(length(values) == length(first(eigenvalues)) for values in eigenvalues) ||
        error("QE ks_energies blocks have inconsistent band counts")
    rows = Matrix{Float64}(permutedims(reduce(hcat, eigenvalues)))
    soc && return QeEigenvalues(kcart, permutedims(rows), nothing)
    up, down = split_lsda(rows, band_structure)
    return QeEigenvalues(kcart, permutedims(up), permutedims(down))
end

end
