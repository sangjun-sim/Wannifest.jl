module VaspXmlIO

using EzXML
using LinearAlgebra: norm

using ..ReferenceData: ReferenceBandData

export read_vasp_band_data

const KPOINT_TOL = 1e-10

_node_text(node) = strip(EzXML.nodecontent(node))
_node_floats(node) = parse.(Float64, split(_node_text(node)))
_attr(node, name::AbstractString) = haskey(node, name) ? node[name] : ""

function _children(node, name::AbstractString)
    return [child for child in EzXML.eachelement(node) if EzXML.nodename(child) == name]
end

function _descendants(node, name::AbstractString, result=EzXML.Node[])
    for child in EzXML.eachelement(node)
        EzXML.nodename(child) == name && push!(result, child)
        _descendants(child, name, result)
    end
    return result
end

function _named_child(node, tag::AbstractString, name::AbstractString)
    for child in EzXML.eachelement(node)
        EzXML.nodename(child) == tag && _attr(child, "name") == name && return child
    end
    return nothing
end

function _named_descendant(node, tag::AbstractString, name::AbstractString)
    matches = [n for n in _descendants(node, tag) if _attr(n, "name") == name]
    isempty(matches) && return nothing
    return first(matches)
end

function _vasp_bool(text::AbstractString)
    value = lowercase(strip(text))
    value in ("t", "true", ".true.") && return true
    value in ("f", "false", ".false.") && return false
    error("Malformed VASP logical value: $text")
end

function _last_calculation(root)
    calculations = _children(root, "calculation")
    isempty(calculations) && error("vasprun.xml contains no <calculation> block; the run did not finish")
    return last(calculations)
end

function _varray_rows(varray, context::AbstractString)
    rows = [_node_floats(v) for v in _children(varray, "v")]
    isempty(rows) && error("vasprun.xml $context is empty")
    return rows
end

"""
    _real_lattice(root)

Lattice vectors of the cell the bands were computed in, as matrix columns in
angstrom. The `finalpos` structure is preferred; `initialpos` is the same cell
for a single-point run and serves as a fallback.
"""
function _real_lattice(root)
    for name in ("finalpos", "initialpos")
        structure = _named_child(root, "structure", name)
        isnothing(structure) && continue
        basis = _named_descendant(structure, "varray", "basis")
        isnothing(basis) && error("vasprun.xml structure $name has no basis varray")
        rows = _varray_rows(basis, "structure/$name/basis")
        length(rows) == 3 && all(length(r) == 3 for r in rows) ||
            error("vasprun.xml structure $name basis must hold three 3-vectors")
        return hcat(rows...)
    end
    error("vasprun.xml contains neither finalpos nor initialpos structure")
end

function _kpoints(root)
    kpoints = only(_children(root, "kpoints"))
    list = _named_child(kpoints, "varray", "kpointlist")
    isnothing(list) && error("vasprun.xml has no kpoints/kpointlist")
    rows = _varray_rows(list, "kpointlist")
    all(length(r) == 3 for r in rows) || error("vasprun.xml kpointlist entries must have three coordinates")
    return rows
end

function _noncollinear(root)
    parameters = _children(root, "parameters")
    isempty(parameters) && return false
    node = _named_descendant(only(parameters), "i", "LNONCOLLINEAR")
    isnothing(node) && return false
    return _vasp_bool(_node_text(node))
end

# Each spin set holds one <set comment="kpoint n"> per k-point whose <r> rows are
# "eigenvalue occupation". Only the first column is kept.
function _spin_eigenvalues(spin_set, nkpoints::Integer)
    kpoint_sets = _children(spin_set, "set")
    length(kpoint_sets) == nkpoints || error(
        "vasprun.xml eigenvalues hold $(length(kpoint_sets)) k-point sets but kpointlist has $nkpoints entries",
    )
    rows = Vector{Vector{Float64}}(undef, nkpoints)
    for (ik, kset) in enumerate(kpoint_sets)
        values = [first(_node_floats(r)) for r in _children(kset, "r")]
        isempty(values) && error("vasprun.xml eigenvalues k-point $ik has no bands")
        rows[ik] = values
    end
    nbands = length(first(rows))
    all(length(r) == nbands for r in rows) || error("vasprun.xml eigenvalues have inconsistent band counts")
    return Matrix{Float64}(permutedims(reduce(hcat, rows)))
end

function _eigenvalues(calculation, nkpoints::Integer)
    blocks = _children(calculation, "eigenvalues")
    isempty(blocks) && error("vasprun.xml calculation has no <eigenvalues> block")
    array = only(_children(first(blocks), "array"))
    outer = only(_children(array, "set"))
    spin_sets = _children(outer, "set")
    isempty(spin_sets) && error("vasprun.xml eigenvalues contain no spin sets")
    length(spin_sets) <= 2 || error("vasprun.xml eigenvalues contain $(length(spin_sets)) spin sets; expected 1 or 2")
    up = _spin_eigenvalues(spin_sets[1], nkpoints)
    down = length(spin_sets) == 2 ? _spin_eigenvalues(spin_sets[2], nkpoints) : nothing
    return up, down
end

function _fermi_energy(path::AbstractString)
    root = EzXML.root(EzXML.readxml(path))
    calculation = _last_calculation(root)
    node = _named_descendant(calculation, "i", "efermi")
    isnothing(node) && error("vasprun.xml has no efermi: $path")
    return parse(Float64, _node_text(node))
end

"""
    _merge_repeated_kpoints(kpoints, up, down)

VASP line-mode `KPOINTS` repeat the junction of two segments (`... M`, `M ...`),
so the explicit list carries the same k-point twice with identical eigenvalues.
Keep the first copy so the path is a plain corner sequence, as a QE
`crystal_b` path is.
"""
function _merge_repeated_kpoints(kpoints, up, down)
    keep = Int[1]
    for index in 2:length(kpoints)
        norm(kpoints[index] .- kpoints[index - 1]) <= KPOINT_TOL && continue
        push!(keep, index)
    end
    length(keep) == length(kpoints) && return kpoints, up, down
    return kpoints[keep], up[keep, :], isnothing(down) ? nothing : down[keep, :]
end

"""
    read_vasp_band_data(bands_xml, scf_xml)

Read a VASP band-structure `vasprun.xml`. K-points come from `kpointlist` in
reduced coordinates and eigenvalues are already in eV. The Fermi energy is read
from `scf_xml`, which should be the self-consistent run: a line-mode band run
prints an `efermi` too, but VASP recomputes it on the path k-points and it is
not the SCF value.
"""
function read_vasp_band_data(bands_xml::AbstractString, scf_xml::AbstractString)::ReferenceBandData
    root = EzXML.root(EzXML.readxml(bands_xml))
    EzXML.nodename(root) == "modeling" || error("not a vasprun.xml file: $bands_xml")
    kpoints = _kpoints(root)
    up, down = _eigenvalues(_last_calculation(root), length(kpoints))
    kpoints, up, down = _merge_repeated_kpoints(kpoints, up, down)
    return ReferenceBandData(
        kpoints,
        up,
        down,
        _fermi_energy(scf_xml),
        _real_lattice(root),
        !isnothing(down),
        _noncollinear(root),
    )
end

end
