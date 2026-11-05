module QeXmlIO

using EzXML

using ..QeXmlSchema: HARTREE_TO_EV, node_text, node_floats, descendants, only_descendant, xml_bool
using ..QeXmlSchema: qe_real_lattice, qe_reciprocal_lattice, split_lsda
using ..ReferenceData: ReferenceBandData, assert_lattice_match

export QeBandData, HARTREE_TO_EV, read_qe_band_data, assert_lattice_match

# QE data is stored in the code-independent reference container.
const QeBandData = ReferenceBandData

function _fermi_energy(path::AbstractString)
    root = EzXML.root(EzXML.readxml(path))
    output = only_descendant(root, "output"; path="output")
    fermi = descendants(output, "fermi_energy")
    length(fermi) == 1 && return HARTREE_TO_EV * parse(Float64, node_text(only(fermi)))
    isempty(fermi) || error("QE XML must contain exactly one output/band_structure/fermi_energy, found $(length(fermi))")
    !isempty(descendants(output, "two_fermi_energies")) && error(
        "QE XML에 fermi_energy가 없고 two_fermi_energies(ef_up, ef_dw)만 있습니다. " *
        "고정 자화의 채널별 Fermi 준위에는 공통 shift를 정의할 수 없어 이 참조는 지원하지 않습니다.",
    )
    error("QE XML must contain output/band_structure/fermi_energy")
end

function read_qe_band_data(bands_xml::AbstractString, scf_xml::AbstractString)::QeBandData
    root = EzXML.root(EzXML.readxml(bands_xml))
    output = only_descendant(root, "output"; path="output")
    band_structure = only_descendant(output, "band_structure"; path="output/band_structure")
    lsda = xml_bool(only_descendant(band_structure, "lsda"))
    noncolin = xml_bool(only_descendant(band_structure, "noncolin"))

    b_qe = qe_reciprocal_lattice(output)
    blocks = descendants(band_structure, "ks_energies")
    isempty(blocks) && error("QE XML contains no band_structure/ks_energies blocks")
    kpoints = Vector{Float64}[]
    eigenvalues = Vector{Vector{Float64}}()
    for block in blocks
        k_cart = node_floats(only_descendant(block, "k_point"; path="ks_energies/k_point"))
        length(k_cart) == 3 || error("QE k_point must have three entries")
        push!(kpoints, vec(b_qe \ k_cart))
        push!(eigenvalues, HARTREE_TO_EV .* node_floats(only_descendant(block, "eigenvalues"; path="ks_energies/eigenvalues")))
    end
    nbnd = length(first(eigenvalues))
    all(length(values) == nbnd for values in eigenvalues) ||
        error("QE ks_energies blocks have inconsistent band counts")

    rows = Matrix{Float64}(permutedims(reduce(hcat, eigenvalues)))
    eigenvalues_up, eigenvalues_down = lsda ? split_lsda(rows, band_structure) : (rows, nothing)
    return QeBandData(
        kpoints,
        eigenvalues_up,
        eigenvalues_down,
        _fermi_energy(scf_xml),
        qe_real_lattice(output),
        lsda,
        noncolin,
    )
end

end
