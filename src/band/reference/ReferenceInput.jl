module ReferenceInput

using ..InputParsing: bool_value, positive_float
using ..InputParsing: optional_string, required_string, required_string_vector
using ..InputParsing: parse_float_pair, reject_unknown_keys, resolve_path
using ..Model: ReferenceConfig, disabled_reference_config

export parse_reference_config

const ALLOWED_REFERENCE_KEYS = Set((
    "enabled",
    "source",
    "bands_dir",
    "scf_dir",
    "k_labels",
    "align",
    "window",
    "color",
    "color_down",
    "linewidth",
))

const REFERENCE_SOURCES = Dict(
    :qe_xml => "data-file-schema.xml",
    :vasp => "vasprun.xml",
)

"""
    _schema_path(base_dir, raw, key, source)

Resolve `reference.bands_dir` / `reference.scf_dir` to the band-data file of
`source`: a directory is joined with the source's file name (QE `prefix.save`
holds `data-file-schema.xml`, a VASP run directory holds `vasprun.xml`), and a
path that already names a file is used as is.
"""
function _schema_path(base_dir::AbstractString, raw::AbstractString, key::AbstractString, source::Symbol)
    filename = REFERENCE_SOURCES[source]
    target = resolve_path(base_dir, raw; empty_value="")
    isempty(target) && error("reference.$key cannot be empty")
    isfile(target) && return target
    isdir(target) || error("reference.$key points to missing directory: $target")
    path = joinpath(target, filename)
    isfile(path) || error("reference.$key has no $filename: $path")
    return path
end

function parse_reference_config(tbl, base_dir::AbstractString)::ReferenceConfig
    isempty(tbl) && return disabled_reference_config()
    reject_unknown_keys(tbl, ALLOWED_REFERENCE_KEYS, "reference")
    bool_value(tbl, "enabled", false; context="reference") || return disabled_reference_config()

    source = Symbol(lowercase(required_string(tbl, "source"; context="reference")))
    haskey(REFERENCE_SOURCES, source) || error("reference.source must be \"qe_xml\" or \"vasp\"")
    bands_dir = required_string(tbl, "bands_dir"; context="reference")
    scf_dir = optional_string(tbl, "scf_dir"; default=bands_dir, context="reference")
    labels = required_string_vector(tbl, "k_labels"; context="reference")
    align = Symbol(lowercase(optional_string(tbl, "align"; default="fermi", context="reference")))
    align in (:fermi, :vbm, :none) || error("reference.align must be fermi, vbm, or none")
    window = haskey(tbl, "window") ? parse_float_pair(tbl["window"], "reference.window") : nothing
    !isnothing(window) && window[1] >= window[2] &&
        error("reference.window must have min < max")

    color = optional_string(tbl, "color"; default="gray60", context="reference")
    return ReferenceConfig(
        true,
        source,
        _schema_path(base_dir, bands_dir, "bands_dir", source),
        _schema_path(base_dir, scf_dir, "scf_dir", source),
        labels,
        align,
        window,
        color,
        optional_string(tbl, "color_down"; default=color, context="reference"),
        positive_float(tbl, "linewidth", 0.8; context="reference"),
    )
end

end
