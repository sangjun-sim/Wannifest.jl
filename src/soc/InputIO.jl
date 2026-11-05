module InputIO

using TOML

using ..InputParsing: reject_unknown_keys, reject_unknown_sibling_tables
using ..InputParsing: required_table, optional_table, required_string, optional_string, required_bool, resolve_path
using ..Model: SocQeConfig, SocWannierConfig, SocRunConfig, SocConfig
using ..WannierTypes: RKey

export read_input

const ALLOWED_ROOT_TABLES = Set(("qe", "wannier", "run", "output"))
const ALLOWED_QE_KEYS = Set(("sp_outdir", "sp_prefix", "sp_input", "soc_outdir", "soc_prefix", "soc_input"))
const ALLOWED_WANNIER_KEYS = Set(("dir_up", "seed_up", "dir_down", "seed_down", "dir_soc", "seed_soc"))
const ALLOWED_RUN_KEYS = Set(("atol", "degeneracy_tol", "window_tol", "eig_tol", "write_vmat", "write_soc", "tmat_R"))
const ALLOWED_OUTPUT_KEYS = Set(("dir",))
const DEFAULT_ATOL = 1e-8
const DEFAULT_DEGENERACY_TOL = 1e-4
const DEFAULT_WINDOW_TOL = 1e-3
const DEFAULT_EIG_TOL = 1e-5
const DEFAULT_OUTPUT_DIR = "soc_out"

function _path(base_dir, tbl, key; context)
    raw = required_string(tbl, key; context)
    occursin('\0', raw) && error("$context.$key must not contain NUL characters")
    return resolve_path(base_dir, raw)
end

"""A QE prefix or Wannier90 seedname: a file stem, never a path."""
function _stem(tbl, key; context)
    name = required_string(tbl, key; context)
    !occursin('\0', name) && basename(name) == name && !(name in (".", "..")) ||
        error("$context.$key must be a file stem without a directory, got \"$name\"")
    return name
end

function _tolerance(tbl, key, default; context)
    haskey(tbl, key) || return default
    value = tbl[key]
    value isa Real && !(value isa Bool) && isfinite(value) && value >= 0 ||
        error("$context.$key must be a finite nonnegative number")
    return Float64(value)
end

function _parse_tmat_R(raw)::Union{Symbol, Vector{RKey}}
    if raw isa AbstractString
        lowercase(strip(raw)) == "hr" || error("soc.run.tmat_R must be \"hr\" or an array of [R1, R2, R3]")
        return :hr
    end
    raw isa AbstractVector && !isempty(raw) || error("soc.run.tmat_R must be \"hr\" or a nonempty array of [R1, R2, R3]")
    Rs = RKey[]
    for (i, entry) in enumerate(raw)
        entry isa AbstractVector && length(entry) == 3 && all(x -> x isa Integer && !(x isa Bool), entry) ||
            error("soc.run.tmat_R entry $i must hold three integer lattice coordinates")
        push!(Rs, (Int(entry[1]), Int(entry[2]), Int(entry[3])))
    end
    allunique(Rs) || error("soc.run.tmat_R lists the same R more than once")
    return Rs
end

function _parse_qe(tbl, base_dir)
    reject_unknown_keys(tbl, ALLOWED_QE_KEYS, "soc.qe")
    path(key) = _path(base_dir, tbl, key; context="soc.qe")
    stem(key) = _stem(tbl, key; context="soc.qe")
    return SocQeConfig(path("sp_outdir"), stem("sp_prefix"), path("sp_input"),
                       path("soc_outdir"), stem("soc_prefix"), path("soc_input"))
end

function _parse_wannier(tbl, base_dir)
    reject_unknown_keys(tbl, ALLOWED_WANNIER_KEYS, "soc.wannier")
    path(key) = _path(base_dir, tbl, key; context="soc.wannier")
    stem(key) = _stem(tbl, key; context="soc.wannier")
    return SocWannierConfig(path("dir_up"), stem("seed_up"), path("dir_down"), stem("seed_down"),
                            path("dir_soc"), stem("seed_soc"))
end

function _parse_run(tbl)
    reject_unknown_keys(tbl, ALLOWED_RUN_KEYS, "soc.run")
    tolerance(key, default) = _tolerance(tbl, key, default; context="soc.run")
    write_vmat = haskey(tbl, "write_vmat") ? required_bool(tbl, "write_vmat"; context="soc.run") : true
    write_soc = haskey(tbl, "write_soc") ? required_bool(tbl, "write_soc"; context="soc.run") : false
    tmat_R = haskey(tbl, "tmat_R") ? _parse_tmat_R(tbl["tmat_R"]) : RKey[(0, 0, 0)]
    return SocRunConfig(tolerance("atol", DEFAULT_ATOL), tolerance("degeneracy_tol", DEFAULT_DEGENERACY_TOL),
                        tolerance("window_tol", DEFAULT_WINDOW_TOL), tolerance("eig_tol", DEFAULT_EIG_TOL),
                        write_vmat, write_soc, tmat_R)
end

"""
    read_input(path) -> SocConfig

Read the `[soc.qe]`, `[soc.wannier]`, `[soc.run]` and `[soc.output]` tables.
Relative paths are resolved against the directory of the TOML file; prefixes and
seednames are kept as given. Calculation files are not opened here.
"""
function read_input(path::AbstractString)::SocConfig
    isfile(path) || error("soc input file does not exist: $path")
    cfg = TOML.parsefile(path)
    root = required_table(cfg, "soc")
    reject_unknown_sibling_tables(root, ALLOWED_ROOT_TABLES, "soc")
    reject_unknown_keys(root, ALLOWED_ROOT_TABLES, "soc")
    base_dir = dirname(abspath(path))
    qe = _parse_qe(required_table(root, "qe"; context="soc"), base_dir)
    wannier = _parse_wannier(required_table(root, "wannier"; context="soc"), base_dir)
    run_tbl = something(optional_table(root, "run"; context="soc"), Dict{String, Any}())
    output_tbl = something(optional_table(root, "output"; context="soc"), Dict{String, Any}())
    reject_unknown_keys(output_tbl, ALLOWED_OUTPUT_KEYS, "soc.output")
    output_dir = resolve_path(base_dir, optional_string(output_tbl, "dir"; default=DEFAULT_OUTPUT_DIR, context="soc.output"))
    return SocConfig(abspath(path), qe, wannier, _parse_run(run_tbl), output_dir)
end

end
