include(joinpath(@__DIR__, "..", "shared", "InputParsing.jl"))

module InputIO

using TOML

using ..InputParsing: namespaced_root, required_table, required_string, optional_string
using ..InputParsing: optional_existing_path
using ..InputParsing: resolve_path, optional_int, parse_vec3_float, parse_matrix3_int
using ..InputParsing: reject_unknown_keys
using ..Model: OrbitalSpec, RunConfig
using ..SpinLayout

export read_input, resolve_geometry_mode, validate_manual_orbitals, warn_if_manual_centers_are_surrogate

const ALLOWED_GEOMETRY_MODES = Set((:none, :wsvec, :manual_centers, :atomic_assumption))
const ALLOWED_SPIN_KEYS = Set(("layout",))

function _parse_geometry_mode(raw)::Symbol
    raw isa AbstractString || error("geometry.mode must be a string")
    mode = Symbol(strip(String(raw)))
    mode in ALLOWED_GEOMETRY_MODES || error("Unsupported geometry mode: $mode")
    return mode
end

function _parse_orbitals(geometry_tbl)::Vector{OrbitalSpec}
    raw = get(geometry_tbl, "orbitals", Any[])
    raw isa AbstractVector || error("geometry.orbitals must be an array of tables")
    specs = OrbitalSpec[]
    for entry in raw
        entry isa AbstractDict || error("Each geometry.orbitals entry must be a table")
        haskey(entry, "label") || error("Each orbital entry needs a label")
        label = String(entry["label"])
        center_frac = parse_vec3_float(get(entry, "center_frac", nothing), "center_frac")
        atom_index = haskey(entry, "atom_index") ? Int(entry["atom_index"]) : nothing
        push!(specs, OrbitalSpec(label, center_frac, atom_index))
    end
    return specs
end

function _parse_spin_layout(spin_tbl)::Symbol
    reject_unknown_keys(spin_tbl, ALLOWED_SPIN_KEYS, "superham.spin")
    return SpinLayout.parse_layout(get(spin_tbl, "layout", nothing); context="superham.spin.layout")
end

function read_input(path::AbstractString)::RunConfig
    cfg = TOML.parsefile(path)
    base_dir = dirname(abspath(path))
    root = namespaced_root(cfg, "superham")

    files_tbl = required_table(root, "files")

    geometry_tbl = get(root, "geometry", Dict{String, Any}())
    geometry_tbl isa AbstractDict || error("[geometry] must be a table")

    spin_tbl = get(root, "spin", Dict{String, Any}())
    spin_tbl isa AbstractDict || error("[spin] must be a table")

    supercell_tbl = get(root, "supercell", Dict{String, Any}())
    supercell_tbl isa AbstractDict || error("[supercell] must be a table")

    hr_path = resolve_path(base_dir, required_string(files_tbl, "hr"))
    structure_path = resolve_path(base_dir, required_string(files_tbl, "structure"))
    win_path = optional_existing_path(base_dir, files_tbl, "win"; context="superham.files")
    wsvec_path = optional_existing_path(base_dir, files_tbl, "wsvec"; context="superham.files")
    centres_path = optional_existing_path(base_dir, files_tbl, "centres"; context="superham.files")
    output_hr = resolve_path(base_dir, optional_string(files_tbl, "output_hr"))
    spin_layout = _parse_spin_layout(spin_tbl)

    strict_geometry = Bool(get(geometry_tbl, "strict", false))
    geometry_mode = _parse_geometry_mode(get(geometry_tbl, "mode", "none"))
    manual_num_wann = optional_int(geometry_tbl, "manual_num_wann")
    orbital_specs = _parse_orbitals(geometry_tbl)
    supercell_matrix = parse_matrix3_int(get(supercell_tbl, "matrix", [1, 0, 0, 0, 1, 0, 0, 0, 1]), "supercell.matrix")

    return RunConfig(
        hr_path,
        structure_path,
        win_path,
        wsvec_path,
        centres_path,
        output_hr,
        spin_layout,
        strict_geometry,
        geometry_mode,
        supercell_matrix,
        manual_num_wann,
        orbital_specs,
    )
end

resolve_geometry_mode(cfg::RunConfig) = cfg.geometry_mode

function validate_manual_orbitals(cfg::RunConfig, num_wann::Int)
    if !isnothing(cfg.manual_num_wann) && cfg.manual_num_wann != num_wann
        error("manual_num_wann=$(cfg.manual_num_wann) does not match hr.dat num_wann=$num_wann")
    end

    if cfg.geometry_mode in (:manual_centers, :atomic_assumption)
        isempty(cfg.orbital_specs) && error("geometry mode $(cfg.geometry_mode) requires geometry.orbitals in input.toml")
        length(cfg.orbital_specs) == num_wann || error("geometry.orbitals length $(length(cfg.orbital_specs)) does not match num_wann=$num_wann")
    end
    return nothing
end

function warn_if_manual_centers_are_surrogate(cfg::RunConfig)
    if cfg.geometry_mode == :manual_centers
        println(stderr, "Warning: geometry mode $(cfg.geometry_mode) uses user-provided orbital centers; exact Wannier90 geometry is only reproduced if these match the converged Wannier centres.")
    end
    return nothing
end

end
