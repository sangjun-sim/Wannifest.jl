module InputIO

using TOML

using ..InputParsing: namespaced_root, optional_table, required_table, required_string, optional_string
using ..InputParsing: required_string_vector, required_bool
using ..InputParsing: parse_vec3_float, resolve_path
using ..InputParsing: reject_unknown_keys, reject_unknown_sibling_tables
using ..Model: BandPlotConfig, CombinedPlotConfig, DosConfig, DosPlotConfig
using ..Model: EnergyConfig, InputFiles, OutputFiles, RunConfig, SpinConfig
using ..ProjectionInput
using ..ObservablesInput
using ..PlotInput
using ..SpinLayout

export read_input

const ALLOWED_MODES = Set((:bands, :dos, :all))
const DEFAULT_HR = "wannier90_hr.dat"
const DEFAULT_KPOINTS = "KPOINTS"
const DEFAULT_STRUCTURE = "POSCAR"
const DEFAULT_BANDS_DATA = "outputs/data/bands.dat"
const DEFAULT_DOS_DATA = "outputs/data/dos.dat"
const DEFAULT_BANDS_PLOT = "outputs/plots/bands.png"
const DEFAULT_DOS_PLOT = "outputs/plots/dos.png"
const DEFAULT_COMBINED_PLOT = "outputs/plots/band_dos.png"
const ALLOWED_RUN_KEYS = Set(("mode", "hr", "wsvec", "kpoints", "structure", "verbose", "hermiticity_tol"))
const ALLOWED_ROOT_TABLES = Set(("run", "output", "dos", "energy", "plot", "band_plot", "dos_plot", "combined_plot", "spin", "projection", "oam", "sam"))
const ALLOWED_SPIN_KEYS = Set(("enabled", "layout", "colors"))

function _parse_mode(raw)::Symbol
    raw isa AbstractString || error("run.mode must be a string")
    mode = Symbol(strip(String(raw)))
    mode in ALLOWED_MODES || error("Unsupported run.mode: $mode")
    return mode
end

function _parse_mesh(raw)::NTuple{3, Int}
    raw isa AbstractString || error("dos.mesh must be a string like \"16x16x16\"")
    parts = split(lowercase(strip(String(raw))), 'x')
    length(parts) == 3 || error("dos.mesh must have form n1xn2xn3")
    vals = ntuple(i -> parse(Int, parts[i]), 3)
    any(v -> v <= 0, vals) && error("dos.mesh entries must be positive")
    return vals
end

function _numeric_value(tbl, key::AbstractString, default::Real; context::AbstractString)
    value = get(tbl, key, default)
    value isa Real || error("$context.$key must be numeric")
    return Float64(value)
end

function _integer_value(tbl, key::AbstractString, default::Integer; context::AbstractString)
    value = get(tbl, key, default)
    value isa Integer || error("$context.$key must be an integer")
    return Int(value)
end

function _bool_value(tbl, key::AbstractString, default::Bool; context::AbstractString)
    return haskey(tbl, key) ? required_bool(tbl, key; context=context) : default
end

function _string_value(tbl, key::AbstractString, default::AbstractString; context::AbstractString)
    return optional_string(tbl, key; default=default, context=context)
end

function _optional_dos_bound(tbl, key::AbstractString)
    haskey(tbl, key) || return nothing
    value = tbl[key]
    value isa Real || error("dos.$key must be numeric when provided")
    return Float64(value)
end

function _table_or_empty(root, key::AbstractString)
    tbl = optional_table(root, key)
    return isnothing(tbl) ? Dict{String, Any}() : tbl
end

function _resolve_existing_file(
    base_dir::AbstractString,
    tbl,
    key::AbstractString,
    default::AbstractString;
    context::AbstractString,
)
    raw = optional_string(tbl, key; default=default, context=context)
    path = resolve_path(base_dir, raw; empty_value="")
    isempty(path) && error("$context.$key cannot be empty")
    isfile(path) || error("$context.$key points to missing file: $path")
    return path
end

function _resolve_optional_existing_file(
    base_dir::AbstractString,
    tbl,
    key::AbstractString,
    ;
    context::AbstractString,
)
    haskey(tbl, key) || return nothing
    raw = optional_string(tbl, key; default="", context=context)
    path = resolve_path(base_dir, raw; empty_value=nothing)
    isnothing(path) && return nothing
    isfile(path) || error("$context.$key points to missing file: $path")
    return path
end

function _positive_float(tbl, key::AbstractString, default::Real; context::AbstractString)
    value = _numeric_value(tbl, key, default; context=context)
    value > 0 || error("$context.$key must be positive")
    return value
end

function _resolve_kpoints(base_dir::AbstractString, run_tbl, mode::Symbol)
    if mode in (:bands, :all)
        return _resolve_existing_file(base_dir, run_tbl, "kpoints", DEFAULT_KPOINTS; context="run")
    end
    if haskey(run_tbl, "kpoints")
        raw = optional_string(run_tbl, "kpoints"; default="", context="run")
        path = resolve_path(base_dir, raw; empty_value="")
        !isempty(path) && !isfile(path) && error("run.kpoints points to missing file: $path")
        return path
    end
    return ""
end

function _resolve_structure(base_dir::AbstractString, run_tbl, mode::Symbol)
    if haskey(run_tbl, "structure")
        raw = optional_string(run_tbl, "structure"; default="", context="run")
        path = resolve_path(base_dir, raw; empty_value="")
        !isempty(path) && !isfile(path) && error("run.structure points to missing file: $path")
        return path
    end

    default_path = resolve_path(base_dir, DEFAULT_STRUCTURE; empty_value="")
    if mode in (:bands, :all)
        isfile(default_path) || error(
            "run.structure is omitted and default file is missing: $default_path. " *
            "Set structure = \"\" to use reduced-coordinate distances without a lattice.",
        )
    end
    return isfile(default_path) ? default_path : ""
end

function _resolve_output(base_dir::AbstractString, output_tbl, key::AbstractString, default::AbstractString)
    raw = optional_string(output_tbl, key; default=default, context="output")
    return resolve_path(base_dir, raw; empty_value="")
end

function _parse_dos_config(dos_tbl)::DosConfig
    mesh = _parse_mesh(get(dos_tbl, "mesh", "16x16x16"))
    shift = parse_vec3_float(get(dos_tbl, "shift", [0.0, 0.0, 0.0]), "dos.shift")
    all(x -> 0.0 <= x < 1.0, shift) || error("dos.shift entries must satisfy 0.0 <= shift < 1.0")
    sigma = _numeric_value(dos_tbl, "sigma", 0.05; context="dos")
    sigma > 0 || error("dos.sigma must be positive")
    dos_npts = _integer_value(dos_tbl, "npts", 1001; context="dos")
    dos_npts >= 2 || error("dos.npts must be >= 2")
    dos_emin = _optional_dos_bound(dos_tbl, "emin")
    dos_emax = _optional_dos_bound(dos_tbl, "emax")
    if !isnothing(dos_emin) && !isnothing(dos_emax)
        dos_emin < dos_emax || error("dos.emin must be smaller than dos.emax")
    end
    return DosConfig(mesh, shift, sigma, dos_npts, dos_emin, dos_emax)
end

function _parse_band_plot_config(band_plot_tbl)::BandPlotConfig
    band_linewidth = _numeric_value(band_plot_tbl, "linewidth", 1.2; context="band_plot")
    band_linewidth > 0 || error("band_plot.linewidth must be positive")
    band_colors = if haskey(band_plot_tbl, "colors")
        required_string_vector(band_plot_tbl, "colors"; context="band_plot")
    elseif haskey(band_plot_tbl, "color")
        [required_string(band_plot_tbl, "color"; context="band_plot")]
    else
        ["#1f77b4"]
    end
    band_ylabel = _string_value(band_plot_tbl, "ylabel", "Energy (eV)"; context="band_plot")
    return BandPlotConfig(band_linewidth, band_colors, band_ylabel)
end

function _parse_dos_plot_config(dos_plot_tbl)::DosPlotConfig
    dos_linewidth = _numeric_value(dos_plot_tbl, "linewidth", 2.0; context="dos_plot")
    dos_linewidth > 0 || error("dos_plot.linewidth must be positive")
    dos_color = _string_value(dos_plot_tbl, "color", "black"; context="dos_plot")
    dos_xlabel = _string_value(dos_plot_tbl, "xlabel", "Energy (eV)"; context="dos_plot")
    dos_ylabel = _string_value(dos_plot_tbl, "ylabel", "DOS (states/eV)"; context="dos_plot")
    return DosPlotConfig(dos_linewidth, dos_color, dos_xlabel, dos_ylabel)
end

function _parse_spin_layout(spin_tbl, enabled::Bool)::Symbol
    if !haskey(spin_tbl, "layout")
        enabled && error("spin.layout is required when spin.enabled=true (allowed: vasp544, qe)")
        return SpinLayout.DEFAULT_LAYOUT
    end
    raw = _string_value(spin_tbl, "layout", ""; context="spin")
    return SpinLayout.parse_layout(raw; context="spin.layout")
end

function _parse_spin_config(spin_tbl)::SpinConfig
    reject_unknown_keys(spin_tbl, ALLOWED_SPIN_KEYS, "spin")
    enabled = _bool_value(spin_tbl, "enabled", false; context="spin")
    layout = _parse_spin_layout(spin_tbl, enabled)
    if haskey(spin_tbl, "colors")
        raw = required_string_vector(spin_tbl, "colors"; context="spin")
        length(raw) == 2 || error("spin.colors must have exactly 2 entries (up, down), got $(length(raw))")
        return SpinConfig(enabled, layout, (raw[1], raw[2]))
    end
    return SpinConfig(enabled, layout, ("#1f77b4", "#d62728"))
end

function _parse_combined_plot_config(combined_tbl)::CombinedPlotConfig
    dos_width_ratio = _numeric_value(combined_tbl, "dos_width_ratio", 0.30; context="combined_plot")
    dos_width_ratio > 0 || error("combined_plot.dos_width_ratio must be positive")
    dos_xlabel = _string_value(combined_tbl, "dos_xlabel", "DOS (states/eV)"; context="combined_plot")
    dos_ylabel = _string_value(combined_tbl, "dos_ylabel", ""; context="combined_plot")
    return CombinedPlotConfig(dos_width_ratio, dos_xlabel, dos_ylabel)
end

function read_input(path::AbstractString)::RunConfig
    cfg = TOML.parsefile(path)
    base_dir = dirname(abspath(path))
    root = namespaced_root(cfg, "band")
    reject_unknown_sibling_tables(root, ALLOWED_ROOT_TABLES, "band")

    run_tbl = required_table(root, "run")
    reject_unknown_keys(run_tbl, ALLOWED_RUN_KEYS, "run")
    output_tbl = _table_or_empty(root, "output")
    dos_tbl = _table_or_empty(root, "dos")
    energy_tbl = get(root, "energy", Dict{String, Any}())
    plot_tbl = _table_or_empty(root, "plot")
    band_plot_tbl = _table_or_empty(root, "band_plot")
    dos_plot_tbl = _table_or_empty(root, "dos_plot")
    combined_tbl = _table_or_empty(root, "combined_plot")
    spin_tbl = _table_or_empty(root, "spin")
    projection_tbl = _table_or_empty(root, "projection")
    oam_tbl = _table_or_empty(root, "oam")
    sam_tbl = _table_or_empty(root, "sam")

    mode = _parse_mode(get(run_tbl, "mode", "all"))
    verbose = _bool_value(run_tbl, "verbose", true; context="run")
    hermiticity_tol = _positive_float(run_tbl, "hermiticity_tol", 1e-8; context="run")
    files = InputFiles(
        _resolve_existing_file(base_dir, run_tbl, "hr", DEFAULT_HR; context="run"),
        _resolve_optional_existing_file(base_dir, run_tbl, "wsvec"; context="run"),
        _resolve_kpoints(base_dir, run_tbl, mode),
        _resolve_structure(base_dir, run_tbl, mode),
    )

    output = OutputFiles(
        _resolve_output(base_dir, output_tbl, "bands_data", DEFAULT_BANDS_DATA),
        _resolve_output(base_dir, output_tbl, "dos_data", DEFAULT_DOS_DATA),
        _resolve_output(base_dir, output_tbl, "bands_plot", DEFAULT_BANDS_PLOT),
        _resolve_output(base_dir, output_tbl, "dos_plot", DEFAULT_DOS_PLOT),
        _resolve_output(base_dir, output_tbl, "combined_plot", DEFAULT_COMBINED_PLOT),
    )

    energy_shift = haskey(energy_tbl, "shift") ? _numeric_value(energy_tbl, "shift", 0.0; context="energy") : 0.0
    if !haskey(energy_tbl, "shift") && verbose
        @warn "Energy shift is 0.0. Band/DOS energies are raw eigenvalues with no Fermi level correction."
    end
    plot_config = PlotInput.parse_plot_config(plot_tbl)
    spin_config = _parse_spin_config(spin_tbl)
    projection_config = ProjectionInput.parse_projection_config(projection_tbl, base_dir)
    oam_config = ObservablesInput.parse_oam_config(oam_tbl)
    sam_config = ObservablesInput.parse_sam_config(sam_tbl)
    if oam_config.enabled
        mode in (:bands, :all) || error("band.oam requires run.mode=\"bands\" or \"all\"")
        projection_config.enabled || error("band.oam requires band.projection metadata")
        projection_config.mode == :win_groups ||
            error("band.oam requires band.projection mode=\"win_groups\"")
    end
    if sam_config.enabled
        mode in (:bands, :all) || error("band.sam requires run.mode=\"bands\" or \"all\"")
        projection_config.enabled || error("band.sam requires band.projection metadata")
        projection_config.mode == :win_groups ||
            error("band.sam requires band.projection mode=\"win_groups\"")
        haskey(spin_tbl, "layout") ||
            error("band.sam requires explicit band.spin.layout (allowed: vasp544, qe)")
    end
    PlotInput.validate_plot_targets(plot_config.targets, mode, projection_config, oam_config, sam_config)

    return RunConfig(
        mode,
        files,
        output,
        _parse_dos_config(dos_tbl),
        EnergyConfig(energy_shift),
        plot_config,
        _parse_band_plot_config(band_plot_tbl),
        _parse_dos_plot_config(dos_plot_tbl),
        _parse_combined_plot_config(combined_tbl),
        spin_config,
        projection_config,
        oam_config,
        sam_config,
        hermiticity_tol,
        verbose,
    )
end

end
