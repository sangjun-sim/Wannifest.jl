module InputIO

using TOML

using ..InputParsing: numeric_value, integer_value, bool_value, positive_float, table_or_empty
using ..InputParsing: namespaced_root, required_table, required_string, optional_string
using ..InputParsing: required_string_vector
using ..InputParsing: parse_vec3_float, resolve_path, optional_existing_path
using ..InputParsing: reject_unknown_keys, reject_unknown_sibling_tables
using ..Model: BandPlotConfig, CombinedPlotConfig, DosConfig, DosPlotConfig
using ..Model: EnergyConfig, InputFiles, OutputFiles, RunConfig, SpinConfig
using ..ProjectionInput
using ..ObservablesInput
using ..PlotInput
using ..ReferenceInput
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
const DEFAULT_COMPARE_DATA = "outputs/data/band_compare.dat"
const DEFAULT_COMPARE_PLOT = "outputs/plots/band_compare.png"
const ALLOWED_RUN_KEYS = Set(("mode", "hr", "hr_down", "wsvec", "wsvec_down", "kpoints", "structure", "verbose", "hermiticity_tol"))
const ALLOWED_ROOT_TABLES = Set(("run", "output", "dos", "energy", "plot", "band_plot", "dos_plot", "combined_plot", "spin", "projection", "oam", "sam", "reference"))
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

function _optional_dos_bound(tbl, key::AbstractString)
    haskey(tbl, key) || return nothing
    value = tbl[key]
    value isa Real || error("dos.$key must be numeric when provided")
    return Float64(value)
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

function _resolve_kpoints(base_dir::AbstractString, run_tbl, mode::Symbol, reference_enabled::Bool)
    if mode in (:bands, :all) && !reference_enabled
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
    sigma = numeric_value(dos_tbl, "sigma", 0.05; context="dos")
    sigma > 0 || error("dos.sigma must be positive")
    dos_npts = integer_value(dos_tbl, "npts", 1001; context="dos")
    dos_npts >= 2 || error("dos.npts must be >= 2")
    dos_emin = _optional_dos_bound(dos_tbl, "emin")
    dos_emax = _optional_dos_bound(dos_tbl, "emax")
    if !isnothing(dos_emin) && !isnothing(dos_emax)
        dos_emin < dos_emax || error("dos.emin must be smaller than dos.emax")
    end
    return DosConfig(mesh, shift, sigma, dos_npts, dos_emin, dos_emax)
end

function _parse_band_plot_config(band_plot_tbl, collinear::Bool)::BandPlotConfig
    band_linewidth = numeric_value(band_plot_tbl, "linewidth", 1.2; context="band_plot")
    band_linewidth > 0 || error("band_plot.linewidth must be positive")
    band_colors = if haskey(band_plot_tbl, "colors")
        required_string_vector(band_plot_tbl, "colors"; context="band_plot")
    elseif haskey(band_plot_tbl, "color")
        [required_string(band_plot_tbl, "color"; context="band_plot")]
    else
        ["#1f77b4"]
    end
    band_ylabel = optional_string(band_plot_tbl, "ylabel"; default="Energy (eV)", context="band_plot")
    legend = bool_value(band_plot_tbl, "legend", collinear; context="band_plot")
    return BandPlotConfig(band_linewidth, band_colors, band_ylabel, legend)
end

function _warn_channel_naming(hr_path::AbstractString, hr_down_path::AbstractString, verbose::Bool)
    verbose || return nothing
    up_name = lowercase(basename(hr_path))
    down_name = lowercase(basename(hr_down_path))
    # `up`/`dn` tokens, or VASP's LWANNIER90 numbering (`wannier90.1_hr.dat` is up,
    # `wannier90.2_hr.dat` is down).
    up_ok = occursin(r"(^|[_\-.])up([_\-.]|$)", up_name) || occursin(r"\.1_hr\.dat$", up_name)
    down_ok = occursin(r"(^|[_\-.])(dn|down)([_\-.]|$)", down_name) || occursin(r"\.2_hr\.dat$", down_name)
    up_ok && down_ok && return nothing
    @warn "run.hr / run.hr_down 파일명에서 up/dn 채널을 확인할 수 없습니다. " *
          "순서가 바뀌면 스핀 라벨이 통째로 뒤집힙니다" hr=up_name hr_down=down_name
    return nothing
end

function _parse_dos_plot_config(dos_plot_tbl)::DosPlotConfig
    dos_linewidth = numeric_value(dos_plot_tbl, "linewidth", 2.0; context="dos_plot")
    dos_linewidth > 0 || error("dos_plot.linewidth must be positive")
    dos_color = optional_string(dos_plot_tbl, "color"; default="black", context="dos_plot")
    dos_xlabel = optional_string(dos_plot_tbl, "xlabel"; default="Energy (eV)", context="dos_plot")
    dos_ylabel = optional_string(dos_plot_tbl, "ylabel"; default="DOS (states/eV)", context="dos_plot")
    return DosPlotConfig(dos_linewidth, dos_color, dos_xlabel, dos_ylabel)
end

function _parse_spin_layout(spin_tbl, enabled::Bool)::Symbol
    if !haskey(spin_tbl, "layout")
        enabled && error("spin.layout is required when spin.enabled=true (allowed: vasp544, qe)")
        return SpinLayout.DEFAULT_LAYOUT
    end
    raw = optional_string(spin_tbl, "layout"; default="", context="spin")
    return SpinLayout.parse_layout(raw; context="spin.layout")
end

function _parse_spin_config(spin_tbl)::SpinConfig
    reject_unknown_keys(spin_tbl, ALLOWED_SPIN_KEYS, "spin")
    enabled = bool_value(spin_tbl, "enabled", false; context="spin")
    layout = _parse_spin_layout(spin_tbl, enabled)
    if haskey(spin_tbl, "colors")
        raw = required_string_vector(spin_tbl, "colors"; context="spin")
        length(raw) == 2 || error("spin.colors must have exactly 2 entries (up, down), got $(length(raw))")
        return SpinConfig(enabled, layout, (raw[1], raw[2]))
    end
    return SpinConfig(enabled, layout, ("#1f77b4", "#d62728"))
end

function _parse_combined_plot_config(combined_tbl)::CombinedPlotConfig
    dos_width_ratio = numeric_value(combined_tbl, "dos_width_ratio", 0.30; context="combined_plot")
    dos_width_ratio > 0 || error("combined_plot.dos_width_ratio must be positive")
    dos_xlabel = optional_string(combined_tbl, "dos_xlabel"; default="DOS (states/eV)", context="combined_plot")
    dos_ylabel = optional_string(combined_tbl, "dos_ylabel"; default="", context="combined_plot")
    return CombinedPlotConfig(dos_width_ratio, dos_xlabel, dos_ylabel)
end

function read_input(path::AbstractString)::RunConfig
    cfg = TOML.parsefile(path)
    base_dir = dirname(abspath(path))
    root = namespaced_root(cfg, "band")
    reject_unknown_sibling_tables(root, ALLOWED_ROOT_TABLES, "band")

    run_tbl = required_table(root, "run")
    reject_unknown_keys(run_tbl, ALLOWED_RUN_KEYS, "run")
    output_tbl = table_or_empty(root, "output")
    dos_tbl = table_or_empty(root, "dos")
    energy_tbl = get(root, "energy", Dict{String, Any}())
    plot_tbl = table_or_empty(root, "plot")
    band_plot_tbl = table_or_empty(root, "band_plot")
    dos_plot_tbl = table_or_empty(root, "dos_plot")
    combined_tbl = table_or_empty(root, "combined_plot")
    spin_tbl = table_or_empty(root, "spin")
    projection_tbl = table_or_empty(root, "projection")
    oam_tbl = table_or_empty(root, "oam")
    sam_tbl = table_or_empty(root, "sam")
    reference_tbl = table_or_empty(root, "reference")

    mode = _parse_mode(get(run_tbl, "mode", "all"))
    reference = ReferenceInput.parse_reference_config(reference_tbl, base_dir)
    reference.enabled && !(mode in (:bands, :all)) &&
        error("band.reference requires run.mode=\"bands\" or \"all\"")
    verbose = bool_value(run_tbl, "verbose", true; context="run")
    hermiticity_tol = positive_float(run_tbl, "hermiticity_tol", 1e-8; context="run")
    hr_path = _resolve_existing_file(base_dir, run_tbl, "hr", DEFAULT_HR; context="run")
    hr_down_path = optional_existing_path(base_dir, run_tbl, "hr_down"; context="run")
    wsvec_path = optional_existing_path(base_dir, run_tbl, "wsvec"; context="run")
    wsvec_down_path = optional_existing_path(base_dir, run_tbl, "wsvec_down"; context="run")
    files = InputFiles(
        hr_path,
        hr_down_path,
        wsvec_path,
        wsvec_down_path,
        _resolve_kpoints(base_dir, run_tbl, mode, reference.enabled),
        _resolve_structure(base_dir, run_tbl, mode),
    )

    output = OutputFiles(
        _resolve_output(base_dir, output_tbl, "bands_data", DEFAULT_BANDS_DATA),
        _resolve_output(base_dir, output_tbl, "dos_data", DEFAULT_DOS_DATA),
        _resolve_output(base_dir, output_tbl, "bands_plot", DEFAULT_BANDS_PLOT),
        _resolve_output(base_dir, output_tbl, "dos_plot", DEFAULT_DOS_PLOT),
        _resolve_output(base_dir, output_tbl, "combined_plot", DEFAULT_COMBINED_PLOT),
        _resolve_output(base_dir, output_tbl, "compare_data", DEFAULT_COMPARE_DATA),
        _resolve_output(base_dir, output_tbl, "compare_plot", DEFAULT_COMPARE_PLOT),
    )

    energy_shift = haskey(energy_tbl, "shift") ? numeric_value(energy_tbl, "shift", 0.0; context="energy") : 0.0
    reference.enabled && reference.align != :none && energy_shift != 0.0 &&
        error("reference.align and energy.shift cannot both be set; set band.energy.shift = 0.0")
    if !haskey(energy_tbl, "shift") && verbose
        @warn "Energy shift is 0.0. Band/DOS energies are raw eigenvalues with no Fermi level correction."
    end
    plot_config = PlotInput.parse_plot_config(plot_tbl)
    isnothing(hr_down_path) && !isnothing(wsvec_down_path) &&
        error("run.wsvec_down requires run.hr_down")
    if !isnothing(hr_down_path)
        mode == :bands || error("run.hr_down (collinear 두 채널)은 run.mode=\"bands\"만 지원합니다")
        haskey(spin_tbl, "enabled") &&
            error("run.hr_down와 spin.enabled는 함께 쓸 수 없습니다 (spin.enabled는 단일 spinful hr 전용)")
        haskey(spin_tbl, "layout") &&
            error("run.hr_down 모드에서는 spin.layout을 지정할 수 없습니다 (각 파일이 단일 스핀 채널)")
        isnothing(wsvec_down_path) == isnothing(wsvec_path) ||
            error("run.wsvec와 run.wsvec_down은 함께 지정해야 합니다")
        _warn_channel_naming(hr_path, hr_down_path, verbose)
    end
    spin_config = _parse_spin_config(spin_tbl)
    projection_config = ProjectionInput.parse_projection_config(projection_tbl, base_dir)
    oam_config = ObservablesInput.parse_oam_config(oam_tbl)
    sam_config = ObservablesInput.parse_sam_config(sam_tbl)
    if !isnothing(hr_down_path)
        projection_config.enabled && error("run.hr_down 모드는 band.projection을 아직 지원하지 않습니다")
        oam_config.enabled && error("run.hr_down 모드는 band.oam을 아직 지원하지 않습니다")
        sam_config.enabled && error("run.hr_down 모드는 band.sam을 아직 지원하지 않습니다")
        invalid = setdiff(plot_config.targets, (:band, :compare))
        isempty(invalid) ||
            error("run.hr_down 모드의 plot.targets는 band, compare만 허용합니다: $invalid")
    end
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
    PlotInput.validate_plot_targets(
        plot_config.targets,
        mode,
        projection_config,
        oam_config,
        sam_config;
        reference_enabled=reference.enabled,
    )

    return RunConfig(
        mode,
        files,
        output,
        _parse_dos_config(dos_tbl),
        EnergyConfig(energy_shift),
        plot_config,
        _parse_band_plot_config(band_plot_tbl, !isnothing(hr_down_path)),
        _parse_dos_plot_config(dos_plot_tbl),
        _parse_combined_plot_config(combined_tbl),
        spin_config,
        projection_config,
        oam_config,
        sam_config,
        reference,
        hermiticity_tol,
        verbose,
    )
end

end
