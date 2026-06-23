module Wannifest

export run_wannifest

const WANNIFEST_DIR = normpath(joinpath(@__DIR__, ".."))
const SRC_DIR = @__DIR__
const STARTUP_FILE_DISABLED = 2

import TOML

include(joinpath(SRC_DIR, "band", "CLI.jl"))
include(joinpath(SRC_DIR, "contour", "CLI.jl"))
include(joinpath(SRC_DIR, "superham", "CLI.jl"))
include(joinpath(SRC_DIR, "supercell", "CLI.jl"))
include(joinpath(SRC_DIR, "flux", "CLI.jl"))

const MODULE_COMMANDS = (
    "band",
    "contour",
    "superham",
    "supercell",
    "flux",
)

const TOML_BACKED_MODULE_COMMANDS = (
    "band",
    "contour",
    "superham",
    "flux",
)

function print_wannifest_usage(io::IO=stdout)
    println(io, """
Usage:
  julia main.jl <command> [options...]
  julia main.jl --input path/to/input.toml [options...]
  julia main.jl path/to/input.toml [options...]

Commands:
  band        Band/DOS calculation
  contour     2D k-plane energy surface and contour plots
  superham    Build a supercell Hamiltonian
  supercell   Build a structural POSCAR supercell
  flux        Add directed flux terms to a Wannier hr.dat file

Examples:
  julia main.jl band --input examples/flake/graphene/input.band.toml --no-plot
  julia main.jl contour --input examples/contour/graphene/input.toml --no-plot
  julia main.jl superham --input examples/flake/graphene/input.superham.toml
  julia main.jl supercell 2 0 0, 0 2 0, 0 0 1
  julia main.jl flux --input examples/flux/graphene/input.toml --no-html --no-diagnostic
""")
end

_is_command(command::AbstractString)::Bool = String(command) in MODULE_COMMANDS

function normalize_input_args(args::Vector{String})
    isempty(args) && return args
    first = args[1]
    if _is_command(first) || first in ("-h", "--help", "help")
        return args
    end
    if startswith(first, "-") || !endswith(lowercase(first), ".toml")
        return args
    end
    return ["--input", first, args[2:end]...]
end

function extract_input_path(args::Vector{String})::Union{Nothing, String}
    i = 1
    while i <= length(args)
        if args[i] == "--input"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --input"))
            return String(args[i])
        end
        i += 1
    end
    return nothing
end

function _single_module_key(cfg)::Union{Nothing, String}
    hits = String[]
    for key in keys(cfg)
        key_text = String(key)
        key_text in TOML_BACKED_MODULE_COMMANDS && push!(hits, key_text)
    end
    unique!(hits)
    if length(hits) > 1
        error("Ambiguous input file: found multiple Wannifest namespaces ($(join(hits, ", ")))")
    end
    return isempty(hits) ? nothing : only(hits)
end

function infer_command_from_input(path::AbstractString)::String
    parent = basename(dirname(abspath(path)))
    parent in TOML_BACKED_MODULE_COMMANDS && return parent

    cfg = TOML.parsefile(path)
    module_key = _single_module_key(cfg)
    !isnothing(module_key) && return module_key

    if haskey(cfg, "supercell") && haskey(cfg, "geometry") && haskey(cfg, "files")
        @warn "Inferring Wannifest command from legacy non-namespaced TOML; prefer [superham.*] tables or pass the command explicitly." command="superham"
        return "superham"
    elseif haskey(cfg, "run") && haskey(cfg, "output") && haskey(cfg, "dos") &&
           haskey(cfg, "band_plot")
        @warn "Inferring Wannifest command from legacy non-namespaced TOML; prefer [band.*] tables or pass the command explicitly." command="band"
        return "band"
    end

    error("Could not infer Wannifest command from input file: $(abspath(path))")
end

function resolve_command(args::Vector{String})
    normalized = normalize_input_args(args)
    isempty(normalized) && return nothing, normalized

    first = normalized[1]
    if _is_command(first)
        return String(first), normalized[2:end]
    elseif first in ("-h", "--help", "help")
        return first, String[]
    end

    input_path = extract_input_path(normalized)
    isnothing(input_path) && return nothing, normalized
    return infer_command_from_input(input_path), normalized
end

function _invoke_cli(module_name::Symbol, args::Vector{String})
    cli_module = Base.invokelatest(getfield, @__MODULE__, module_name)
    run_main_fn = Base.invokelatest(getfield, cli_module, :run_main)
    return Base.invokelatest(run_main_fn, args)
end

function _print_command_help(command::AbstractString)
    if command == "band"
        println("Usage: julia main.jl band [--input path/to/input.toml] [--no-plot]")
    elseif command == "contour"
        println("Usage: julia main.jl contour [--input path/to/input.toml] [--output-dir dir] [--no-plot]")
    elseif command == "superham"
        println("Usage: julia main.jl superham [--input path/to/input.toml] [--kpoint kx ky kz] [--output-hr path]")
    elseif command == "supercell"
        println("""
Usage:
  julia main.jl supercell nx ny nz
  julia main.jl supercell a b c, d e f, g h i
""")
    elseif command == "flux"
        println("Usage: julia main.jl flux [--input path/to/input.toml] [--output path] [--html path] [--no-html] [--diagnostic path] [--no-diagnostic] [--validate-roundtrip]")
    else
        print_wannifest_usage()
    end
    return 0
end

_wants_help(args::Vector{String})::Bool = any(arg -> arg in ("-h", "--help", "help"), args)

const RUNNERS = Dict{String, Function}(
    "band" => args -> _invoke_cli(:BandCLI, args),
    "contour" => args -> _invoke_cli(:ContourCLI, args),
    "superham" => args -> _invoke_cli(:SuperhamCLI, args),
    "supercell" => args -> _invoke_cli(:SupercellCLI, args),
    "flux" => args -> _invoke_cli(:FluxCLI, args),
)

function run_command(command::AbstractString, args::Vector{String})::Int
    if command in ("-h", "--help", "help")
        print_wannifest_usage()
        return 0
    end
    _wants_help(args) && return _print_command_help(String(command))
    runner = get(RUNNERS, String(command), nothing)
    isnothing(runner) && throw(ArgumentError("Unknown Wannifest command: $command"))
    return runner(args)
end

function run_wannifest(args::Vector{String})::Int
    command, rest = resolve_command(args)
    isnothing(command) && (print_wannifest_usage(); return 1)
    return run_command(command, rest)
end

end
