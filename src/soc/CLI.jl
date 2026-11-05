module SocCLI

include(joinpath(@__DIR__, "Core.jl"))

const InputIO = SocCore.InputIO
const Service = SocCore.Service

export run_main, parse_args, print_usage

function print_usage()
    println("""
Usage:
  wannifest soc path/to/input.toml [--no-vmat] [--output-dir dir]
""")
end

function parse_args(args::Vector{String})
    input_path = joinpath(@__DIR__, "input.toml")
    output_dir = nothing
    write_vmat = nothing
    show_help = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--input"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --input"))
            input_path = args[i]
        elseif arg == "--output-dir"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --output-dir"))
            output_dir = args[i]
        elseif arg == "--no-vmat"
            write_vmat = false
        elseif arg in ("-h", "--help")
            show_help = true
        elseif i == 1 && !startswith(arg, "-")
            input_path = arg
        else
            throw(ArgumentError("Unknown argument: $arg"))
        end
        i += 1
    end
    return (
        input_path=String(input_path),
        output_dir=output_dir,
        write_vmat=write_vmat,
        show_help=show_help,
    )
end

_format(value) = isnothing(value) ? "n/a" : string(value)

function _print_summary(result)
    println("SOC run complete.")
    println("Output directory: ", abspath(result.output_dir))
    println("V(k): ", isnothing(result.vmat_path) ? "not written (--no-vmat / write_vmat = false)" : basename(result.vmat_path))
    println("T(k): ", basename(result.tmat_k_path))
    println("T(R): ", basename(result.tmat_r_path))
    println("Metadata: ", basename(result.meta_path))
    for label in ("up", "down", "soc")
        haskey(result.eig_max_diff, label) &&
            println(".eig vs QE XML max|ΔE| ($label): ", result.eig_max_diff[label], " eV")
    end
    println("max_k ‖T(k)T(k)† − I‖: ", _format(result.unitarity))
    println("weight0 = ‖T(0)‖² / Σ_R ‖T(R)‖²: ", result.onsite.weight0)
    println("unitarity0 = ‖T(0)T(0)† − I‖: ", _format(result.onsite.unitarity0))
    println("perm0 = ‖|T(0)| − P‖_max: ", _format(result.onsite.perm0))
    correction = result.correction
    if !isnothing(correction)
        println("T(0)† H_sp(R) T(0): ", basename(correction.rotated_path))
        println("ΔH(R) = H_soc(R) − T(0)† H_sp(R) T(0): ", basename(correction.delta_path))
        println("‖ΔH(0)‖_max: ", correction.onsite_max_abs, " eV; max_{R≠0} ‖ΔH(R)‖_max: ", correction.offsite_max_abs, " eV")
        println("ΔH Hermiticity error: ", correction.hermiticity, " eV")
    end
    return nothing
end

function run_main(args::Vector{String})::Int
    opts = parse_args(args)
    if opts.show_help
        print_usage()
        return 0
    end
    cfg = InputIO.read_input(opts.input_path)
    result = Service.run(
        cfg;
        write_vmat=isnothing(opts.write_vmat) ? cfg.run.write_vmat : opts.write_vmat,
        output_dir=isnothing(opts.output_dir) ? cfg.output_dir : abspath(opts.output_dir),
        progress=true,
    )
    _print_summary(result)
    return 0
end

end
