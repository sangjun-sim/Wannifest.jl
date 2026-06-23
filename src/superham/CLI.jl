module SuperhamCLI

include(joinpath(@__DIR__, "Core.jl"))

const InputIO = SuperhamCore.InputIO
const Service = SuperhamCore.Service

export run_main, parse_args, print_usage

function print_usage()
    println("""
Usage:
  julia main.jl superham [--input path/to/input.toml] [--kpoint kx ky kz] [--output-hr path]
""")
end

function parse_args(args::Vector{String})
    input_path = joinpath(@__DIR__, "input.toml")
    kpoint = Float64[0.0, 0.0, 0.0]
    output_hr = nothing
    show_help = false
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--input"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --input"))
            input_path = args[i]
        elseif arg == "--kpoint"
            i + 3 <= length(args) || throw(ArgumentError("Expected three numbers after --kpoint"))
            kpoint = parse.(Float64, args[i + 1:i + 3])
            i += 3
        elseif arg == "--output-hr"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --output-hr"))
            output_hr = args[i]
        elseif arg in ("-h", "--help")
            show_help = true
        else
            throw(ArgumentError("Unknown argument: $arg"))
        end
        i += 1
    end
    return (
        input_path = String(input_path),
        kpoint = kpoint,
        output_hr = output_hr,
        show_help = show_help,
    )
end

function _render_summary(result, input_path::AbstractString)
    println("Input file: ", abspath(input_path))
    println("Geometry mode: ", result.config.geometry_mode)
    println("Geometry source: ", result.geometry_source)
    println("Supercell wsvec output policy: ", result.build_report.wsvec_output_policy)
    println("Supercell center output policy: ", result.build_report.center_output_policy)
    println("Primitive Hermiticity error: ", result.primitive_hermiticity_error)
    println("Supercell multiplicity: ", result.multiplicity)
    println("Primitive orbitals: ", result.model.num_wann)
    println("Supercell orbitals: ", result.super_model.num_wann)
    println("Folded-spectrum check: PASS (max diff = ", result.folded_spectrum_diff, ")")
    println("Primitive eigenvalues at k: ", join(string.(result.primitive_eigenvalues), ", "))
    if !isnothing(result.primitive_wsvec_eigenvalues)
        println("Primitive eigenvalues at k (wsvec-aware phase): ", join(string.(result.primitive_wsvec_eigenvalues), ", "))
    end
    println("Supercell eigenvalues at K: ", join(string.(result.supercell_eigenvalues), ", "))

    if !isnothing(result.output_hr)
        println("Exported supercell hr: ", abspath(String(result.output_hr)))
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
    result = Service.run(cfg; kpoint=opts.kpoint, output_hr=opts.output_hr)
    _render_summary(result, opts.input_path)
    return 0
end

end
