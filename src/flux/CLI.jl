module FluxCLI

include(joinpath(@__DIR__, "Core.jl"))

const InputIO = FluxCore.InputIO
const Service = FluxCore.Service

export run_main, parse_args, print_usage

function print_usage()
    println("""
Usage:
  julia main.jl flux [--input path/to/input.toml] [--output path] [--html path] [--no-html] [--diagnostic path] [--no-diagnostic] [--validate-roundtrip]
""")
end

function parse_args(args::Vector{String})
    input_path = joinpath(@__DIR__, "input.toml")
    output_hr = nothing
    html_path = nothing
    diagnostic_path = nothing
    make_html = true
    make_diagnostic = nothing
    validate_roundtrip = false
    show_help = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--input"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --input"))
            input_path = args[i]
        elseif arg == "--output"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --output"))
            output_hr = args[i]
        elseif arg == "--html"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --html"))
            html_path = args[i]
            make_html = true
        elseif arg == "--no-html"
            make_html = false
        elseif arg == "--diagnostic"
            i += 1
            i <= length(args) || throw(ArgumentError("Missing value for --diagnostic"))
            diagnostic_path = args[i]
            make_diagnostic = true
        elseif arg == "--no-diagnostic"
            make_diagnostic = false
        elseif arg == "--validate-roundtrip"
            validate_roundtrip = true
        elseif arg in ("-h", "--help")
            show_help = true
        else
            throw(ArgumentError("Unknown argument: $arg"))
        end
        i += 1
    end

    return (
        input_path=String(input_path),
        output_hr=output_hr,
        html_path=html_path,
        diagnostic_path=diagnostic_path,
        make_html=make_html,
        make_diagnostic=make_diagnostic,
        validate_roundtrip=validate_roundtrip,
        show_help=show_help,
    )
end

function _print_summary(result)
    println("Flux run complete.")
    println("Output hr: ", abspath(result.output_hr))
    if !isnothing(result.html_path)
        println("Flux HTML: ", abspath(result.html_path))
    end
    if !isnothing(result.diagnostic_path)
        println("Flux diagnostic: ", abspath(result.diagnostic_path))
    end
    println("Flux seed edges: ", length(result.edges))
    if !isnothing(result.diagnostic)
        println("Plaquette diagnostics: ", length(result.diagnostic.plaquettes))
        if !isempty(result.diagnostic.site_flows)
            status = result.diagnostic.continuity_passed ? "PASS" : "FAIL"
            println("Continuity diagnostic: ", status)
        end
    end
    result.roundtrip_validated && println("Round-trip validation: PASS")
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
        output_hr=opts.output_hr,
        html_path=opts.html_path,
        diagnostic_path=opts.diagnostic_path,
        make_html=opts.make_html,
        make_diagnostic=isnothing(opts.make_diagnostic) ? cfg.diagnostic.enabled : opts.make_diagnostic,
        validate_roundtrip=opts.validate_roundtrip,
    )
    _print_summary(result)
    return 0
end

end
