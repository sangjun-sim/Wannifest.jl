module ContourCLI

const CONTOUR_DIR = @__DIR__

include(joinpath(CONTOUR_DIR, "Core.jl"))

const InputIO = ContourCore.InputIO
const Service = ContourCore.Service

export run_main, parse_args

function parse_args(args::Vector{String})
    input_path = ""
    output_dir = nothing
    no_plot = false
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
        elseif arg == "--no-plot"
            no_plot = true
        elseif arg in ("-h", "--help")
            show_help = true
        else
            if i == 1 && !startswith(arg, "-")
                input_path = arg
            else
                throw(ArgumentError("Unknown argument: $arg"))
            end
        end
        i += 1
    end
    if !show_help && isempty(input_path)
        throw(ArgumentError("Missing required --input path/to/input.toml"))
    end
    return (
        input_path=String(input_path),
        output_dir=output_dir,
        no_plot=no_plot,
        show_help=show_help,
    )
end

function run_main(args::Vector{String})::Int
    opts = parse_args(args)
    if opts.show_help
        println("Usage: julia main.jl contour [--input path/to/input.toml] [--output-dir dir] [--no-plot]")
        return 0
    end

    config = InputIO.read_input(opts.input_path; output_dir_override=opts.output_dir)
    make_plot = !opts.no_plot
    result = Service.run(config; make_plot=make_plot)
    Service.print_summary(result; make_plot=make_plot)
    if make_plot && result.config.plot.interactive
        println("Interactive contour plot open. Press Enter to close...")
        readline()
    end
    return 0
end

end
