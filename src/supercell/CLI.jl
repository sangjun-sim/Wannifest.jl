module SupercellCLI

using LinearAlgebra

include(joinpath(@__DIR__, "Core.jl"))

const InputIO = SupercellCore.InputIO
const Symmetry = SupercellCore.Symmetry
const Service = SupercellCore.Service

export run_main, print_help, parse_args

function print_help()
    println("""
Usage:
  julia main.jl supercell nx ny nz
  julia main.jl supercell a b c, d e f, g h i

Examples:
  julia main.jl supercell 2 2 1
  julia main.jl supercell 2 0 0, 0 2 0, 0 0 1

Defaults:
  input POSCAR
  output POSCAR.supercell
  validate=true, use_symmetry=true, symprec=1.0e-5, angle_tolerance=-1.0, digits=12
""")
end

function print_symmetry_summary(title::String, summary::Symmetry.SymmetrySummary)
    println(title)
    println("  SG: $(summary.international_symbol) ($(summary.spacegroup_number))")
    println("  Hall: $(summary.hall_number)")
    println("  choice: $(summary.choice)")
    println("  point group: $(summary.pointgroup_symbol)")
end

function parse_args(args::Vector{String})
    show_help = any(arg -> arg in ("-h", "--help"), args)

    if show_help
        return (
            config = nothing,
            show_help = show_help,
        )
    end

    matrix_rows = InputIO.parse_matrix_args(args)
    config = InputIO.default_config(matrix_rows)

    return (
        config = config,
        show_help = show_help,
    )
end

function _render_summary(result::Service.SupercellRunResult)
    if result.use_symmetry
        print_symmetry_summary("Input symmetry:", result.input_summary)
        print_symmetry_summary("Output symmetry:", result.output_summary)
        println()
    end

    println("Supercell checks:")
    println("  use_symmetry: $(result.use_symmetry)")
    println("  basis: $(result.basis)")
    println("  det(M): $(round(Int, det(Matrix{Float64}(result.matrix))))")
    println("  expected atoms: $(result.expected_atoms)")
    println("  actual atoms: $(result.actual_atoms)")
    println("  expected volume ratio: $(result.expected_ratio)")
    println("  actual volume ratio: $(round(result.actual_ratio; digits=8))")
    println("  output POSCAR: $(abspath(result.output_path))")

    if !isnothing(result.validation)
        validation = result.validation
        println()
        println("Symmetry mapping:")
        println("  n_operations: $(validation.n_operations)")
        println("  n_checks: $(validation.n_checks)")
        println("  max mismatch: $(validation.max_mismatch) A")
        println("  mean mismatch: $(validation.mean_mismatch) A")
        println("  tolerance: $(validation.tolerance) A")
        println("  mismatch table: $(abspath(validation.mismatch_path))")
    end

    return nothing
end

function run_main(args::Vector{String})::Int
    opts = parse_args(args)
    if opts.show_help
        print_help()
        return 0
    end

    result = Service.run(opts.config)
    _render_summary(result)
    return 0
end

end
