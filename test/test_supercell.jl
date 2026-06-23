const SupercellCLI = Wannifest.SupercellCLI
const SupercellInputIO = Wannifest.SupercellCLI.InputIO
const SupercellService = Wannifest.SupercellCLI.Service
const SupercellTransform = Wannifest.SupercellCLI.SupercellCore.Transform
const SupercellGeometry = Wannifest.SupercellCLI.SupercellCore.SupercellGeometry
const SupercellPoscarIO = Wannifest.SupercellCLI.SupercellCore.PoscarIO
const SupercellValidate = Wannifest.SupercellCLI.SupercellCore.Validate
const SAMPLE_SUPERCELL_DIR = joinpath(WANNIFEST_DIR, "test", "sample", "supercell")

function supercell_run_config(
    input_path::AbstractString,
    output_path::AbstractString;
    matrix_rows::Matrix{Int}=Matrix{Int}(I, 3, 3),
    use_symmetry::Bool=true,
    validate::Bool=false,
)
    return SupercellInputIO.RunConfig(
        String(input_path),
        String(output_path),
        nothing,
        :input,
        use_symmetry,
        nothing,
        1.0e-5,
        -1.0,
        12,
        validate,
        matrix_rows,
    )
end

@testset "supercell CLI input and structure readers" begin
    toml_error = error_message(() -> SupercellService.load_input_cell(joinpath(SAMPLE_SUPERCELL_DIR, "input.toml")))
    @test occursin("TOML structure files are not supported", toml_error)

    poscar_cell = SupercellService.load_input_cell(joinpath(SAMPLE_SUPERCELL_DIR, "POSCAR_rhl_sample"))
    @test size(poscar_cell.lattice) == (3, 3)
    @test SupercellPoscarIO.natoms(poscar_cell) == 18

    diag_cfg = SupercellCLI.parse_args(["2", "2", "1"]).config
    @test diag_cfg.input_path == "POSCAR"
    @test diag_cfg.output_path == "POSCAR.supercell"
    @test diag_cfg.matrix_rows == [2 0 0; 0 2 0; 0 0 1]
    @test diag_cfg.basis == :input
    @test diag_cfg.use_symmetry
    @test diag_cfg.validate
    @test diag_cfg.symprec == 1.0e-5
    @test diag_cfg.angle_tolerance == -1.0
    @test diag_cfg.digits == 12

    @test SupercellInputIO.parse_matrix_args(["1", "0", "0,", "0", "1", "0,", "0", "0", "1"]) == Matrix{Int}(I, 3, 3)

    @test occursin("accepts only matrix arguments", error_message(() -> SupercellCLI.parse_args(["2", "2", "1", "--input", "POSCAR"])))
    @test occursin("Use commas, not semicolons", error_message(() -> SupercellInputIO.parse_matrix_args(["1 0 0; 0 1 0; 0 0 1"])))
    @test occursin("Diagonal supercell form", error_message(() -> SupercellInputIO.parse_matrix_args(["1", "0", "0", "0", "1", "0", "0", "0", "1"])))
    @test occursin("exactly nine entries", error_message(() -> SupercellInputIO.parse_matrix_args(["1 0 0, 0 1 0, 0 0 1"])))
    @test occursin("full-rank", error_message(() -> SupercellInputIO.parse_matrix_args(["1", "0", "0"])))
end

@testset "supercell symmetry labels" begin
    @test SupercellValidate._schonflies_label(Matrix{Int}(I, 3, 3)) == "E"
    @test SupercellValidate._schonflies_label([1 0 0; 0 -1 0; 0 0 -1]) == "C2"
    @test SupercellValidate._schonflies_label(-Matrix{Int}(I, 3, 3)) == "i"
    @test SupercellValidate._schonflies_label([1 0 0; 0 1 0; 0 0 -1]) == "\u03c3"

    label = SupercellValidate._operation_label(Matrix{Int}(I, 3, 3), [-0.0, -0.0, -0.0])
    @test "op1: " * label == "op1: {E | [-0.000000,-0.000000,-0.000000]}"
end

@testset "supercell transform core" begin
    cell = SupercellPoscarIO.StructureCell(
        "toy",
        Matrix{Float64}(I, 3, 3),
        reshape([0.0, 0.0, 0.0], 3, 1),
        ["X"],
        [1],
        "inline",
    )
    matrix = [2 0 0; 0 1 0; 0 0 1]
    super = SupercellTransform.build_supercell(cell, matrix)
    @test SupercellPoscarIO.natoms(super) == 2
    @test super.lattice ≈ [2.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test SupercellTransform.supercell_multiplicity(matrix) == 2

    shear_rows = [1 1 0; 0 1 0; 0 0 1]
    geom = SupercellGeometry.from_user_matrix(shear_rows)
    @test geom.input_matrix == shear_rows
    @test geom.lattice_matrix == [1 0 0; 1 1 0; 0 0 1]
    @test geom.multiplicity == 1

    shear = SupercellTransform.build_supercell(cell, shear_rows)
    @test SupercellPoscarIO.natoms(shear) == 1
    @test shear.lattice ≈ [1.0 0.0 0.0; 1.0 1.0 0.0; 0.0 0.0 1.0]
end

@testset "supercell service preserves input basis" begin
    mktempdir() do dir
        poscar = joinpath(dir, "POSCAR")
        write(poscar, """
Input basis cell
1.0
1.0 0.0 0.0
0.4 1.2 0.0
0.1 0.2 2.0
X
1
Direct
0.0 0.0 0.0
""")
        output_path = joinpath(dir, "POSCAR.out")
        cfg = supercell_run_config(poscar, output_path; use_symmetry=false, validate=false)
        result = SupercellService.run(cfg)
        output_cell = SupercellPoscarIO.read_poscar(output_path)
        @test result.basis == :input
        @test !result.use_symmetry
        @test output_cell.lattice ≈ [1.0 0.4 0.1; 0.0 1.2 0.2; 0.0 0.0 2.0]

        output_with_summary = joinpath(dir, "POSCAR.summary")
        summary_cfg = supercell_run_config(poscar, output_with_summary; use_symmetry=true, validate=false)
        summary_result = SupercellService.run(summary_cfg)
        summary_cell = SupercellPoscarIO.read_poscar(output_with_summary)
        @test summary_result.use_symmetry
        @test summary_result.basis == :input
        @test summary_cell.lattice ≈ [1.0 0.4 0.1; 0.0 1.2 0.2; 0.0 0.0 2.0]
    end
end
