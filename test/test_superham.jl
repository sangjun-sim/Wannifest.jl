@testset "superham atomic center warnings and wrapping" begin
    mktempdir() do dir
        poscar = joinpath(dir, "POSCAR")
        write(poscar, """
Two atom cell
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Si
2
Direct
0.0 0.0 0.0
0.5 0.0 0.0
""")
        cell = SuperhamCore.PoscarIO.read_poscar(poscar)
        specs = [
            SuperhamCore.Model.OrbitalSpec("Si1_s_up", (1.0, 0.0, 0.0), 1),
            SuperhamCore.Model.OrbitalSpec("Si1_s_dn", (0.0, 0.0, 0.0), 1),
            SuperhamCore.Model.OrbitalSpec("Si2_s_up", (0.5, 0.0, 0.0), 2),
            SuperhamCore.Model.OrbitalSpec("Si2_s_dn", (0.5, 0.0, 0.0), 2),
        ]
        centers = @test_logs (:warn, r"atomic-centered approximation.*POSCAR.*wannier90\.win.*converged Wannier centers") CenterIO.build_atomic_centers(
            poscar,
            Matrix{Float64}(cell.lattice),
            specs;
            cell=cell,
        )
        @test centers.centers_frac[:, 1] ≈ centers.centers_frac[:, 2]
        @test centers.centers_frac[:, 3] ≈ centers.centers_frac[:, 4]
        @test centers.centers_frac[:, 1] ≈ [0.0, 0.0, 0.0]

        block_specs = [
            SuperhamCore.Model.OrbitalSpec("Si1_s_up", (0.0, 0.0, 0.0), 1),
            SuperhamCore.Model.OrbitalSpec("Si2_s_up", (0.5, 0.0, 0.0), 2),
            SuperhamCore.Model.OrbitalSpec("Si1_s_dn", (0.0, 0.0, 0.0), 1),
            SuperhamCore.Model.OrbitalSpec("Si2_s_dn", (0.5, 0.0, 0.0), 2),
        ]
        block_centers = @test_logs (:warn, r"atomic-centered approximation.*POSCAR.*wannier90\.win.*converged Wannier centers") CenterIO.build_atomic_centers(
            poscar,
            Matrix{Float64}(cell.lattice),
            block_specs;
            cell=cell,
        )
        @test block_centers.centers_frac[:, 1] ≈ block_centers.centers_frac[:, 3]
        @test block_centers.centers_frac[:, 2] ≈ block_centers.centers_frac[:, 4]

        canonical = CenterIO.canonicalize_centers(block_centers, :vasp544)
        @test canonical.labels == ["Si1_s_up", "Si1_s_dn", "Si2_s_up", "Si2_s_dn"]
        @test canonical.centers_frac[:, 1] ≈ canonical.centers_frac[:, 2]
        @test canonical.centers_frac[:, 3] ≈ canonical.centers_frac[:, 4]

        centers_frac6 = zeros(Float64, 3, 6)
        centers_frac6[1, :] .= (1:6) ./ 10
        centers6 = SuperhamCore.Model.CenterTable(
            centers_frac6,
            copy(centers_frac6),
            ["A_up", "B_up", "C_up", "A_dn", "B_dn", "C_dn"],
            "inline",
            :manual_centers,
        )
        canonical6 = CenterIO.canonicalize_centers(centers6, :vasp544)
        @test canonical6.labels == ["A_up", "A_dn", "B_up", "B_dn", "C_up", "C_dn"]
        @test canonical6.labels != ["A_up", "B_dn", "A_dn", "C_up", "B_up", "C_dn"]
        @test canonical6.centers_frac[1, :] ≈ [0.1, 0.4, 0.2, 0.5, 0.3, 0.6]

        edge_centers = SuperhamCore.Model.CenterTable(
            reshape([-1.0e-15, 1.0 - 1.0e-15, 1.0 + 1.0e-15], 3, 1),
            zeros(Float64, 3, 1),
            ["edge"],
            "inline",
            :manual_centers,
        )
        edge_model = SuperhamCore.Model.HrModel(
            "edge wrap",
            Matrix{Float64}(I, 3, 3),
            2π .* Matrix{Float64}(I, 3, 3),
            1,
            Dict{SuperhamCore.Model.RKey, Matrix{ComplexF64}}((0, 0, 0) => ComplexF64[0.0;;]),
            Dict{SuperhamCore.Model.RKey, Int}((0, 0, 0) => 1),
            nothing,
            edge_centers,
        )
        propagated = SuperhamCore.SupercellHam.propagate_centers(
            edge_model,
            reshape([0, 0, 0], 3, 1),
            Matrix{Float64}(I, 3, 3),
        )
        @test propagated.centers_frac[:, 1] ≈ [0.0, 0.0, 0.0]
    end
end

@testset "superham k-space eigensystem rejects non-Hermitian matrices" begin
    lattice = Matrix{Float64}(I, 3, 3)
    model = SuperhamCore.Model.HrModel(
        "bad onsite",
        lattice,
        2π .* inv(lattice)',
        2,
        Dict{SuperhamCore.Model.RKey, Matrix{ComplexF64}}(
            (0, 0, 0) => ComplexF64[0.0 1.0; 0.0 0.0],
        ),
        Dict{SuperhamCore.Model.RKey, Int}((0, 0, 0) => 1),
        nothing,
        nothing,
    )
    @test occursin(
        "H(k) is not Hermitian",
        error_message(() -> SuperhamCore.Kspace.eigenvalues_k(model, [0.0, 0.0, 0.0])),
    )
end

@testset "superham build report" begin
    lattice = Matrix{Float64}(I, 3, 3)
    reciprocal = 2π .* inv(lattice)'
    hoppings = Dict{SuperhamCore.Model.RKey, Matrix{ComplexF64}}(
        (0, 0, 0) => ComplexF64[0.0 + 0.0im;;],
    )
    ndegen = Dict{SuperhamCore.Model.RKey, Int}((0, 0, 0) => 1)
    centers = SuperhamCore.Model.CenterTable(
        reshape([0.0, 0.0, 0.0], 3, 1),
        reshape([0.0, 0.0, 0.0], 3, 1),
        ["X"],
        "inline",
        :manual_centers,
    )
    model = SuperhamCore.Model.HrModel(
        "onsite",
        lattice,
        reciprocal,
        1,
        hoppings,
        ndegen,
        nothing,
        centers,
    )
    geom = SuperhamCore.SupercellGeometry.from_user_matrix([2 0 0; 0 1 0; 0 0 1])
    result = SuperhamCore.SupercellHam.build_supercell(model, geom)

    @test :blocks ∉ propertynames(model)
    @test :blocks in propertynames(model, true)
    @test SuperhamCore.Model.hr_normalization(model) == :raw
    @test result.report.wsvec_input == false
    @test result.report.wsvec_output_policy == :none
    @test result.report.center_output_policy == :propagated
    @test result.model.num_wann == 2
    @test isnothing(result.model.wsvec)
    @test !isnothing(result.model.centers)
end

@testset "superham spin layout input" begin
    mktempdir() do dir
        input_path = joinpath(dir, "input.superham.toml")
        write(input_path, """
[superham.files]
hr = "model_hr.dat"
structure = "POSCAR"

[superham.spin]
layout = "vasp544"

[superham.supercell]
matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
""")
        cfg = SuperhamCore.InputIO.read_input(input_path)
        @test cfg.spin_layout == :vasp544

        bad_input = replace(read(input_path, String), "layout = \"vasp544\"" => "layout = \"interleaved\"")
        bad_path = joinpath(dir, "bad.superham.toml")
        write(bad_path, bad_input)
        @test occursin("Unsupported superham.spin.layout", error_message(() -> SuperhamCore.InputIO.read_input(bad_path)))

        missing_wsvec = replace(read(input_path, String), "structure = \"POSCAR\"" => "structure = \"POSCAR\"\nwsvec = \"missing_wsvec.dat\"")
        missing_wsvec_path = joinpath(dir, "missing-wsvec.superham.toml")
        write(missing_wsvec_path, missing_wsvec)
        @test occursin("superham.files.wsvec points to missing file", error_message(() -> SuperhamCore.InputIO.read_input(missing_wsvec_path)))
    end
end
