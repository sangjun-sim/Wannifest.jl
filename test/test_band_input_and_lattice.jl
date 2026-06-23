@testset "band input defaults and lattice policy" begin
    example_input = joinpath(GRAPHENE_EXAMPLE_DIR, "input.band.toml")
    cfg = InputIO.read_input(example_input)
    @test cfg.mode == :all
    @test basename(cfg.files.hr_path) == "graphene_hr.dat"
    @test isnothing(cfg.files.wsvec_path)
    @test basename(cfg.files.kpoints_path) == "KPOINTS"
    @test basename(cfg.files.structure_path) == "POSCAR"
    @test basename(cfg.output.bands_data) == "bands.dat"
    @test !occursin("[band.solver]", read(example_input, String))

    mktempdir() do dir
        cp(joinpath(GRAPHENE_EXAMPLE_DIR, "graphene_hr.dat"), joinpath(dir, "wannier90_hr.dat"))
        cp(joinpath(GRAPHENE_EXAMPLE_DIR, "KPOINTS"), joinpath(dir, "KPOINTS"))
        cp(joinpath(GRAPHENE_EXAMPLE_DIR, "POSCAR"), joinpath(dir, "POSCAR"))
        input_path = joinpath(dir, "minimal.toml")
        write(input_path, """
[band.run]
mode = "bands"

[band.energy]
shift = 0.0
""")
        minimal = InputIO.read_input(input_path)
        @test basename(minimal.files.hr_path) == "wannier90_hr.dat"
        @test isnothing(minimal.files.wsvec_path)
        @test basename(minimal.output.bands_data) == "bands.dat"
        @test minimal.dos.mesh == (16, 16, 16)
        @test minimal.hermiticity_tol == 1e-8
        @test minimal.verbose
        @test minimal.combined_plot.dos_xlabel == "DOS (states/eV)"
        @test minimal.plot.targets == Symbol[]

        cp(joinpath(SAMPLE_SPECTRA_DIR, "wsvec_test.dat"), joinpath(dir, "wannier90_wsvec.dat"))
        auto_wsvec_cfg = InputIO.read_input(input_path)
        @test isnothing(auto_wsvec_cfg.files.wsvec_path)

        explicit_wsvec_input = joinpath(dir, "explicit_wsvec.toml")
        write(explicit_wsvec_input, """
[band.run]
mode = "bands"
wsvec = "wannier90_wsvec.dat"
verbose = false
hermiticity_tol = 1.0e-6

[band.energy]
shift = 0.0

[band.combined_plot]
dos_xlabel = "Projected DOS"
dos_ylabel = "Energy"
""")
        explicit_wsvec_cfg = InputIO.read_input(explicit_wsvec_input)
        @test basename(explicit_wsvec_cfg.files.wsvec_path) == "wannier90_wsvec.dat"
        @test !explicit_wsvec_cfg.verbose
        @test explicit_wsvec_cfg.hermiticity_tol == 1.0e-6
        @test explicit_wsvec_cfg.combined_plot.dos_xlabel == "Projected DOS"
        @test explicit_wsvec_cfg.combined_plot.dos_ylabel == "Energy"

        plot_targets_input = joinpath(dir, "plot_targets.toml")
        write(plot_targets_input, """
[band.run]
mode = "all"

[band.energy]
shift = 0.0

[band.plot]
targets = ["band", "dos", "combined"]
""")
        plot_targets_cfg = InputIO.read_input(plot_targets_input)
        @test plot_targets_cfg.plot.targets == [:band, :dos, :combined]

        bad_plot_target_input = joinpath(dir, "bad_plot_target.toml")
        write(bad_plot_target_input, replace(read(plot_targets_input, String), "\"combined\"" => "\"bad\""))
        @test occursin(
            "plot.targets entries must be one of",
            error_message(() -> InputIO.read_input(bad_plot_target_input)),
        )

        duplicate_plot_target_input = joinpath(dir, "duplicate_plot_target.toml")
        write(duplicate_plot_target_input, replace(read(plot_targets_input, String), "\"combined\"" => "\"band\""))
        @test occursin(
            "plot.targets contains duplicate target",
            error_message(() -> InputIO.read_input(duplicate_plot_target_input)),
        )

        invalid_plot_mode_input = joinpath(dir, "invalid_plot_mode.toml")
        write(invalid_plot_mode_input, replace(read(plot_targets_input, String), "mode = \"all\"" => "mode = \"bands\""))
        @test occursin(
            "plot target \"dos\" requires",
            error_message(() -> InputIO.read_input(invalid_plot_mode_input)),
        )

        projection_plot_without_projection_input = joinpath(dir, "projection_plot_without_projection.toml")
        write(
            projection_plot_without_projection_input,
            replace(read(plot_targets_input, String), "[\"band\", \"dos\", \"combined\"]" => "[\"fatband\"]"),
        )
        @test occursin(
            "plot target \"fatband\" requires band.projection.enabled=true",
            error_message(() -> InputIO.read_input(projection_plot_without_projection_input)),
        )

        explicit_spin_input = joinpath(dir, "explicit_spin.toml")
        write(explicit_spin_input, """
[band.run]
mode = "bands"

[band.energy]
shift = 0.0

[band.spin]
enabled = true
layout = "vasp544"
colors = ["up", "down"]
""")
        explicit_spin_cfg = InputIO.read_input(explicit_spin_input)
        @test explicit_spin_cfg.spin.enabled
        @test explicit_spin_cfg.spin.layout == :vasp544
        @test explicit_spin_cfg.spin.colors == ("up", "down")

        missing_spin_layout_input = joinpath(dir, "missing_spin_layout.toml")
        write(missing_spin_layout_input, """
[band.run]
mode = "bands"

[band.energy]
shift = 0.0

[band.spin]
enabled = true
""")
        @test occursin(
            "spin.layout is required",
            error_message(() -> InputIO.read_input(missing_spin_layout_input)),
        )

        legacy_spin_layout_input = joinpath(dir, "legacy_spin_layout.toml")
        write(legacy_spin_layout_input, """
[band.run]
mode = "bands"

[band.energy]
shift = 0.0

[band.spin]
enabled = true
layout = "interleaved"
""")
        @test occursin(
            "Unsupported spin.layout",
            error_message(() -> InputIO.read_input(legacy_spin_layout_input)),
        )

        unknown_band_table_input = joinpath(dir, "unknown_band_table.toml")
        write(unknown_band_table_input, """
[band.run]
mode = "bands"

[band.energy]
shift = 0.0

[band.solver]
threads = 1
""")
        @test occursin("Unsupported band table", error_message(() -> InputIO.read_input(unknown_band_table_input)))
    end

    mktempdir() do dir
        write(joinpath(dir, "wannier90_hr.dat"), "")
        unknown_run_input = joinpath(dir, "unknown_run.toml")
        write(unknown_run_input, """
[band.run]
mode = "bands"
legacy_option = true

[band.energy]
shift = 0.0
""")
        @test occursin("Unsupported run option", error_message(() -> InputIO.read_input(unknown_run_input)))

        bad_mesh_input = joinpath(dir, "bad_mesh.toml")
        write(bad_mesh_input, """
[band.run]
mode = "dos"

[band.energy]
shift = 0.0

[band.dos]
mesh = 16
""")
        @test occursin("dos.mesh must be a string", error_message(() -> InputIO.read_input(bad_mesh_input)))
    end
end

@testset "band structure lattice readers" begin
    toml_path = joinpath(GRAPHENE_EXAMPLE_DIR, "input.band.toml")
    toml_error = error_message(() -> LatticeIO.read_lattice(toml_path))
    @test occursin("TOML structure files are not supported", toml_error)

    poscar_lattice = LatticeIO.read_lattice(joinpath(GRAPHENE_EXAMPLE_DIR, "POSCAR"))
    @test size(poscar_lattice.real_lattice) == (3, 3)
    @test poscar_lattice.source == abspath(joinpath(GRAPHENE_EXAMPLE_DIR, "POSCAR"))

    mktempdir() do dir
        win_path = joinpath(dir, "wannier90.win")
        write(win_path, """
begin unit_cell_cart bohr
  1.8897261245650618  0.0                 0.0
  0.0                 1.8897261245650618  0.0
  0.0                 0.0                 1.8897261245650618
end unit_cell_cart
begin atoms_cart bohr
  C  0.0  0.0  0.0
end atoms_cart
""")
        win_lattice = LatticeIO.read_wannier_win(win_path)
        @test size(win_lattice.real_lattice) == (3, 3)
        @test occursin("unit_cell_cart", win_lattice.source)

        no_atoms_path = joinpath(dir, "no_atoms.win")
        write(no_atoms_path, """
begin unit_cell_cart
  1.0  0.0  0.0
  0.0  1.0  0.0
  0.0  0.0  1.0
end unit_cell_cart
""")
        no_atoms_lattice = LatticeIO.read_wannier_win(no_atoms_path)
        @test no_atoms_lattice.real_lattice ≈ Matrix{Float64}(I, 3, 3)
    end
end

@testset "band k-path conventions" begin
    lattice = LatticeIO.read_lattice(joinpath(SAMPLE_SPECTRA_DIR, "POSCAR_cubic"))
    cartesian = KPath.parse_kpoints(joinpath(SAMPLE_SPECTRA_DIR, "KPOINTS_cartesian"); lattice=lattice)
    cartesian_result = KPath.generate_kpath(cartesian; lattice=lattice)
    @test cartesian_result isa KPath.KPathResult
    @test cartesian_result.is_physical_distance

    mktempdir() do dir
        poscar = joinpath(dir, "POSCAR")
        write(poscar, """
Nonorthogonal cell
1.0
1.0 0.0 0.0
0.5 0.8660254037844386 0.0
0.0 0.0 1.0
Si
1
Direct
0.0 0.0 0.0
""")
        nonorth_lattice = LatticeIO.read_lattice(poscar)
        target = [0.25, 0.5, 0.0]
        kcart = nonorth_lattice.reciprocal_lattice * target
        kpoints_path = joinpath(dir, "KPOINTS")
        write(kpoints_path, """
Nonorthogonal Cartesian path
4
L
Cartesian
$(kcart[1]) $(kcart[2]) $(kcart[3]) ! A
0.0 0.0 0.0 ! G
""")
        parsed = KPath.parse_kpoints(kpoints_path; lattice=nonorth_lattice)
        @test parsed.segments[1].k_start ≈ target
    end

    reduced = KPath.parse_kpoints(joinpath(GRAPHENE_EXAMPLE_DIR, "KPOINTS"))
    reduced_result = KPath.generate_kpath(reduced)
    @test !reduced_result.is_physical_distance

    cfg = InputIO.read_input(joinpath(GRAPHENE_EXAMPLE_DIR, "input.band.toml"))
    hr = WannierHrIO.read_hr(joinpath(GRAPHENE_EXAMPLE_DIR, "graphene_hr.dat"))
    @test_logs (:info, r"No lattice provided") Bands.compute_bands(
        hr,
        reduced_result,
        cfg,
    )

    evals, _ = WannierEigensystem.solve_kpoint(hr, [1 / 3, 1 / 3, 0.0])
    @test maximum(abs.(evals)) < 1e-10

    merged = KPath.KPathData(
        [
            KPath.KSegment([0.0, 0.0, 0.0], [0.5, 0.0, 0.0], "", "X"),
            KPath.KSegment([0.5, 0.0, 0.0], [0.5, 0.5, 0.0], "", ""),
        ],
        2,
    )
    @test KPath.generate_kpath(merged).tick_labels == ["", "X (reduced)", ""]

    mktempdir() do dir
        zero_segment = joinpath(dir, "KPOINTS")
        write(zero_segment, """
Zero-length with blanks
3
Line-mode
Reciprocal

0.0 0.0 0.0 ! G

0.0 0.0 0.0 ! G
""")
        zero_error = error_message(() -> KPath.parse_kpoints(zero_segment))
        @test occursin("line 6", zero_error)
    end
end

@testset "hr read normalization policy" begin
    mktempdir() do dir
        path = joinpath(dir, "ndegen_two_hr.dat")
        write(path, """
Raw ndegen policy
  1
  1
  2
  0  0  0  1  1   4.000000   0.000000
""")

        raw = WannierHrIO.read_hr(path)

        @test raw.hoppings[(0, 0, 0)][1, 1] == 4.0 + 0.0im
        @test raw.normalization == :raw
        @test WannierKspace.hamiltonian_k_plain(raw, [0.0, 0.0, 0.0])[1, 1] ≈ 2.0 + 0.0im
        @test WannierHrIO.normalized_hoppings(raw)[(0, 0, 0)][1, 1] ≈ 2.0 + 0.0im
    end
end

@testset "hr read spin layout canonicalization" begin
    mktempdir() do dir
        path = joinpath(dir, "vasp_block_hr.dat")
        H = ComplexF64[
            11.0  13.0   0.0   0.0
            31.0  33.0   0.0   0.0
             0.0   0.0  22.0  24.0
             0.0   0.0  42.0  44.0
        ]
        open(path, "w") do io
            println(io, "VASP block spin layout")
            println(io, "  4")
            println(io, "  1")
            println(io, "  1")
            for m in 1:4, n in 1:4
                println(io, "  0  0  0  $m  $n   $(real(H[m, n]))   $(imag(H[m, n]))")
            end
        end

        qe = WannierHrIO.read_hr(path; spin_layout=:qe)
        vasp = WannierHrIO.read_hr(path; spin_layout=:vasp544)
        index_map = SpinLayout.source_to_canonical_indices(4, :vasp544)
        expected = zeros(ComplexF64, 4, 4)
        for m in 1:4, n in 1:4
            expected[index_map[m], index_map[n]] = H[m, n]
        end

        @test qe.hoppings[(0, 0, 0)] == H
        @test vasp.hoppings[(0, 0, 0)] == expected
        @test index_map == [1, 3, 2, 4]

        path6 = joinpath(dir, "vasp_block_hr_6.dat")
        H6 = Matrix{ComplexF64}(Diagonal(ComplexF64.(1:6)))
        open(path6, "w") do io
            println(io, "VASP block spin layout 6")
            println(io, "  6")
            println(io, "  1")
            println(io, "  1")
            for m in 1:6, n in 1:6
                println(io, "  0  0  0  $m  $n   $(real(H6[m, n]))   $(imag(H6[m, n]))")
            end
        end

        index_map6 = SpinLayout.source_to_canonical_indices(6, :vasp544)
        double_map6 = index_map6[index_map6]
        expected_once = zeros(ComplexF64, 6, 6)
        expected_twice = zeros(ComplexF64, 6, 6)
        for m in 1:6, n in 1:6
            expected_once[index_map6[m], index_map6[n]] = H6[m, n]
            expected_twice[double_map6[m], double_map6[n]] = H6[m, n]
        end

        vasp6 = WannierHrIO.read_hr(path6; spin_layout=:vasp544)
        @test diag(vasp6.hoppings[(0, 0, 0)]) == ComplexF64[1, 4, 2, 5, 3, 6]
        @test vasp6.hoppings[(0, 0, 0)] == expected_once
        @test vasp6.hoppings[(0, 0, 0)] != expected_twice
    end
end
