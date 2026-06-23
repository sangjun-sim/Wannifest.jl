@testset "atomic OAM operators" begin
    Sel = ObservablesModel.OamOrbitalSelection
    Entry = LocalAxisRotation.LocalBasisEntry

    p_entries = Entry[
        Entry(1, "A1", "pz", "unpolarized"),
        Entry(2, "A1", "px", "unpolarized"),
        Entry(3, "A1", "py", "unpolarized"),
    ]
    p_ops = AtomicOam.build_l_operators(3, p_entries, [Sel("A1", "p")])
    @test p_ops.Lx ≈ p_ops.Lx'
    @test p_ops.Ly ≈ p_ops.Ly'
    @test p_ops.Lz ≈ p_ops.Lz'
    @test p_ops.Lx * p_ops.Ly - p_ops.Ly * p_ops.Lx ≈ im * p_ops.Lz
    @test p_ops.L2 ≈ 2I
    p_plus = ComplexF64[0, 1, im] ./ sqrt(2)
    @test p_ops.Lz * p_plus ≈ p_plus
    @test AtomicOam.oam_expectations(p_ops, reshape(p_plus, 3, 1))[1, 3] ≈ 1.0

    d_entries = Entry[
        Entry(1, "D1", "dz2", "unpolarized"),
        Entry(2, "D1", "dxz", "unpolarized"),
        Entry(3, "D1", "dyz", "unpolarized"),
        Entry(4, "D1", "dx2-y2", "unpolarized"),
        Entry(5, "D1", "dxy", "unpolarized"),
    ]
    d_ops = AtomicOam.build_l_operators(5, d_entries, [Sel("D1", "d")])
    @test d_ops.Lx * d_ops.Ly - d_ops.Ly * d_ops.Lx ≈ im * d_ops.Lz
    @test d_ops.L2 ≈ 6I
    @test sort(round.(eigvals(Hermitian(d_ops.Lz)); digits=10)) ≈ [-2, -1, 0, 1, 2]

    t2g_entries = Entry[
        Entry(1, "Ru1", "dxy", "unpolarized"),
        Entry(2, "Ru1", "dxz", "unpolarized"),
        Entry(3, "Ru1", "dyz", "unpolarized"),
    ]
    t2g_ops = AtomicOam.build_l_operators(3, t2g_entries, [Sel("Ru1", "t2g")])
    @test t2g_ops.Lx * t2g_ops.Ly - t2g_ops.Ly * t2g_ops.Lx ≈ -im * t2g_ops.Lz
    @test t2g_ops.L2 ≈ 2I

    spinful_entries = Entry[
        Entry(1, "A1", "pz", "up"),
        Entry(2, "A1", "pz", "dn"),
        Entry(3, "A1", "px", "up"),
        Entry(4, "A1", "px", "dn"),
        Entry(5, "A1", "py", "up"),
        Entry(6, "A1", "py", "dn"),
    ]
    spinful_ops = AtomicOam.build_l_operators(6, spinful_entries, [Sel("A1", "p")])
    @test spinful_ops.L2[1:2:5, 1:2:5] ≈ 2I
    @test spinful_ops.L2[2:2:6, 2:2:6] ≈ 2I

    @test occursin(
        "lacks OAM orbital",
        error_message(() -> AtomicOam.build_l_operators(2, p_entries[1:2], [Sel("A1", "p")])),
    )
    @test occursin(
        "unknown site label",
        error_message(() -> AtomicOam.build_l_operators(3, p_entries, [Sel("Missing", "p")])),
    )
    @test_logs (:warn, r"eg-only OAM block is zero") AtomicOam.build_l_operators(
        2,
        Entry[
            Entry(1, "E1", "dz2", "unpolarized"),
            Entry(2, "E1", "dx2-y2", "unpolarized"),
        ],
        [Sel("E1", "eg")],
    )
end

@testset "band OAM output" begin
    Sel = ObservablesModel.OamOrbitalSelection
    Entry = LocalAxisRotation.LocalBasisEntry

    function write_hr(path::AbstractString, H::AbstractMatrix{<:Complex})
        nw = size(H, 1)
        open(path, "w") do io
            println(io, "toy p-shell Lz")
            println(io, nw)
            println(io, 1)
            println(io, "  1")
            for i in 1:nw, j in 1:nw
                value = H[i, j]
                println(io, "  0  0  0  $i  $j  $(real(value))  $(imag(value))")
            end
        end
    end

    mktempdir() do dir
        entries = Entry[
            Entry(1, "A1", "pz", "unpolarized"),
            Entry(2, "A1", "px", "unpolarized"),
            Entry(3, "A1", "py", "unpolarized"),
        ]
        ops = AtomicOam.build_l_operators(3, entries, [Sel("A1", "p")])
        write_hr(joinpath(dir, "wannier90_hr.dat"), ops.Lz)
        write(joinpath(dir, "KPOINTS"), """
toy path
2
Line-mode
Reciprocal
0.0 0.0 0.0 ! G
0.5 0.0 0.0 ! X
""")
        write(joinpath(dir, "basis.win"), """
num_wann = 3

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
A 0.0 0.0 0.0
end atoms_frac

begin projections
A : p
end projections
""")
        input_path = joinpath(dir, "oam.toml")
        write(input_path, """
[band.run]
mode = "bands"
structure = ""
verbose = false

[band.energy]
shift = 0.0

[band.projection]
enabled = true
mode = "win_groups"
win = "basis.win"
groups = [
  ["A_p", ["A"], ["p"], "#1f77b4"],
]

[band.oam]
enabled = true
orbitals = [
  ["A1", "p"],
	]
	degeneracy_tol = 1.0e-8
	plot_components = ["Lz", "L_norm"]
	""")
	        cfg = InputIO.read_input(input_path)
	        @test cfg.oam.enabled
	        @test cfg.oam.selections == [Sel("A1", "p")]
	        @test cfg.oam.plot_components == [:Lz, :L_norm]
	        result = Service.run(cfg; make_plot=false)
	        @test !isnothing(result.band_result.oam)
	        @test result.band_result.oam[1, :, 3] ≈ [-1.0, 0.0, 1.0]
        @test result.band_result.oam[1, :, 1] ≈ [0.0, 0.0, 0.0] atol=1e-12
        @test result.band_result.oam[1, :, 2] ≈ [0.0, 0.0, 0.0] atol=1e-12
        @test result.band_result.oam[1, :, 5] ≈ [2.0, 2.0, 2.0]
        oam_path = ObservablesOutput.oam_data_path(cfg.output.bands_data)
	        oam_text = read(oam_path, String)
	        @test occursin("# observable = projected_atomic_oam", oam_text)
	        @test occursin("# orbitals = A1:p", oam_text)
	        @test occursin("# L_norm = sqrt(<Lx>^2 + <Ly>^2 + <Lz>^2), not sqrt(<L2>).", oam_text)
	        @test occursin("1.0000000000  1.0000000000  2.0000000000", oam_text)

	        PlotService._observables_plot_module_ref[] = nothing
	        PlotService.maybe_plot(cfg, result.band_result, nothing, false)
	        @test isnothing(PlotService._observables_plot_module_ref[])
	        PlotService.maybe_plot(cfg, result.band_result, nothing, true)
	        @test isnothing(PlotService._observables_plot_module_ref[])
	        @test !isfile(joinpath(dir, "outputs", "plots", "bands_oam_lz.png"))
	        @test !isfile(joinpath(dir, "outputs", "plots", "bands_oam_l_norm.png"))
	        observables_plot_module = PlotService.observables_plot_module()
	        oam_plot_path_fn = Base.invokelatest(getfield, observables_plot_module, :oam_plot_path)
	        plot_oam_bands_fn = Base.invokelatest(getfield, observables_plot_module, :plot_oam_bands)
	        @test Base.invokelatest(
	            oam_plot_path_fn,
	            cfg.output.bands_plot,
	            :Lz,
	        ) == joinpath(dir, "outputs", "plots", "bands_oam_lz.png")
	        oam_plot = Base.invokelatest(
	            plot_oam_bands_fn,
	            result.band_result,
	            cfg;
	            component=:Lz,
	            save=false,
	        )
	        @test !isnothing(oam_plot)

	        oam_plot_input = joinpath(dir, "oam_plot.toml")
	        write(
	            oam_plot_input,
	            replace(read(input_path, String), "[band.projection]" => "[band.plot]\ntargets = [\"oam\"]\n\n[band.projection]"),
	        )
	        oam_plot_cfg = InputIO.read_input(oam_plot_input)
	        PlotService.maybe_plot(oam_plot_cfg, result.band_result, nothing, true)
	        @test isfile(joinpath(dir, "outputs", "plots", "bands_oam_lz.png"))
	        @test isfile(joinpath(dir, "outputs", "plots", "bands_oam_l_norm.png"))

	        bad_key_input = joinpath(dir, "bad_key.toml")
	        write(bad_key_input, replace(read(input_path, String), "degeneracy_tol = 1.0e-8" => "basis = \"basis.toml\""))
	        @test occursin("Unsupported oam option", error_message(() -> InputIO.read_input(bad_key_input)))

	        bad_component_input = joinpath(dir, "bad_component.toml")
	        write(bad_component_input, replace(read(input_path, String), "[\"Lz\", \"L_norm\"]" => "[\"Sz\"]"))
	        @test occursin("oam.plot_components entries must be one of", error_message(() -> InputIO.read_input(bad_component_input)))

	        bad_shell_input = joinpath(dir, "bad_shell.toml")
	        write(bad_shell_input, replace(read(input_path, String), "[\"A1\", \"p\"]" => "[\"A1\", \"px\"]"))
	        @test occursin("orbital_shell must be one of", error_message(() -> InputIO.read_input(bad_shell_input)))

        dos_oam_input = joinpath(dir, "dos_oam.toml")
        write(dos_oam_input, replace(read(input_path, String), "mode = \"bands\"" => "mode = \"dos\""))
        @test occursin("band.oam requires run.mode", error_message(() -> InputIO.read_input(dos_oam_input)))

        index_input = joinpath(dir, "index_oam.toml")
        write(index_input, replace(read(input_path, String), "mode = \"win_groups\"\nwin = \"basis.win\"\ngroups = [\n  [\"A_p\", [\"A\"], [\"p\"], \"#1f77b4\"],\n]" => "mode = \"index_groups\"\ngroups = [[\"A\", [1], \"blue\"]]"))
        @test occursin("band.oam requires band.projection mode", error_message(() -> InputIO.read_input(index_input)))
    end
end
