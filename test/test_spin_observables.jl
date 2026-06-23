@testset "atomic spin operators" begin
    Entry = LocalAxisRotation.LocalBasisEntry

    entries = Entry[
        Entry(1, "A1", "s", "up"),
        Entry(2, "A1", "s", "dn"),
    ]
    ops = AtomicSpin.build_spin_operators(2, entries)
    @test ops.Sx ≈ ops.Sx'
    @test ops.Sy ≈ ops.Sy'
    @test ops.Sz ≈ ops.Sz'
    @test ops.Sx * ops.Sy - ops.Sy * ops.Sx ≈ im * ops.Sz
    @test ops.S2 ≈ 0.75I

    up = ComplexF64[1, 0]
    dn = ComplexF64[0, 1]
    sx_plus = ComplexF64[1, 1] ./ sqrt(2)
    @test AtomicSpin.spin_expectations(ops, reshape(up, 2, 1))[1, 3] ≈ 0.5
    @test AtomicSpin.spin_expectations(ops, reshape(dn, 2, 1))[1, 3] ≈ -0.5
    @test AtomicSpin.spin_expectations(ops, reshape(sx_plus, 2, 1))[1, 1] ≈ 0.5

    @test occursin(
        "spinful up/dn",
        error_message(() -> AtomicSpin.build_spin_operators(1, Entry[Entry(1, "A1", "s", "unpolarized")])),
    )
    @test occursin(
        "spinful up/dn",
        error_message(() -> AtomicSpin.build_spin_operators(3, Entry[
            Entry(1, "A1", "s", "up"),
            Entry(2, "A1", "s", "dn"),
            Entry(3, "B1", "s", "unpolarized"),
        ])),
    )
    @test occursin(
        "complete up/dn pair",
        error_message(() -> AtomicSpin.build_spin_operators(1, Entry[Entry(1, "A1", "s", "up")])),
    )
    @test occursin(
        "duplicate SAM basis entry",
        error_message(() -> AtomicSpin.build_spin_operators(2, Entry[
            Entry(1, "A1", "s", "up"),
            Entry(2, "A1", "s", "up"),
        ])),
    )
end

@testset "SAM basis ordering" begin
    Entry = LocalAxisRotation.LocalBasisEntry

    mktempdir() do dir
        win_path = joinpath(dir, "spinors.win")
        write(win_path, """
num_wann = 4
spinors = .true.

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
X 0.0 0.0 0.0
Y 0.5 0.0 0.0
end atoms_frac

begin projections
X:s
Y:s
end projections
""")
        qe_basis = Win90Basis.read_win_basis(win_path; spin_layout=:qe)
        vasp_basis = Win90Basis.read_win_basis(win_path; spin_layout=:vasp544)
        @test [orb.spin for orb in qe_basis.orbitals] == [:up, :dn, :up, :dn]
        @test [orb.spin for orb in vasp_basis.orbitals] == [:up, :up, :dn, :dn]
        vasp_entries = LocalAxisRotation.basis_entries_from_win(vasp_basis)
        vasp_canonical = Projection._canonicalize_basis_entries(vasp_entries, 4, :vasp544)
        ops = AtomicSpin.build_spin_operators(4, vasp_canonical)
        @test ops.Sx[1, 2] ≈ 0.5
        @test ops.Sx[3, 4] ≈ 0.5
        @test diag(ops.Sz) ≈ ComplexF64[0.5, -0.5, 0.5, -0.5]
        @test ops.S2 ≈ 0.75I
    end
end

@testset "band SAM output" begin
    function write_hr(path::AbstractString, H::AbstractMatrix{<:Complex})
        nw = size(H, 1)
        open(path, "w") do io
            println(io, "toy spin Sz")
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
        write_hr(joinpath(dir, "wannier90_hr.dat"), ComplexF64[0.5 0.0; 0.0 -0.5])
        write(joinpath(dir, "KPOINTS"), """
toy path
2
Line-mode
Reciprocal
0.0 0.0 0.0 ! G
0.5 0.0 0.0 ! X
""")
        write(joinpath(dir, "basis.win"), """
num_wann = 2
spinors = .true.

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
A 0.0 0.0 0.0
end atoms_frac

begin projections
A:s
end projections
""")
        input_path = joinpath(dir, "sam.toml")
        write(input_path, """
[band.run]
mode = "bands"
structure = ""
verbose = false

[band.energy]
shift = 0.0

[band.spin]
layout = "qe"

[band.projection]
enabled = true
mode = "win_groups"
win = "basis.win"
groups = [
  ["A_s", ["A"], ["s"], "blue"],
]

[band.sam]
enabled = true
degeneracy_tol = 1.0e-8
""")
        cfg = InputIO.read_input(input_path)
        @test cfg.sam.enabled
        @test cfg.sam.plot_components == [:Sz]
        @test !cfg.spin.enabled
        result = Service.run(cfg; make_plot=false)
        @test !isnothing(result.band_result.sam)
        @test result.band_result.sam[1, :, 3] ≈ [-0.5, 0.5]
        @test result.band_result.sam[1, :, 5] ≈ [0.75, 0.75]
        sam_path = ObservablesOutput.sam_data_path(cfg.output.bands_data)
        sam_text = read(sam_path, String)
        @test occursin("# observable = spin_angular_momentum", sam_text)
        @test occursin("# spin_layout = qe", sam_text)
        @test occursin("# S_norm = sqrt(<Sx>^2 + <Sy>^2 + <Sz>^2), not sqrt(<S2>).", sam_text)

        PlotService._observables_plot_module_ref[] = nothing
        PlotService.maybe_plot(cfg, result.band_result, nothing, false)
        @test isnothing(PlotService._observables_plot_module_ref[])
        PlotService.maybe_plot(cfg, result.band_result, nothing, true)
        @test isnothing(PlotService._observables_plot_module_ref[])
        @test !isfile(joinpath(dir, "outputs", "plots", "bands_sam_sz.png"))
        observables_plot_module = PlotService.observables_plot_module()
        sam_plot_path_fn = Base.invokelatest(getfield, observables_plot_module, :sam_plot_path)
        plot_sam_bands_fn = Base.invokelatest(getfield, observables_plot_module, :plot_sam_bands)
        @test Base.invokelatest(
            sam_plot_path_fn,
            cfg.output.bands_plot,
            :Sz,
        ) == joinpath(dir, "outputs", "plots", "bands_sam_sz.png")
        sam_plot = Base.invokelatest(
            plot_sam_bands_fn,
            result.band_result,
            cfg;
            component=:Sz,
            save=false,
        )
        @test !isnothing(sam_plot)

        sam_plot_input = joinpath(dir, "sam_plot.toml")
        write(
            sam_plot_input,
            replace(
                read(input_path, String),
                "[band.projection]" => "[band.plot]\ntargets = [\"sam\"]\n\n[band.projection]",
                "degeneracy_tol = 1.0e-8" => "degeneracy_tol = 1.0e-8\nplot_components = [\"Sz\", \"S_norm\"]",
            ),
        )
        sam_plot_cfg = InputIO.read_input(sam_plot_input)
        @test sam_plot_cfg.sam.plot_components == [:Sz, :S_norm]
        PlotService.maybe_plot(sam_plot_cfg, result.band_result, nothing, true)
        @test isfile(joinpath(dir, "outputs", "plots", "bands_sam_sz.png"))
        @test isfile(joinpath(dir, "outputs", "plots", "bands_sam_s_norm.png"))

        bad_key_input = joinpath(dir, "bad_key.toml")
        write(bad_key_input, replace(read(input_path, String), "degeneracy_tol = 1.0e-8" => "basis = \"basis.toml\""))
        @test occursin("Unsupported sam option", error_message(() -> InputIO.read_input(bad_key_input)))

        bad_component_input = joinpath(dir, "bad_component.toml")
        write(bad_component_input, replace(
            read(input_path, String),
            "degeneracy_tol = 1.0e-8" => "degeneracy_tol = 1.0e-8\nplot_components = [\"Lz\"]",
        ))
        @test occursin("sam.plot_components entries must be one of", error_message(() -> InputIO.read_input(bad_component_input)))

        missing_sam_target_input = joinpath(dir, "missing_sam_target.toml")
        write(
            missing_sam_target_input,
            replace(
                read(input_path, String),
                "[band.projection]" => "[band.plot]\ntargets = [\"sam\"]\n\n[band.projection]",
                "[band.sam]\nenabled = true\ndegeneracy_tol = 1.0e-8\n" => "",
            ),
        )
        @test occursin(
            "plot target \"sam\" requires band.sam.enabled=true",
            error_message(() -> InputIO.read_input(missing_sam_target_input)),
        )

        dos_sam_input = joinpath(dir, "dos_sam.toml")
        write(dos_sam_input, replace(read(input_path, String), "mode = \"bands\"" => "mode = \"dos\""))
        @test occursin("band.sam requires run.mode", error_message(() -> InputIO.read_input(dos_sam_input)))

        no_projection_input = joinpath(dir, "no_projection.toml")
        write(no_projection_input, replace(
            read(input_path, String),
            """[band.projection]
enabled = true
mode = "win_groups"
win = "basis.win"
groups = [
  ["A_s", ["A"], ["s"], "blue"],
]

""" => "",
        ))
        @test occursin("band.sam requires band.projection metadata", error_message(() -> InputIO.read_input(no_projection_input)))

        no_layout_input = joinpath(dir, "no_layout.toml")
        write(no_layout_input, replace(read(input_path, String), """
[band.spin]
layout = "qe"

""" => ""))
        @test occursin("band.sam requires explicit band.spin.layout", error_message(() -> InputIO.read_input(no_layout_input)))
    end
end
