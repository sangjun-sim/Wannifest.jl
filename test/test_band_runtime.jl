@testset "band DOS mesh and metadata" begin
    shifted = Dos.generate_kmesh(2, 2, 1; shift=(0.5, 0.5, 0.0))
    @test shifted[1] == [0.25, 0.25, 0.0]

    cfg = InputIO.read_input(joinpath(GRAPHENE_EXAMPLE_DIR, "input.band.toml"))
    result = Dos.compute_dos([0.0 1.0; 0.5 1.5], cfg)
    @test length(result.energies) == cfg.dos.npts
    @test result.num_bands == 2
    @test result.integral >= 0
    @test result.center_of_mass ≈ Dos.density_center_of_mass(result.energies, result.dos)
    @test isfinite(result.center_of_mass)
    @test !result.window_is_auto

    function polynomial_dos_center(npts::Int)
        energies = collect(range(0.0, 1.0; length=npts))
        density = 1.0 .+ energies
        return Dos.density_center_of_mass(energies, density)
    end
    exact_polynomial_center = 5.0 / 9.0
    coarse_center = polynomial_dos_center(11)
    fine_center = polynomial_dos_center(101)
    @test abs(fine_center - exact_polynomial_center) < abs(coarse_center - exact_polynomial_center) / 50
    @test fine_center ≈ exact_polynomial_center atol=2e-5

    PlotService._plotbands_module_ref[] = nothing
    PlotService.maybe_plot(cfg, nothing, nothing, false)
    @test isnothing(PlotService._plotbands_module_ref[])
    plotbands_module = PlotService._plotbands_module()
    @test PlotService._plotbands_module() === plotbands_module
    widths = Base.invokelatest(getfield(plotbands_module, :_combined_widths), cfg)
    @test widths[2] / widths[1] ≈ cfg.combined_plot.dos_width_ratio
    hidden_ticks = Base.invokelatest(getfield(plotbands_module, :_hidden_ticks))
    @test hidden_ticks == (Float64[], String[])
    center_marker_count(plt) = count(
        series -> get(series.plotattributes, :linestyle, nothing) == :dot,
        plt.series_list,
    )
    @test Base.invokelatest(getfield(plotbands_module, :_center_visible), result.center_of_mass, cfg)
    @test !Base.invokelatest(getfield(plotbands_module, :_center_visible), Inf, cfg)
    dos_plot = Base.invokelatest(getfield(plotbands_module, :plot_dos), result, cfg; save=false)
    @test center_marker_count(dos_plot) == 1

    band_result = Model.BandResult(
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
        [0.0, 1.0],
        [0.0 1.0; 0.5 1.5],
        [0.0, 1.0],
        ["G", "X"],
        UnitRange{Int}[1:2],
        true,
    )
    mktempdir() do dir
        combined_cfg = Model.RunConfig(
            cfg.mode,
            cfg.files,
            Model.OutputFiles(
                joinpath(dir, "bands.dat"),
                joinpath(dir, "dos.dat"),
                joinpath(dir, "bands.png"),
                joinpath(dir, "dos.png"),
                joinpath(dir, "combined.png"),
            ),
            cfg.dos,
            cfg.energy,
            cfg.plot,
            cfg.band_plot,
            cfg.dos_plot,
            cfg.combined_plot,
            cfg.spin,
            cfg.projection,
            cfg.oam,
            cfg.sam,
            cfg.hermiticity_tol,
            cfg.verbose,
        )
        combined_plot = Base.invokelatest(
            getfield(plotbands_module, :plot_combined),
            band_result,
            result,
            combined_cfg,
        )
        @test center_marker_count(combined_plot) == 1
    end

    spin_cfg = BandCore.Model.RunConfig(
        cfg.mode,
        cfg.files,
        cfg.output,
        cfg.dos,
        cfg.energy,
        cfg.plot,
        cfg.band_plot,
        cfg.dos_plot,
        cfg.combined_plot,
        BandCore.Model.SpinConfig(true, :qe, ("up", "down")),
        cfg.hermiticity_tol,
        cfg.verbose,
    )
    @test Base.invokelatest(getfield(plotbands_module, :_band_color), spin_cfg, 1, 48) == "up"
    @test Base.invokelatest(getfield(plotbands_module, :_band_color), spin_cfg, 2, 48) == "down"
    @test Base.invokelatest(getfield(plotbands_module, :_band_color), spin_cfg, 29, 48) == "up"
    @test Base.invokelatest(getfield(plotbands_module, :_band_plot_order), spin_cfg, 6) == collect(1:6)
    @test Base.invokelatest(getfield(plotbands_module, :_band_plot_order), cfg, 6) == collect(1:6)

    spin_dos = Dos.compute_dos([-1.0 1.0 -2.0 2.0], spin_cfg)
    @test !isnothing(spin_dos.dos_down)
    @test spin_dos.energies[argmax(spin_dos.dos)] < 0
    @test spin_dos.energies[argmax(spin_dos.dos_down)] > 0
    @test spin_dos.center_of_mass ≈
        Dos.density_center_of_mass(spin_dos.energies, spin_dos.dos .+ spin_dos.dos_down)

    @test SpinLayout.band_indices(6, :qe) == ([1, 3, 5], [2, 4, 6])
    @test SpinLayout.band_indices(6, :vasp544) == ([1, 3, 5], [2, 4, 6])
    vasp_map6 = SpinLayout.source_to_canonical_indices(6, :vasp544)
    @test vasp_map6 == [1, 3, 5, 2, 4, 6]
    @test vasp_map6[vasp_map6] == [1, 5, 4, 3, 2, 6]
    @test vasp_map6[vasp_map6] != collect(1:6)

    vasp_spin_cfg = BandCore.Model.RunConfig(
        cfg.mode,
        cfg.files,
        cfg.output,
        cfg.dos,
        cfg.energy,
        cfg.plot,
        cfg.band_plot,
        cfg.dos_plot,
        cfg.combined_plot,
        BandCore.Model.SpinConfig(true, :vasp544, ("up", "down")),
        cfg.hermiticity_tol,
        cfg.verbose,
    )
    @test Base.invokelatest(getfield(plotbands_module, :_band_color), vasp_spin_cfg, 1, 6) == "up"
    @test Base.invokelatest(getfield(plotbands_module, :_band_color), vasp_spin_cfg, 2, 6) == "down"
    @test Base.invokelatest(getfield(plotbands_module, :_band_color), vasp_spin_cfg, 3, 6) == "up"
    @test Base.invokelatest(getfield(plotbands_module, :_band_color), vasp_spin_cfg, 4, 6) == "down"
    @test Base.invokelatest(getfield(plotbands_module, :_band_plot_order), vasp_spin_cfg, 6) == collect(1:6)
    vasp_spin_dos = Dos.compute_dos([-1.0 1.0 -2.0 2.0], vasp_spin_cfg)
    @test !isnothing(vasp_spin_dos.dos_down)
    @test vasp_spin_dos.energies[argmax(vasp_spin_dos.dos)] < 0
    @test vasp_spin_dos.energies[argmax(vasp_spin_dos.dos_down)] > 0
end

@testset "band execution policy" begin
    old_threads = BLAS.get_num_threads()
    target_threads = old_threads == 1 ? 2 : 1
    try
        observed_threads = Execution.with_blas_threads(target_threads) do
            BLAS.get_num_threads()
        end
        @test observed_threads == target_threads
        @test BLAS.get_num_threads() == old_threads
        @test_throws ErrorException Execution.with_blas_threads(target_threads) do
            error("forced BLAS restoration check")
        end
        @test BLAS.get_num_threads() == old_threads
    finally
        BLAS.set_num_threads(old_threads)
    end
end

@testset "band output writers" begin
    cfg = InputIO.read_input(joinpath(GRAPHENE_EXAMPLE_DIR, "input.band.toml"))
    projection = Model.BandProjectionResult(
        ["g1", "g2"],
        ["red", "blue"],
        reshape([0.1, 0.9, 0.2, 0.8], 1, 2, 2),
        true,
        true,
    )
    band_result = Model.BandResult(
        [[0.0, 0.0, 0.0]],
        [0.0],
        [0.0 1.0],
        [0.0],
        ["G"],
        UnitRange{Int}[1:1],
        true,
        projection,
    )
    dos_result = Model.DosResult(
        [0.0, 1.0],
        [1.0, 2.0],
        [0.5, 0.25],
        1.5,
        2,
        false,
        Model.ProjectedDosResult(["g1"], ["red"], reshape([0.4, 0.6], 2, 1), true, false),
    )

    mktempdir() do dir
        bands_path = joinpath(dir, "data", "bands.dat")
        Output.write_bands_data(bands_path, band_result; config=cfg)
        bands_text = read(bands_path, String)
        @test occursin("# columns = segment  kx  ky  kz  k_distance  band_1  band_2", bands_text)
        @test occursin("# energy_shift = 0.0", bands_text)

        weights_path = joinpath(dir, "data", "weights.dat")
        Output.write_projection_weights_data(weights_path, band_result)
        weights_text = read(weights_path, String)
        @test occursin("# labels = g1, g2", weights_text)
        @test occursin("# columns: ik ib label weight", weights_text)

        dos_path = joinpath(dir, "data", "dos.dat")
        Output.write_dos_data(dos_path, dos_result)
        @test occursin("# energy  dos_up  dos_down", read(dos_path, String))

        pdos_path = joinpath(dir, "data", "pdos.dat")
        Output.write_pdos_data(pdos_path, dos_result)
        pdos_text = read(pdos_path, String)
        @test occursin("# energy  g1", pdos_text)
        @test occursin("# atom_counts = 1", pdos_text)
        @test occursin("# covers_all = false", pdos_text)
    end
end
