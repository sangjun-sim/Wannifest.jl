@testset "band service smoke" begin
    cfg = InputIO.read_input(joinpath(GRAPHENE_EXAMPLE_DIR, "input.band.toml"))
    result = Service.run(cfg; make_plot=false)
    @test result.hermiticity_ok
    @test !isnothing(result.band_result)
    @test !isempty(result.band_result.kpoints_frac)
    @test result.band_result.is_physical_distance
    bands_data = read(cfg.output.bands_data, String)
    @test occursin("# is_physical_distance = true", bands_data)
    @test occursin("# energy_shift = 0.0", bands_data)
    @test occursin("# hr = ", bands_data)

    mktempdir() do dir
        cp(joinpath(SAMPLE_SPECTRA_DIR, "chain_1d_hr.dat"), joinpath(dir, "wannier90_hr.dat"))
        write(joinpath(dir, "KPOINTS"), """
1D chain path
4
Line-mode
Reciprocal
0.0 0.0 0.0 ! G
0.5 0.0 0.0 ! X
""")
        write(joinpath(dir, "wannier90_wsvec.dat"), """
incomplete wsvec
 0  0  0  1  1
 1
 0  0  0
""")
        bands_input = joinpath(dir, "bands.toml")
        write(bands_input, """
[band.run]
mode = "bands"
structure = ""
verbose = false

[band.energy]
shift = 0.0
""")
        bands_cfg = InputIO.read_input(bands_input)
        @test isnothing(bands_cfg.files.wsvec_path)
        bands_result = Service.run(bands_cfg; make_plot=false)
        @test size(bands_result.band_result.eigenvalues) == (4, 1)
        PlotService._plotbands_module_ref[] = nothing
        Service.run(bands_cfg; make_plot=true)
        @test !isfile(bands_cfg.output.bands_plot)
        @test isnothing(PlotService._plotbands_module_ref[])

        band_plot_input = joinpath(dir, "bands_plot.toml")
        write(band_plot_input, """
[band.run]
mode = "bands"
structure = ""
verbose = false

[band.output]
bands_plot = "outputs/plots/explicit_bands.png"

[band.energy]
shift = 0.0

[band.plot]
targets = ["band"]
""")
        band_plot_cfg = InputIO.read_input(band_plot_input)
        Service.run(band_plot_cfg; make_plot=true)
        @test isfile(band_plot_cfg.output.bands_plot)

        no_plot_input = joinpath(dir, "no_plot.toml")
        write(no_plot_input, replace(read(band_plot_input, String), "explicit_bands.png" => "no_plot.png"))
        no_plot_cfg = InputIO.read_input(no_plot_input)
        Service.run(no_plot_cfg; make_plot=false)
        @test !isfile(no_plot_cfg.output.bands_plot)

        explicit_bad_wsvec = joinpath(dir, "explicit_bad_wsvec.toml")
        write(explicit_bad_wsvec, """
[band.run]
mode = "bands"
wsvec = "wannier90_wsvec.dat"
structure = ""

[band.energy]
shift = 0.0
""")
        explicit_bad_cfg = InputIO.read_input(explicit_bad_wsvec)
        @test occursin("wsvec is missing", error_message(() -> Service.run(explicit_bad_cfg; make_plot=false)))

        write(joinpath(dir, "wannier90_hr.dat"), """
wsvec overrides hr ndegen test
  1
  1
  2
  0  0  0  1  1   2.000000   0.000000
""")
        write(joinpath(dir, "complete_wsvec.dat"), """
complete wsvec
 0  0  0  1  1
 1
 0  0  0
""")
        write(joinpath(dir, "wsvec_overrides_ndegen.toml"), """
[band.run]
mode = "bands"
wsvec = "complete_wsvec.dat"
structure = ""
verbose = false

[band.energy]
shift = 0.0
""")
        wsvec_overrides_cfg = InputIO.read_input(joinpath(dir, "wsvec_overrides_ndegen.toml"))
        wsvec_overrides_result = Service.run(wsvec_overrides_cfg; make_plot=false)
        @test wsvec_overrides_result.band_result.eigenvalues[1, 1] ≈ 2.0

        write(joinpath(dir, "basis.win"), """
num_wann = 1

begin unit_cell_cart
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
Ru 0.0 0.0 0.0
end atoms_frac

begin projections
Ru : s
end projections
""")
        win_input = joinpath(dir, "win_projection.toml")
        write(win_input, """
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
weights_data = "win_weights.dat"

[[band.projection.groups]]
label = "Ru_s"
species = ["Ru"]
orbitals = ["s"]
""")
        win_result = Service.run(InputIO.read_input(win_input); make_plot=false)
        @test win_result.band_result.projection.labels == ["Ru_s"]
        @test all(win_result.band_result.projection.weights[:, 1, 1] .≈ 1.0)
    end

    mktempdir() do dir
        cp(joinpath(GRAPHENE_EXAMPLE_DIR, "graphene_hr.dat"), joinpath(dir, "graphene_hr.dat"))
        cp(joinpath(GRAPHENE_EXAMPLE_DIR, "KPOINTS"), joinpath(dir, "KPOINTS"))
        cp(joinpath(GRAPHENE_EXAMPLE_DIR, "POSCAR"), joinpath(dir, "POSCAR"))
        projected_input = joinpath(dir, "projected.toml")
        write(projected_input, """
[band.run]
mode = "all"
hr = "graphene_hr.dat"
kpoints = "KPOINTS"
structure = "POSCAR"
verbose = false

[band.dos]
mesh = "1x1x1"
sigma = 0.2
npts = 101

[band.energy]
shift = 0.0

[band.projection]
enabled = true
mode = "index_groups"
color_group = ["orb1"]
plot_style = "empty_circle"
colorbar_colormap = "plasma"
circle_max_size = 12.5
circle_stroke_width = 1.4
weights_data = "projection_weights.dat"
pdos_data = "pdos.dat"

[[band.projection.groups]]
label = "orb1"
indices = [1]
color = "#1f77b4"

[[band.projection.groups]]
label = "orb2"
indices = [2]
color = "#ff7f0e"
""")
        projected_cfg = InputIO.read_input(projected_input)
        @test projected_cfg.projection.enabled
        @test projected_cfg.projection.color_group == ["orb1"]
        @test projected_cfg.projection.plot_style == :empty_circle
        @test projected_cfg.projection.colorbar_colormap == "plasma"
        @test projected_cfg.projection.circle_max_size == 12.5
        @test projected_cfg.projection.circle_stroke_width == 1.4
        projected_hr = WannierHrIO.read_hr(projected_cfg.files.hr_path)
        projected_spec = Projection.build_projection_spec(projected_cfg, projected_hr.num_wann)
        direct_dos = Dos.run_dos(projected_hr, projected_cfg; projection_spec=projected_spec)
        @test !isnothing(direct_dos.projected)
        direct_pdos_sum = vec(sum(direct_dos.projected.pdos; dims=2))
        @test maximum(abs.(direct_pdos_sum .- direct_dos.dos)) < 1e-10
        @test length(direct_dos.projected.centers_of_mass) == 2
        @test all(isfinite, direct_dos.projected.centers_of_mass)
        for ig in eachindex(direct_dos.projected.labels)
            @test direct_dos.projected.centers_of_mass[ig] ≈
                Dos.density_center_of_mass(direct_dos.energies, view(direct_dos.projected.pdos, :, ig))
        end
        plotbands_module = PlotService._plotbands_module()
        marker_sizes = Base.invokelatest(
            getfield(plotbands_module, :_projection_marker_sizes),
            [0.0, 0.5, 1.0],
            projected_cfg,
        )
        @test marker_sizes ≈ [0.0, 6.25, 12.5]
        projected_result = Service.run(projected_cfg; make_plot=false)
        @test !isnothing(projected_result.band_result.projection)
        PlotService._plotbands_module_ref[] = nothing
        Service.run(projected_cfg; make_plot=true)
        @test !isfile(projected_cfg.output.bands_plot)
        @test !isfile(projected_cfg.output.dos_plot)
        @test !isfile(projected_cfg.output.combined_plot)
        @test !isfile(projected_cfg.projection.projected_bands_plot)
        @test !isfile(projected_cfg.projection.pdos_plot)
        @test !isfile(projected_cfg.projection.projected_combined_plot)
        @test isnothing(PlotService._plotbands_module_ref[])

        projected_plot_input = joinpath(dir, "projected_plot.toml")
        write(
            projected_plot_input,
            replace(read(projected_input, String), "[band.projection]" => "[band.plot]\ntargets = [\"fatband\", \"pdos\"]\n\n[band.projection]"),
        )
        projected_plot_cfg = InputIO.read_input(projected_plot_input)
        Service.run(projected_plot_cfg; make_plot=true)
        @test isfile(projected_plot_cfg.projection.projected_bands_plot)
        @test isfile(projected_plot_cfg.projection.pdos_plot)
        @test !isfile(projected_plot_cfg.output.bands_plot)
        @test !isfile(projected_plot_cfg.output.dos_plot)
        @test !isfile(projected_plot_cfg.output.combined_plot)
        @test !isfile(projected_plot_cfg.projection.projected_combined_plot)

        combined_plot_input = joinpath(dir, "combined_plot.toml")
        write(
            combined_plot_input,
            replace(
                read(projected_input, String),
                "[band.projection]" => "[band.plot]\ntargets = [\"combined\", \"fatband_pdos\"]\n\n[band.projection]",
            ),
        )
        combined_plot_cfg = InputIO.read_input(combined_plot_input)
        Service.run(combined_plot_cfg; make_plot=true)
        @test isfile(combined_plot_cfg.output.combined_plot)
        @test isfile(combined_plot_cfg.projection.projected_combined_plot)
        empty_circle_plot = Base.invokelatest(
            getfield(plotbands_module, :plot_projected_bands),
            projected_result.band_result,
            projected_cfg;
            save=false,
        )
        scatter_alphas = [
            get(series.plotattributes, :markeralpha, nothing)
            for series in empty_circle_plot.series_list
            if series.plotattributes[:seriestype] == :scatter
        ]
        @test !isempty(scatter_alphas)
        @test all(alpha -> isnothing(alpha) || alpha > 0, scatter_alphas)
        single_group_scatter_count = length(scatter_alphas)

        projected_string_input = joinpath(dir, "projected_string.toml")
        write(
            projected_string_input,
            replace(read(projected_input, String), "color_group = [\"orb1\"]" => "color_group = \"orb1\""),
        )
        @test occursin(
            "projection.color_group must be an array of strings",
            error_message(() -> InputIO.read_input(projected_string_input)),
        )

        projected_multi_input = joinpath(dir, "projected_multi.toml")
        write(
            projected_multi_input,
            replace(read(projected_input, String), "color_group = [\"orb1\"]" => "color_group = [\"orb1\", \"orb2\"]"),
        )
        projected_multi_cfg = InputIO.read_input(projected_multi_input)
        @test projected_multi_cfg.projection.color_group == ["orb1", "orb2"]
        group_indices = Base.invokelatest(
            getfield(plotbands_module, :_projection_group_indices),
            projected_result.band_result.projection.labels,
            projected_multi_cfg.projection.color_group,
        )
        @test group_indices == [1, 2]
        multi_circle_plot = Base.invokelatest(
            getfield(plotbands_module, :plot_projected_bands),
            projected_result.band_result,
            projected_multi_cfg;
            save=false,
        )
        multi_scatter_count = count(
            series -> series.plotattributes[:seriestype] == :scatter,
            multi_circle_plot.series_list,
        )
        @test multi_scatter_count == 2 * single_group_scatter_count
        @test !isnothing(projected_result.dos_result.projected)
        pdos_plot = Base.invokelatest(
            getfield(plotbands_module, :plot_projected_dos),
            projected_result.dos_result,
            projected_cfg;
            save=false,
        )
        @test pdos_plot.subplots[1][:legend_position] == :none
        center_marker_count(plt) = count(
            series -> get(series.plotattributes, :linestyle, nothing) == :dot,
            plt.series_list,
        )
        @test center_marker_count(pdos_plot) >= length(projected_result.dos_result.projected.labels)
        projected_combined_plot = Base.invokelatest(
            getfield(plotbands_module, :plot_projected_combined),
            projected_result.band_result,
            projected_result.dos_result,
            projected_cfg,
        )
        @test projected_combined_plot.subplots[2][:legend_position] == :none
        @test center_marker_count(projected_combined_plot) >= length(projected_result.dos_result.projected.labels)
        @test isfile(projected_cfg.projection.weights_data)
        @test isfile(projected_cfg.projection.pdos_data)
        weights = projected_result.band_result.projection.weights
        @test all(abs.(dropdims(sum(weights; dims=3); dims=3) .- 1.0) .< 1e-10)
        pdos_sum = vec(sum(projected_result.dos_result.projected.pdos; dims=2))
        @test maximum(abs.(pdos_sum .- projected_result.dos_result.dos)) < 1e-10

        summary_io = IOBuffer()
        Service.print_summary(projected_result; make_plot=true, io=summary_io)
        summary = String(take!(summary_io))
        @test occursin("Band run complete.", summary)
        @test occursin("DOS integral:", summary)
        @test occursin("DOS center of mass:", summary)
        @test occursin("PDOS center of mass:", summary)
        @test occursin("orb1:", summary)
        @test !occursin("Bands data:", summary)
        @test !occursin("Band projection weights:", summary)
        @test !occursin("Projected bands plot:", summary)
        @test !occursin("DOS data:", summary)
        @test !occursin("PDOS data:", summary)
        @test !occursin("PDOS plot:", summary)
        @test !occursin("Combined plot:", summary)
        @test !occursin("Projected combined plot:", summary)
    end

    io = IOBuffer()
    Service.print_summary(BandCore.Model.BandRunResult(cfg, nothing, nothing, false); make_plot=false, io=io)
    @test occursin("Hermiticity check: failed", String(take!(io)))

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
        BandCore.Model.SpinConfig(true, :vasp544, ("up", "down")),
        cfg.hermiticity_tol,
        cfg.verbose,
    )
    spin_io = IOBuffer()
    Service.print_summary(
        BandCore.Model.BandRunResult(spin_cfg, nothing, nothing, true);
        make_plot=false,
        io=spin_io,
    )
    spin_summary = String(take!(spin_io))
    @test occursin("Hermiticity check: ok\nSpin layout: vasp544", spin_summary)
    @test !occursin("Spin layout check:", spin_summary)
    @test !occursin("WARNING:", spin_summary)
end
