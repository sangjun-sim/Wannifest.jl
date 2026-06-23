function write_square_hr(path::AbstractString; t::Float64=1.0)
    write(path, """
square lattice test
  1
  4
  1 1 1 1
  1  0  0  1  1  $( -t )  0.0
 -1  0  0  1  1  $( -t )  0.0
  0  1  0  1  1  $( -t )  0.0
  0 -1  0  1  1  $( -t )  0.0
""")
end

function write_two_band_hr(path::AbstractString)
    write(path, """
two-band spin layout test
  2
  1
  1
  0  0  0  1  1  -1.0  0.0
  0  0  0  2  1   0.0  0.0
  0  0  0  1  2   0.0  0.0
  0  0  0  2  2   1.0  0.0
""")
end

function write_multi_band_hr(path::AbstractString)
    write(path, """
three-band plot test
  3
  1
  1
  0  0  0  1  1  -1.0  0.0
  0  0  0  2  1   0.0  0.0
  0  0  0  3  1   0.0  0.0
  0  0  0  1  2   0.0  0.0
  0  0  0  2  2   0.0  0.0
  0  0  0  3  2   0.0  0.0
  0  0  0  1  3   0.0  0.0
  0  0  0  2  3   0.0  0.0
  0  0  0  3  3   1.0  0.0
""")
end

function write_contour_input(path::AbstractString; hr::AbstractString="wannier90_hr.dat", extra::AbstractString="")
    write(path, """
[contour.run]
hr = "$hr"
structure = ""
verbose = false

[contour.plane]
axes = ["kx", "ky"]
fixed_axis = "kz"
fixed_value = 0.25
range_x = [-0.25, 0.25]
range_y = [-0.25, 0.25]
mesh = "3x2"

[contour.energy]
shift = 0.5
bands = [1]

$extra
""")
end

@testset "contour input and mesh" begin
    mktempdir() do dir
        write_square_hr(joinpath(dir, "wannier90_hr.dat"))
        input_path = joinpath(dir, "input.contour.toml")
        write_contour_input(input_path)

        cfg = ContourInputIO.read_input(input_path)
        @test cfg.files.hr_path == joinpath(dir, "wannier90_hr.dat")
        @test cfg.files.structure_path == ""
        @test cfg.output.output_dir == joinpath(dir, "outputs")
        @test cfg.spin_layout == :qe
        @test cfg.energy.bands == [1]

        custom_cfg = ContourInputIO.read_input(input_path; output_dir_override="custom")
        @test custom_cfg.output.output_dir == joinpath(dir, "custom")

        mesh = ContourMesh.generate_plane_mesh(cfg.plane)
        @test mesh.nx == 3
        @test mesh.ny == 2
        @test mesh.kpoints[1] == [-0.25, -0.25, 0.25]
        @test mesh.kpoints[2] == [0.0, -0.25, 0.25]
        @test mesh.kpoints[4] == [-0.25, 0.25, 0.25]
        @test ContourMesh.grid_index(mesh, 1) == (iy=1, ix=1)
        @test ContourMesh.grid_index(mesh, 4) == (iy=2, ix=1)
        @test occursin("outside", error_message(() -> ContourMesh.grid_index(mesh, 7)))

        opts = Wannifest.ContourCLI.parse_args(["--input", input_path, "--output-dir", "custom", "--no-plot"])
        @test opts.input_path == input_path
        @test opts.output_dir == "custom"
        @test opts.no_plot

        bad_axes = joinpath(dir, "bad_axes.toml")
        write(bad_axes, """
[contour.run]
hr = "wannier90_hr.dat"
structure = ""
verbose = false

[contour.plane]
axes = ["kx", "kx"]
""")
        @test occursin("distinct", error_message(() -> ContourInputIO.read_input(bad_axes)))

        bad_mesh = joinpath(dir, "bad_mesh.toml")
        write(bad_mesh, """
[contour.run]
hr = "wannier90_hr.dat"
structure = ""
verbose = false

[contour.plane]
mesh = "3"
""")
        @test occursin("form NxM", error_message(() -> ContourInputIO.read_input(bad_mesh)))

        bad_spin = joinpath(dir, "bad_spin.toml")
        write_contour_input(bad_spin; extra="""
[contour.spin]
layout = "blocked"
""")
        @test occursin("Unsupported contour.spin.layout", error_message(() -> ContourInputIO.read_input(bad_spin)))

        output_key = joinpath(dir, "output_key.toml")
        write_contour_input(output_key; extra="""
[contour.output]
data = "not_allowed.dat"
""")
        @test occursin("Unsupported contour table", error_message(() -> ContourInputIO.read_input(output_key)))
    end
end

@testset "contour energy surface and output" begin
    mktempdir() do dir
        write_square_hr(joinpath(dir, "wannier90_hr.dat"))
        input_path = joinpath(dir, "input.contour.toml")
        write_contour_input(input_path)
        cfg = ContourInputIO.read_input(input_path)
        hr = ContourCore.WannierHrIO.read_hr(cfg.files.hr_path; spin_layout=cfg.spin_layout)
        mesh = ContourMesh.generate_plane_mesh(cfg.plane)
        surface = ContourSurface.compute_energy_surface(hr, mesh, cfg)

        @test size(surface.energies) == (2, 3, 1)
        for iy in eachindex(surface.y_axis), ix in eachindex(surface.x_axis)
            kx = surface.x_axis[ix]
            ky = surface.y_axis[iy]
            expected = -2 * (cos(2π * kx) + cos(2π * ky)) - cfg.energy.shift
            @test surface.energies[iy, ix, 1] ≈ expected
        end

        ContourPlotService._plotcontour_module_ref[] = nothing
        result = ContourService.run(cfg; make_plot=false)
        @test result.hermiticity_ok
        @test isempty(result.plot_handles)
        @test isnothing(ContourPlotService._plotcontour_module_ref[])
        @test isfile(result.data_path)
        data_text = read(result.data_path, String)
        @test occursin("# columns = ix iy kx ky kz band energy", data_text)
        @test count(line -> !startswith(line, "#") && !isempty(strip(line)), split(data_text, '\n')) == 6

        plot_input = joinpath(dir, "plot.contour.toml")
        write_contour_input(plot_input; extra="""
[contour.plot]
mode = "heatmap"
interactive = false
size = [320, 240]
""")
        plot_result = ContourService.run(ContourInputIO.read_input(plot_input); make_plot=true)
        @test length(plot_result.plot_handles) == 1
        @test isfile(joinpath(dir, "outputs", "plots", "contour_heatmap.png"))
        plot_module = ContourPlotService._plotcontour_module()
        @test Base.invokelatest(isdefined, plot_module, :open_interactive)
        opened = String[]
        open_fn = Base.invokelatest(getfield, plot_module, :open_interactive)
        html_paths = Base.invokelatest(
            open_fn,
            plot_result.surface,
            plot_result.config;
            opener=path -> push!(opened, path),
        )
        @test html_paths == [joinpath(dir, "outputs", "plots", "contour_heatmap.html")]
        @test opened == html_paths
        @test occursin("scrollZoom: true", read(only(html_paths), String))
    end
end

@testset "contour plots selected bands together" begin
    mktempdir() do dir
        write_multi_band_hr(joinpath(dir, "wannier90_hr.dat"))
        input_path = joinpath(dir, "multi.contour.toml")
        write(input_path, """
[contour.run]
hr = "wannier90_hr.dat"
structure = ""
verbose = false

[contour.plane]
mesh = "2x2"

[contour.energy]
bands = [1, 2, 3]

[contour.plot]
mode = "both"
interactive = false
size = [320, 240]
""")
        result = ContourService.run(ContourInputIO.read_input(input_path); make_plot=true)
        @test length(result.plot_handles) == 2
        plot_dir = joinpath(dir, "outputs", "plots")
        @test isfile(joinpath(plot_dir, "contour_surface.png"))
        @test isfile(joinpath(plot_dir, "contour_contour.png"))
        @test !any(contains("band"), readdir(plot_dir))

        plot_module = ContourPlotService._plotcontour_module()
        open_fn = Base.invokelatest(getfield, plot_module, :open_interactive)
        html_paths = Base.invokelatest(open_fn, result.surface, result.config; open_browser=false)
        @test sort(basename.(html_paths)) == ["contour_contour.html", "contour_surface.html"]
        @test all(isfile, html_paths)
    end
end

@testset "contour spinful input uses band index only" begin
    mktempdir() do dir
        write_two_band_hr(joinpath(dir, "wannier90_hr.dat"))
        input_path = joinpath(dir, "spinful.contour.toml")
        write(input_path, """
[contour.run]
hr = "wannier90_hr.dat"
structure = ""
verbose = false

[contour.plane]
mesh = "1x1"

[contour.energy]
bands = [1, 2]

[contour.spin]
layout = "vasp544"
""")
        cfg = ContourInputIO.read_input(input_path)
        result = ContourService.run(cfg; make_plot=false)
        @test cfg.spin_layout == :vasp544
        @test result.surface.bands == [1, 2]
        @test result.surface.energies[1, 1, 1] ≈ -1.0
        @test result.surface.energies[1, 1, 2] ≈ 1.0
        @test !occursin("spin_z", read(result.data_path, String))
    end
end
