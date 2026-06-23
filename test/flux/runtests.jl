using Test
using Wannifest

const FluxCLI = Wannifest.FluxCLI
const FluxCore = FluxCLI.FluxCore

function write_minimal_win(path::AbstractString)
    write(path, """
num_wann = 2

begin unit_cell_cart
ang
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
C 0.0 0.0 0.0
C 0.5 0.0 0.0
end atoms_frac

begin projections
C: s
end projections
""")
    return path
end

function write_minimal_hr(path::AbstractString)
    hops = Dict{Tuple{Int, Int, Int}, Matrix{ComplexF64}}(
        (0, 0, 0) => zeros(ComplexF64, 2, 2),
    )
    FluxCore.HrFormat.write_hr_blocks_normalized(path, "minimal flux fixture", 2, hops)
    return path
end

function write_partner_hr(path::AbstractString)
    h0 = zeros(ComplexF64, 2, 2)
    hm = zeros(ComplexF64, 2, 2)
    hp = zeros(ComplexF64, 2, 2)
    hm[1, 2] = ComplexF64(0.5, 0.0)
    hp[2, 1] = ComplexF64(0.5, 0.0)
    hops = Dict{Tuple{Int, Int, Int}, Matrix{ComplexF64}}(
        (0, 0, 0) => h0,
        (-1, 0, 0) => hm,
        (1, 0, 0) => hp,
    )
    FluxCore.HrFormat.write_hr_blocks_normalized(path, "partner flux fixture", 2, hops)
    return path
end

function write_flux_input(
    path::AbstractString,
    hr_path::AbstractString,
    win_path::AbstractString,
    term_body::AbstractString;
    plot_body::AbstractString="",
)
    write(path, """
[flux.run]
hr = "$hr_path"
win = "$win_path"

[flux.geometry]
search_bounds = [1, 0, 0]
distance_tol = 1.0e-8

[flux.plot]
interactive = true
$plot_body

	[[flux.terms]]
	term = [
	$term_body
	]
	""")
    return path
end

function fixture_paths()
    dir = mktempdir()
    hr = write_minimal_hr(joinpath(dir, "wannier90_hr.dat"))
    win = write_minimal_win(joinpath(dir, "wannier90.win"))
    return dir, hr, win
end

foreach(name -> include(joinpath(@__DIR__, name)), ("test_poscar_fallback.jl", "test_plot_style.jl", "test_explicit_cell.jl", "test_diagnostics.jl"))

@testset "flux input parser" begin
    dir, hr, win = fixture_paths()
    input = write_flux_input(
        joinpath(dir, "input.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],",
    )
    cfg = FluxCLI.InputIO.read_input(input)
    @test cfg.files.hr_path == hr
    @test cfg.files.win_path == win
    @test isnothing(cfg.files.poscar_path)
    @test isnothing(cfg.plot.cell_bounds)
    @test isempty(cfg.plot.arrow_styles)
    @test cfg.terms[1].rows[1].nn == 1
    @test cfg.terms[1].rows[1].R == (0, 0, 0)
    @test cfg.terms[1].rows[1].value == ComplexF64(0.0, -0.02)

    styled = write_flux_input(
        joinpath(dir, "styled.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],";
        plot_body="""
cell_bounds = [2, 0, 0]
arrow_styles = [
  ["C1:s", 2.5, "#123456"],
]
""",
    )
    styled_cfg = FluxCLI.InputIO.read_input(styled)
    @test styled_cfg.plot.cell_bounds == (2, 0, 0)
    @test length(styled_cfg.plot.arrow_styles) == 1
    @test styled_cfg.plot.arrow_styles[1].selector == "C1:s"
    @test styled_cfg.plot.arrow_styles[1].size == 2.5
    @test styled_cfg.plot.arrow_styles[1].color == "#123456"

    bad = write_flux_input(
        joinpath(dir, "bad.flux.toml"),
        hr,
        win,
        "  [1, [1, \"C2\"], [0, 0, 0], [0.0, -0.02]],",
    )
    @test_throws ErrorException FluxCLI.InputIO.read_input(bad)

    missing_win = joinpath(dir, "missing_win.flux.toml")
    write(missing_win, """
[flux.run]
hr = "$hr"

[[flux.terms]]
term = [
	  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],
	]
	""")
    @test_throws ErrorException FluxCLI.InputIO.read_input(missing_win)

    bad_orientation = replace(read(input, String), "[[flux.terms]]" => "[[flux.terms]]\norientation = \"canonical\"")
    write(joinpath(dir, "bad_orientation.flux.toml"), bad_orientation)
    @test_throws ErrorException FluxCLI.InputIO.read_input(joinpath(dir, "bad_orientation.flux.toml"))

    bad_row = write_flux_input(
        joinpath(dir, "bad_row.flux.toml"),
        hr,
        win,
        "  [1, [1, 2]],",
    )
    @test_throws ErrorException FluxCLI.InputIO.read_input(bad_row)

    bad_style = write_flux_input(
        joinpath(dir, "bad_style.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],";
        plot_body="arrow_styles = [[\"C1:s\", -1.0, \"#123456\"]]",
    )
    @test_throws ErrorException FluxCLI.InputIO.read_input(bad_style)
end

@testset "numeric endpoint flux run" begin
    dir, hr, win = fixture_paths()
    input = write_flux_input(
        joinpath(dir, "numeric.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],\n  [1, [1, 2], [-1, 0, 0], [0.0, -0.02]],",
    )
    out_hr = joinpath(dir, "numeric_flux_hr.dat")
    out_html = joinpath(dir, "numeric_flux.html")
    result = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(input);
        output_hr=out_hr,
        html_path=out_html,
        validate_roundtrip=true,
    )
    @test isfile(out_hr)
    @test isfile(out_html)
    @test length(result.edges) == 2

    parsed = FluxCore.WannierHrIO.read_hr(out_hr)
    hops = FluxCore.WannierHrIO.normalized_hoppings(parsed)
    @test all(==(1), values(parsed.ndegen))
    @test hops[(0, 0, 0)][1, 2] == ComplexF64(0.0, -0.02)
    @test hops[(0, 0, 0)][2, 1] == ComplexF64(0.0, 0.02)
    @test hops[(-1, 0, 0)][1, 2] == ComplexF64(0.0, -0.02)
    @test hops[(1, 0, 0)][2, 1] == ComplexF64(0.0, 0.02)
    html = read(out_html, String)
    @test occursin("Plotly.newPlot", html)
    @test occursin("scatter3d", html)
    @test occursin("cone", html)
    sizeref_match = match(r"\"sizeref\":([0-9.]+)", html)
    @test !isnothing(sizeref_match)
    @test parse(Float64, sizeref_match.captures[1]) > 0
    @test occursin("\"x\":[1.0,0.5]", html)
    @test occursin("\"x\":[0.25],\"y\":[0.0],\"z\":[0.0],\"u\":[0.5]", html)
    @test !occursin("\"x\":[-0.5,0.0]", html)
    @test occursin("supercell sites", html)
    @test occursin("C1 [-1,0,0]", html)
    @test occursin("C2 [1,0,0]", html)
end

@testset "numeric endpoint preserves requested direction" begin
    dir, hr, win = fixture_paths()
    input = write_flux_input(
        joinpath(dir, "reverse_numeric.flux.toml"),
        hr,
        win,
        "  [1, [2, 1], [0, 0, 0], [0.0, -0.04]],\n  [1, [2, 1], [1, 0, 0], [0.0, -0.04]],",
    )
    out_hr = joinpath(dir, "reverse_numeric_flux_hr.dat")
    result = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(input);
        output_hr=out_hr,
        html_path=joinpath(dir, "reverse_numeric_flux.html"),
        validate_roundtrip=true,
    )
    @test length(result.edges) == 2
    @test all(edge -> edge.from_index == 2 && edge.to_index == 1, result.edges)

    parsed = FluxCore.WannierHrIO.read_hr(out_hr)
    hops = FluxCore.WannierHrIO.normalized_hoppings(parsed)
    @test hops[(0, 0, 0)][2, 1] == ComplexF64(0.0, -0.04)
    @test hops[(0, 0, 0)][1, 2] == ComplexF64(0.0, 0.04)
    @test hops[(1, 0, 0)][2, 1] == ComplexF64(0.0, -0.04)
    @test hops[(-1, 0, 0)][1, 2] == ComplexF64(0.0, 0.04)
end

@testset "flux preserves existing Hermitian partner blocks" begin
    dir, _, win = fixture_paths()
    hr = write_partner_hr(joinpath(dir, "partner_hr.dat"))
    input = write_flux_input(
        joinpath(dir, "partner.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],\n  [1, [1, 2], [-1, 0, 0], [0.0, -0.02]],",
    )
    out_hr = joinpath(dir, "partner_flux_hr.dat")
    result = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(input);
        output_hr=out_hr,
        html_path=joinpath(dir, "partner_flux.html"),
        validate_roundtrip=true,
    )
    @test length(result.edges) == 2

    parsed = FluxCore.WannierHrIO.read_hr(out_hr)
    hops = FluxCore.WannierHrIO.normalized_hoppings(parsed)
    @test hops[(-1, 0, 0)][1, 2] == ComplexF64(0.5, -0.02)
    @test hops[(1, 0, 0)][2, 1] == ComplexF64(0.5, 0.02)
end

@testset "string endpoint flux run" begin
    dir, hr, win = fixture_paths()
    input = write_flux_input(
        joinpath(dir, "string.flux.toml"),
        hr,
        win,
        "  [1, [\"C1\", \"C2\"], [0, 0, 0], [0.0, -0.03]],\n  [1, [\"C1\", \"C2\"], [-1, 0, 0], [0.0, -0.03]],",
    )
    result = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(input);
        output_hr=joinpath(dir, "string_flux_hr.dat"),
        html_path=joinpath(dir, "string_flux.html"),
        validate_roundtrip=true,
    )
    @test length(result.edges) == 2
    @test all(edge -> edge.from_label == "C1:s", result.edges)
    @test all(edge -> edge.to_label == "C2:s", result.edges)
end

@testset "CLI smoke" begin
    dir, hr, win = fixture_paths()
    input = write_flux_input(
        joinpath(dir, "cli.flux.toml"),
        hr,
        win,
        "  [1, [\"C1\", \"C2\"], [0, 0, 0], [0.0, -0.01]],",
    )
    out_hr = joinpath(dir, "cli_flux_hr.dat")
    out_html = joinpath(dir, "cli_flux.html")
    @test FluxCLI.run_main([
        "--input", input,
        "--output", out_hr,
        "--html", out_html,
        "--validate-roundtrip",
    ]) == 0
    @test isfile(out_hr)
    @test isfile(out_html)
end
