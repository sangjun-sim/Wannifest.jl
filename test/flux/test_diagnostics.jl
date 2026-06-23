function write_diagnostic_flux_input(
    path::AbstractString,
    hr_path::AbstractString,
    win_path::AbstractString,
    term_body::AbstractString,
    diagnostic_body::AbstractString,
)
    write(path, """
[flux.run]
hr = "$hr_path"
win = "$win_path"

[flux.geometry]
search_bounds = [1, 0, 0]
distance_tol = 1.0e-8

[flux.plot]
interactive = false

[flux.diagnostic]
$diagnostic_body

[[flux.terms]]
term = [
$term_body
]
""")
    return path
end

function write_triangle_win(path::AbstractString)
    write(path, """
num_wann = 3

begin unit_cell_cart
ang
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart

begin atoms_frac
C 0.0 0.0 0.0
C 0.5 0.0 0.0
C 0.0 0.5 0.0
end atoms_frac

begin projections
C: s
end projections
""")
    return path
end

function write_triangle_hr(path::AbstractString)
    hops = Dict{Tuple{Int, Int, Int}, Matrix{ComplexF64}}(
        (0, 0, 0) => zeros(ComplexF64, 3, 3),
    )
    FluxCore.HrFormat.write_hr_blocks_normalized(path, "triangle flux diagnostic fixture", 3, hops)
    return path
end

function triangle_fixture_paths()
    dir = mktempdir()
    hr = write_triangle_hr(joinpath(dir, "wannier90_hr.dat"))
    win = write_triangle_win(joinpath(dir, "wannier90.win"))
    return dir, hr, win
end

const BALANCED_DIAGNOSTIC_TERMS = """
  [1, [1, 2], [0, 0, 0], [0.0, -0.01]],
  [2, [2, 3], [0, 0, 0], [0.0, -0.01]],
  [1, [1, 3], [0, 0, 0], [0.0, 0.01]],
"""

const TRIANGLE_DIAGNOSTIC = """
continuity = true
continuity_tol = 1.0e-12
plaquettes = [
  ["cell-triangle", [
    [1, [0, 0, 0]],
    [2, [0, 0, 0]],
    [3, [0, 0, 0]],
  ]],
]
"""

@testset "flux diagnostic parser" begin
    dir, hr, win = triangle_fixture_paths()
    input = write_diagnostic_flux_input(
        joinpath(dir, "diagnostic.flux.toml"),
        hr,
        win,
        BALANCED_DIAGNOSTIC_TERMS,
        TRIANGLE_DIAGNOSTIC,
    )
    cfg = FluxCLI.InputIO.read_input(input)
    @test cfg.diagnostic.enabled
    @test cfg.diagnostic.continuity
    @test cfg.diagnostic.continuity_tol == 1.0e-12
    @test length(cfg.diagnostic.plaquettes) == 1
    @test length(cfg.diagnostic.plaquettes[1].vertices) == 3

    bad_output = replace(read(input, String), "continuity = true" => "continuity = true\noutput = \"bad.tsv\"")
    bad_path = joinpath(dir, "bad_output.flux.toml")
    write(bad_path, bad_output)
    @test_throws ErrorException FluxCLI.InputIO.read_input(bad_path)

    bad_loop = replace(read(input, String), "[3, [0, 0, 0]]," => "")
    bad_loop_path = joinpath(dir, "bad_loop.flux.toml")
    write(bad_loop_path, bad_loop)
    @test_throws ErrorException FluxCLI.InputIO.read_input(bad_loop_path)
end

@testset "flux plaquette imaginary diagnostic" begin
    dir, hr, win = triangle_fixture_paths()
    input = write_diagnostic_flux_input(
        joinpath(dir, "plaquette.flux.toml"),
        hr,
        win,
        BALANCED_DIAGNOSTIC_TERMS,
        TRIANGLE_DIAGNOSTIC,
    )
    out_diag = joinpath(dir, "diagnostic.tsv")
    result = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(input);
        output_hr=joinpath(dir, "plaquette_hr.dat"),
        diagnostic_path=out_diag,
        make_html=false,
        validate_roundtrip=true,
    )
    @test isfile(out_diag)
    @test !isnothing(result.diagnostic)
    plaquette = only(result.diagnostic.plaquettes)
    @test plaquette.name == "cell-triangle"
    @test plaquette.edge_imags ≈ [-0.01, -0.01, -0.01]
    @test plaquette.imag_sum ≈ -0.03
    @test result.diagnostic.continuity_passed
    @test all(row -> row.passed, result.diagnostic.site_flows)
    text = read(out_diag, String)
    @test occursin("# plaquette_imag", text)
    @test occursin("cell-triangle\t-0.030000000000\t3\t-0.010000000000,-0.010000000000,-0.010000000000", text)
    @test occursin("# continuity", text)
end

@testset "flux diagnostic validation and continuity failure" begin
    dir, hr, win = triangle_fixture_paths()
    missing = write_diagnostic_flux_input(
        joinpath(dir, "missing_edge.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.01]],",
        TRIANGLE_DIAGNOSTIC,
    )
    @test_throws ErrorException FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(missing);
        output_hr=joinpath(dir, "missing_hr.dat"),
        make_html=false,
    )

    broad_selector = write_diagnostic_flux_input(
        joinpath(dir, "broad_selector.flux.toml"),
        hr,
        win,
        BALANCED_DIAGNOSTIC_TERMS,
        replace(TRIANGLE_DIAGNOSTIC, "[1, [0, 0, 0]]," => "[\"C\", [0, 0, 0]],"),
    )
    @test_throws ErrorException FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(broad_selector);
        output_hr=joinpath(dir, "broad_selector_hr.dat"),
        make_html=false,
    )

    flow_input = write_diagnostic_flux_input(
        joinpath(dir, "flow_fail.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],",
        "continuity = true\ncontinuity_tol = 1.0e-12",
    )
    result = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(flow_input);
        output_hr=joinpath(dir, "flow_fail_hr.dat"),
        make_html=false,
    )
    @test !result.diagnostic.continuity_passed
    @test any(row -> !row.passed, result.diagnostic.site_flows)

    zero_flow_input = write_diagnostic_flux_input(
        joinpath(dir, "zero_flow.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [1.0, 0.0]],",
        "continuity = true\ncontinuity_tol = 1.0e-12",
    )
    zero_flow = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(zero_flow_input);
        output_hr=joinpath(dir, "zero_flow_hr.dat"),
        make_html=false,
    )
    @test zero_flow.diagnostic.continuity_passed
    @test all(row -> row.flow_in == 0.0 && row.flow_out == 0.0, zero_flow.diagnostic.site_flows)
end

@testset "flux diagnostic CLI paths" begin
    dir, hr, win = triangle_fixture_paths()
    input = write_diagnostic_flux_input(
        joinpath(dir, "cli_diag.flux.toml"),
        hr,
        win,
        BALANCED_DIAGNOSTIC_TERMS,
        TRIANGLE_DIAGNOSTIC,
    )
    default_path = joinpath(dir, "outputs", "diagnostic", "wannier90_flux_diagnostic.tsv")
    @test FluxCLI.run_main(["--input", input, "--output", joinpath(dir, "cli_hr.dat"), "--no-html"]) == 0
    @test isfile(default_path)

    custom_path = joinpath(dir, "custom.tsv")
    @test FluxCLI.run_main([
        "--input", input,
        "--output", joinpath(dir, "cli_custom_hr.dat"),
        "--no-html",
        "--diagnostic", custom_path,
    ]) == 0
    @test isfile(custom_path)

    suppressed_dir, suppressed_hr, suppressed_win = triangle_fixture_paths()
    suppressed_input = write_diagnostic_flux_input(
        joinpath(suppressed_dir, "suppressed.flux.toml"),
        suppressed_hr,
        suppressed_win,
        BALANCED_DIAGNOSTIC_TERMS,
        TRIANGLE_DIAGNOSTIC,
    )
    suppressed_path = joinpath(suppressed_dir, "outputs", "diagnostic", "wannier90_flux_diagnostic.tsv")
    @test FluxCLI.run_main([
        "--input", suppressed_input,
        "--output", joinpath(suppressed_dir, "suppressed_hr.dat"),
        "--no-html",
        "--no-diagnostic",
    ]) == 0
    @test !isfile(suppressed_path)
end
