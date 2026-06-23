@testset "explicit cell term validation" begin
    dir, hr, win = fixture_paths()

    malformed_R = write_flux_input(
        joinpath(dir, "malformed_R.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0], [0.0, -0.02]],",
    )
    @test_throws ErrorException FluxCLI.InputIO.read_input(malformed_R)

    no_match = write_flux_input(
        joinpath(dir, "no_match.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [1, 0, 0], [0.0, -0.02]],",
    )
    @test_throws ErrorException FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(no_match);
        output_hr=joinpath(dir, "no_match_flux_hr.dat"),
        validate_roundtrip=true,
    )

    broad = write_flux_input(
        joinpath(dir, "broad.flux.toml"),
        hr,
        win,
        "  [1, [\"C\", \"C\"], [0, 0, 0], [0.0, -0.02]],",
    )
    result = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(broad);
        output_hr=joinpath(dir, "broad_flux_hr.dat"),
        validate_roundtrip=true,
    )
    @test length(result.edges) == 1
    @test result.edges[1].R == (0, 0, 0)

    duplicate = write_flux_input(
        joinpath(dir, "duplicate.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],\n  [1, [2, 1], [0, 0, 0], [0.0, 0.02]],",
    )
    @test_throws ErrorException FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(duplicate);
        output_hr=joinpath(dir, "duplicate_flux_hr.dat"),
        validate_roundtrip=true,
    )
end

@testset "phase sign uses Wannier bra-ket convention" begin
    dir, hr, win = fixture_paths()
    input = write_flux_input(
        joinpath(dir, "phase_direction.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],",
    )
    out_html = joinpath(dir, "phase_direction.html")
    FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(input);
        output_hr=joinpath(dir, "phase_direction_hr.dat"),
        html_path=out_html,
        validate_roundtrip=true,
    )
    html = read(out_html, String)
    @test occursin("\"x\":[0.25],\"y\":[0.0],\"z\":[0.0],\"u\":[0.5]", html)
    @test !occursin("\"x\":[0.0],\"y\":[0.0],\"z\":[0.0],\"u\":[-0.5]", html)
end
