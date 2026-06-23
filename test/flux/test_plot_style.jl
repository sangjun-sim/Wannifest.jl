@testset "plot style overrides" begin
    dir, hr, win = fixture_paths()
    input = write_flux_input(
        joinpath(dir, "styled_plot.flux.toml"),
        hr,
        win,
        "  [1, [1, 2], [0, 0, 0], [0.0, -0.02]],\n  [1, [1, 2], [-1, 0, 0], [0.0, -0.02]],";
        plot_body="""
cell_bounds = [2, 0, 0]
arrow_styles = [
  ["C1:s", 2.5, "#123456"],
]
""",
    )
    out_html = joinpath(dir, "styled_plot_flux.html")
    FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(input);
        output_hr=joinpath(dir, "styled_plot_flux_hr.dat"),
        html_path=out_html,
        validate_roundtrip=true,
    )
    html = read(out_html, String)
    @test occursin("C1 [-2,0,0]", html)
    @test !occursin("C1 [-3,0,0]", html)
    @test occursin("\"x\":[-2.0,-2.5]", html)
    @test occursin("\"x\":[2.0,2.5]", html)
    @test occursin("\"sizeref\":2.5", html)
    @test occursin("#123456", html)
end
