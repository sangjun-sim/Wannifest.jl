function write_minimal_poscar(path::AbstractString)
    write(path, """
minimal POSCAR flux fixture
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
C
2
Direct
0.0 0.0 0.0
0.5 0.0 0.0
""")
    return path
end

function write_duplicate_species_poscar(path::AbstractString)
    write(path, """
duplicate species POSCAR flux fixture
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
C C
1 1
Direct
0.0 0.0 0.0
0.5 0.0 0.0
""")
    return path
end

function write_vasp4_poscar(path::AbstractString)
    write(path, """
vasp4 POSCAR flux fixture
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
2
Direct
0.0 0.0 0.0
0.5 0.0 0.0
""")
    return path
end

function write_zero_hr(path::AbstractString, num_wann::Int)
    hops = Dict{Tuple{Int, Int, Int}, Matrix{ComplexF64}}(
        (0, 0, 0) => zeros(ComplexF64, num_wann, num_wann),
    )
    FluxCore.HrFormat.write_hr_blocks_normalized(path, "zero flux fixture", num_wann, hops)
    return path
end

function write_poscar_flux_input(
    path::AbstractString,
    hr_path::AbstractString,
    poscar_path::AbstractString,
    term_body::AbstractString;
    basis_body::Union{Nothing, AbstractString}=nothing,
)
    basis_text = isnothing(basis_body) ? "" : "\n[flux.basis]\n$(basis_body)\n"
    write(path, """
[flux.run]
hr = "$hr_path"
poscar = "$poscar_path"

[flux.geometry]
search_bounds = [1, 0, 0]
distance_tol = 1.0e-8
	$basis_text
	[[flux.terms]]
	term = [
	$term_body
	]
""")
    return path
end

@testset "POSCAR fallback basis" begin
    dir = mktempdir()
    hr = write_zero_hr(joinpath(dir, "wannier90_hr.dat"), 2)
    poscar = write_minimal_poscar(joinpath(dir, "POSCAR"))
    input = write_poscar_flux_input(
        joinpath(dir, "poscar.flux.toml"),
        hr,
        poscar,
        "  [1, [\"C1\", \"C2\"], [0, 0, 0], [0.0, -0.02]],",
    )
    result = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(input);
        output_hr=joinpath(dir, "poscar_flux_hr.dat"),
        html_path=joinpath(dir, "poscar_flux.html"),
        validate_roundtrip=true,
    )
    @test length(result.edges) == 1
    @test all(edge -> edge.from_label == "C1:orb1", result.edges)
    @test all(edge -> edge.to_label == "C2:orb1", result.edges)
end

@testset "POSCAR duplicate species count policy" begin
    dir = mktempdir()
    hr = write_zero_hr(joinpath(dir, "wannier90_hr.dat"), 4)
    poscar = write_duplicate_species_poscar(joinpath(dir, "POSCAR"))

    same_input = write_poscar_flux_input(
        joinpath(dir, "duplicate_same.flux.toml"),
        hr,
        poscar,
        "  [1, [\"C\", \"C\"], [0, 0, 0], [0.0, -0.03]],";
        basis_body="orbitals_per_atom = [[\"C\", 2], [\"C\", 2]]",
    )
    same_cfg = FluxCLI.InputIO.read_input(same_input)
    @test same_cfg.basis.orbitals_per_atom["C"] == 2

    group_input = write_poscar_flux_input(
        joinpath(dir, "duplicate_group.flux.toml"),
        hr,
        poscar,
        "  [1, [\"C1\", \"C2\"], [0, 0, 0], [0.0, -0.03]],";
        basis_body="orbitals_per_atom = [1, 3]",
    )
    group_cfg = FluxCLI.InputIO.read_input(group_input)
    @test group_cfg.basis.orbitals_per_species_group == [1, 3]
    result = FluxCLI.Service.run(
        group_cfg;
        output_hr=joinpath(dir, "duplicate_group_flux_hr.dat"),
        html_path=joinpath(dir, "duplicate_group_flux.html"),
        validate_roundtrip=true,
    )
    @test any(edge -> edge.from_label == "C1:orb1" && edge.to_label == "C2:orb3", result.edges)

    bad_input = write_poscar_flux_input(
        joinpath(dir, "duplicate_conflict.flux.toml"),
        hr,
        poscar,
        "  [1, [\"C\", \"C\"], [0, 0, 0], [0.0, -0.03]],";
        basis_body="orbitals_per_atom = [[\"C\", 1], [\"C\", 3]]",
    )
    @test_throws ErrorException FluxCLI.InputIO.read_input(bad_input)
end

@testset "POSCAR without element names uses Type labels" begin
    dir = mktempdir()
    hr = write_zero_hr(joinpath(dir, "wannier90_hr.dat"), 2)
    poscar = write_vasp4_poscar(joinpath(dir, "POSCAR"))
    input = write_poscar_flux_input(
        joinpath(dir, "vasp4.flux.toml"),
        hr,
        poscar,
        "  [1, [\"Type1\", \"Type1\"], [0, 0, 0], [0.0, -0.02]],";
        basis_body="orbitals_per_atom = [1]",
    )
    result = FluxCLI.Service.run(
        FluxCLI.InputIO.read_input(input);
        output_hr=joinpath(dir, "vasp4_flux_hr.dat"),
        html_path=joinpath(dir, "vasp4_flux.html"),
        validate_roundtrip=true,
    )
    @test any(edge -> edge.from_label == "Type1_1:orb1" && edge.to_label == "Type1_2:orb1", result.edges)
end

@testset "POSCAR species selector with explicit orbital counts" begin
    dir = mktempdir()
    hr = write_zero_hr(joinpath(dir, "wannier90_hr.dat"), 4)
    poscar = write_minimal_poscar(joinpath(dir, "POSCAR"))
    input = write_poscar_flux_input(
        joinpath(dir, "poscar_species.flux.toml"),
        hr,
        poscar,
        "  [1, [\"C\", \"C\"], [0, 0, 0], [0.0, -0.03]],";
        basis_body="orbitals_per_atom = [[\"C\", 2]]",
    )
    cfg = FluxCLI.InputIO.read_input(input)
    @test cfg.basis.orbitals_per_atom["C"] == 2

    result = FluxCLI.Service.run(
        cfg;
        output_hr=joinpath(dir, "poscar_species_flux_hr.dat"),
        html_path=joinpath(dir, "poscar_species_flux.html"),
        validate_roundtrip=true,
    )
    @test !isempty(result.edges)
    @test any(edge -> edge.from_label == "C1:orb1" && edge.to_label == "C2:orb2", result.edges)
    @test all(edge -> startswith(edge.from_label, "C") && startswith(edge.to_label, "C"), result.edges)
end
