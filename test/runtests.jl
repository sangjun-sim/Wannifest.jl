using Test
using LinearAlgebra
using TOML
using Wannifest

const REPO = dirname(@__DIR__)
const Band = Wannifest.BandCLI.BandCore
const Contour = Wannifest.ContourCLI.ContourCore
const Flux = Wannifest.FluxCLI.FluxCore
const Superham = Wannifest.SuperhamCLI.SuperhamCore
const Soc = Wannifest.SocCLI.SocCore

function write_input(path, config)
    open(io -> TOML.print(io, config), path, "w")
    return path
end

@testset "Wannifest" begin
    mktempdir() do tmp
        calc = mkpath(joinpath(tmp, "calculation"))
        elsewhere = mkpath(joinpath(tmp, "elsewhere"))
        hr_path = joinpath(calc, "chain_hr.dat")
        # One orbital: E(k) = 0.4 + 2 cos(2πkx), including raw ndegen weights.
        write(hr_path, """
        Analytic chain
        1
        3
        2 1 2
        -1 0 0 1 1 2.0 0.0
         0 0 0 1 1 0.4 0.0
         1 0 0 1 1 2.0 0.0
        """)
        write(joinpath(calc, "POSCAR"), "Chain\n1\n1 0 0\n0 1 0\n0 0 1\nX\n1\nDirect\n0 0 0\n")
        write(joinpath(calc, "KPOINTS"), "Chain\n3\nLine-mode\nReciprocal\n0 0 0 ! G\n0.5 0 0 ! X\n")

        @testset "Command inference and invalid inputs" begin
            path = joinpath(tmp, "input.toml")
            for command in ("band", "contour", "flux", "superham", "soc")
                write_input(path, Dict(command => Dict("run" => Dict())))
                @test Wannifest.resolve_command([path]) == (command, ["--input", path])
                @test Wannifest.resolve_command(["--input", path])[1] == command
            end
            write_input(path, Dict("band" => Dict(), "soc" => Dict()))
            @test_throws ErrorException Wannifest.resolve_command([path])
            write_input(path, Dict("unrecognized" => Dict()))
            @test_throws ErrorException Wannifest.resolve_command([path])
            @test_throws ArgumentError Wannifest.resolve_command(["--input"])
            @test_throws ArgumentError Wannifest.run_wannifest(["band", "--unrecognized"])
        end

        @testset "Distributed example inputs" begin
            examples = [
                ("band/rucl3", "band"), ("contour/graphene", "contour"),
                ("fatband/1tp-mos2", "band"), ("fatband/csv3sb5", "band"),
                ("flux/graphene", "flux"), ("oam/rucl3", "band"),
                ("sam/bi2se3", "band"), ("soc/fe", "soc"),
            ]
            for (example, command) in examples
                path = joinpath(REPO, "examples", example, "input.toml")
                @test Wannifest.infer_command_from_input(path) == command
                cli = getproperty(Wannifest, Symbol(uppercasefirst(command) * "CLI"))
                @test !isnothing(cli.InputIO.read_input(path))
            end
        end

        @testset "Hamiltonian normalization, wsvec, and Hermiticity" begin
            hr = Band.WannierHrIO.read_hr(hr_path)
            @test hr.num_wann == 1
            @test Band.WannierHrIO.pair_hermiticity_error(hr) == 0
            for (kx, energy) in ((0.0, 2.4), (0.25, 0.4), (0.5, -1.6))
                @test Band.WannierEigensystem.solve_kpoint_values(hr, [kx, 0, 0]) ≈ [energy]
            end
            types = Band.WannierTypes
            wsvec = types.WsvecTable(Dict(
                ((-1, 0, 0), 1, 1) => types.WsvecEntry(2, [0 -1; 0 0; 0 0]),
                ((0, 0, 0), 1, 1) => types.WsvecEntry(1, zeros(Int, 3, 1)),
                ((1, 0, 0), 1, 1) => types.WsvecEntry(2, [0 1; 0 0; 0 0]),
            ))
            # The extra shifts give E(k) = 0.4 + cos(2πkx) + cos(4πkx).
            for (kx, energy) in ((0.0, 2.4), (0.25, -0.6), (0.5, 0.4))
                @test Band.WannierEigensystem.solve_kpoint_values(hr, [kx, 0, 0]; wsvec) ≈ [energy]
            end
            hr.hoppings[(1, 0, 0)][1, 1] += im
            @test_throws ErrorException Band.WannierEigensystem.solve_kpoint_values(hr, zeros(3))
            broken = joinpath(tmp, "broken_hr.dat")
            write(broken, "Truncated\n1\n1\n1\n")
            @test_throws ErrorException Band.WannierHrIO.read_hr(broken)
        end

        @testset "Band and contour from another working directory" begin
            band_path = write_input(joinpath(calc, "band.toml"), Dict("band" => Dict(
                "run" => Dict("mode" => "bands", "hr" => "chain_hr.dat", "kpoints" => "KPOINTS", "structure" => "POSCAR", "verbose" => false),
                "energy" => Dict("shift" => 0.4),
                "plot" => Dict("targets" => ["band"]),
            )))
            contour_path = write_input(joinpath(calc, "contour.toml"), Dict("contour" => Dict(
                "run" => Dict("hr" => "chain_hr.dat", "structure" => ""),
                "plane" => Dict("mesh" => "3x2", "range_x" => [0.0, 0.5], "range_y" => [0.0, 0.5]),
                "energy" => Dict("shift" => 0.4),
            )))
            cd(elsewhere) do
                cfg = Band.InputIO.read_input(band_path)
                @test cfg.files.hr_path == hr_path
                result = Band.Service.run(cfg; make_plot=false)
                @test vec(result.band_result.eigenvalues) ≈ [2.0, 0.0, -2.0] atol=1e-12
                @test isfile(cfg.output.bands_data)
                cfg2 = Contour.InputIO.read_input(contour_path)
                result2 = Contour.Service.run(cfg2; make_plot=false)
                @test result2.surface.energies[:, :, 1] ≈ [2.0 0.0 -2.0; 2.0 0.0 -2.0] atol=1e-12
                @test isfile(result2.data_path)
                @test !ispath(joinpath(elsewhere, "outputs"))
            end
            cmd = `$(Base.julia_cmd()) --startup-file=no $(joinpath(REPO, "main.jl")) $contour_path --no-plot`
            # Start the CLI with its normal load path, outside Pkg.test's sandbox.
            cmd = addenv(Cmd(cmd; dir=elsewhere), "JULIA_LOAD_PATH" => nothing)
            @test success(pipeline(cmd, stdout=devnull, stderr=stderr))
        end

        @testset "Flux export and continuity" begin
            cfg = Flux.InputIO.read_input(joinpath(REPO, "examples/flux/graphene/input.toml"))
            result = Flux.Service.run(cfg; output_hr=joinpath(tmp, "flux_hr.dat"),
                diagnostic_path=joinpath(tmp, "flux.tsv"), make_html=false, validate_roundtrip=true)
            @test result.diagnostic.continuity_passed
            @test length(result.edges) == 6
            @test result.roundtrip_validated
            hr = Band.WannierHrIO.read_hr(result.output_hr)
            @test Band.WannierHrIO.pair_hermiticity_error(hr) < 1e-12
            @test any(H -> any(x -> !iszero(imag(x)), H), values(hr.hoppings))
            @test isfile(result.diagnostic_path)
        end

        @testset "Supercell band folding" begin
            path = write_input(joinpath(calc, "superham.toml"), Dict("superham" => Dict(
                "files" => Dict("hr" => "chain_hr.dat", "structure" => "POSCAR", "output_hr" => "super_hr.dat"),
                "supercell" => Dict("matrix" => [[2, 0, 0], [0, 1, 0], [0, 0, 1]]),
            )))
            result = Superham.Service.run(Superham.InputIO.read_input(path))
            @test result.multiplicity == 2
            @test result.supercell_eigenvalues ≈ [-1.6, 2.4]
            @test result.folded_spectrum_diff < 1e-12
            exported = Band.WannierHrIO.read_hr(result.output_hr)
            @test Band.WannierEigensystem.solve_kpoint_values(exported, zeros(3)) ≈ [-1.6, 2.4]
        end

        @testset "Orbital and spin angular momentum" begin
            for (l, operators) in ((1, Band.AtomicAngularMomentum.p_operators()), (2, Band.AtomicAngularMomentum.d_operators()))
                x, y, z = operators
                @test x * y - y * x ≈ im * z
                @test x * x + y * y + z * z ≈ l * (l + 1) * I
                @test all(ishermitian, operators)
            end
            for spins in (("up", "dn"), ("dn", "up"))
                entries = [(index=i, site="X1", orbital="s", spin=spins[i]) for i in 1:2]
                ops = Band.AtomicSpin.build_spin_operators(2, entries)
                psi = zeros(ComplexF64, 2, 1)
                psi[findfirst(==("up"), spins), 1] = 1
                @test vec(Band.AtomicSpin.spin_expectations(ops, psi)) ≈ [0, 0, 0.5, 0.5, 0.75]
                @test_throws ErrorException Band.AtomicSpin.build_spin_operators(1, entries[1:1])
            end
        end

        @testset "SOC k matching, Fourier sign, and binary output" begin
            kpoints = [0.0 0.25; 0.0 0.0; 0.0 0.0]
            @test Soc.KpointMatch.match_kpoints(kpoints, kpoints[:, [2, 1]]; atol=1e-10, label="test") == [2, 1]
            @test_throws ErrorException Soc.KpointMatch.match_kpoints(kpoints, zeros(3, 2); atol=1e-10, label="duplicates")
            @test_throws ErrorException Soc.KpointMatch.match_kpoints(kpoints, kpoints .+ 1; atol=1e-10, label="missing")
            Tk = reshape(ComplexF64[1, 2im], 1, 1, 2)
            Rs = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]
            @test vec(Soc.Tmat.tmat_r(Tk, kpoints, Rs)) ≈ [0.5 + im, -0.5, 1.5]
            for labels in (kpoints, [0 1; 0 0; 0 0])
                path = joinpath(tmp, "tensor.bin")
                Soc.BinaryIO.write_tensor(path, "tmat_k", Tk, labels)
                tensor = Soc.BinaryIO.read_tensor(path)
                @test tensor.name == "tmat_k"
                @test tensor.data == Tk
                @test tensor.labels == labels
                open(io -> write(io, UInt8(0)), path, "a")
                @test_throws ErrorException Soc.BinaryIO.read_tensor(path)
            end
            path = joinpath(tmp, "nscf.in")
            write(path, "&SYSTEM\n nosym=.true., noinv=.true., ! noinv=.false. in a comment\n/\n")
            @test isnothing(Soc.QeInputCheck.check_no_symmetry(path))
            write(path, "&SYSTEM\n nosym=.true., ! noinv=.true. only in a comment\n/\n")
            @test_throws ErrorException Soc.QeInputCheck.check_no_symmetry(path)
        end
    end
end
