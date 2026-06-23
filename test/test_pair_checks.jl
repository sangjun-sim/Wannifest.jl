@testset "shared Hermiticity pair checks" begin
    H = ComplexF64[1.0 2.0 + 3.0im; 4.0 5.0]
    Hm = H'
    @test Hermiticity.pair_hermiticity_error(H, Hm) == 0.0

    bad = copy(Matrix(Hm))
    bad[1, 2] += 0.25
    @test Hermiticity.pair_hermiticity_error(H, bad) ≈ 0.25

    hops = Dict(
        (0, 0, 0) => ComplexF64[1.0 0.0; 0.0 2.0],
        (1, 0, 0) => H,
        (-1, 0, 0) => Matrix(Hm),
    )
    @test PairChecks.pair_dict_max_error(hops) == 0.0
    @test PairChecks.check_hr_pair_symmetry(hops) == 0.0

    missing = Dict((1, 0, 0) => H)
    @test PairChecks.pair_dict_max_error(missing) == Inf
    @test occursin(
        "Missing Hermitian partner",
        error_message(() -> PairChecks.check_hr_pair_symmetry(missing)),
    )

    bad_hops = copy(hops)
    bad_hops[(-1, 0, 0)] = bad
    @test occursin(
        "Hermiticity check failed",
        error_message(() -> PairChecks.check_hr_pair_symmetry(bad_hops; atol=1e-12)),
    )
end
