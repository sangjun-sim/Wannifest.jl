module FluxCore

include(joinpath(@__DIR__, "..", "shared", "CrystalCell.jl"))
include(joinpath(@__DIR__, "..", "shared", "PoscarIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "LatticeIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "SpinLayout.jl"))
include(joinpath(@__DIR__, "..", "shared", "Hermiticity.jl"))
include(joinpath(@__DIR__, "..", "shared", "PairChecks.jl"))
include(joinpath(@__DIR__, "..", "shared", "WannierTypes.jl"))
include(joinpath(@__DIR__, "..", "shared", "HrHermiticity.jl"))
include(joinpath(@__DIR__, "..", "shared", "HrFormat.jl"))
include(joinpath(@__DIR__, "..", "shared", "WannierHrIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "InputParsing.jl"))
include(joinpath(@__DIR__, "..", "shared", "Win90OrbitalTokens.jl"))
include(joinpath(@__DIR__, "..", "shared", "Win90ProjectionIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "OrbitalProjection.jl"))
include(joinpath(@__DIR__, "..", "shared", "Win90Basis.jl"))

include(joinpath(@__DIR__, "Model.jl"))
include(joinpath(@__DIR__, "InputIO.jl"))
include(joinpath(@__DIR__, "BasisSource.jl"))
include(joinpath(@__DIR__, "BondGeometry.jl"))
include(joinpath(@__DIR__, "FluxCompiler.jl"))
include(joinpath(@__DIR__, "Diagnostics.jl"))
include(joinpath(@__DIR__, "PlotFlux3D.jl"))
include(joinpath(@__DIR__, "Service.jl"))

end
