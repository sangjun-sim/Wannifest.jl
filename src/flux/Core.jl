module FluxCore

include(joinpath(@__DIR__, "..", "shared", "geometry", "CrystalCell.jl"))
include(joinpath(@__DIR__, "..", "shared", "io", "PoscarIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "basis", "Win90Keywords.jl"))
include(joinpath(@__DIR__, "..", "shared", "io", "LatticeIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "basis", "SpinLayout.jl"))
include(joinpath(@__DIR__, "..", "shared", "validation", "Hermiticity.jl"))
include(joinpath(@__DIR__, "..", "shared", "validation", "PairChecks.jl"))
include(joinpath(@__DIR__, "..", "shared", "wannier", "WannierTypes.jl"))
include(joinpath(@__DIR__, "..", "shared", "validation", "HrHermiticity.jl"))
include(joinpath(@__DIR__, "..", "shared", "io", "HrFormat.jl"))
include(joinpath(@__DIR__, "..", "shared", "io", "WannierHrIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "utils", "InputParsing.jl"))
include(joinpath(@__DIR__, "..", "shared", "basis", "Win90OrbitalTokens.jl"))
include(joinpath(@__DIR__, "..", "shared", "basis", "Win90ProjectionIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "operators", "OrbitalProjection.jl"))
include(joinpath(@__DIR__, "..", "shared", "basis", "Win90Basis.jl"))

include(joinpath(@__DIR__, "Model.jl"))
include(joinpath(@__DIR__, "InputIO.jl"))
include(joinpath(@__DIR__, "BasisSource.jl"))
include(joinpath(@__DIR__, "BondGeometry.jl"))
include(joinpath(@__DIR__, "FluxCompiler.jl"))
include(joinpath(@__DIR__, "Diagnostics.jl"))
include(joinpath(@__DIR__, "PlotFlux3D.jl"))
include(joinpath(@__DIR__, "Service.jl"))

end
