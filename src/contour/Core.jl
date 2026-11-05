module ContourCore

include(joinpath(@__DIR__, "..", "shared", "geometry", "CrystalCell.jl"))
include(joinpath(@__DIR__, "..", "shared", "basis", "SpinLayout.jl"))
include(joinpath(@__DIR__, "..", "shared", "validation", "Hermiticity.jl"))
include(joinpath(@__DIR__, "..", "shared", "validation", "PairChecks.jl"))
include(joinpath(@__DIR__, "..", "shared", "wannier", "WannierTypes.jl"))
include(joinpath(@__DIR__, "..", "shared", "io", "WannierHrIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "io", "WannierWsvecIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "wannier", "WannierKspace.jl"))
include(joinpath(@__DIR__, "..", "shared", "wannier", "WannierEigensystem.jl"))
include(joinpath(@__DIR__, "..", "shared", "utils", "InputParsing.jl"))
include(joinpath(@__DIR__, "..", "band", "Execution.jl"))
include(joinpath(@__DIR__, "Model.jl"))
include(joinpath(@__DIR__, "InputIO.jl"))
include(joinpath(@__DIR__, "Mesh.jl"))
include(joinpath(@__DIR__, "Surface.jl"))
include(joinpath(@__DIR__, "Output.jl"))
include(joinpath(@__DIR__, "PlotInteractive.jl"))
include(joinpath(@__DIR__, "PlotService.jl"))
include(joinpath(@__DIR__, "Service.jl"))

end
