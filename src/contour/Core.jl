module ContourCore

include(joinpath(@__DIR__, "..", "shared", "CrystalCell.jl"))
include(joinpath(@__DIR__, "..", "shared", "SpinLayout.jl"))
include(joinpath(@__DIR__, "..", "shared", "Hermiticity.jl"))
include(joinpath(@__DIR__, "..", "shared", "PairChecks.jl"))
include(joinpath(@__DIR__, "..", "shared", "WannierTypes.jl"))
include(joinpath(@__DIR__, "..", "shared", "WannierHrIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "WannierWsvecIO.jl"))
include(joinpath(@__DIR__, "..", "shared", "WannierWsvecGenerate.jl"))
include(joinpath(@__DIR__, "..", "shared", "WannierKspace.jl"))
include(joinpath(@__DIR__, "..", "shared", "WannierEigensystem.jl"))
include(joinpath(@__DIR__, "..", "shared", "InputParsing.jl"))
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
