using Test
using LinearAlgebra
using Wannifest

const WANNIFEST_DIR = normpath(joinpath(@__DIR__, ".."))
const EXAMPLES_DIR = joinpath(WANNIFEST_DIR, "examples")
_example_dir(candidates...) = begin
    for candidate in candidates
        isdir(candidate) && return candidate
    end
    return first(candidates)
end
const GRAPHENE_EXAMPLE_DIR = _example_dir(
    joinpath(EXAMPLES_DIR, "graphene"),
    joinpath(EXAMPLES_DIR, "flake", "graphene"),
)
const RUCL3_EXAMPLE_DIR = _example_dir(
    joinpath(EXAMPLES_DIR, "rucl3"),
    joinpath(EXAMPLES_DIR, "band", "rucl3"),
)
const MOS2_EXAMPLE_DIR = joinpath(EXAMPLES_DIR, "1tp-mos2")
const SAMPLE_SPECTRA_DIR = joinpath(WANNIFEST_DIR, "test", "sample", "spectra")

const BandCore = Wannifest.BandCLI.BandCore
const ContourCore = Wannifest.ContourCLI.ContourCore
const SuperhamCore = Wannifest.SuperhamCLI.SuperhamCore
const Bands = BandCore.Bands
const CenterIO = SuperhamCore.CenterIO
const Dos = BandCore.Dos
const Execution = BandCore.Execution
const InputIO = Wannifest.BandCLI.InputIO
const KPath = BandCore.KPath
const LatticeIO = BandCore.LatticeIO
const LocalAxisRotation = BandCore.LocalAxisRotation
const Model = BandCore.Model
const AtomicOam = BandCore.AtomicOam
const AtomicSpin = BandCore.AtomicSpin
const OrbitalProjection = BandCore.OrbitalProjection
const Output = BandCore.Output
const ObservablesModel = BandCore.ObservablesModel
const ObservablesOutput = BandCore.ObservablesOutput
const PlotService = BandCore.PlotService
const Projection = BandCore.Projection
const Service = Wannifest.BandCLI.Service
const SpinLayout = BandCore.SpinLayout
const Hermiticity = BandCore.Hermiticity
const PairChecks = BandCore.PairChecks
const WannierEigensystem = BandCore.WannierEigensystem
const WannierHrIO = BandCore.WannierHrIO
const WannierKspace = BandCore.WannierKspace
const WannierTypes = BandCore.WannierTypes
const WannierWsvecGenerate = BandCore.WannierWsvecGenerate
const Win90Basis = BandCore.Win90Basis
const WannierWsvecIO = BandCore.WannierWsvecIO
const ContourInputIO = Wannifest.ContourCLI.InputIO
const ContourMesh = ContourCore.Mesh
const ContourModel = ContourCore.Model
const ContourOutput = ContourCore.Output
const ContourPlotService = ContourCore.PlotService
const ContourService = Wannifest.ContourCLI.Service
const ContourSurface = ContourCore.Surface

function error_message(f)
    try
        f()
        return ""
    catch err
        return sprint(showerror, err)
    end
end
