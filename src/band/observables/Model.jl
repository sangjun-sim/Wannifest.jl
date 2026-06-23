module ObservablesModel

export OamOrbitalSelection, OamConfig, SamConfig
export disabled_oam_config, disabled_sam_config

struct OamOrbitalSelection
    site::String
    orbital_shell::String
end

struct OamConfig
    enabled::Bool
    selections::Vector{OamOrbitalSelection}
    degeneracy_tol::Float64
    plot_components::Vector{Symbol}
end

disabled_oam_config() = OamConfig(false, OamOrbitalSelection[], 1.0e-4, Symbol[])

struct SamConfig
    enabled::Bool
    degeneracy_tol::Float64
    plot_components::Vector{Symbol}
end

disabled_sam_config() = SamConfig(false, 1.0e-4, Symbol[])

end
