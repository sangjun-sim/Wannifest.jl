module BasisLabelNormalize

export canonical_orbital, canonical_spin

function canonical_orbital(label::AbstractString)
    clean = lowercase(strip(String(label)))
    clean == "dx2y2" && return "dx2-y2"
    return clean
end

function canonical_spin(spin)::String
    clean = lowercase(strip(String(spin)))
    clean == "down" && return "dn"
    clean == "none" && return "unpolarized"
    return clean
end

end
