"""
TODO: DOCSTRING
"""
struct PowerLawSpectrum{T} <: AbstractCoronalSpectrum
    Γ::T
end

function coronal_spectrum(spectrum::PowerLawSpectrum, g)
    g^(-spectrum.Γ)
end


export PowerLawSpectrum