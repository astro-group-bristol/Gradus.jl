"""
struct PowerLawSpectrum{T} <: AbstractCoronalSpectrum
    Γ::T
end

Spectrum for the power law simply returns the photon index `Γ`. `I` stores the 'coronal_spectrum' which is a function 
of `spectrum` and `g` since it is defined as `g^Γ`.

"""
struct PowerLawSpectrum{T} <: AbstractCoronalSpectrum
    Γ::T
end

"""
    function coronal_spectrum(spectrum::PowerLawSpectrum, g)

Function takes the PowerLawSpectrum type spectrum (photon index `Γ`) and `g` as arguments and returns the coronal 
spectrum as `g` to the power of `Γ`. Since the number of photons travelling along any given ray must be conserved 
for different values of the photon index, it serves as a method to accurately account for
change in the number of photons in each bin of the emissivity profile, with the default value of `Γ` = 2 [as
per Gonzalez et. al (2017)]
"""
function coronal_spectrum(spectrum::PowerLawSpectrum, g)
    g^(-spectrum.Γ)
end


export PowerLawSpectrum
