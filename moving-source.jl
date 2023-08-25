using Gradus

"""
struct BeamedPointSource{T} <: Gradus.AbstractCoronaModel{T}
    r::T
    β::T
end

Point source corona moving away from the black hole with specified starting height above the disc `r` and source speed `β`.
Point sources are the most plausible example of source that would support beaming (ref: Gonzalez et al 2017).

"""
struct BeamedPointSource{T} <: Gradus.AbstractCoronaModel{T}
    r::T
    β::T
end

# these are specific to BeamedPointSource, so we'll scope them in a module
# so that they don't pollute the namespace of Gradus

module __BeamedPointSource

import ..Gradus: AbstractMetric, metric_components, SVector

drdt(g, β) = β * √(-g[1] / g[2])

drdt(m::AbstractMetric, x, β) = drdt(metric_components(m, SVector(x[2], x[3])), β)
end # module

function Gradus.sample_position_velocity(m::AbstractMetric, model::BeamedPointSource)
    x = SVector{4}(0, model.r, 1e-4, 0)
    g = metric_components(m, SVector(x[2], x[3]))
    v̄ = SVector(1, __BeamedPointSource.drdt(g, model.β), 0, 0)
    v = Gradus.constrain_normalize(m, x, v̄; μ = 1)
    x, v
end

# since we have axis-symmetry, exploit this method for much faster computation
Gradus.emissivity_profile(
    ::Nothing,
    m::AbstractMetric,
    d::AbstractAccretionDisc,
    model::BeamedPointSource,
    spec::Gradus.AbstractCoronalSpectrum;
    kwargs...,
) = Gradus._point_source_symmetric_emissivity_profile(m, d, model, spec; kwargs...)

