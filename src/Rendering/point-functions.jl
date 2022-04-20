abstract type AbstractPointFunction end

struct PointFunction{F} <: AbstractPointFunction
    f::F
end

struct FilterPointFunction{F,T} <: AbstractPointFunction
    f::F
    default::T
end

@inline function (pf::AbstractPointFunction)(
    m::AbstractMetricParams{T},
    gp::GradusBase.AbstractGeodesicPoint{T},
    max_time;
    kwargs...,
)::T where {T}
    convert(T, pf.f(m, gp, max_time; kwargs...))
end

@inline function apply(
    pf::AbstractPointFunction,
    rc::SolutionRenderCache{M,T,G};
    kwargs...,
) where {M,T,G}
    map(sol -> pf.f(rc.m, get_endpoint(m, sol), rc.max_time; kwargs...), rc.geodesics)
end

@inline function apply(
    pf::AbstractPointFunction,
    rc::EndpointRenderCache{M,T,P};
    kwargs...,
) where {M,T,P}
    map(gp -> pf.f(rc.m, gp, rc.max_time; kwargs...), rc.points)
end

@inline function Base.:∘(pf1::AbstractPointFunction, pf2::AbstractPointFunction)
    PointFunction((m, gp, max_time; kwargs...) -> pf1.f(pf2.f(m, gp, max_time; kwargs...)))
end

@inline function Base.:∘(pf1::AbstractPointFunction, pf2::FilterPointFunction{F}) where {F}
    PointFunction(
        (m, gp, max_time; kwargs...) -> begin
            pass_on = pf2.f(m, gp, max_time; kwargs...)
            if pass_on
                pf1.f(m, gp, max_time; kwargs...)
            else
                pf2.default
            end
        end,
    )
end

module ConstPointFunctions
import ..Rendering: PointFunction, FilterPointFunction

const filter_early_term =
    FilterPointFunction((m, gp, max_time; kwargs...) -> gp.t < max_time, NaN)

const filter_intersected =
    FilterPointFunction((m, gp, max_time; kwargs...) -> gp.retcode == :Intersected, NaN)

const affine_time = PointFunction((m, gp, max_time; kwargs...) -> gp.t)

const shadow = affine_time ∘ filter_early_term
end # module
