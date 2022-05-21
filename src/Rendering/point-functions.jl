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
    ThreadsX.map(
        sol -> pf.f(rc.m, getgeodesicpoint(m, sol), rc.max_time; kwargs...),
        rc.geodesics,
    )
end

@inline function apply(
    pf::AbstractPointFunction,
    rc::EndpointRenderCache{M,T,P};
    kwargs...,
) where {M,T,P}
    ThreadsX.map(gp -> pf.f(rc.m, gp, rc.max_time; kwargs...), rc.points)
end

@inline function Base.:∘(pf1::AbstractPointFunction, pf2::AbstractPointFunction)
    let f1 = pf1.f, f2 = pf2.f
        PointFunction((m, gp, max_time; kwargs...) -> f1(f2(m, gp, max_time; kwargs...)))
    end
end

@inline function Base.:∘(pf1::AbstractPointFunction, pf2::FilterPointFunction{F}) where {F}
    let f1 = pf1.f, f2 = pf2.f
        PointFunction(
            (m, gp, max_time; kwargs...) -> begin
                pass_on = f2(m, gp, max_time; kwargs...)
                if pass_on
                    f1(m, gp, max_time; kwargs...)
                else
                    pf2.default
                end
            end,
        )
    end
end
