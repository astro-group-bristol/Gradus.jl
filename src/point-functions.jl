"""
    struct PointFunction <: AbstractPointFunction
    PointFunction(func)

$(FIELDS)

Point functions are functions that are used to calculate physical parameters from geodesic
integrations, and to compose more complex models. A number of default and utility `PointFunction`
are defined in [`Gradus.ConstPointFunctions`](@ref).

Principally, point functions return a single value per geodesic, and are used to fill rendered
images with values, e.g. colouring redshift.

Point functions may be instantiated by wrapping a function with the following signature
```julia
function func(m::AbstractMetric{T}, gp::AbstractGeodesicPoint, max_time::T; kwargs...)::T where {T}
    # ...
end

pf = PointFunction(func)
```
- The `AbstractMetric` argument may be used to dispatch for different metrics.
- `gp` is an [`GradusBase.AbstractGeodesicPoint`](@ref) corresponding to a given geodesic.
- The `max_time` parameter is the maximum integration time used to integrate the geodesics. This may be useful when trying to determine whether a geodesic terminated early or not.

They may be invoked by invoking the instance
```julia
result = pf(m, gp, max_time)
```

!!! note
    As of version 0.1.0, the `kwargs` parameter is reserved only for passing optional results
    when chaining multiple point functions (see below). This is subject to revision and breaking changes
    in future versions.

Multiple [`AbstractPointFunction`](@ref) may be chained together using the `∘` operator, and
are evaluated from right to left
```julia
pf3 = pf2 ∘ pf1
```

This may be useful for constructing filters using [`FilterPointFunction`](@ref). When used with
two `PointFunction` objects, the output of the previous `PointFunction` is passed to the next
via the `value` keyword argument.
"""
struct PointFunction{F} <: AbstractPointFunction
    "Wrapped function."
    f::F
end

"""
    struct FilterPointFunction <: AbstractPointFunction
    FilterPointFunction(func, default_value)

$(FIELDS)

Point functions used to filter geodesics. They may be constructed with
```julia
function func(m::AbstractMetric{T}, gp::AbstractGeodesicPoint, max_time::T; kwargs...)::Bool where {T}
    # ... return Bool
end

fpf = FilterPointFunction(func, NaN64)
```

The second argument to the constructor is the default value, given to the pixel if the boolean
condition of `func` is `false`.

# Example

A filter for geodesics within a certain radius, used to only calculate redshift within
10 ``\\r_\\text{g}``
```julia
func(m, gp, max_time) = gp.u[2] < 10.0
pf = ConstPointFunctions.redshift(m, u) ∘ FilterPointFunction(func, NaN64)
```
"""
struct FilterPointFunction{F,T} <: AbstractPointFunction
    "Wrapped function."
    f::F
    "Default return value if condition is false."
    default::T
end

# todo: this should type specialize better
FilterPointFunction(f) = FilterPointFunction(f, NaN)

@inline function (pf::AbstractPointFunction)(
    m::AbstractMetric{T},
    gp::GradusBase.AbstractGeodesicPoint{T},
    max_time;
    kwargs...,
)::T where {T}
    convert(T, pf.f(m, gp, max_time; kwargs...))
end

@inline function apply(pf::AbstractPointFunction, rc::SolutionRenderCache; kwargs...)
    _threaded_map(
        sol -> pf.f(rc.m, process_solution(m, sol), rc.max_time; kwargs...),
        rc.geodesics,
    )
end

@inline function apply(pf::AbstractPointFunction, rc::EndpointRenderCache; kwargs...)
    _threaded_map(gp -> pf.f(rc.m, gp, rc.max_time; kwargs...), rc.points)
end

@inline function Base.:∘(pf1::AbstractPointFunction, pf2::AbstractPointFunction)
    let f1 = pf1.f, f2 = pf2.f
        PointFunction(
            (m, gp, max_time; kwargs...) ->
                f1(m, gp, max_time; value = f2(m, gp, max_time; kwargs...)),
        )
    end
end

@inline function Base.:∘(pf1::AbstractPointFunction, pf2::FilterPointFunction)
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


# some utility methods
FilterStatusCode(code::StatusCodes.T, default = NaN) =
    FilterPointFunction((m, gp, λ) -> gp.status == code, default)

export apply, PointFunction, FilterPointFunction, FilterStatusCode
