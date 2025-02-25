abstract type AbstractRootAlgorithm end
struct RootsAlg <: AbstractRootAlgorithm end
Base.@kwdef struct NonLinearAlg{A} <: AbstractRootAlgorithm
    alg::A = SimpleNonlinearSolve.SimpleBroyden()
end

DEFAULT_ROOT_SOLVER() = NonLinearAlg()

"""
    root_solve(f_objective, initial_value, args)

Wrapper to different root solving backends to make root solve fast and efficient
"""
function root_solve(
    f_objective,
    initial_value::T,
    args;
    kwargs...,
) where {T<:Union{<:Number,<:SVector{1}}}
    root_solve(DEFAULT_ROOT_SOLVER(), f_objective, initial_value, args; kwargs...)
end
function root_solve(
    ::RootsAlg,
    f_objective,
    initial_value::T,
    args;
    abstol = 1e-9,
    kwargs...,
) where {T<:Union{<:Number,<:SVector{1}}}
    x0 = Roots.find_zero(
        r -> f_objective(r, args),
        initial_value,
        Roots.Order0();
        atol = abstol,
    )
    resid = f_objective(x0, args)
    x0, resid
end
function root_solve(
    alg::NonLinearAlg,
    f_objective,
    initial_value::T,
    args;
    abstol = 1e-9,
    kwargs...,
) where {T<:Union{<:Number,<:SVector{1}}}
    x0, f = if T <: Number
        function _obj_wrapper(x::SVector, p)
            @inbounds SVector{1,eltype(x)}(f_objective(x[1], p))
        end
        SVector{1}(initial_value), _obj_wrapper
    else
        initial_value, f_objective
    end
    prob = SimpleNonlinearSolve.NonlinearProblem{false}(f, x0, args)
    sol = solve(prob, alg.alg, abstol = abstol, reltol = abstol, maxiters = 500, kwargs...)
    sol.u[1], sol.resid[1]
end


@inline function _symmetric_matrix(comps::AbstractVector{T})::SMatrix{4,4,T} where {T}
    @SMatrix [
        comps[1] 0 0 comps[5]
        0 comps[2] 0 0
        0 0 comps[3] 0
        comps[5] 0 0 comps[4]
    ]
end

@inline function _threaded_map(f, itr::T) where {T}
    N = length(itr)
    items = !(T <: AbstractArray) ? collect(itr) : itr
    output = Vector{Core.Compiler.return_type(f, Tuple{eltype(items)})}(undef, N)
    Threads.@threads for i = 1:N
        @inbounds output[i] = f(items[i])
    end
    output
end

function spherical_to_cartesian(v)
    x = v[1] * cos(v[3]) * sin(v[2])
    y = v[1] * sin(v[3]) * sin(v[2])
    z = v[1] * cos(v[2])
    SVector{3}(x, y, z)
end

# specialisation for four-vector
spherical_to_cartesian(v::SVector{4}) = spherical_to_cartesian(@views(v[2:end]))

function cartesian_squared_distance(::AbstractMetric, x1, x2)
    # all metrics are currently in boyer lindquist coords
    y1 = spherical_to_cartesian(x1)
    y2 = spherical_to_cartesian(x2)
    diff = @. (y2 - y1)^2
    sum(diff)
end

cartesian_distance(m::AbstractMetric, x1, x2) = √(cartesian_squared_distance(m, x1, x2))

# both energy and angular momentum
# assume time only coupled to radial coordinate
# need to think of a nice way to keep this efficient
# whilst allowing metrics with other couplings

"""
    E(m::AbstractMatrix{T}, v)
    E(m::AbstractMetric{T}, u, v)

Compute the energy for a numerically evaluated metric, and some velocity four vector `v`,
```math
E = - p_t = - g_{t\\nu} p^\\nu.
```

For null geodesics, the velocity is the momentum ``v^\\nu = p^\\nu``. For massive geodesics,
the mass ``\\mu`` needs to be known to compute ``\\mu v^\\nu = p^\\nu``.
"""
function E(metric::AbstractMatrix{T}, v) where {T}
    T(@inbounds -(metric[1, 1] * v[1] + metric[1, 4] * v[4]))
end
E(m::AbstractMetric, u, v) = E(metric(m, u), v)
E(m::AbstractMetric, gp::AbstractGeodesicPoint) = E(m, gp.x, gp.v)


"""
    Lz(m::AbstractMatrix{T}, v)
    Lz(m::AbstractMetric{T}, u, v)

Compute the angular momentum for a numerically evaluated metric, and some velocity four vector `v`.
```math
L_z = p_\\phi = - g_{\\phi\\nu} p^\\nu.
```
"""
function Lz(metric::AbstractMatrix{T}, v) where {T}
    T(@inbounds metric[4, 4] * v[4] + metric[1, 4] * v[1])
end
Lz(m::AbstractMetric, u, v) = Lz(metric(m, u), v)
Lz(m::AbstractMetric, gp::AbstractGeodesicPoint) = Lz(m, gp.x, gp.v)

@inline function _optional_abs(value::T, signed::Bool)::T where {T}
    if signed
        value
    else
        abs(value)
    end
end

_equatorial_project(r, θ; signed = false) = r * _optional_abs(sin(θ), signed)
_equatorial_project(x::SVector; signed = false) =
    _equatorial_project(x[2], x[3], signed = signed)

_spinaxis_project(r, θ; signed = false) = r * _optional_abs(cos(θ), signed)
_spinaxis_project(x::SVector; signed = false) =
    _spinaxis_project(x[2], x[3], signed = signed)

_rotate_about_spinaxis(n::SVector{3}, ϕ) = SVector(n[1] * cos(ϕ), n[1] * sin(ϕ), n[3])

_zero_if_nan(x::T) where {T} = isnan(x) ? zero(T) : x

@inline function _smooth_interpolate(
    x::T,
    x₀;
    δx = T(2.5),
    smoothing_offset = T(1e4),
) where {T}
    if x ≤ x₀
        one(T)
    elseif x₀ ≤ x ≤ x₀ + δx
        t = (x - x₀) / δx
        atan(smoothing_offset * t) * 2 / π
    else
        zero(T)
    end
end

@generated function _unroll_for(f, ::Val{N}, items) where {N}
    exprs = [:(f(items[$i])) for i = 1:N]
    quote
        $(exprs...)
    end
end

function quadrature_domain(a, b, X::AbstractVector)
    q = (b - a) / 2
    ((xi + 1) * q + a for xi in X)
end

function quadrature_sum(values::AbstractVector, weights::AbstractVector)
    sum(i * j for (i, j) in zip(values, weighs))
end

export cartesian_squared_distance, cartesian_distance, spherical_to_cartesian
