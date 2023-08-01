function _make_interpolation(x, y)
    DataInterpolations.LinearInterpolation(y, x)
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

export cartesian_squared_distance, cartesian_distance, spherical_to_cartesian
