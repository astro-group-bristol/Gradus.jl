function split_options(::AbstractMetric, opts)
    opts, (;)
end

"""
    local_momentum(r_obs, α, β)

Calculate the local momentum vector ``\\bar{p}_{(\\mu)}`` corresponding to the impact parameters 
``\\alpha`` and ``\\beta`` for an observer located at ``r_\\text{obs}``.

Note, assumes ``\\mu = 0``.
"""
function local_momentum(r_obs, α, β)
    b = β / r_obs
    a = α / r_obs
    prpt = -inv(sqrt(1 + a^2 + b^2))
    pθpt = b * prpt
    pϕpt = a * prpt
    SVector(1, prpt, pθpt, pϕpt)
end

"""
    lnr_momentum_to_global_velocity_transform(m, x)

Returns a function that performs the transformation from LNRF momentum to global velocity.
That is, this function retuns the map

```math
\\bar{p}_{(\\nu)} \\mapsto v^\\mu = g^{\\mu \\sigma} \\mathbf{e}_{\\sigma}^{(\\nu)} \\bar{p}_{(\\nu)}.
```
"""
function lnr_momentum_to_global_velocity_transform(m::AbstractMetric, x)
    g = metric(m, x)
    Tx = hcat(lnrbasis(g)...)
    ginv = inv(g)
    # transform to global frame, then raise indices
    function _transform(p̄)
        ginv * (Tx * p̄)
    end
end

"""
    map_impact_parameters(m::AbstractMetric{T}, x, α, β)

Map impact parameters `α` and `β` to an (unconstrained) velocity vector at some position `x` in the given metric `m`.

Assumes stationary observer in a locally non-rotating frame (LNRF), calculating the momenta

```math
p_{\\mu} = e^{(\\nu)}_{\\phantom{(\\nu)} \\mu} p_{(\\nu)},
```
where ``e^{(\\nu)}_{\\phantom{(\\nu)} \\mu}`` is the LNRF. The metric is then used to find the 
contravariant velocity.

This functions uses [`lnr_momentum_to_global_velocity_transform`](@ref) to calculate the LNRF to
global frame transformation, and determines the local momenta with [`local_momentum`](@ref).

The impact parameters are interpreted as follows:

- if the geodesic were a straight line path, the impact parameter in a given dimension is the distance to the origin from 
the closest point along the geodesic in ``r_\\text{g}``.

For example, for the Schwarzschild metric with ``M = 1``, the impact parameters ``(\\alpha = 2, \\beta = 0)`` would travel 
tangential to the event horizon (``r_\\text{s} = 2M``) if space were flat. 
"""
@inline function map_impact_parameters(m::AbstractMetric, x, α::Number, β::Number)
    _map_impact_parameters(m, x, α, β)
end
function map_impact_parameters(m::AbstractMetric, x, α::AbstractVector, β::AbstractVector)
    _map_impact_parameters(m, x, ((a, b) for (a, b) in zip(α, β))) |> collect
end
function map_impact_parameters(m::AbstractMetric, x, α::AbstractVector, β::Number)
    _map_impact_parameters(m, x, ((a, β) for a in α)) |> collect
end
function map_impact_parameters(m::AbstractMetric, x, α::Number, β::AbstractVector)
    _map_impact_parameters(m, x, ((α, b) for b in β)) |> collect
end

# implementation specifics
function _map_impact_parameters(m::AbstractMetric, x, iterator)
    xfm = lnr_momentum_to_global_velocity_transform(m, x)
    (xfm(local_momentum(x[2], α, β)) for (α, β) in iterator)
end
function _map_impact_parameters(m::AbstractMetric, x, α, β)
    xfm = lnr_momentum_to_global_velocity_transform(m, x)
    xfm(local_momentum(x[2], α, β))
end

function faraday_tensor(m::AbstractMetric, x)
    ST = SVector{4,eltype(x)}
    dA = ForwardDiff.jacobian(t -> electromagnetic_potential(m, t), SVector(x[2], x[3]))
    ∂A = hcat(zeros(ST), dA, zeros(ST))
    g = inv(metric(m, x))
    # raise first index: F^μ_κ
    # @tullio F[μ, κ] := g[μ, σ] * (∂A[σ, κ] - ∂A[κ, σ])

    # faster
    g * (∂A - ∂A')
end

export map_impact_parameters
