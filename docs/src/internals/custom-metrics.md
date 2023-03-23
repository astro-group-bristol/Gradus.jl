# Implementing a new metric

```@meta
CurrentModule = Gradus
```

Gradus.jl is able to integrate any 3+1 dimensional metric. A new metric may be defined by implementing one of the abstract types with a concrete type, and defining a number of methods. Depending on what you want to be able to do with a metric, different functions need to be implemented. 

Gradus also provides a few derivative abstract types to implement to ensure the most efficient code is executed for a given metric (see [Metric parameter types](@ref) below).

## Example: Schwarzschild

As a minimal example, here is how the Schwarzschild metric may be implemented. First, we must define what the _metric parameters_ for this metric are. These are effectively constants of the spacetime, representing physical quantities that appear in the metric expression. For the Schwarzschild metric, this is only the black hole mass ``M``, but e.g. the Kerr metric also has the black hole spin ``a``.

We can choose the integration strategy by sub-typing an abstract type representing different classes of spacetimes. For the Schwarzschild metric, we will use the static, axis-symmetric class, with the automatic differentiation (AD) backend. With AD, we only need to specify the non-zero components of the metric as Julia functions, and the rest is done for us.

For ease, we choose the [Eddington-Finkelstein coordinates](https://en.wikipedia.org/wiki/Eddington%E2%80%93Finkelstein_coordinates) of the Schwarzschild solution, which may be written

```math
\text{d}s^2 =
    - \left( 1 - \frac{2 M}{r} \right) \text{d}t^2
    + \left( 1 - \frac{2 M}{r} \right)^{-1} \text{d}r^2
    + r^2 \text{d}\theta^2
    + r^2 \sin^2(\theta) \text{d}\phi^2.
```

Here is a possible implementation for Gradus.jl:
```julia
using Gradus

@with_kw struct EddingtonFinkelsteinAD{T} <: AbstractStaticAxisSymmetric{T}
    M = 1.0
end

function GradusBase.metric_components(m::EddingtonFinkelsteinAD{T}, rθ) where {T}
    (r, θ) = rθ
    M = m.M

    tt = -(1 - (2M / r))
    rr = -inv(tt)
    θθ = r^2
    ϕϕ = r^2 * sin(θ)^2

    (tt, rr, θθ, ϕϕ, T(0.0))
end

GradusBase.inner_radius(m::EddingtonFinkelsteinAD) = 2 * m.M
```
A few notes:
- We use `@with_kw` from [Parameters.jl](https://github.com/mauro3/Parameters.jl) to define various utility constructors for us.
- [`GradusBase.metric_components`](@ref) must return five elements for [`AbstractStaticAxisSymmetric`](@ref), where the last element is the off-axis ``g_{t \phi}`` matrix element, which in this case is always 0.
- The [`GradusBase.inner_radius`](@ref) function defines the inner-radius of integration chart. This defines where the integration should terminate to avoid running indefinitely, and is, in this case, set to the event-horizon of our metric.

That's all we need! This metric is now ready to be traced in the usual way.

!!! note
    For more examples of how to implement different metrics, click on the "source" button of a metric in [Implemented Metrics](@ref). Alternatively, view the source code directly [here](https://github.com/astro-group-bristol/Gradus.jl/tree/main/src/metrics).


## Metric parameter types

The following types may be implemented to add new metrics. Each type has different requirements for its interface.

### First-Order

```@docs
AbstractFirstOrderMetric
Gradus.four_velocity
Gradus.calc_lq
Gradus.Vr
Gradus.Vθ
```

### Second-Order

```@docs
AbstractMetric2ndOrder
metric_components
AbstractStaticAxisSymmetric
```