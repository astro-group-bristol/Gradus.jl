# Available metrics

```@meta
CurrentModule = Gradus
```

Gradus.jl implements a library of metrics ready to use for integrations and rendering.

```@index
Pages = ["metrics.md"]
Modules = [Gradus]
Order = [:type]
```


!!! note
    To implement your own custom metrics, please see [Implementing a new metric](@ref).

## Catalogue of spacetimes

```@autodocs
Modules = [Gradus]
Filter = t -> typeof(t) === UnionAll && t <: Gradus.AbstractMetricParams
```