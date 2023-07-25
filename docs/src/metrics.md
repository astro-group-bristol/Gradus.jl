# Catalogue of metrics

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
    To implement your own custom metrics, please see [Implementing a new metric](@ref). If you have a complex metric, please open an issue requesting for it to be added.

## Currently available metrics

```@autodocs
Modules = [Gradus]
Filter = t -> typeof(t) === UnionAll && t <: Gradus.AbstractMetric
```
