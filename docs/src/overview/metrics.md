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

## First-Order
```@docs
BoyerLindquistFO
```

## Second-Order
```@docs
BoyerLindquistAD
JohannsenAD
KerrRefractiveAD
DilatonAxionAD
MorrisThorneAD
```