# Defining a new metric

```@meta
CurrentModule = Gradus
```

Gradus.jl is able to integrate any 3+1 dimensional metric, however provides a few derivative abstract types to implement to ensure the most efficient code is executed for a given metric. Currently, Gradus differentiates the following specializations:


## First-Order

```@docs
FirstOrderMethods.AbstractFirstOrderMetricParams
FirstOrderMethods.four_velocity
FirstOrderMethods.calc_lq
FirstOrderMethods.Vr
FirstOrderMethods.Vθ
```

## Second-Order

```@docs
GeodesicTracer.AbstractAutoDiffMetricParams
GeodesicTracer.metric_components
GeodesicTracer.AbstractAutoDiffStaticAxisSymmetricParams
```