## Getting started

```@meta
CurrentModule = Gradus
```

Unlike conventional [ray-tracing](https://en.wikipedia.org/wiki/Ray_tracing_(graphics)), ray-tracing in general relativity (GR) has the added complication that the trajectory of light is altered by the curvature of space. In particular, the spacetime around compact singularities, such as black holes, may be significantly curved in weird and wonderful ways, depending on the nature of the object being studied. When attempting to visualise or calculate observational signatures related to these objects, is important to account for so-called _GR effects_; these effects not only alter how things _look_, but also the _energetics_ of the system itself.

This short getting-started guide should hopefully illustrate some of the key ideas that need to be kept in mind when studying accretion processes and simulating spectra of black holes.

## The trajectory of photons

```@docs
geodesic_equation
```

## Where to go from here?

The documentation is a rich resource for information related to using Gradus, and tailoring the toolkit for your needs. Take a look at [Examples](@ref) for a number of quick recipes, or try [Implementing a new metric](@ref) and study a different spacetime.