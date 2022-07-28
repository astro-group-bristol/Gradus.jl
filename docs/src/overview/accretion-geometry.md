# Accretion geometry

```@meta
CurrentModule = Gradus.AccretionGeometry
```

Gradus.jl supports the ability to implement custom accretion geometry, or even load in mesh files in any standard format using [MeshIO.jl](https://github.com/JuliaIO/MeshIO.jl). Geometry may be standard spherically symmetric accretion discs, or any other custom type.

!!! note
    Currently geometry is optically thick _always_. Radiative transfer will be added soon.

```@docs
AbstractAccretionGeometry
in_nearby_region
has_intersect
```

## Accretion discs

```@docs
AbstractAccretionDisc
GeometricThinDisc
```

## Meshes

```@docs

```