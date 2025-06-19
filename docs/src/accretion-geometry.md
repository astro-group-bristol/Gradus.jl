# Accretion geometry

```@meta
CurrentModule = Gradus
```

Gradus.jl supports the ability to implement custom accretion geometry, or even load in mesh files in any standard format using [MeshIO.jl](https://github.com/JuliaIO/MeshIO.jl). Geometry may be standard spherically symmetric accretion discs, or any other custom type.

!!! note
    Currently geometry is optically thick _always_. Radiative transfer will be added soon.

```@docs
AbstractAccretionGeometry
Gradus.in_nearby_region
Gradus.has_intersect
```

## Adding new accretion geometry

## Accretion discs

```@docs
AbstractAccretionDisc
distance_to_disc
AbstractThickAccretionDisc
cross_section
```

### Available accretion geometry

```@docs
ThinDisc
ThickDisc
ShakuraSunyaev
```

## Meshes

```@docs
MeshAccretionGeometry
```