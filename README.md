<p align="center">
  <img width="30%" alt="BRImage" src="docs/src/assets/logo.png">
</p>

# Gradus.jl

<a href="https://codecov.io/gh/astro-group-bristol/Gradus.jl">
    <img src="https://codecov.io/gh/astro-group-bristol/Gradus.jl/branch/main/graph/badge.svg?token=A91E22KZR5"/>
</a>
<a href="https://github.com/astro-group-bristol/Gradus.jl/actions/workflows/test.yml">
    <img src="https://github.com/astro-group-bristol/Gradus.jl/actions/workflows/test.yml/badge.svg"/>
</a>
<a href="https://github.com/astro-group-bristol/Gradus.jl/actions/workflows/docs.yml">
    <img alt="Docs" src="https://github.com/astro-group-bristol/Gradus.jl/actions/workflows/docs.yml/badge.svg"/>
</a>
<a href="https://doi.org/10.5281/zenodo.6471796">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.6471796.svg" alt="DOI">
</a> 
<a href="https://astro-group-bristol.github.io/Gradus.jl/dev/">
    <img alt="Docs" src="https://img.shields.io/badge/docs-dev-blue.svg"/>
</a>
<a href="https://github.com/JuliaTesting/Aqua.jl">
    <img alt="Docs" src="https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg"/>
</a>

Spacetime generic, general relativistic ray tracing (GRRT) toolkit in Julia.

<p align="center"> <i> This package is in development and subject to rapid breaking changes, with documentation updates lagging behind features.</i> </p>

## About

Gradus.jl is a suite of tools related to tracing geodesics, requiring only a specification of the non-zero metric components of a chosen spacetime. Algorithms for calculating various physical quantities are implemented generically for different classes of spacetime, and can be used without extension on new metrics.

Currently, Gradus.jl can be used for any static, axis-symmetric spacetime to calculate:

- geodesic orbits and special radii (event horizon shapes, ISCO radii)
- black hole shadows
- redshift images
- Cunningham transfer functions
- line profiles and spectra
- reverberation transfer functions
- time-lags from different coronal models

See [Examples](https://astro-group-bristol.github.io/Gradus.jl/dev/examples/examples/).

All of the above may be calculated for different accretion disc geometries, including geometrically thin and thick discs, with radiative transport for optically thin material planned. Gradus.jl can also use mesh files to specify non-symmetric geometry.

Gradus.jl uses [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) and [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) as the backend for integrating and solving the geodesic equation for arbitrary metrics, and vendors the DifferentialEquations.jl solver and callback system, making Gradus.jl easy to extend for new problems. Gradus.jl currently supports multi-CPU integration and analysis, with GPU support on the horizon.


## Usage

We assume you already have Julia >1.6.

- Recommended:

All non-General dependencies for Gradus.jl are in the [AstroRegistry](https://github.com/astro-group-bristol/AstroRegistry) which can be added to Julia with:

```julia
julia>] add registry https://github.com/astro-group-bristol/AstroRegistry
```

Gradus.jl can then be fetched easily:
```julia
julia>] add Gradus
julia> using Gradus
```

- Alternate:

Gradus.jl depends on [Buckets.jl](https://github.com/fjebaker/Buckets.jl), which can be installed directly from GitHub:

```julia
import Pkg;
Pkg.add(url = "https://github.com/fjebaker/Buckets.jl")
```
After this, Gradus.jl can be installed, also directly from GitHub:
```julia
Pkg.add(url = "https://github.com/astro-group-bristol/Gradus.jl")

using Gradus
```

## See also 

- [The Julia programming language](https://github.com/JuliaLang/Julia)
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)

<hr>

<p align="center"> Astrophysics Group Bristol </p>