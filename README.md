<p align="center">
  <img width="30%" alt="BRImage" src="docs/src/assets/logo.png">
</p>

# Gradus.jl

<a href="https://ascl.net/2503.035">
    <img src="https://img.shields.io/badge/ascl-2503.035-blue.svg?colorB=262255" alt="ascl:2503.035" />
</a>
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

Spacetime generic, general relativistic ray-tracing (GRRT) in Julia. See our [documentation](https://astro-group-bristol.github.io/Gradus.jl/dev/) for more details and usage.

<p align="center"> <i> This package is in development and subject to rapid breaking changes, with documentation updates lagging behind features.</i> </p>

## About

Gradus.jl is a suite of tools related to tracing geodesics and calculating observational signatures. Gradus.jl requires only a specification of the non-zero metric components of a chosen spacetime in order to solve the [geodesic equation](https://en.wikipedia.org/wiki/Solving_the_geodesic_equations) and compute a wide variety of trajectories and orbits. Various algorithms for calculating physical quantities are implemented generically, so they may be used with different classes of spacetime with minimal implementation.

Currently, Gradus.jl can be used for any static, axis-symmetric spacetime to calculate:

- geodesic orbits and special radii (event horizon shapes, ISCO radii, etc.)
- null / time / space like trajectories including for charged particles
- black hole shadows
- redshift images
- Cunningham transfer functions
- line profiles and spectra
- reverberation transfer functions
- time-lags from different coronal models
- emissivity profiles on the accretion disc
- covariant radiative transfer
- various toy accretion models (thin disc, $\alpha$-discs, rotationally-supported polish doughnut, etc)
- non-symmetric disc geometries
- mesh file geometry

The library is written to make adding new features as effortless as possible. See [Examples](https://astro-group-bristol.github.io/Gradus.jl/dev/examples/) for more. Many new features are currently being developed as our research advances.

Gradus.jl uses [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) and [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) as the backend for integrating and solving the geodesic equation for arbitrary metrics, and vendors the DifferentialEquations.jl solver and callback system, making Gradus.jl easy to extend for new problems. Gradus.jl currently supports multi-CPU integration and analysis, with GPU support on the horizon.

## Usage

We assume you already have Julia >1.6.

All non-General dependencies for Gradus.jl are in the [AstroRegistry](https://github.com/astro-group-bristol/AstroRegistry) which can be added to Julia with:

```julia
julia>] registry add https://github.com/astro-group-bristol/AstroRegistry
```

NB: the [Julia General Registry](https://github.com/JuliaRegistries/General) is required to install Gradus, however this should be configured by default with any Julia installation.

Gradus.jl can then be fetched easily:
```julia
julia>] add Gradus
julia> using Gradus
```

See [GettingStarted](https://astro-group-bristol.github.io/Gradus.jl/dev/getting-started/) for setting up your first traces.

## See also

- [The Julia programming language](https://github.com/JuliaLang/Julia)
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
- [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)

## Citing

If you use Gradus.jl in your work, please cite [Baker & Young (2025)](./CITATION.bib).

There are various other ways to reference Gradus.jl:
- [ascl:2503.035](https://ascl.net/2503.035)
- [10.5281/zenodo.6471796](https://zenodo.org/records/6471796)

<hr>

<p align="center"> Astrophysics Group Bristol </p>
