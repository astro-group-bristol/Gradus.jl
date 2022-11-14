
```@raw html
<p align="center" pa="0" ma="0">
<img width="30%" src="assets/uob-astro-grey.png">
</p>
```

# Gradus.jl Documentation

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


## Setup

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

## About

Gradus.jl is a research tool for calculating geodesic paths in arbitrary space-times. It is currently work-in-progress, and breaking changes are frequent, as the interface is redesigned to match changing use-cases.

It is part of a larger developing eco-system of *Strong Gravity Codes*, created by members of the University of Bristol Astrophysics Group

- Fergus Baker (PhD Student)
- Dr. Andrew Young (Associate Professor)

For more University of Bristol Astrophysics Group codes, see [our GitHub organisation](https://github.com/astro-group-bristol).