# Gradus.jl Documentation

A pure Julia geodesic integration system built on [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) using automatic differentiation (AD) and computer algebra systems (CAS) to efficiently compute the geodesic equation. This package requires only a specification of the non-zero metric components in order to solve the 2nd order geodesic system. Alternatively, an implementation of the four velocity components may be specified to integrate a regular 1st order system.

The motivation behind this package began with an interest in studying reverberation lags around accreting black holes, however the scope has since expanded to facilitate the exploration of generic metrics through time-like, space-like, and null geodesics. 

Our aim is to make testing modified Kerr metrics and alternative gravity theories _fast_.

```@raw html
<p align="center" pa="0" ma="0">
<img width="30%" src="assets/uob-astro-grey.png">
</p>
```

Gradus.jl allows for drastically different relativistic simulations to be computed with a composable and reusable API, permitting an end user to simply and expressively calculate physical formulae, create observational signatures, and interface with other popular astrophysics tools. Gradus.jl implements a number of high level abstractions, on the path towards a fully parallelized, high performance numerical relativity ecosystem, scalable from personal computers to super computers.

## Setup

Requires Julia >v1.6. First, install [Buckets.jl](https://github.com/fjebaker/Buckets.jl) and then add Gradus.jl:

```
import Pkg
Pkg.add(url="https://github.com/fjebaker/Buckets.jl")
Plg.add(url="https://github.com/astro-group-bristol/Gradus.jl")
```

## About

Gradus.jl is a research tool for calculating geodesic paths in arbitrary space-times. It is currently work-in-progress, and breaking changes are frequent, as the interface is redesigned to match changing use-cases.

It is part of a larger developing eco-system of *Strong Gravity Codes*, created by members of the University of Bristol Astrophysics Group

- Fergus Baker (PhD Student)
- Dr. Andrew Young (Associate Professor)

For more University of Bristol Astrophysics Group codes, see [our GitHub organisation](https://github.com/astro-group-bristol).