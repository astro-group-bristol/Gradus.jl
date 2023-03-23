## Getting started

```@meta
CurrentModule = Gradus
```

Unlike conventional [ray-tracing](https://en.wikipedia.org/wiki/Ray_tracing_(graphics)), ray-tracing in general relativity (GR) has the added complication that the trajectory of light is altered by the curvature of space. In particular, the spacetime around compact singularities, such as black holes, may be significantly curved in weird and wonderful ways, depending on the nature of the object being studied. When attempting to visualise or calculate observational signatures related to these objects, is important to account for so-called _GR effects_: these effects not only alter how things _look_, but also the _energetics_ of the system itself.

This short getting-started guide should hopefully tour you through some of the features of Gradus.jl, and how it can be used to study accretion processes, and different spacetimes.

## 1 Defining a spacetime

```@docs
geodesic_equation
```

The workhorse of Gradus.jl is [`tracegeodesics`](@ref). This function is responsible to setting up ordinary differential equation and solving them, either sequentially, in parallel, or on different hardware. To get started, we must minimally choose a spacetime to trace in, an initial position and an initial velocity for our test photon. 

To demonstrate of the features in the library, we will choose the simplest [Schwarzschild spacetime](https://en.wikipedia.org/wiki/Schwarzschild_metric), which describes a spherically symmetric black hole with mass $M$. We will implement this metric ourselves.

!!! note
    Many metrics have already been implemented in Gradus.jl; for a comprehensive list, see [Catalogue of metrics](@ref).

```julia
using Gradus

struct Schwarzschild{T} <: AbstractStaticAxisSymmetric{T}
    M::T
end

function Gradus.metric_components(m::Schwarzschild, x)
    r, θ = x
    M = m.M

    dt2 = -(1 - (2M / r))
    dr2 = -inv(dt2)
    dθ2 = r^2
    dϕ2 = r^2 * sin(θ)^2
    dtdϕ = zero(r)

    SVector(dt2, dr2, dθ2, dϕ2, dtdϕ)
end

Gradus.inner_radius(m::Schwarzschild) = 2 * m.M
```

Going through this line by line:

```julia
struct Schwarzschild{T} <: AbstractStaticAxisSymmetric{T}
    M::T
end
```

- First we define a `struct` that will parameterise our spacetime. In this case, the mass $M$. We declare our struct to be a _subtype_ of [`AbstractStaticAxisSymmetric`](@ref) since we know our metric will be static (no time dependence) and axis-symmetric (no $\phi$ dependence). This describes the more general [Petrov type D](https://en.wikipedia.org/wiki/Petrov_classification) class of spacetimes, and allows Gradus.jl to make a number of simplifying assumptions under the hood about how this spacetime will behave. 
- The `T` parameter is the number type of this metric, and dictates the precision of all numerics in the trace. Therefore, if `M` is a `Float32`, Gradus.jl will raise errors if you attempt 64-bit floating point operations when tracing. This is _by design_, since many GPU architectures prefer `Float32` for speed, especially when precision is less important, and throwing errors is preferable to debugging type coercions. 

```julia
function Gradus.metric_components(m::Schwarzschild, x)
    r, θ = x
    M = m.M

    dt2 = -(1 - (2M / r))
    dr2 = -inv(dt2)
    dθ2 = r^2
    dϕ2 = r^2 * sin(θ)^2
    dtdϕ = zero(r)

    SVector(dt2, dr2, dθ2, dϕ2, dtdϕ)
end
```

- Here we have given the actual implementation of our metric. Since the metric is static, axis-symmetric, the position vector `x` only contains the radial and poloidal coordinates, and expects the `metric_components` function to return the five matrix elements of the metric. For the Schwarzschild metric, the $\text{d}t \text{d}\phi$ component is zero everywhere. We set this to `zero(r)`, which is a Julia function that returns 0 but of the same type as `r`.

!!! note
    Currently, Gradus.jl uses exclusively Boyer-Lindquist coordinates for its metrics. However, new coordinates can be implemented, and documentation for this will come soon.

```julia
GradusBase.inner_radius(m::KerrMetric) = 2 * m.M
```

- Finally, we specify some inner radius for the integration. This is the cutoff around the origin at which radius the geodesic integration will stop to avoid numerical errors. Here, it is just the Schwarzschild radius, or the outer event horizon. Gradus.jl can calculate different horizons from the metric automatically, which can be useful if you don't know the solution ahead of time, or if the solution is non-symmetric in $\theta$. But if we know it, we can benefit from a small performance boost by implementing it directly.

!!! note
    For a full description of implementing a metric, see [Implementing a new metric](@ref).

If you're familiar with other GRRT softwares, you might be wondering "where do we define the Christoffel symbols?", or "do I not need a prescription for Carter's constant?". Thanks to [automatic differentiation (AD)](https://en.wikipedia.org/wiki/Automatic_differentiation), we can calculate the Christoffel symbols _on the fly_! We determine the metric Jacobian with respect to coordinates of interest, and then sparsely compute the Christoffel symbols for the given spacetime class. For full details, see [Geodesic integration strategies](@ref).

## 2 Photon trajectories

The metric parameters and position can be easily and arbitrarily chosen, however the velocity has a precondition which must be satisfied.

```@docs
constrain
```

The [`constrain`](@ref) function is automatically invoked by [`tracegeodesics`](@ref) to normalize velocity vectors appropriately.

We'll position ourselves at a great distance from the singularity at the origin, around $1000 \, r_\text{g}$ away. We'll furthermore setup our spacetime with an arbitrary choice of mass `M = 1.0`, which acts to rescale our system (since all units in Gradus.jl are in standard GR units).

The initial velocity vector we will somewhat arbitrarily set to be directed towards the black hole ($v^r = -1$), with a small $v^\phi$ component so it grazes past the singularity.

The trajectory is calculated with a call to [`tracegeodesics`](@ref):

```julia
m = Schwarzschild(1.0)
x = SVector(0.0, 1000.0, π/2, 0.0)
v = SVector(0.0, -1.0, 0.0, -8e-6)

# maximum affine time ~ 2 * x[2]
λ_max = 2000.0
sol = tracegeodesics(m, x, v, λ_max)
```

The trajectory can be visualized with the use of [Plots.jl](https://docs.juliaplots.org/latest/):

```julia
using Plots

# plot solution trajectory
plot_paths(sol)
# plot 
plot_horizon!(m)
```

![](./figs/getting-started-1-trajectories.svg)

Choosing the initial velocity in this manner lacks interpretation. We can instead use so-called _impact parameters_ $(\alpha, \beta)$. These may be thought of as follows: 

```@docs
impact_parameters_to_three_velocity
```

Finally, if you imagine a two dimensional image plane, where $x$ is the horizontal and $y$ the vertical coordinate, the $\alpha$ impact parameter corresponds to that closest approach along the $x$ axis, and $\beta$ along the $y$ axis.

Hopefully that makes sense. With this, we can more easily setup a handful of geodesics to trace and know that they will roughly travel close to the central singularity.

```julia
# grid of impact parameters in horizontal direction
# keeping β fixed at 0
α = range(-10.0, 10.0, 30)
vs = map_impact_parameters(m, x, α, 0.0)

# need a position for each velocity vector
xs = fill(x, size(vs))

sols = tracegeodesics(m, xs, vs, λ_max)

# plot
p = plot_paths(sols, legend=false)
plot_horizon!(m, color = :black)
```

![](./figs/getting-started-2-multi-trajectories.svg)

When we invoke [`tracegeodesics`](@ref) in this way, Gradus.jl will automatically distribute the workload onto as many threads as Julia was started with. For example, starting julia with

```bash
julia -t6
```

will spawn 6 worker threads for Gradus.jl to use. Passing `-tauto` will allow Julia to use as many threads as your hardware supports.

!!! note
    For more about parallelism in Gradus.jl, see [Parallelism and ensembles](@ref). 

## 3 Rendering an image

A common task we'll want to do is render an image; that is, assign some $(\alpha, \beta)$ to each pixel in a 2-dimensional plane, located at the position `x`. Each pair of impact parameters is then traced along it, and its corresponding pixel coloured according to some function of the geodesic endpoint.

The simplest non-trivial thing we can do is colour the pixel by the time component of the final position. Since we are not interested in what happens to the geodesic along the path, only the start and end points, we can pass `save_on = false` to [`tracegeodesics`](@ref). This tells the solvers to not save intermediary points along the solution, and thereby avoid the overhead of allocating a memory we have no wish to use.

```julia
# set up our image parameters
α = range(-10.0, 10.0, 100)
β = range(-10.0, 10.0, 100)

# this will set up a 100x100 matrix of velocity vectors
# so we use `vec` to flatten the structure
vs = vec([map_impact_parameters(m, x, a, b) for a in α, b in β])
xs = fill(x, size(vs))

# trace in parallel
sols = tracegeodesics(m, xs, vs, λ_max, save_on = false)
```

To help us process the solutions, Gradus.jl exports a number of utility functions. There is one in specific we will want to use:

```@docs
process_solution
``` 
 
The [`GeodesicPoint`](@ref) struct contains everything we might want to know about the start and endpoint of a geodesic solution, including four-velocities and the nature of the termination (fell into the black hole, went to infinity, collided with some geometry, etc.). 

We can easily filter those geodesic that fell into the black hole, and extract their final coordinate times $x^t({\lambda_\text{final}})$:

```julia
points = process_solution.(sols.u)
# reshape into the same dimensions as the image
points = reshape(points, (100, 100))

times = map(points) do gp
    # check if went off the integration chart on the inner boundary
    if gp.status == StatusCodes.WithinInnerBoundary
        # get the time coordinate
        gp.x[1]
    else
        NaN
    end
end

heatmap(α, β, times, aspect_ratio = 1)
```

![](./figs/getting-started-3-basic-shadow.png)

This is the so-called _shadow_ of a black hole.

## 4 Defining and using `PointFunction`s

The above is a little verbose, when all we really wanted to do trace a given grid of $(\alpha, \beta)$, and then compute some physical quantity at each pixel. Having to manually write out the for loop and remember to reshape the solutions array is a trifle unnecessary and error-prone. Furthermore, to someone reading our code, it may not always be obvious what physical quantity it is that we are calculating from just the `for` loop.

Gradus.jl instead uses a callback function, which in its own parlance is dubbed the [`PointFunction`](@ref). These functions are isolated small physical quantities, that allow us to compose a more complex observable. Many point functions have already been implemented ready for use.

As a motivating example, say we wanted to write the above as a [`PointFunction`](@ref):

```julia
time_coord = PointFunction((m, gp, λ) -> gp.x[1])
```

Point functions always receive the metric parameters `m`, a geodesic point `gp`, and the final time of the integration `λ`. To then filter those geodesics which fell into the event horizon, we can use a [`FilterPointFunction`](@ref) and compose them. Here, we use one of the utility methods [`FilterStatusCode`](@ref).

```julia
filter_event_horizon = FilterStatusCode(StatusCodes.WithinInnerBoundary)
# compose in reverse order
pf = time_coord ∘ filter_event_horizon
```

We can then apply our point function on the geodesic points:

```julia
times = pf.(m, points, λ_max)
```

Point functions can also be used in other contexts. For example, [`rendergeodesics`](@ref) is a utility method to help render images, and one of the keywords we can pass is `pf`, so that each pixel value is coloured by the point function we gave. We can create a higher resolution render of the above easily using [`rendergeodesics`](@ref):

```julia
# this function returns the impact parameter axes
α, β, image = rendergeodesics(
    m, 
    x,
    # no longer need to specify the velocities
    # these are automatically calculated
    λ_max, 
    pf = pf, 
    # image parameters
    image_width = 800, 
    image_height = 800,
    # the "zoom" -- field of view scale
    fov = 52,
    verbose = true
)

heatmap(α, β, image, aspect_ratio = 1)
```

![](./figs/getting-started-4-hr-shadow.png)

This is good, but let's make our render a bit more interesting.

### Short aside: `rendergeodesics` and `tracegeodesics`

The two functions [`rendergeodesics`](@ref) and [`tracegeodesics`](@ref) should be regarded as the fundamental tools that Gradus.jl provides. The former is used to create visualisations, used to _see_ what is going on with your system in images. The latter is the entry point for modelling physical processes, being much more versatile than [`rendergeodesics`](@ref), but also requiring preparatory work.

## 5 Adding geometry

The simplest thing we can do is put a disc around our black hole and visualize that system. Gradus.jl implements many different accretion disc types, but some kind of torus would be good to start with.

!!! note
    There are many different disc types already implemented, see [Available accretion geometry](@ref). For adding your own geometry, see [Adding new accretion geometry](@ref).

We may be imaginative and specify our own cross-section for the disc. We need to specify the function in the positive `z` axis, and this will be mirrored in the $x-y$ plane and rotated around the black hole.

The cross section need not be physical. Choosing some arbitrary shape, we can preview what our cross section will look like over a sample range:

```julia
function cross_section(x)
    # centered circle on 8 rg
    center = 8
    radius = 3

    if (x < center - radius) || (radius + center < x)
        zero(x)
    else
        r = x - center
        sqrt(radius^2 - r^2) + (0.5sin(3x))
    end
end

# preview the cross section over a sample range
sample = collect(range(0.0, 20.0, 300))
y = cross_section.(sample)

plot(sample, y, xlabel = "r", ylabel = "height", aspect_ratio = 1)
```

![](./figs/getting-started-5-cross-section.svg)

We then wrap our cross section function as a [`ThickDisc`](@ref) type:

```julia
d = ThickDisc(x -> cross_section(x[2]))
```

The thick disc callback receives the full four-position, so we forward only the radial component. 

We now need to update our point function so that it filters those geodesics which intersected with the geometry instead of those that fell into the black hole. This is a standard function already implemented in [`ConstPointFunctions`](@ref); only a composition is needed:

```julia
pf_geometry = time_coord ∘ ConstPointFunctions.filter_intersected
```

We then make a handful of small changes to make our image more interesting, and render just as before, passing the disc in to the [`rendergeodesics`](@ref) function:

```julia
# change inclination
x = SVector(0.0, 1000.0, deg2rad(70), 0.0)

α, β, image = rendergeodesics(
    m, 
    x,
    # add the disc argument
    d,
    λ_max, 
    # new point function
    pf = pf_geometry, 
    # slightly wider image
    image_width = 1200, 
    image_height = 800,
    # zoom out a little
    fov = 22,
    verbose = true
)

heatmap(α, β, image, aspect_ratio = 1)
```

![](./figs/getting-started-6-weird-disc-render.png)

## 6 Calculating physical quantities

A common quantity to look at when ray tracing is the _redshift_ of a photon; that is, the ratio of the energy where the photon was emitted to where it was observed. Mathematically, this is the quantity

```math
g = \frac{\left. k_\mu u^\mu \right\rvert_\text{obs}}{\left. k_\nu u^\nu \right\rvert_\text{em}},
```

where the subscript denote the observer and emitters position respectively. Here, $k_\mu$ is the covariant momentum of the photon, and $u^\mu$ is the velocity of the disc patch, or of the observer respectively. 

We can choose any velocity profile we like, but for simplicity we use the velocity of the stable circular orbit at the corresponding radius where the photon hit the disc. The above formula for the redshift $g$ is already implemented with this velocity profile for us -- we need only specify which spacetime we are in and where our observer is positioned:

```julia
redshift = ConstPointFunctions.redshift(m, x)
# compose to filter those that intersected with the geometry
redshift_geometry = redshift ∘ ConstPointFunctions.filter_intersected
```

This is another [`PointFunction`](@ref), and is used in the same way. Rendering as before:

```julia
α, β, image = rendergeodesics(
    m, 
    x,
    d,
    λ_max, 
    # new point function
    pf = redshift_geometry,
    image_width = 1200, 
    image_height = 800,
    fov = 22,
    verbose = true
)

heatmap(α, β, image, aspect_ratio = 1)
```

![](./figs/getting-started-7-weird-redshift.png)

## 7 Changing metric

To change the metric, we need only pass a new metric to [`ConstPointFunctions.redshift`](@ref) to update how the redshift is calculated, and to [`rendergeodesics`](@ref) to update how the geodesic equation is integrated. To switch to e.g. the Johannsen (2013)[^1] metric, the following modifications are needed:

```julia
j_m = JohannsenMetric(M=1.0, a = 0.7, α13 = 2.0, ϵ3 = 1.0)

# pass the new metric
j_redshift = ConstPointFunctions.redshift(j_m, x)
j_redshift_geometry = j_redshift ∘ ConstPointFunctions.filter_intersected

α, β, image = rendergeodesics(
    # pass the new metric
    j_m, 
    x,
    d,
    λ_max, 
    # and the new point function
    pf = j_redshift_geometry,
    image_width = 1200, 
    image_height = 800,
    fov = 22,
    verbose = true
)

heatmap(α, β, image, aspect_ratio = 1)
```

![](./figs/getting-started-8-jm-weird-redshift.png)

## 8 Calculating line profiles

As a final step, we can calculate the line profile of emission from our system with the Johannsen metric. We can compare this to the case where the disc is geometrically thin in the equatorial plane, and furthermore compare this to the Schwarzschild spacetime.

Line profiles are calculated with [`lineprofile`](@ref), a function that accepts much the same arguments as the previous tracing and rendering functions. We will also limit the domain of the integration to in the upper hemisphere only -- thereby avoiding any false images in our line profile calculations. This choice is physically motiviated, as the gaps in the inner regions of the disc are often assumed to be opaque due to extreme ionization of the matter.

Domain limiting can be done by adding a callback to the integrator.

!!! note
    For more on using callbacks and finer control of the integrator, see [Using callbacks](@ref).

```julia
# define custom bins for g
bins = collect(range(0.1, 1.4, 200))

# define the plane to perform the binning over
plane = PolarPlane(GeometricGrid(); Nr = 1000, Nθ = 1000, r_max = 50.0)
```

In the line above we created an explicit [`PolarPlane`](@ref), since we no longer wish to integrate on a rectangular grid. We want to primarily sample the region close to the event horizon where all of the interesting physics is taking place, and as such we scale the radial coordinate geometrically. Adjusting `Nr` and `Nθ` lets us control the "resolution" of our render, and will smoothen the line profiles. The values chosen here are to balance resolution and computational time (~ 30 seconds on a 2021 M1 Mac laptop).

We can preview what the grid will look like (though at lower resolution to avoid unnecessary noise):

```julia
plot(
    PolarPlane(GeometricGrid(); Nr = 10, Nθ = 20, r_max = 50.0)
)
```

![](./figs/getting-started-9-polar-plane.svg)

Each point on this plane represent a photon which will be traced, and the intensity scaled according to the area the point covers on the image plane.

With that, we are ready to calculate the line profiles. To avoid having to reuse large parts of our code, we can write a short function that wraps [`lineprofile`](@ref):

```julia
function calculate_line_profile(m, x, d, bins, plane)
    _, f = lineprofile(
        m, 
        x, 
        d, 
        algorithm = BinnedLineProfile(), 
        # no false images
        callback = domain_upper_hemisphere(),
        verbose = true,
        bins = bins,
        plane = plane,
    )
    return f
end
```

Note that [`lineprofile`](@ref) returns both the redshift $g$ (bins) and flux at each $g$. Since we specified the binning, we can ignore the first return value, and keep only the flux.

```julia
d_j_thin = GeometricThinDisc(Gradus.isco(j_m), 200.0, π / 2)
# and for the schwarzschild metric
d_s_thin = GeometricThinDisc(Gradus.isco(m), 200.0, π / 2)

f_j_thick_disc = calculate_line_profile(j_m, x, d, bins, plane)
f_s_thick_disc = calculate_line_profile(m, x, d, bins, plane)
f_j_thin_disc = calculate_line_profile(j_m, x, d_j_thin, bins, plane)
f_s_thin_disc = calculate_line_profile(m, x, d_s_thin, bins, plane)

plot(bins, f_j_thick_disc, label = "Johannsen[thick]")
plot!(bins, f_s_thick_disc, label = "Schwarzschild[thick]")
plot!(bins, f_j_thin_disc, label = "Johannsen[thin]")
plot!(bins, f_s_thin_disc, label = "Schwarzschild[thin]")
```

![](./figs/getting-started-10-line-profiles.svg)

!!! note
    For how to use these line profiles and other observables in fitting programs, see [Exporting data products](@ref).

## Where to go from here?

The documentation is a rich resource for information related to using Gradus, and tailoring the toolkit for your needs. Take a look at [Examples](@ref) for a number of quick recipes, or try [Implementing a new metric](@ref) and study a different spacetime.


[^1]: Johannsen, Tim. ‘Regular Black Hole Metric with Three Constants of Motion’. Physical Review D 88, no. 4 (1 August 2013): 044002. https://doi.org/10.1103/PhysRevD.88.044002.
