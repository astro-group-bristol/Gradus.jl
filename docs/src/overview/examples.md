# Examples

```@contents
Pages = ["examples.md"]
Depth = 3
```

```@setup env
using Gradus
using StaticArrays
using Plots
ENV["GKSwstype"]="nul"
Plots.default(show=false)
```

## Tracing geodesic paths

```@example env
using Gradus
using Plots
using StaticArrays

m = JohannsenPsaltisAD(M=1.0, a=0.6, ϵ3=2.0)
# observer position
u = @SVector [0.0, 1000.0, π/2, 0.0]

# set up impact parameter space
α = collect(range(-10.0, 10.0, 20))
β = [0.0 for _ in α]

# build initial velocity and position vectors
vs = map_impact_parameters(m, u, α, β)
us = [u for _ in vs]

sols = tracegeodesics(
    m, us, vs, (0.0, 2000.0);
    abstol = 1e-12, reltol = 1e-12
)

# only use the subset of the solution we're plotting
trange = range(990, 1035, 5000)

p = plot(projection = :polar, legend = false, range = (0, 10))
for s in sols
    r = [s(t)[2] for t in trange]
    ϕ = [s(t)[4] for t in trange]
    plot!(p, ϕ, r)
end

# plot event horizon 
r0 = Gradus.GradusBase.inner_radius(m)
plot!(p, collect(range(0, 2π, 200)), [r0 for _ in 1:200], color = :black, linewidth = 2)

p
```

## Redshift image

!!! note
    The [`Gradus.ConstPointFunctions.redshift`](@ref) function is an analytic solution for redshift, which may not be implemented for every type of metric or disc geometry. See [Interpolating redshifts](@ref) for a more flexible numeric alternative.

```@example env
using Gradus
using StaticArrays
using Plots

# metric and metric parameters
m = BoyerLindquistAD(M=1.0, a=1.0)
# observer position
u = @SVector [0.0, 1000.0, deg2rad(60), 0.0]
# accretion disc
d = GeometricThinDisc(1.0, 50.0, deg2rad(90))

# define point function which filters geodesics that intersected the accretion disc
# and use those to calculate redshift
pf = Gradus.ConstPointFunctions.redshift ∘ Gradus.ConstPointFunctions.filter_intersected

img = rendergeodesics(
    m,
    u,
    d,
    # maximum integration time
    2000.0,
    fov_factor = 6.0,
    image_width = 700,
    image_height = 240,
    verbose = true,
    pf = pf
)

heatmap(img)
```

## Line profile

Using the redshift example, we can bin a line-profile using [StatsBase.jl](https://juliastats.org/StatsBase.jl/stable/empirical/#StatsBase.Histogram). We'll calculate the iron line profile, with a delta-emission at 6.4 keV.

```@example env
using StatsBase

# remove nans and flatten the redshift image
redshift_data = filter(!isnan, vec(img))

# transpose to iron-line
data = redshift_data .* 6.4

x_bins = range(0.0, 10.0, 100) 
lineprof = fit(Histogram, data, x_bins)

plot(x_bins[1:end-1], lineprof.weights, seriestype = :steppre)
```

## Interpolating redshifts

In cases where no analytic redshift solution is known, we can instead interpolate a numeric approximation. For example, interpolating the plunging region velocities and using the analytic solution for general static, axis symmetric metrics outside of the ISCO can be achieved with:

```@example env
using Gradus
using StaticArrays
using Plots

# metric and metric parameters
m = BoyerLindquistAD(M=1.0, a=0.4)
# observer's initial position
u = @SVector [0.0, 1000.0, deg2rad(85), 0.0]
# accretion disc
d = GeometricThinDisc(1.0, 50.0, deg2rad(90))

pl_int = interpolate_plunging_velocities(m)

redshift = interpolate_redshift(pl_int, u)

pf = redshift ∘ Gradus.ConstPointFunctions.filter_intersected

img = rendergeodesics(
    m,
    u,
    d,
    # maximum integration time
    2000.0,
    fov_factor = 6.0,
    image_width = 700,
    image_height = 240,
    verbose = true,
    pf = pf
)

heatmap(img)
```

For more complex disc geometry: TODO

## Circular orbits

Simple equatorial circular orbits are straight forward to calculate with Gradus.jl:

```@example env
using Gradus
using Plots
using StaticArrays

m = BoyerLindquistAD(M=1.0, a=0.8)

p = plot(aspect_ratio=1)

for r in [3.0, 4.0, 5.0, 6.0]
    v = Gradus.CircularOrbits.fourvelocity(m, r)
    # trace the circular orbit
    path = tracegeodesics(m, @SVector([0.0, r, π/2, 0.0]), v, (0.0, 300.0), μ = 1.0)
    r = [path(t)[2] for t in range(0.0, 100, 200)]
    ϕ = [path(t)[4] for t in range(0.0, 100, 200)]

    x = @. r * cos(ϕ)
    y = @. r * sin(ϕ)

    plot!(p, x, y, label = false)
end

p
```

## ISCO

The [Gradus.isco](@ref) may be calculated with a simple convenience function, as may the energy associated with a given stable circular orbit.

```@example env
using Gradus
using Plots

# prepare plot
p = plot(legend=:bottomright, ylabel = "E", xlabel = "r", xscale = :log10)

# choice of spin to plot energy curves for
for a in [0.0, 0.4, 0.6]
    m = BoyerLindquistAD(M = 1.0, a = a)

    rs = range(Gradus.isco(m), 100.0, 500)
    energy = map(rs) do r
        Gradus.CircularOrbits.energy(m, r)
    end

    plot!(rs, energy, label = "a=$a")
end

# calculate the ISCO as a function of spin
data = map(range(-1.0, 0.8, 100)) do a
    m = BoyerLindquistAD(M = 1.0, a = a)
    r = Gradus.isco(m)
    Gradus.CircularOrbits.energy(m, r), r
end

# overlay onto plot
plot!(last.(data), first.(data), color=:black, linestyle=:dash, label="ISCO")

p
```

## Event horizons and naked singularities

Here is an example of how to use [`Gradus.event_horizon`](@ref) to plot the shape of an event horizon in two dimensions. In the case of a naked singularity, as with the certain parameters combinations in the [`JohannsenPsaltisAD`](@ref) metric, we see a disconnected region in the plot.

```@example env
using Gradus
using Plots

function draw_horizon(p, m)
    rs, θs = Gradus.event_horizon(m, resolution = 200)
    radius = rs

    x = @. radius * sin(θs)
    y = @. radius * cos(θs)
    plot!(p, x, y, label = "a = $(m.a)")
end

p = plot(aspect_ratio = 1)
for a in [0.0, 0.5, 0.6, 0.7, 0.8]
    m = JohannsenPsaltisAD(M = 1.0, a = a, ϵ3 = 2.0)
    draw_horizon(p, m)
end
p
```

We can also calculate parameter combinations that lead to naked singularities, and plot the parameter space domains to show exclusion zones:

```@example env
function calc_exclusion(as, ϵs)
    regions = [
        Gradus.is_naked_singularity(JohannsenPsaltisAD(M = 1.0, a = a, ϵ3 = ϵ))
        for a in as, ϵ in ϵs
    ]

    map(i -> i ? 1.0 : NaN, regions)
end

# define ranges (small in this example as a little computationally intense)
as = range(0, 1.0, 40)
ϵs = range(-10, 10, 40)

img = calc_exclusion(as, ϵs)
heatmap(
    as, 
    ϵs, 
    img', 
    color = :black, 
    colorbar = false, 
    xlabel = "a", 
    ylabel = "ϵ"
)
```

## Cunningham transfer functions

Recreating Fig. 1 and 2 from [Bambi et al. (2017)](https://iopscience.iop.org/article/10.3847/1538-4357/aa74c0) for the transfer functions of a Kerr black hole

```@example env
using Gradus
using StaticArrays
using Plots

p = plot(legend = false)
for angle in [3, 35, 50, 65, 74, 85]
    @show angle
    u = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
    ctf = @time Gradus.cunningham_transfer_function(
        m, u, d, 4.0, 2000.0; finite_diff_order = 5
    )
    mask = @. (ctf.gstar > 0.001) & (ctf.gstar < 0.999)
    @views plot!(p, ctf.gstar[mask], ctf.f[mask])
end
p
```

And Fig. 2:

```@example env
p = plot(legend = false)
for a in [0.0, 0.25, 0.5, 0.75, 0.9, 0.998]
    @show a
    m = BoyerLindquistAD(1.0, a)
    ctf = @time Gradus.cunningham_transfer_function(
        m, u, d, 7.0, 2000.0; finite_diff_order = 5
    )
    mask = @. (ctf.gstar > 0.001) & (ctf.gstar < 0.999)
    @views plot!(p, ctf.gstar[mask], ctf.f[mask])
end
p
```