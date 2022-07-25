# Examples

## Redshift image

```julia
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

![redshift](../assets/examples/redshift.png)

## Line profile

Using the redshift example, we can bin a line-profile using [StatsBase.jl](https://juliastats.org/StatsBase.jl/stable/empirical/#StatsBase.Histogram). We'll calculate the iron line profile, with a delta-emission at 6.4 keV.

```julia
using StatsBase

# remove nans and flatten the redshift image
redshift_data = filter(!isnan, vec(img))

# transpose to iron-line
data = redshift_data .* 6.4

x_bins = range(0.0, 10.0, 100) 
lineprof = fit(Histogram, data, x_bins)

plot(x_bins[1:end-1], lineprof.weights, seriestype = :steppost)
```

![iron-line-profile](../assets/examples/line-profile.svg)