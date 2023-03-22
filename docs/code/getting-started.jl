using Gradus

struct Schwarzschild{T} <: AbstractStaticAxisSymmetricParameters{T}
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

# ---------------------------------------------------------

m = Schwarzschild(1.0)

x = SVector(0.0, 1000.0, π / 2, 0.0)
v = SVector(0.0, -1.0, 0.0, -8e-6)

sol = tracegeodesics(m, x, v, 2000.0)

using Plots

p = plot_paths(sol)
plot_horizon!(m, color = :black)

# ---------------------------------------------------------

# grid of impact parameters in horizontal direction
# keeping β fixed at 0
α = range(-10.0, 10.0, 30)
vs = map_impact_parameters(m, x, α, 0.0)

# need a position for each velocity vector
xs = fill(x, size(vs))

sols = tracegeodesics(m, xs, vs, 2000.0)

# plot
p = plot_paths(sols, legend = false)
plot_horizon!(m, color = :black)

# ---------------------------------------------------------

# set up our image parameters
α = range(-10.0, 10.0, 100)
β = range(-10.0, 10.0, 100)

# this will set up a 100x100 matrix of velocity vectors
# so we use `vec` to flatten the structure
vs = vec([map_impact_parameters(m, x, a, b) for a in α, b in β])
xs = fill(x, size(vs))

# trace in parallel
sols = tracegeodesics(m, xs, vs, 2000.0, save_on = false)

# ---------------------------------------------------------

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

p = heatmap(α, β, times, aspect_ratio = 1)

# ---------------------------------------------------------

time_coord = PointFunction((m, gp, λ) -> gp.x[1])

# ---------------------------------------------------------

filter_event_horizon = FilterStatusCode(StatusCodes.WithinInnerBoundary)
# compose in reverse order
pf = time_coord ∘ filter_event_horizon

# ---------------------------------------------------------

times = pf.(m, points, 2000.0)

# ---------------------------------------------------------

α, β, image = rendergeodesics(
    m,
    x,
    # no longer need to specify the velocities
    # these are automatically calculated
    2000.0,
    pf = pf,
    # image parameters
    image_width = 800,
    image_height = 800,
    # the "zoom" -- field of view scale
    fov = 52,
    verbose = true,
)

heatmap(α, β, image, aspect_ratio = 1)

# ---------------------------------------------------------

function cross_section(x)
    # centered circle on 12 rg
    center = 15
    radius = 4

    if (x < center - radius) || (radius + center < x)
        zero(x)
    else
        r = x - center
        sqrt(radius^2 - r^2) + (0.5sin(3x))
    end
end

# preview the cross section over a sample range
sample = collect(range(0.0, 30.0, 300))
y = cross_section.(sample)

p = plot(sample, y, xlabel = "r", ylabel = "height", aspect_ratio = 1)


# ---------------------------------------------------------

d = ThickDisc(x -> cross_section(x[2]))

# ---------------------------------------------------------

pf_geometry = time_coord ∘ ConstPointFunctions.filter_intersected

# ---------------------------------------------------------

# change inclination
x = SVector(0.0, 1000.0, deg2rad(70), 0.0)

α, β, image = rendergeodesics(
    m,
    x,
    # add the disc argument
    d,
    2000.0,
    # new point function
    pf = pf_geometry,
    # slightly wider image
    image_width = 1200,
    image_height = 800,
    # zoom out a little
    fov = 22,
    verbose = true,
)

heatmap(α, β, image, aspect_ratio = 1)
# ---------------------------------------------------------

redshift = ConstPointFunctions.redshift(m, x)
# compose to filter those that intersected with the geometry
redshift_geometry = redshift ∘ ConstPointFunctions.filter_intersected

# ---------------------------------------------------------

α, β, image = rendergeodesics(
    m,
    x,
    d,
    2000.0,
    # new point function
    pf = redshift_geometry,
    image_width = 1200,
    image_height = 800,
    fov = 22,
    verbose = true,
)

heatmap(α, β, image, aspect_ratio = 1)
# --------------------------------------

j_m = JohannsenMetric(M = 1.0, a = 0.7, α13 = 2.0, ϵ3 = 1.0)
is_naked_singularity(j_m)

# pass the new metric
j_redshift = ConstPointFunctions.redshift(j_m, x)
j_redshift_geometry = j_redshift ∘ ConstPointFunctions.filter_intersected

α, β, image = rendergeodesics(
    # pass the new metric
    j_m,
    x,
    d,
    2000.0,
    # and the new point function
    pf = j_redshift_geometry,
    image_width = 1200,
    image_height = 800,
    fov = 22,
    verbose = true,
)

p = heatmap(α, β, image, aspect_ratio = 1)

# ---------------------------------------------------------

# define custom bins for g
bins = collect(range(0.1, 1.4, 200))

# define the plane to perform the binning over
plane = PolarPlane(GeometricGrid(); Nr = 700, Nθ = 1300, r_max = 50.0)

# ---------------------------------------------------------

p = plot(PolarPlane(GeometricGrid(); Nr = 10, Nθ = 20, r_max = 50.0))

# ---------------------------------------------------------
_, f = lineprofile(
    j_m,
    x,
    d,
    algorithm = BinnedLineProfile(),
    # no false images
    callback = domain_upper_hemisphere(),
    verbose = true,
    bins = bins,
    plane = plane,
)

# geometric thin disc in the equitorial plane
d_thin = GeometricThinDisc(Gradus.isco(j_m), 200.0, π / 2)

_, f_thin = lineprofile(
    j_m,
    x,
    # new disc
    d_thin,
    algorithm = BinnedLineProfile(),
    callback = domain_upper_hemisphere(),
    verbose = true,
    bins = bins,
    plane = plane,
)

p = plot(bins, f, label = "thick")
plot!(bins, f_thin, label = "thin")


d_thin2 = GeometricThinDisc(0.0, 400.0, π / 2)
g, kk = lineprofile(m, x, d_thin2, verbose = true)
g2, kk2 = lineprofile(m, x, d_thin2, verbose = true, algorithm = BinnedLineProfile())
plot(g, kk)
plot!(g2, kk2)

# ---------------------------------------------------------
