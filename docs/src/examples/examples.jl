using Gradus
using StaticArrays
using Plots
ENV["GKSwstype"] = "nul"
Plots.default(show = false)

using StatsBase

function ex_tracing()
    m = JohannsenPsaltisMetric(M = 1.0, a = 0.6, ϵ3 = 2.0)
    # observer position
    u = @SVector [0.0, 1000.0, π / 2, 0.0]

    # set up impact parameter space
    α = collect(range(-10.0, 10.0, 20))
    β = [0.0 for _ in α]

    # build initial velocity and position vectors
    vs = map_impact_parameters(m, u, α, β)
    us = [u for _ in vs]

    sols = tracegeodesics(m, us, vs, (0.0, 2000.0); abstol = 1e-12, reltol = 1e-12)

    # only use the subset of the solution we're plotting
    trange = range(990, 1035, 5000)

    p = plot(projection = :polar, legend = false, range = (0, 10))
    for s in sols
        r = [s(t)[2] for t in trange]
        ϕ = [s(t)[4] for t in trange]
        plot!(p, ϕ, r)
    end

    # plot event horizon 
    r0 = inner_radius(m)
    plot!(p, collect(range(0, 2π, 200)), [r0 for _ = 1:200], color = :black, linewidth = 2)

    savefig(p, "example-tracing.svg")
end

function ex_redshift()
    # metric and metric parameters
    m = KerrMetric(M = 1.0, a = 1.0)
    # observer position
    u = @SVector [0.0, 1000.0, deg2rad(60), 0.0]
    # accretion disc
    d = GeometricThinDisc(1.0, 50.0, deg2rad(90))

    # define point function which filters geodesics that intersected the accretion disc
    # and use those to calculate redshift
    pf = ConstPointFunctions.redshift(m, u) ∘ ConstPointFunctions.filter_intersected

    α, β, img = rendergeodesics(
        m,
        u,
        d,
        # maximum integration time
        2000.0,
        fov_factor = 6.0,
        image_width = 700,
        image_height = 240,
        verbose = true,
        pf = pf,
    )

    p = heatmap(α, β, img)
    savefig(p, "example-redshift.png")


    # remove nans and flatten the redshift image
    redshift_data = filter(!isnan, vec(img))

    # transpose to iron-line
    data = redshift_data .* 6.4

    x_bins = range(0.0, 10.0, 100)
    lineprof = fit(Histogram, data, x_bins)

    p = plot(x_bins[1:end-1], lineprof.weights, seriestype = :steppre)
    savefig(p, "example-redshift-histogram.svg")
end

function ex_interpolating()
    # metric and metric parameters
    m = KerrMetric(M = 1.0, a = 0.4)
    # observer's initial position
    u = @SVector [0.0, 1000.0, deg2rad(85), 0.0]
    # accretion disc
    d = GeometricThinDisc(1.0, 50.0, deg2rad(90))

    pl_int = interpolate_plunging_velocities(m)

    redshift = interpolate_redshift(pl_int, u)

    pf = redshift ∘ ConstPointFunctions.filter_intersected

    α, β, img = rendergeodesics(
        m,
        u,
        d,
        # maximum integration time
        2000.0,
        fov_factor = 6.0,
        image_width = 700,
        image_height = 240,
        verbose = true,
        pf = pf,
    )

    p = heatmap(α, β, img)
    savefig(p, "example-interpolated.png")
end

function ex_circular_orbits()
    m = KerrMetric(M = 1.0, a = 0.8)

    p = plot(aspect_ratio = 1)

    for r in [3.0, 4.0, 5.0, 6.0]
        v = CircularOrbits.fourvelocity(m, r)
        # trace the circular orbit
        path = tracegeodesics(m, @SVector([0.0, r, π / 2, 0.0]), v, (0.0, 300.0), μ = 1.0)
        r = [path(t)[2] for t in range(0.0, 100, 200)]
        ϕ = [path(t)[4] for t in range(0.0, 100, 200)]

        x = @. r * cos(ϕ)
        y = @. r * sin(ϕ)

        plot!(p, x, y, label = false)
    end

    p

    savefig(p, "example-circular-orbits.svg")
end

function ex_isco()
    # prepare plot
    p = plot(legend = :bottomright, ylabel = "E", xlabel = "r", xscale = :log10)

    # choice of spin to plot energy curves for
    for a in [0.0, 0.4, 0.6]
        m = KerrMetric(M = 1.0, a = a)

        rs = range(Gradus.isco(m), 100.0, 500)
        energy = map(rs) do r
            CircularOrbits.energy(m, r)
        end

        plot!(rs, energy, label = "a=$a")
    end

    # calculate the ISCO as a function of spin
    data = map(range(-1.0, 0.8, 100)) do a
        m = KerrMetric(M = 1.0, a = a)
        r = Gradus.isco(m)
        CircularOrbits.energy(m, r), r
    end

    # overlay onto plot
    plot!(last.(data), first.(data), color = :black, linestyle = :dash, label = "ISCO")

    savefig(p, "example-isco.svg")
end

function ex_horizon()
    function draw_horizon(p, m)
        rs, θs = event_horizon(m, resolution = 200)
        radius = rs

        x = @. radius * cos(θs)
        y = @. radius * sin(θs)
        plot!(p, x, y, label = "a = $(m.a)")
    end

    p = plot(aspect_ratio = 1)
    for a in [0.0, 0.5, 0.6, 0.7, 0.8]
        m = JohannsenPsaltisMetric(M = 1.0, a = a, ϵ3 = 2.0)
        draw_horizon(p, m)
    end
    p

    savefig(p, "example-horizon.svg")


    function calc_exclusion(as, ϵs)
        regions = [
            is_naked_singularity(JohannsenPsaltisMetric(M = 1.0, a = a, ϵ3 = ϵ)) for
            a in as, ϵ in ϵs
        ]

        map(i -> i ? 1.0 : NaN, regions)
    end

    # define ranges (small in this example as a little computationally intense)
    as = range(0, 1.0, 40)
    ϵs = range(-10, 10, 40)

    img = calc_exclusion(as, ϵs)
    p = heatmap(as, ϵs, img', color = :black, colorbar = false, xlabel = "a", ylabel = "ϵ")

    savefig(p, "example-exclusion.png")
end

function ex_transfer_functions()
    m = KerrMetric(M = 1.0, a = 0.998)
    d = GeometricThinDisc(0.0, 200.0, π / 2)

    p = plot(legend = false)
    for angle in [3, 35, 50, 65, 74, 85]
        u = @SVector [0.0, 1000.0, deg2rad(angle), 0.0]
        ctf = cunningham_transfer_function(m, u, d, 4.0, 2000.0)
        mask = @. (ctf.gstar > 0.001) & (ctf.gstar < 0.999)
        @views plot!(p, ctf.gstar[mask], ctf.f[mask])
    end
    p

    savefig(p, "example-bambi-fig1.svg")

    # new position vector
    u = @SVector [0.0, 1000.0, deg2rad(30), 0.0]

    p = plot(legend = false)
    for a in [0.0, 0.25, 0.5, 0.75, 0.9, 0.998]
        m = KerrMetric(1.0, a)
        ctf = cunningham_transfer_function(m, u, d, 7.0, 2000.0)
        mask = @. (ctf.gstar > 0.001) & (ctf.gstar < 0.999)
        @views plot!(p, ctf.gstar[mask], ctf.f[mask])
    end
    p

    savefig(p, "example-bambi-fig2.svg")
end

function ex_concentric_rings()
    # their papers has a=-a
    m = KerrMetric(M = 1.0, a = -0.4)
    u = @SVector [0.0, 1000, acos(0.25), 0.0]
    d = GeometricThinDisc(0.0, 100.0, π / 2)

    radii = 2.6:1.0:7.6

    p = plot(aspect_ratio = 1, legend = false)

    # crosshair on origin
    hline!(p, [0.0], color = :black, linestyle = :dash)
    vline!(p, [0.0], color = :black, linestyle = :dash)

    for r in radii
        α, β = impact_parameters_for_radius(m, u, d, r)
        plot!(p, α, β)
    end

    p

    savefig(p, "example-concentric.svg")
end

function ex_doughnut()
    m = KerrMetric(1.0, 0.0)
    u = @SVector [0.0, 1000.0, deg2rad(85), 0.0]

    # define the disc shape -- return a negative number 
    # where the disc should not be intersected, else the cross 
    # sectional height
    d = ThickDisc() do u
        r = u[2]
        if r < 9.0 || r > 11.0
            return -1.0
        else
            x = r - 10.0
            sqrt(1 - x^2)
        end
    end

    # and then render as usual
    α, β, img = rendergeodesics(
        m,
        u,
        d,
        2000.0,
        fov_factor = 18.0,
        image_width = 700,
        image_height = 350,
        verbose = true,
    )

    p = heatmap(α, β, img, aspect_ratio = 1)
    savefig(p, "example-thick-disc-doughtnut.png")
end

function ex_lineprofile()
    d = GeometricThinDisc(0.0, 400.0, π / 2)
    u = @SVector [0.0, 1000.0, deg2rad(40), 0.0]
    m = KerrMetric(1.0, 0.998)

    # maximal integration radius
    maxrₑ = 50.0

    # emissivity function
    ε(r) = r^(-3)

    # g grid to do flux integration over
    gs = range(0.0, 1.2, 500)
    _, flux = lineprofile(gs, ε, m, u, d, maxrₑ = maxrₑ)

    # transform to observed energy
    energy = gs .* 6.4

    p = plot(energy, flux, legend = false)
    savefig(p, "example-line-profile.svg")
end

ex_tracing()
ex_redshift()
ex_interpolating()
ex_circular_orbits()
ex_isco()
ex_horizon()
ex_transfer_functions()
ex_concentric_rings()
ex_doughnut()
ex_lineprofile()
