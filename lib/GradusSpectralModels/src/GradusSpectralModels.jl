module GradusSpectralModels

using Gradus, SpectralFitting

struct LineProfile{D,T} <: AbstractTableModel{T,Additive}
    table::D
    K::T
    "Spin"
    a::T
    "Observer inclination (degrees off of the spin axis)."
    θ::T
    "Inner radius of the accretion disc."
    rin::T
    "Outer radius of the accretion disc."
    rout::T
    "Central emission line energy (keV)."
    E₀::T
end

function Base.copy(m::LineProfile)
    table = Gradus.CunninghamTransferTable(m.table.table.params, m.table.table.grids)
    setup = Gradus.IntegrationSetup(
        m.table.setup.h,
        m.table.setup.profile,
        m.table.setup.integrand,
        m.table.setup.quadrature_rule,
        deepcopy(m.table.setup.index_cache),
        m.table.setup.g_grid_upscale,
        m.table.setup.n_radii,
        m.table.setup.t0,
    )
    typeof(m)(
        (; setup = setup, table = table),
        (copy(getproperty(m, f)) for f in fieldnames(typeof(m))[2:end])...,
    )
end

function LineProfile(
    profile,
    table::Gradus.CunninghamTransferTable;
    K = FitParam(1.0),
    a = FitParam(0.998),
    θ = FitParam(45.0),
    rin = FitParam(1.0),
    rout = FitParam(100.0, upper_limit = 100.0),
    E₀ = FitParam(1.0),
    kwargs...,
)
    setup = integration_setup(profile, table((get_value(a), get_value(θ))); kwargs...)
    LineProfile((; setup = setup, table = table), K, a, θ, rin, rout, E₀)
end

function SpectralFitting.invoke!(output, domain, model::LineProfile)
    grid = model.table.table((model.a, model.θ))
    rmin = if model.rin < grid.r_grid[1]
        grid.r_grid[1]
    else
        model.rin
    end
    rout = if model.rout < rmin
        rmin
    else
        model.rout
    end
    Gradus.integrate_lineprofile!(
        output,
        model.table.setup,
        grid,
        domain;
        rmin = rmin,
        rmax = rout,
        g_scale = model.E₀,
    )
    output
end

export LineProfile

end # module GradusSpectralModels
