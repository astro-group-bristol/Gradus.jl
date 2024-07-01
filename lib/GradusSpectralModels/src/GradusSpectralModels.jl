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
        table.setup.h,
        table.setup.time,
        table.setup.integrand,
        table.setup.pure_radial,
        table.setup.quadrature_rule,
        deepcopy(table.index_cache),
        table.setup.g_grid_upscale,
        table.setup.n_radii,
    )
    typeof(m)(
        (; setup = setup, table = table),
        (copy(getproperty(m, f)) for f in fieldnames(typeof(f))[2:end])...,
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
    Gradus.integrate_lineprofile!(
        output,
        model.table.setup,
        grid,
        domain;
        rmin = rmin,
        rmax = model.rout,
        g_scale = model.E₀,
    )
    output
end

export LineProfile

end # module GradusSpectralModels
