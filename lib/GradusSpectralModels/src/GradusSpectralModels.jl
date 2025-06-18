module GradusSpectralModels

using Gradus, SpectralFitting

struct LineProfile{D,T} <: AbstractTableModel{T,Additive}
    table::D
    K::T
    "Spin"
    a::T
    "Observer inclination (degrees off of the spin axis)."
    θ_obs::T
    "Inner radius of the accretion disc."
    inner_r::T
    "Outer radius of the accretion disc."
    outer_r::T
    "Central emission line energy (keV)."
    lineE::T
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
    a = FitParam(0.998, upper_limit = 1.0),
    θ_obs = FitParam(45.0, upper_limit = 90.0),
    inner_r = FitParam(1.0),
    outer_r = FitParam(100.0, upper_limit = 1000.0),
    lineE = FitParam(6.4),
    kwargs...,
)
    setup = integration_setup(profile, table((get_value(a), get_value(θ_obs))); kwargs...)
    LineProfile((; setup = setup, table = table), K, a, θ_obs, inner_r, outer_r, lineE)
end

function SpectralFitting.invoke!(output, domain, model::LineProfile)
    grid = model.table.table((model.a, model.θ_obs))
    rmin = (model.inner_r < grid.r_grid[1]) ? grid.r_grid[1] : model.inner_r
    outer_r = (model.outer_r < rmin) ? rmin : model.outer_r
    Gradus.integrate_lineprofile!(
        output,
        model.table.setup,
        grid,
        domain;
        rmin = rmin,
        rmax = outer_r,
        g_scale = model.lineE,
    )
    output
end

export LineProfile

end # module GradusSpectralModels
