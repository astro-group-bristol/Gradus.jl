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

function LineProfile(
    profile,
    table::Gradus.CunninghamTransferTable;
    K = FitParam(1.0),
    a = FitParam(0.998),
    θ = FitParam(45.0),
    rin = FitParam(1.0),
    rout = FitParam(100.0, upper_limit = 100.0),
    E₀ = FitParam(1.0),
)
    setup = integration_setup(profile, table((get_value(θ), get_value(a))))
    LineProfile((; setup = setup, table = table), K, a, θ, rin, rout, E₀)
end

function SpectralFitting.invoke!(output, domain, model::LineProfile)
    grid = model.table.table((model.θ, model.a))
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
    )
    output
end

export LineProfile

end # module GradusSpectralModels
