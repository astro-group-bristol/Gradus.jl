module GradusDiffEqGPUExt

using Gradus
using Gradus: EnsembleProblem, TracingConfiguration, solve

using DiffEqGPU

const EnsembleGPU = Union{<:EnsembleGPUArray,<:EnsembleGPUKernel}

function Gradus.ensemble_solve_tracing_problem(
    ensemble::EnsembleGPU,
    problem::EnsembleProblem,
    config::TracingConfiguration;
    progress_bar = nothing,
    save_on = false, # capture to avoid passing
    solver_opts...,
)
    if !isnothing(progress_bar)
        @warn "Progress meter is currently not supported for GPU ensemble algorithms."
    end
    solve(
        problem,
        config.solver,
        ensemble,
        ;
        trajectories = config.trajectories,
        abstol = config.abstol,
        reltol = config.reltol,
        solver_opts...,
    )
end

end # module
