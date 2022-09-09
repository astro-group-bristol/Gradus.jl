# single position and single velocity
function __tracegeodesics(
    m::AbstractMetricParams{T},
    init_pos::AbstractVector{T},
    init_vel::AbstractVector{T},
    time_domain::Tuple{T,T},
    solver;
    kwargs...,
) where {T<:Number}
    prob = integrator_problem(m, init_pos, init_vel, time_domain)
    solve_geodesic(solver, prob; kwargs...)
end

# single position and single velocity
# ensembled with an indexable function
function __tracegeodesics(
    m::AbstractMetricParams{T},
    init_pos::AbstractVector{T},
    vel_func::Function,
    time_domain::Tuple{T,T},
    solver;
    trajectories::Int,
    ensemble = EnsembleThreads(),
    kwargs...,
) where {T<:Number}
    prob = integrator_problem(m, init_pos, vel_func(1), time_domain)
    ens_prob = EnsembleProblem(
        prob,
        prob_func = (prob, i, repeat) ->
            integrator_problem(m, init_pos, vel_func(i), time_domain),
        safetycopy = false,
    )
    solve_geodesic(solver, ens_prob, ensemble; trajectories = trajectories, kwargs...)
end

# indexables
function __tracegeodesics(
    m::AbstractMetricParams{T},
    init_positions,
    init_velocities,
    time_domain::Tuple{T,T},
    solver;
    ensemble = EnsembleThreads(),
    kwargs...,
) where {T}
    prob = integrator_problem(m, init_positions[1], init_velocities[1], time_domain)
    ens_prob = EnsembleProblem(
        prob,
        prob_func = (prob, i, repeat) ->
            integrator_problem(m, init_positions[i], init_velocities[i], time_domain),
        safetycopy = false,
    )

    solve_geodesic(
        solver,
        ens_prob,
        ensemble;
        trajectories = length(init_velocities),
        kwargs...,
    )
end

function solve_geodesic(solver, prob, ensemble; solver_opts...)
    solve(prob, solver, ensemble; solver_opts..., kwargshandle = KeywordArgError)
end

function solve_geodesic(solver, prob; solver_opts...)
    solve(prob, solver, ; solver_opts..., kwargshandle = KeywordArgError)
end
