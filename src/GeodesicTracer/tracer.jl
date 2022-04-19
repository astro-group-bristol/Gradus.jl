# single position and single velocity
function __tracegeodesics(
    m::AbstractMetricParams{T},
    init_pos::AbstractVector{T},
    init_vel::AbstractVector{T},
    time_domain::Tuple{T,T},
    solver;
    callback,
    kwargs...
) where {T}
    prob = integrator_problem(m, init_pos, init_vel, time_domain)
    cbs = create_callback_set(m, callback)
    solve_geodesic(solver, prob; callback = cbs, kwargs...)
end

# columnar of positions and velocity
function __tracegeodesics(
    m::AbstractMetricParams{T},
    init_positions::AbstractArray{T,2},
    init_velocities::AbstractArray{T,2},
    time_domain::Tuple{T,T},
    solver;
    ensemble = EnsembleThreads(),
    callback,
    kwargs...
) where {T}
    prob = integrator_problem(
        m,
        @view(init_positions[:, 1]),
        @view(init_velocities[:, 1]),
        time_domain
    )
    ens_prob = EnsembleProblem(
        prob,
        prob_func = (prob, i, repeat) -> integrator_problem(
            m,
            @view(init_positions[:, i]),
            @view(init_velocities[:, i]),
            time_domain
        ),
        safetycopy = false
    )

    cbs = create_callback_set(m, callback)
    solve_geodesic(
        solver,
        ens_prob,
        ensemble;
        trajectories = size(init_velocities, 2),
        callback = cbs,
        kwargs...
    )
end

# indexables
function __tracegeodesics(
    m::AbstractMetricParams{T},
    init_positions,
    init_velocities,
    time_domain::Tuple{T,T},
    solver;
    ensemble = EnsembleThreads(),
    callback,
    kwargs...
) where {T}
    prob = integrator_problem(m, init_positions[1], init_velocities[1], time_domain)
    ens_prob = EnsembleProblem(
        prob,
        prob_func = (prob, i, repeat) ->
            integrator_problem(m, init_positions[i], init_velocities[i], time_domain),
        safetycopy = false
    )

    cbs = create_callback_set(m, callback)
    solve_geodesic(
        solver,
        ens_prob,
        ensemble;
        trajectories = length(init_velocities),
        callback = cbs,
        kwargs...
    )
end

function solve_geodesic(solver, prob, ensemble; solver_opts...)
    solve(prob, solver, ensemble; solver_opts...)
end

function solve_geodesic(solver, prob; solver_opts...)
    solve(prob, solver, ; solver_opts...)
end
