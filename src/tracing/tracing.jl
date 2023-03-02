struct EnsembleEndpointThreads end

"""
    tracegeodesics(
        m::AbstractMetricParams,
        position,
        velocity::V,
        args...;
        solver = Tsit5(),
        ensemble = EnsembleThreads(),
        trajectories = nothing,
        μ = 0.0,
        chart = chart_for_metric(m),
        callback = nothing,
        solver_opts...,
    )

Integrate a geodesic for metric parameterised by `m`, for some initial positions and velocities.
The positions and velocities may be

  - a single position and velocity in the form of a vector of numbers,
  - a collection of positions and velocities, as either a vector of vectors, or as a matrix,
  - a single position and a function with signature ``vel_func(i)`` which returns a four-velocity.

The matrix specification reads each corresponding column as the initial position and velocity. When a collection of
positions and velocities is supplied, this method dispatched `EnsembleProblem`, offering `ensemble` as a `solver_opts`,
specifying the ensemble method to use.

`solver_opts` are the common solver options in DifferentialEquations.
"""
function tracegeodesics(
    m::AbstractMetricParams,
    position,
    velocity::V,
    args...;
    solver = Tsit5(),
    ensemble = EnsembleThreads(),
    trajectories = nothing,
    μ = 0.0,
    q = 0.0,
    chart = chart_for_metric(m),
    callback = nothing,
    solver_opts...,
) where {V}
    problem = geodesic_problem(
        m,
        position,
        velocity,
        args...;
        μ = μ,
        q = q,
        chart = chart,
        callback = callback,
    )
    # how many trajectories
    if V <: Function
        if isnothing(trajectories)
            error("When velocity is a function, trajectories must be defined")
        end
        solve_geodesic_problem(
            problem,
            solver,
            restric_ensemble(m, ensemble),
            trajectories;
            solver_opts...,
        )
    elseif eltype(velocity) <: SVector
        solve_geodesic_problem(
            problem,
            solver,
            restric_ensemble(m, ensemble),
            length(velocity);
            solver_opts...,
        )
    else
        if !isnothing(trajectories)
            error(
                "Trajectories should be `nothing` when solving only a single geodesic problem.",
            )
        end
        solve_geodesic_problem(problem, solver; solver_opts...)
    end
end

# ensemble dispatch
@inline function solve_geodesic_problem(
    ens_prob::EnsembleProblem,
    solver,
    ensemble,
    trajectories::Int;
    abstol = 1e-9,
    reltol = 1e-9,
    progress_bar = nothing,
    solver_opts...,
)
    prob = if !isnothing(progress_bar)
        function output_bump_progress(sol, i)
            ProgressMeter.next!(progress_bar)
            (sol, false)
        end
        # bootstrap incrementing progress bar
        EnsembleProblem(
            ens_prob.prob;
            prob_func = ens_prob.prob_func,
            output_func = output_bump_progress,
            safetycopy = ens_prob.safetycopy,
        )
    else
        ens_prob
    end
    ensol = solve(
        prob,
        solver,
        ensemble;
        trajectories = trajectories,
        abstol = abstol,
        reltol = reltol,
        solver_opts...,
        kwargshandle = KeywordArgError,
    )
    ensol
end

# non-ensemble dispatch
@inline function solve_geodesic_problem(
    prob::ODEProblem,
    solver;
    abstol = 1e-9,
    reltol = 1e-9,
    progress_bar = nothing,
    solver_opts...,
)
    if !isnothing(progress_bar)
        @warn "Ignoring progress bar for single geodesic"
    end
    sol = solve(
        prob,
        solver,
        ;
        abstol = abstol,
        reltol = reltol,
        solver_opts...,
        kwargshandle = KeywordArgError,
    )
    sol
end

# thread reusing
@inline function solve_geodesic_problem(
    prob::EnsembleProblem{<:ODEProblem{S}},
    solver,
    ensemble::EnsembleEndpointThreads,
    trajectories::Int;
    save_on = false,
    progress_bar = nothing,
    solver_opts...,
) where {S}
    if save_on
        error("Cannot use `EnsembleEndpointThreads` with `save_on`")
    end
    pf = prob.prob_func
    # init one integrator for each thread
    integrators = map(
        i -> _init_integrator(
            pf(prob.prob, 1, 0),
            ;
            solver = solver,
            save_on = save_on,
            solver_opts...,
        ),
        1:Threads.nthreads(),
    )
    # pre-allocate all of the returns
    output =
        Vector{Core.Compiler.return_type(_solve_reinit!, Tuple{eltype(integrators),S})}(
            undef,
            trajectories,
        )

    # solve
    Threads.@threads for i = 1:trajectories
        integ = integrators[Threads.threadid()]
        p = pf(prob.prob, i, 0)
        output[i] = _solve_reinit!(integ, p.u0, p.p)
        # update progress bar 
        if !isnothing(progress_bar)
            ProgressMeter.next!(progress_bar)
        end
    end
    output
end

@inline function _init_integrator(
    problem;
    solver = Tsit5(),
    save_on = true,
    abstol = 1e-9,
    reltol = 1e-9,
    solver_opts...,
)
    SciMLBase.init(
        problem,
        solver;
        save_on = save_on,
        abstol = abstol,
        reltol = reltol,
        solver_opts...,
        kwargshandle = KeywordArgError,
    )
end

@inline function _init_integrator(
    m,
    u,
    v,
    args...;
    chart = chart_for_metric(m),
    solver_opts...,
)
    prob = geodesic_problem(m, u, v, args...; chart = chart)
    _init_integrator(prob; solver_opts...)
end

@inline function _solve_reinit!(integrator, u0, p = nothing)
    reinit!(
        integrator,
        u0,
        reset_dt = true,
        reinit_callbacks = true,
        erase_sol = true,
        reinit_cache = true,
    )
    auto_dt_reset!(integrator)
    # if new parameters have been passed, update them
    if !isnothing(p)
        integrator.p = update_integration_parameters!(integrator.sol.prob.p, p)
    end
    process_solution(solve!(integrator))
end

export tracegeodesics, geodesic_problem
