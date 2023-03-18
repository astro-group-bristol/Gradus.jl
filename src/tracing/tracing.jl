struct TraceGeodesic{T} <: AbstractTraceParameters
    μ::T
    q::T
    TraceGeodesic(; μ = 0.0, q = 0.0) = new{typeof(μ)}(μ, q)
end

struct TraceRadiativeTransfer{T} <: AbstractTraceParameters
    μ::T
    q::T
    ν::T
    TraceRadiativeTransfer(; μ = 0.0, q = 0.0, ν = 1.0) = new{typeof(μ)}(μ, q, ν)
end

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
function tracegeodesics(m::AbstractMetricParams, args...; μ = 0.0, q = 0.0, kwargs...)
    config, solver_opts = tracing_configuration(m, args...; kwargs...)
    problem = assemble_tracing_problem(TraceGeodesic(; μ = μ, q = q), config)
    solve_tracing_problem(problem, config; solver_opts...)
end

# ensemble dispatch
solve_tracing_problem(
    problem::EnsembleProblem,
    config::TracingConfiguration;
    solver_opts...,
) = ensemble_solve_tracing_problem(config.ensemble, problem, config; solver_opts...)
# non-ensemble dispatch
@inline function solve_tracing_problem(
    prob::ODEProblem,
    config::TracingConfiguration;
    progress_bar = nothing,
    solver_opts...,
)
    if !isnothing(progress_bar)
        @warn "Ignoring progress bar for single geodesic"
    end
    sol = solve(
        prob,
        config.solver,
        ;
        abstol = config.abstol,
        reltol = config.reltol,
        solver_opts...,
        # throw error if wrong keyword arguments passed
        kwargshandle = KeywordArgError,
    )
    sol
end

# ensemble dispatch
function ensemble_solve_tracing_problem(
    ensemble,
    problem::EnsembleProblem,
    config::TracingConfiguration;
    progress_bar = nothing,
    solver_opts...,
)
    prob = if !isnothing(progress_bar)
        # bootstrap incrementing progress bar
        function _output_bump_progress(sol, i)
            ProgressMeter.next!(progress_bar)
            (sol, false)
        end
        # rebuild with new output function
        EnsembleProblem(
            problem.prob;
            prob_func = problem.prob_func,
            output_func = _output_bump_progress,
            safetycopy = problem.safetycopy,
        )
    else
        problem
    end
    solve(
        prob,
        config.solver,
        config.ensemble,
        ;
        trajectories = config.trajectories,
        abstol = config.abstol,
        reltol = config.reltol,
        solver_opts...,
        # throw error if wrong keyword arguments passed
        kwargshandle = KeywordArgError,
    )
end
# thread reusing dispatch
@inline function ensemble_solve_tracing_problem(
    ensemble::EnsembleEndpointThreads,
    problem::EnsembleProblem{<:ODEProblem{S}},
    config::TracingConfiguration;
    progress_bar = nothing,
    save_on = false,
    solver_opts...,
) where {S}
    if save_on
        error("Cannot use `EnsembleEndpointThreads` with `save_on`")
    end
    pf = problem.prob_func
    # init one integrator for each thread
    integrators = map(
        i -> _init_integrator(
            pf(problem.prob, 1, 0),
            ;
            solver = config.solver,
            abstol = config.abstol,
            reltol = config.reltol,
            save_on = save_on,
            solver_opts...,
        ),
        1:Threads.nthreads(),
    )
    # pre-allocate all of the returns
    output =
        Vector{Core.Compiler.return_type(_solve_reinit!, Tuple{eltype(integrators),S})}(
            undef,
            config.trajectories,
        )

    # solve
    Threads.@threads for i = 1:config.trajectories
        integ = integrators[Threads.threadid()]
        p = pf(problem.prob, i, 0)
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
    abstol = 1e-9,
    reltol = 1e-9,
    save_on = true,
    solver_opts...,
)
    SciMLBase.init(
        problem,
        solver;
        abstol = abstol,
        reltol = reltol,
        save_on = save_on,
        solver_opts...,
        # throw error if wrong keyword arguments passed
        kwargshandle = KeywordArgError,
    )
end

_init_integrator(m::AbstractMetricParams, args...; trace = TraceGeodesic(), kwargs...) =
    _init_integrator(trace, m, args...; kwargs...)

@inline function _init_integrator(
    trace::AbstractTraceParameters,
    m::AbstractMetricParams,
    args...;
    kwargs...,
)
    config, solver_opts = tracing_configuration(m, args...; kwargs...)
    problem = assemble_tracing_problem(trace, config)
    _init_integrator(problem; solver_opts...)
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

export tracegeodesics
