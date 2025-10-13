struct TraceGeodesic{T} <: AbstractTrace
    μ::T
    q::T
    function TraceGeodesic(; μ = 0, q = 0)
        new{promote_type(typeof(μ), typeof(q))}(μ, q)
    end
end

struct TraceRadiativeTransfer{T} <: AbstractTrace
    μ::T
    q::T
    ν::T
    I₀::T
    function TraceRadiativeTransfer(; μ = 0, q = 0, ν = 1, I₀ = 1.0)
        new{promote_type(typeof(μ), typeof(q), typeof(ν), typeof(I₀))}(μ, q, ν, I₀)
    end
end

"""
    tracegeodesics(
        m::AbstractMetric,
        x, v, [disc], λ_domain;
        kwargs...
    )
    tracegeodesics(
        m::AbstractMetric,
        x, p::AbstractImagePlane, [disc], λ_domain;
        kwargs...
    )
    tracegeodesics(
        m::AbstractMetric,
        model::AbstractCoronaModel, [disc], λ_domain;
        kwargs...
    )

Trace a geodesic in the spacetime given by `m` with initial four-position `x` and four-velocity `v`.
- If both `x` and `v` are arrays of `SVector`, then many geodesics are traced in parallel.
- If `x` is an `SVector` and the third argument is of type [`AbstractImagePlane`](@ref), then every geodesic
of the image plane is traced in parallel.
- If `x` is an SVector and `v` is a function, then `trajectories` must be specified as a `kwarg`, and tracing is performed
in parallel.

The functional form for the velocity must map an index to a velocity

```julia
function velocity_function(i)
    return SVector(...)
end
```

The possible keyword arguments are

```julia
    chart = chart_for_metric(m),
    callback = nothing,
    solver = Tsit5(),
    ensemble = EnsembleThreads(),
    trajectories = nothing,
    abstol = 1e-9,
    reltol = 1e-9,
    solver_opts...,
```

with `solver_opts` being forwarded to the [SciML `solve` function](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/#solver_options).
"""
function tracegeodesics(
    m::AbstractMetric{T},
    args...;
    μ = zero(T),
    q = zero(T),
    trace = TraceGeodesic(; μ = μ, q = q),
    kwargs...,
) where {T}
    tracegeodesics(trace, m, args...; kwargs...)
end
function tracegeodesics(trace::AbstractTrace, m::AbstractMetric, args...; kwargs...)
    config, solver_opts = tracing_configuration(trace, m, args...; kwargs...)
    problem = assemble_tracing_problem(trace, config)
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
        kwargshandle = SciMLBase.KeywordArgError,
        verbose = config.verbose,
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
        ensemble,
        ;
        trajectories = config.trajectories,
        abstol = config.abstol,
        reltol = config.reltol,
        solver_opts...,
        # throw error if wrong keyword arguments passed
        kwargshandle = SciMLBase.KeywordArgError,
        verbose = config.verbose,
    )
end
# thread reusing dispatch
function ensemble_solve_tracing_problem(
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
    n_threads = Threads.nthreads()
    # init one integrator for each thread
    integrators = map(
        _ -> _init_integrator(
            pf(problem.prob, 1, 0),
            ;
            solver = config.solver,
            abstol = config.abstol,
            reltol = config.reltol,
            verbose = config.verbose,
            save_on = save_on,
            solver_opts...,
        ),
        1:n_threads,
    )
    # pre-allocate all of the returns
    output =
        Vector{Core.Compiler.return_type(_solve_reinit!, Tuple{eltype(integrators),S})}(
            undef,
            config.trajectories,
        )

    # solve
    Threads.@threads for i = 1:config.trajectories
        integ = integrators[_thread_id(n_threads)]
        p = pf(problem.prob, i, 0)
        output[i] = _solve_reinit!(integ, p.u0, p.p)
        # update progress bar
        if !isnothing(progress_bar)
            ProgressMeter.next!(progress_bar)
        end
    end
    output
end

function _init_integrator(
    problem;
    solver = Tsit5(),
    abstol = DEFAULT_TOLERANCE,
    reltol = DEFAULT_TOLERANCE,
    save_on = true,
    verbose = true,
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
        kwargshandle = SciMLBase.KeywordArgError,
        verbose = verbose,
    )
end

function _init_integrator(
    m::AbstractMetric{T},
    args...;
    μ = zero(T),
    q = zero(T),
    trace = TraceGeodesic(; μ = μ, q = q),
    kwargs...,
) where {T}
    config, solver_opts = tracing_configuration(trace, m, args...; kwargs...)
    problem = assemble_tracing_problem(trace, config)
    _init_integrator(problem; verbose = config.verbose, solver_opts...)
end


function _solve_reinit!(integrator, u0, p = nothing)
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
    else
        set_status_code!(integrator.p, StatusCodes.NoStatus)
    end
    sol = unpack_solution(solve!(integrator))
    return sol
end

export tracegeodesics
