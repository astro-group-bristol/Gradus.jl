using Gradus
using BenchmarkTools

# includes
include("integrator/benchmark-tracing.jl")

# build benchmark suite
suite = BenchmarkGroup()

suite["tracing"] = tracing_suite

# run and display
results = run(suite, verbose = true)
display(results)
