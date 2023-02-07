
@inline function _threaded_map(f, itr)::Vector{Base.@default_eltype(Base.Generator(f, itr))}
    # see https://github.com/JuliaFolds/Transducers.jl/issues/13
    Transducers.tcollect(itr |> Transducers.Map(f))
end