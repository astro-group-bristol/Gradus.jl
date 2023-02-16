
@inline function _threaded_map(f, itr)
    N = length(itr)
    items = !(typeof(itr) <: AbstractArray) ? collect(itr) : itr
    output = Vector{Core.Compiler.return_type(f, Tuple{eltype(items)})}(undef, N)
    Threads.@threads for i = 1:N
        @inbounds output[i] = f(items[i])
    end
    output
end
