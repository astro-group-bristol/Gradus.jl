function _make_interpolation(x, y)
    DataInterpolations.LinearInterpolation(y, x)
end

@inline function _threaded_map(f, itr)
    N = length(itr)
    items = !(typeof(itr) <: AbstractArray) ? collect(itr) : itr
    output = Vector{Core.Compiler.return_type(f, Tuple{eltype(items)})}(undef, N)
    Threads.@threads for i = 1:N
        @inbounds output[i] = f(items[i])
    end
    output
end

function spherical_to_cartesian(v)
    x = v[1] * cos(v[3]) * sin(v[2])
    y = v[1] * sin(v[3]) * sin(v[2])
    z = v[1] * cos(v[2])
    SVector{3}(x, y, z)
end

# specialisation for four-vector
spherical_to_cartesian(v::SVector{4}) = spherical_to_cartesian(@views(v[2:end]))

function cartesian_squared_distance(::AbstractMetric, x1, x2)
    # all metrics are currently in boyer lindquist coords
    y1 = spherical_to_cartesian(x1)
    y2 = spherical_to_cartesian(x2)
    diff = @. (y2 - y1)^2
    sum(diff)
end

cartesian_distance(m::AbstractMetric, x1, x2) = √(cartesian_squared_distance(m, x1, x2))

export cartesian_squared_distance, cartesian_distance, spherical_to_cartesian
