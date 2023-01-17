module Grids

function inverse_grid(min, max, N)
    Iterators.reverse(inv(x) for x in range(1 / max, 1 / min, N))
end

function geometric_grid(min, max, N)
    K = (max / min)^(1 / N)
    map(i -> min * K^(i - 1), 1:N)
end

end # module
