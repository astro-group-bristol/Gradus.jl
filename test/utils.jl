function count_inner_boundary(m, simsols)
    points = unpack_solution.(m, simsols.u)
    count(i -> i.status == StatusCodes.WithinInnerBoundary, points)
end

function returntype(f, args...)
    T = map(typeof, args)
    info = code_typed(f, T)
    first(info).second
end
