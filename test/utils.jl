function count_inner_boundary(m, simsols)
    points = process_solution.(m, simsols.u)
    count(i -> i.status == StatusCodes.WithinInnerBoundary, points)
end
