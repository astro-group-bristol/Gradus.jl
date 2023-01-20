function count_inner_boundary(m, simsols)
    points = getgeodesicpoint.(m, simsols.u)
    count(i -> i.status == StatusCodes.WithinInnerBoundary, points)
end
