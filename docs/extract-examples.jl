lines = readlines(@__DIR__() * "/src/examples.md")

examples = String[]

itt = Iterators.Stateful(lines)
for line in itt
    if line == "```julia"
        code = String[]
        while (line = popfirst!(itt)) != "```"
            push!(code, line)
        end
        push!(examples, join(code, "\n"))
        # check if we're supposed to output a picture or not
        img = match(r"!\[.*\]\((.*)\)\s*$", peek(itt))
        if !isnothing(img)
            filename = last(splitdir(first(img.captures)))
            savepath = joinpath(@__DIR__() * "/src/example-figures/", filename)
            save_line = "savefig(\"$savepath\")"
            push!(examples, save_line)
        end
    end
end

open(@__DIR__() * "/code/examples.jl", "w") do f
    for code in examples
        write(f, code)
        write(f, "\n", "#"^100, "\n\n")
    end
end
