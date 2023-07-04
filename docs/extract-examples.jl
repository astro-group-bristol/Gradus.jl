lines = readlines(@__DIR__() * "/src/examples/examples.md")

examples = String[]

itt = Iterators.Stateful(lines)
for line in itt
    if line == "```julia"
        code = String[]
        while (line = popfirst!(itt)) != "```"
            push!(code, line)
        end
        push!(examples, join(code, "\n"))
    end
end

open(@__DIR__() * "/code/examples.jl", "w") do f
    for code in examples
        write(f, code)
        write(f, "\n", "#" ^ 100, "\n\n")
    end
end