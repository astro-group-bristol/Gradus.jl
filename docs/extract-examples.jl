function add_julia_code!(examples::Vector{String}, src_file)
    lines = readlines(src_file)
    itt = Iterators.Stateful(enumerate(lines))
    for (lineno, line) in itt
        if line == "```julia"
            code = String["# $(src_file):$(lineno)"]
            while ((_, line) = popfirst!(itt))[2] != "```"
                push!(code, line)
            end
            push!(examples, join(code, "\n"))
            # check if we're supposed to output a picture or not
            (_, next_line) = peek(itt)
            img = match(r"!\[.*\]\((.*)\)\s*$", next_line)
            if !isnothing(img)
                filename = last(splitdir(first(img.captures)))
                savepath = joinpath(@__DIR__() * "/src/example-figures/", filename)
                save_line = "savefig(\"$savepath\")"
                push!(examples, save_line)
            end
        end
    end
end

file1 = @__DIR__() * "/src/examples.md"
file2 = @__DIR__() * "/src/getting-started.md"

examples = String[]
add_julia_code!(examples, file1)
add_julia_code!(examples, file2)

open(@__DIR__() * "/examples-code.jl", "w") do f
    for code in examples
        write(f, code)
        write(f, "\n", "#"^100, "\n\n")
    end
end
