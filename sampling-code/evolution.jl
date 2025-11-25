include("init.jl")

function main()

    individuals = Individual[]
    n = 0

    while n < 200
        W = wsample([-1.0, 0.0, 1.0], Weights([1, 2, 1]), (N, N))
        individual = Individual(ones(N), W)
        population = [deepcopy(individual) for _ in 1:L]

        G = 1000
        fitness = Vector{Float64}[]
        for g in 1:G
            foreach(individual -> mutation!(individual, μ), population)

            Threads.@threads for n in 1:L
                callback_fitness!(population[n], x0, η, t)
            end

            f = sort([individual.fitness for individual in population], rev=true)
            println("$(n) $(g) $(f[1]) $(mean(f))")
            push!(fitness, f)
            probability_selection!(population, ν)
        end
        sort!(population, by=individual -> individual.fitness, rev=true)
        if population[1].fitness > -0.1
            push!(individuals, population[rand(findall(individual -> individual.fitness > -0.1, population))])
            n += 1
        end
    end

    jldsave("individuals_evolution.jld2"; individuals)
end

main()
