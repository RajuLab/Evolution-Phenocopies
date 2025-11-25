include("init.jl")

function main()

    N = 40
    M = 3

    history = load("learned.jld2")
    reference = history["reference"]
    x0 = reference.u[1]
    t = extrema(reference.t)

    η = zeros(N)

    binf, lnf = history["binf"], history["lnf"]

    minf = minimum(binf)
    maxf = maximum(binf)
    stepf = step(binf)

    indexf(f) = floor(Int, (f - minf) / stepf) + 1


    Individuals = [Vector{Individual}() for _ in 1:5]

    Threads.@threads for individuals in Individuals
        W = wsample([-1.0, 0.0, 1.0], Weights([1, 2, 1]), (N, N))
        individual = Individual(ones(N), W)
        callback_fitness!(individual, x0, η, t)

        n = 0
        n_ = 0

        while n < 200

            individual_ = deepcopy(individual)
            mutation_individual!(individual_)
            callback_fitness!(individual_, x0, η, t)

            println("$(n) $(n_)")

            if rand(Uniform(0, 1)) < exp((lnf[indexf(individual_.fitness)] - lnf[indexf(individual.fitness)]))
                individual = individual_
            end

            if individual.fitness > -0.1
                if n_ == 120
                    push!(individuals, deepcopy(individual))
                    n += 1
                    n_ = 0
                end
                n_ += 1
            end
        end
    end
    individuals = dropdims(stack_t(Individuals), dims=2)
    jldsave("individuals_mcmc.jld2"; individuals)
end

main()
