include("init.jl")

function main()

    minf = -sqrt(M)
    maxf = 0.0

    binf = range(start=minf, stop=maxf, length=40)
    stepf = step(binf)

    indexf(f) = floor(Int, (f - minf) / stepf) + 1

    γ = ones(N)
    W = wsample([-1.0, 0.0, 1.0], Weights([1, 2, 1]), (N, N))
    individual = Individual(γ, W)
    callback_fitness!(individual, x0, η, t)

    lnf = zeros(length(binf))
    C = 1.0
    MCS = 2 * 10^6

    histogram = zeros(20, length(binf))
    for n in 1:20
        println(n)
        for mcs in 1:MCS
            individual_ = deepcopy(individual)
            mutation!(individual_, μ)
            callback_fitness!(individual_, x0, η, t)

            if rand(Uniform(0, 1)) < exp((lnf[indexf(individual_.fitness)] - lnf[indexf(individual.fitness)]))
                individual = individual_
            end
            lnf[indexf(individual.fitness)] -= C
            histogram[n, indexf(individual.fitness)] += 1
        end
        C /= 2.0
    end

    jldsave("learned.jld2"; binf, lnf, histogram)

end

main()
