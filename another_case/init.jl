using Random
using Distributions

using StatsBase

using LinearAlgebra
using DifferentialEquations, DiffEqCallbacks

using DelimitedFiles
using JLD2

stack_t(vector_of_vector) = stack(vector_of_vector, dims=1)

mutable struct Individual

    γ::Vector{Float64}
    W::Matrix{Float64}

    fitness::Float64

    Individual(γ, W) = new(γ, W, -Inf)
end

weights = Weights([1, 2, 1])
function mutation!(individual, mutation_size)
    idxs = rand(1:N, (2, mutation_size))
    for (i, j) in eachcol(idxs)
        individual.W[i, j] = sample([-1.0, 0.0, 1.0], weights)
    end
end

function probability_selection!(population, selection_rate)

    population_size = length(population)

    fitness_func = [exp(selection_rate * individual.fitness) for individual in population]
    if !all(fitness_func .== 0)
        fitness_func ./= sum(fitness_func)
    else
        fitness_func .= 1.0 / length(fitness_func)
    end

    selected = sample(1:population_size, Weights(fitness_func), population_size)

    population .= [deepcopy(individual) for individual in population[selected]]
end

function sigmoid(z)
    β = 40.0
    1 / (1 + exp(-β * z))
end

function grn!(dx, x, p, t)
    γ, W, η = p

    Wx = W * x
    @. dx = γ * (sigmoid(Wx / N^(1 / 2) + η) - x)

    nothing
end

function grn(x, p)
    γ, W, η = p

    Wx = W * x
    @. γ * (sigmoid(Wx / N^(1 / 2) + η) - x)
end

function sol_individual(individual::Individual, x0, η, t; kwargs...)
    p = (individual.γ, individual.W, η)

    prob = ODEProblem(grn!, x0, t, p)
    solve(prob, RK4(), saveat=0.05; kwargs...)
    # callback=TerminateSteadyState(1e-5, 1e-3))
end

function distance_t(x1, x2)
    mean(norm(diff_t) for diff_t in eachcol(x1 - x2))
end

function callback_fitness(individual, x0, η, t)
    sol = sol_individual(individual, x0, η, t)
    -distance_t(sol[1:M, :], x̂[1:M, :])
end

function callback_fitness!(individual, x0, η, t)
    individual.fitness = callback_fitness(individual, x0, η, t)
end

N = 40
M = 3

# reference = load("reference1.jld2")["reference"]
# reference = load("reference2.jld2")["reference"]
# reference = load("reference3.jld2")["reference"]
# reference = load("reference4.jld2")["reference"]
reference = load("reference5.jld2")["reference"]

x̂ = stack(reference.u)
if size(x̂, 1) == N
    x0 = x̂[:, 1]
else
    x0 = reference.x0
end
t̂ = reference.t
t = extrema(t̂)

η = zeros(N)

L = 120
μ = 2
ν = 2.0

function callback_perturbation(t′, Δ)
    function condition(x, t, integrator)
        t′ - t + 1e-5
    end

    function affect!(integrator)
        x = integrator.u
        @. x += rand(Uniform(-Δ, Δ))
    end

    ContinuousCallback(condition, affect!, save_positions=(false, false))
end

function environmental_copy(individual, t′, Δ)
    sol_individual(individual, x0, η, t, callback=callback_perturbation(t′, Δ))
end

function environmental_copies(individual, t′, Δ, L′)
    copies = Vector{ODESolution}(undef, L′)
    Threads.@threads for l in 1:L′
        copies[l] = environmental_copy(individual, t′, Δ)
    end
    copies
end

function mutant(individual, μ′)
    individual′ = deepcopy(individual)
    mutation!(individual′, μ′)
    individual′
end

function mutational_copy(individual, μ′)
    sol_individual(mutant(individual, μ′), x0, η, t)
end

function mutational_copies(individual, μ′, L′)
    copies = Vector{ODESolution}(undef, L′)
    Threads.@threads for l in 1:L′
        copies[l] = mutational_copy(individual, μ′)
    end
    copies
end
