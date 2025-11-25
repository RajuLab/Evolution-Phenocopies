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
    return 1 / (1 + exp(-β * z))
end

function grn!(dx, x, p, t)
    γ, W, η = p

    Wx = W * x

    @. dx = γ * (sigmoid(Wx / N^(1 / 2) + η) - x)

    nothing
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

const N = 40
const M = 3

const reference = load("shared/reference.jld2")["reference"]
const x̂ = stack(reference.u)
const x0 = x̂[:, 1]
const t̂ = reference.t
const t = extrema(t̂)

const η = zeros(N)

const L = 120
const μ = 2
const ν = 2.0