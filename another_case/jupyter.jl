include("init.jl")
ENV["LINES"] = 7

using CairoMakie
const cm = CairoMakie
# using GLMakie
# const cm = GLMakie
mt_fonts_dir = joinpath(dirname(pathof(Makie.MathTeXEngine)), "..", "assets", "fonts", "NewComputerModern")
set_theme!(theme_minimal(),
    size=(640, 450),
    fonts=(
        regular=joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
        bold=joinpath(mt_fonts_dir, "NewCM10-Bold.otf")),
    Lines=(
        linewidth=2,
    ),
    Axis=(
        xticklabelsize=24, xlabelsize=28,
        yticklabelsize=24, ylabelsize=28,
    ),
    Axis3=(
        xticklabelsize=24, xlabelsize=28,
        yticklabelsize=24, ylabelsize=28,
        zticklabelsize=24, zlabelsize=28,
    ),
)

function plot(; aspect=nothing, size=nothing, kwargs...)
    if aspect != :equal
        if isnothing(size)
            f = Figure()
        else
            f = Figure(size=size)
        end
        ax = cm.Axis(f[1, 1]; kwargs...)
        f, ax
    else
        if isnothing(size)
            f = Figure(size=(540, 540))
        else
            f = Figure(size=size)
        end
        ax = cm.Axis(f[1, 1], aspect=AxisAspect(1); kwargs...)
        f, ax
    end
end

function plot3(; aspect=nothing, size=nothing, kwargs...)
    if aspect != :equal
        if isnothing(size)
            f = Figure()
        else
            f = Figure(size=size)
        end
        ax = cm.Axis3(f[1, 1]; kwargs...)
        f, ax
    else
        if isnothing(size)
            f = Figure(size=(540, 540))
        else
            f = Figure(size=size)
        end
        ax = cm.Axis3(f[1, 1], aspect=(1, 1, 1); kwargs...)
        f, ax
    end
end

function show(f)
    display(f)
    nothing
end

using ColorSchemes: tab10
using ColorTypes: RGB
ucolors = [
    RGB(0 / 255, 90 / 255, 255 / 255)
    RGB(255 / 255, 75 / 255, 0 / 255)
    RGB(255 / 255, 241 / 255, 0 / 255)
    RGB(3 / 255, 175 / 255, 122 / 255)
    RGB(77 / 255, 196 / 255, 255 / 255)
    RGB(255 / 255, 128 / 255, 130 / 255)
    RGB(246 / 255, 170 / 255, 0 / 255)
    RGB(153 / 255, 0 / 255, 153 / 255)
    RGB(128 / 255, 64 / 255, 0 / 255)
    RGB(0 / 255, 0 / 255, 0 / 255)
]

using LaTeXStrings
using MultivariateStats

function grn(x, individual, η)
    γ, W = individual.γ, individual.W
    Wx = W * x
    @. γ * (sigmoid(Wx / N^(1 / 2) + η) - x)
end

set_theme!(theme_minimal(),
    size=(640, 450),
    fonts=(
        regular=joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
        bold=joinpath(mt_fonts_dir, "NewCM10-Bold.otf")),
    Lines=(
        linewidth=2,
    ),
    Axis=(
        xticklabelsize=24, xlabelsize=28,
        yticklabelsize=24, ylabelsize=28,
    ),
    Axis3=(
        xticklabelsize=24, xlabelsize=28,
        yticklabelsize=24, ylabelsize=28,
        zticklabelsize=24, zlabelsize=28,
    ),
)

function covariance_ellipse(data1, data2; coeff=1.0)
    X = stack([data1, data2])
    covariance_X = cov(X)
    eigen_covariance_X = eigen(covariance_X)

    a1, a2 = coeff * sqrt.(eigen_covariance_X.values)
    vectors = eigen_covariance_X.vectors

    z = stack([vectors * [a1 * cos(t), a2 * sin(t)] for t in 0:0.01:2π])

    z[1, :] .+ mean(data1), z[2, :] .+ mean(data2)
end

return nothing

# dats = load("shared/mds_ingredients.jld2")
# D, tr = dats["D"], dats["trajectories"]
# K1, K2, K3, K4 = dats["K1"], dats["K2"], dats["K3"], dats["K4"]
# mds = fit(MDS, D; distances=true, maxoutdim=20)
# X = predict(mds)
# tr1 = tr[K1+230]
# tr2 = tr[73]
# tr3 = tr[K1+5]
# tr4 = tr[K3+203]
# jldsave("mds1.jld2"; evolution=(X[1:2, 1:K2], K1), mcmc=(X[1:2, K2+1:K4], K3 - K2), tr1, tr2, tr3, tr4)

# D = pairwise(Euclidean(), vec([sol[:, end] for sol in sols]))
# h = hclust(D)
# ho = h.order
# basins_ = cutree(h, k=2)
# basins = reshape(basins_, length(z1s), length(z2s))