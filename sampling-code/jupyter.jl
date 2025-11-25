include("init.jl")
ENV["LINES"] = 7

using CairoMakie
cm = CairoMakie
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
        xticklabelsize=20, xlabelsize=24,
        yticklabelsize=20, ylabelsize=24,
    ),
    Axis3=(
        xticklabelsize=20, xlabelsize=24,
        yticklabelsize=20, ylabelsize=24,
        zticklabelsize=20, zlabelsize=24,
    ),
)

function plot(; aspect=nothing)
    if aspect != :equal
        f = Figure()
        ax = cm.Axis(f[1, 1])
        f, ax
    else
        f = Figure(size=(540, 540))
        ax = cm.Axis(f[1, 1], aspect=AxisAspect(1))
        f, ax
    end
end

function plot3d(; aspect=nothing)
    if aspect != :equal
        f = Figure()
        ax = cm.Axis3(f[1, 1])
        f, ax
    else
        f = Figure(size=(540, 540))
        ax = cm.Axis3(f[1, 1], aspect=AxisAspect(1))
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

return nothing

# include("init.jl")

# using Plots
# default(fontfamily="Computer Modern", titlefontsize=16, guidefontsize=18, tickfontsize=14, legendfontsize=12, lw=2.5, size=(640, 400), xticks=4, yticks=4, zticks=4, margin=3Plots.mm, legend=false, grid=false)
# ENV["LINES"] = 7

# ucolors = [
#     RGB(0 / 255, 90 / 255, 255 / 255)
#     RGB(255 / 255, 75 / 255, 0 / 255)
#     RGB(255 / 255, 241 / 255, 0 / 255)
#     RGB(3 / 255, 175 / 255, 122 / 255)
#     RGB(77 / 255, 196 / 255, 255 / 255)
#     RGB(255 / 255, 128 / 255, 130 / 255)
#     RGB(246 / 255, 170 / 255, 0 / 255)
#     RGB(153 / 255, 0 / 255, 153 / 255)
#     RGB(128 / 255, 64 / 255, 0 / 255)
#     RGB(0 / 255, 0 / 255, 0 / 255)
# ]

# function pshow(args...; kwargs...)
#     display(plot!(args...; kwargs...))
# end

# return nothing
