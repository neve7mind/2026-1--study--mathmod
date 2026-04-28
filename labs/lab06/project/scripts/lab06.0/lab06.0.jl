using Pkg
Pkg.activate("../project")
using DrWatson
using DifferentialEquations
using Plots

gr()
default(fmt = :png, size = (800, 500), titlefont = font(10))

N = 12000            # Общая численность популяции
I0 = 120             # Число заболевших в начальный момент
R0 = 20              # Здоровые с иммунитетом
S0 = N - I0 - R0     # Восприимчивые особи
u0 = [S0, I0, R0]

alpha = 0.01         # Коэффициент заболеваемости
beta = 0.02          # Коэффициент выздоровления
tspan = (0.0, 100.0)

function epidemic_case1!(du, u, p, t)
    S, I, R = u
    du[1] = 0
    du[2] = -beta * I
    du[3] = beta * I
end

function epidemic_case2!(du, u, p, t)
    S, I, R = u
    du[1] = -alpha * S
    du[2] = alpha * S - beta * I
    du[3] = beta * I
end

prob1 = ODEProblem(epidemic_case1!, u0, tspan)
sol1 = solve(prob1, Tsit5(), saveat=0.1)

prob2 = ODEProblem(epidemic_case2!, u0, tspan)
sol2 = solve(prob2, Tsit5(), saveat=0.1)

p1 = plot(sol1, title="Эпидемия: случай I(t) <= I*",
          xlabel="t", ylabel="Численность",
          label=["S(t)" "I(t)" "R(t)"], lw=2)

p2 = plot(sol2, title="Эпидемия: случай I(t) > I*",
          xlabel="t", ylabel="Численность",
          label=["S(t)" "I(t)" "R(t)"], lw=2)

layout = plot(p1, p2, layout=(2,1))
display(layout)

mkpath(plotsdir("lab06"))
savefig(layout, plotsdir("lab06", "epidemic_results.png"))
