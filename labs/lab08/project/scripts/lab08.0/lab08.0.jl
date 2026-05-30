using Pkg
Pkg.activate("../project")
using DrWatson
@quickactivate "project"

using DifferentialEquations, Plots

const p_cr = 20.0   # Критическая стоимость продукта
const τ1 = 10.0     # Длительность производственного цикла фирмы 1
const p1 = 9.0      # Себестоимость продукта у фирмы 1
const τ2 = 16.0     # Длительность производственного цикла фирмы 2
const p2 = 7.0      # Себестоимость продукта у фирмы 2
const N = 10.0      # Число потребителей (в тыс. единиц)
const q = 1.0       # Максимальная потребность одного человека

const a1 = p_cr / (τ1^2 * p1^2 * N * q)
const a2 = p_cr / (τ2^2 * p2^2 * N * q)
const b  = p_cr / (τ1 * τ2 * p1 * p2 * N * q)
const c1 = (p_cr - p1) / (τ1 * p1)
const c2 = (p_cr - p2) / (τ2 * p2)

const M0 = [2.0, 1.0]
const tspan = (0.0, 30.0)

function competition_case1!(dM, M, p, t)
    dM[1] = (c1 / c1) * M[1] - (a1 / c1) * M[1]^2 - (b / c1) * M[1] * M[2]
    dM[2] = (c2 / c1) * M[2] - (a2 / c1) * M[2]^2 - (b / c1) * M[1] * M[2]
end

prob1 = ODEProblem(competition_case1!, M0, tspan)

sol1 = solve(prob1, RK4(), reltol=1e-8, abstol=1e-8)

plt1 = plot(sol1, vars=(0, 1), label="Фирма 1", color=:blue, lw=2)
plot!(plt1, sol1, vars=(0, 2), label="Фирма 2", color=:green, lw=2,
      title="Случай 1: Рыночная конкуренция", xlabel="Безразмерное время θ", ylabel="Оборотные средства M")

function competition_case2!(dM, M, p, t)
    dM[1] = (c1 / c1) * M[1] - (a1 / c1) * M[1]^2 - ((b / c1) + 0.002) * M[1] * M[2]
    dM[2] = (c2 / c1) * M[2] - (a2 / c1) * M[2]^2 - (b / c1) * M[1] * M[2]
end

prob2 = ODEProblem(competition_case2!, M0, tspan)

sol2 = solve(prob2, RK4(), reltol=1e-8, abstol=1e-8)

plt2 = plot(sol2, vars=(0, 1), label="Фирма 1 (с лобби)", color=:blue, lw=2)
plot!(plt2, sol2, vars=(0, 2), label="Фирма 2", color=:green, lw=2,
      title="Случай 2: Социально-психологический фактор", xlabel="Безразмерное время θ", ylabel="Оборотные средства M")

p_final = plot(plt1, plt2, layout=(2, 1), size=(800, 600))
display(p_final)

savefig(p_final, plotsdir("firms_competition.png"))
