using Pkg
Pkg.activate("../project")
using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots

const N = 1030      # Максимальное количество людей, которых может заинтересовать товар
const n0 = 11       # Количество людей, знающих о товаре в начальный момент времени
const tspan = (0.0, 30.0) # Временной промежуток моделирования

u0 = [Float64(n0)]

function advertisement_case1!(du, u, p, t)
    alpha1, alpha2, N_pop = p
    n = u[1]
    du[1] = (alpha1 + alpha2 * n) * (N_pop - n)
end

function advertisement_case2!(du, u, p, t)
    alpha1, alpha2, N_pop = p
    n = u[1]
    du[1] = (alpha1 + alpha2 * n) * (N_pop - n)
end

function advertisement_case3!(du, u, p, t)

    _, _, N_pop = p
    n = u[1]
    alpha1_t = 0.55 * sin(t)
    alpha2_t = 0.55 * cos(t)
    du[1] = (alpha1_t + alpha2_t * n) * (N_pop - n)
end

p1 = (0.55, 0.00005, N)
prob1 = ODEProblem(advertisement_case1!, u0, tspan, p1)
sol1 = solve(prob1, Tsit5(), saveat=0.1)

p2 = (0.00005, 0.55, N)
prob2 = ODEProblem(advertisement_case2!, u0, tspan, p2)

sol2 = solve(prob2, Rosenbrock23(), saveat=0.01)

p3 = (0.0, 0.0, N) # alpha1 и alpha2 вычисляются внутри, передаем только N
prob3 = ODEProblem(advertisement_case3!, u0, (0.0, 0.1), p3) # Короткий интервал из-за экстремального роста
sol3 = solve(prob3, Rosenbrock23(), saveat=0.0005)

default(titlefont=font(12, "Computer Modern"),
        guidefont=font(10, "Computer Modern"),
        tickfont=font(8, "Computer Modern"),
        legendfont=font(9, "Computer Modern"))

plt1 = plot(sol1, xlims=tspan, ylims=(0, N*1.05),
            title="Эффективность рекламы (Случай 1: a1 >> a2)",
            xlabel="Время t", ylabel="Число информированных клиентов n(t)",
            label="n(t)", color=:blue, lw=2, legend=:bottomright)
savefig(plt1, plotsdir("advertising_case1.png"))

t_max_speed = 0.0
for r in sol2.u
    if r[1] >= N/2
        idx = findfirst(x -> x == r, sol2.u)
        t_max_speed = sol2.t[idx]
        break
    end
end

plt2 = plot(sol2, xlims=(0.0, 0.1), ylims=(0, N*1.05),
            title="Эффективность рекламы (Случай 2: a1 << a2)",
            xlabel="Время t", ylabel="Число информированных клиентов n(t)",
            label="n(t)", color=:red, lw=2, legend=:bottomright)
vline!([t_max_speed], label="Max скорость (t ≈ $(round(t_max_speed, digits=4)))", color=:green, linestyle=:dash)
savefig(plt2, plotsdir("advertising_case2.png"))

plt3 = plot(sol3, ylims=(0, N*1.05),
            title="Эффективность рекламы (Случай 3: динамические коэф.)",
            xlabel="Время t", ylabel="Число информированных клиентов n(t)",
            label="n(t)", color=:purple, lw=2, legend=:bottomright)
savefig(plt3, plotsdir("advertising_case3.png"))

plot(plt1, plt2, plt3, layout=(3, 1), size=(800, 1000))
