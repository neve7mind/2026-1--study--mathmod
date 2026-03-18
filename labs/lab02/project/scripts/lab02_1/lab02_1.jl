using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots
default(fmt = :png)

k = 10.0
n = 3.4
fi = 3*pi/4 # Заданный угол направления лодки

script_name = "lab02_task"
mkpath(plotsdir(script_name))

# Функция ДУ
function f(r, p, theta)
    return r / sqrt(n^2 - 1)
end

r0_1 = k / (n + 1)
tspan1 = (0.0, 2*pi)
prob1 = ODEProblem(f, r0_1, tspan1)
sol1 = solve(prob1, Tsit5(), saveat=0.01)

r0_2 = k / (n - 1)
tspan2 = (-pi, pi)
prob2 = ODEProblem(f, r0_2, tspan2)
sol2 = solve(prob2, Tsit5(), saveat=0.01)

# Данные для траектории лодки
theta_boat = [fi, fi]
r_boat = [0, 20]

# График 1
p1 = plot(sol1.t, sol1.u, proj=:polar, lims=(0,20), title="Сценарий 1", label="Катер", lw=2)
plot!(p1, theta_boat, r_boat, label="Лодка", linestyle=:dash, color=:red)

# График 2
p2 = plot(sol2.t, sol2.u, proj=:polar, lims=(0,20), title="Сценарий 2", label="Катер", lw=2)
plot!(p2, theta_boat, r_boat, label="Лодка", linestyle=:dash, color=:red)

final_plot = plot(p1, p2, layout=(1,2), size=(1000, 500))

r_meet1 = sol1(fi)
r_meet2 = sol2(fi)

println("Точка пересечения 1: r = ", round(r_meet1, digits=3), " км, θ = ", round(fi, digits=3))
println("Точка пересечения 2: r = ", round(r_meet2, digits=3), " км, θ = ", round(fi, digits=3))

savefig(final_plot, plotsdir(script_name, "pursuit.png"))
final_plot
