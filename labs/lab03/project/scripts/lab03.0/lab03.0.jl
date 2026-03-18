using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots
default(fmt = :png)

a = 0.35   # потери X, не связанные с боем
h = 0.45   # потери Y, не связанные с боем
b = 0.75   # эффективность войск Y
c = 0.55   # эффективность войск X

x0 = 35000
y0 = 24000
u0 = [x0, y0]
tspan = (0.0, 1.5)

P(t) = sin(t) + 1.5
Q(t) = cos(t) + 1.2

function combat_regular!(du, u, p, t)
    x, y = u
    du[1] = -a * x - b * y + P(t)
    du[2] = -c * x - h * y + Q(t)
end

function combat_mixed!(du, u, p, t)
    x, y = u
    du[1] = -a * x - b * y + P(t)
    du[2] = -c * x * y - h * y + Q(t) # партизаны Y более скрытны
end

prob1 = ODEProblem(combat_regular!, u0, tspan)
sol1 = solve(prob1, Tsit5(), saveat=0.01)

prob2 = ODEProblem(combat_mixed!, u0, tspan)
sol2 = solve(prob2, Tsit5(), saveat=0.01)

script_name = "lab03"
mkpath(plotsdir(script_name))

p1 = plot(sol1, title="Регулярные войска vs Регулярные",
          xlabel="Время", ylabel="Численность", label=["Армия X" "Армия Y"], lw=2)

p2 = plot(sol2, title="Регулярные войска vs Партизаны",
          xlabel="Время", ylabel="Численность", label=["Армия X" "Армия Y"], lw=2)

final_plot = plot(p1, p2, layout=(1,2), size=(1000, 450))
savefig(final_plot, plotsdir(script_name, "combat_models.png"))

final_plot
