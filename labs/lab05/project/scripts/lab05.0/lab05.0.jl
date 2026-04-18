using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots

gr()
default(fmt = :png, size = (800, 500))

a = 0.2    # Смертность хищников
c = 0.5    # Прирост жертв
b = 0.05   # Прирост хищников (эффективность охоты)
d = 0.02   # Смертность жертв (интенсивность охоты)

u0 = [5.0, 10.0]
tspan = (0.0, 400.0)

function lotka_volterra!(du, u, p, t)
    x, y = u
    du[1] = -a*x + b*x*y  # Хищники
    du[2] = c*y - d*x*y   # Жертвы
end

prob = ODEProblem(lotka_volterra!, u0, tspan)
sol = solve(prob, Tsit5(), saveat=0.1)

p1 = plot(sol,
          title = "Динамика популяций (Задание 1)",
          xlabel = "Время t",
          ylabel = "Численность",
          label = ["Хищники (x)" "Жертвы (y)"],
          lw = 1.5)

p2 = plot(sol, vars=(2, 1),
          title = "Фазовый портрет (Задание 1)",
          xlabel = "Численность жертв (y)",
          ylabel = "Численность хищников (x)",
          label = "Траектория",
          lw = 1.5)

scatter!([a/b], [c/d], color=:red, label="Стац. состояние A(25, 4)")

layout = plot(p1, p2, layout = (2, 1))
display(layout)

mkpath(plotsdir("lab05"))
savefig(layout, plotsdir("lab05", "predator_prey_results.png"))
