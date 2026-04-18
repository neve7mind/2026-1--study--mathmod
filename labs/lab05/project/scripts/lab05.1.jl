# ## Лабораторная работа № 5. Модель хищник-жертва
using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots

# Настройки графики для отчета
gr()
default(fmt = :png, size = (800, 500), titlefont = font(10))

# ## 1. Параметры варианта 19
a = 0.19
b = 0.039
c = 0.39
d = 0.019

u0 = [7.0, 11.0] # Начальные численности [хищники, жертвы]
tspan = (0.0, 100.0)

# ## 2. Определение системы
function lotka_volterra!(du, u, p, t)
    x, y = u
    du[1] = -a*x + b*x*y  # Изменение числа хищников
    du[2] = c*y - d*x*y   # Изменение числа жертв
end

# ## 3. Решение
prob = ODEProblem(lotka_volterra!, u0, tspan)
sol = solve(prob, Tsit5(), saveat=0.1)

# ## 4. Визуализация
# График изменения численности со временем
p1 = plot(sol,
          title = "Динамика популяций (Вариант 19)",
          xlabel = "Время t",
          ylabel = "Численность",
          label = ["Хищники (x)" "Жертвы (y)"],
          lw = 2)

# Фазовый портрет
p2 = plot(sol, vars=(2, 1),
          title = "Фазовый портрет",
          xlabel = "Численность жертв (y)",
          ylabel = "Численность хищников (x)",
          label = "Траектория",
          lw = 2)

# Отметка стационарного состояния
scatter!([a/b], [c/d], color=:red, label="Стац. точка A(20.5, 4.9)")

layout = plot(p1, p2, layout = (2, 1))
display(layout)

# Сохранение для Quarto/Beamer
mkpath(plotsdir("lab05"))
savefig(layout, plotsdir("lab05", "variant19_results.png"))
