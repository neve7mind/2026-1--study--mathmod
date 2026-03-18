# # Лабораторная работа №3. Модель боевых действий
# **Вариант:** 19 (на основе предоставленных параметров)
#
# ## Инициализация проекта
using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots
default(fmt = :png)

# ## 1. Параметры модели
# Коэффициенты потерь и эффективности (на основе примера из методички):
# - a, h: потери, не связанные с боевыми действиями
# - b, c: эффективность боевых действий сторон
a = 0.4   # потери армии X
h = 0.7   # потери армии Y
b = 0.8   # эффективность армии
c = 0.5   # эффективность армии X

# Начальные условия:
x0 = 20000
y0 = 9000
u0 = [x0, y0]
tspan = (0.0, 1.0) # Временной интервал

# Функции подкрепления P(t) и Q(t)
P(t) = sin(t) + 1
Q(t) = cos(t) + 1

# Названия для сохранения
script_name = "lab03"
mkpath(plotsdir(script_name))

# ## 2. Определение моделей

# ### Модель №1: Регулярные войска против регулярных
# $\frac{dx}{dt} = -ax(t) - by(t) + P(t)$
# $\frac{dy}{dt} = -cx(t) - hy(t) + Q(t)$
function combat_regular!(du, u, p, t)
    x, y = u
    du[1] = -a * x - b * y + P(t)
    du[2] = -c * x - h * y + Q(t)
end

# ### Модель №2: Регулярные войска против партизан
# $\frac{dx}{dt} = -ax(t) - by(t) + P(t)$
# $\frac{dy}{dt} = -cx(t)y(t) - hy(t) + Q(t)$
function combat_mixed!(du, u, p, t)
    x, y = u
    du[1] = -a * x - b * y + P(t)
    du[2] = -c * x * y - h * y + Q(t)
end

# ## 3. Решение и визуализация

# ### Расчет для случая №1
prob1 = ODEProblem(combat_regular!, u0, tspan)
sol1 = solve(prob1, Tsit5(), saveat=0.01)

p1 = plot(sol1, title="Модель №1 (Регулярные vs Регулярные)",
          xlabel="Время", ylabel="Численность",
          label=["Армия X" "Армия Y"], lw=2)

# ### Расчет для случая №2
prob2 = ODEProblem(combat_mixed!, u0, tspan)
sol2 = solve(prob2, Tsit5(), saveat=0.01)

p2 = plot(sol2, title="Модель №2 (Регулярные vs Партизаны)",
          xlabel="Время", ylabel="Численность",
          label=["Армия X" "Армия Y"], lw=2)

# ## 4. Итоговые графики
final_plot = plot(p1, p2, layout=(1,2), size=(1000, 450))

# Сохранение результатов
savefig(final_plot, plotsdir(script_name, "combat_results.png"))
final_plot
