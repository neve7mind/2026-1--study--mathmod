# ## Инициализация проекта
using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots
default(fmt = :png)

gr()
default(fmt = :png, size = (800, 450))

# ## 1. Параметры модели
w0 = 2.0
gamma_val = 0.4 # Характеризует потери энергии
x0 = [0.0, 1.0]
tspan = (0.0, 15.0)

# ## 2. Определение системы с затуханием
function oscillator_damped!(du, u, p, t)
    x, y = u
    du[1] = y
    du[2] = -w0^2 * x - 2 * gamma_val * y
end

# ## 3. Решение и визуализация
prob2 = ODEProblem(oscillator_damped!, x0, tspan)
sol2 = solve(prob2, Tsit5(), saveat=0.05)

p2_ts = plot(sol2, title="Решение: С затуханием", xlabel="t", label=["x(t)" "y(t)"])
p2_ph = plot(sol2, vars=(1, 2), title="Фазовый портрет", xlabel="x", ylabel="y", label="Траектория")
layout2 = plot(p2_ts, p2_ph, layout=(1,2))

display(layout2)

# Сохранение
mkpath(plotsdir("lab04.1"))
savefig(layout2, plotsdir("lab04.1", "oscillator_damped.png"))
