# ## Инициализация проекта
using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots

gr()
default(fmt = :png, size = (800, 450))

# ## 1. Параметры модели
w0 = 2.0
gamma_val = 0.4
x0 = [0.0, 1.0]
tspan = (0.0, 15.0)

# Внешняя сила f(t)
f_ext(t) = 2.0 * sin(2.0 * t)

# ## 2. Определение системы с внешней силой
function oscillator_forced!(du, u, p, t)
    x, y = u
    du[1] = y
    du[2] = -w0^2 * x - 2 * gamma_val * y + f_ext(t)
end

# ## 3. Решение и визуализация
prob3 = ODEProblem(oscillator_forced!, x0, tspan)
sol3 = solve(prob3, Tsit5(), saveat=0.05)

p3_ts = plot(sol3, title="Решение: Внешняя сила", xlabel="t", label=["x(t)" "y(t)"])
p3_ph = plot(sol3, vars=(1, 2), title="Фазовый портрет", xlabel="x", ylabel="y", label="Траектория")
layout3 = plot(p3_ts, p3_ph, layout=(1,2))

display(layout3)

# Сохранение
mkpath(plotsdir("lab04.2"))
savefig(layout3, plotsdir("lab04.2", "oscillator_forced.png"))
