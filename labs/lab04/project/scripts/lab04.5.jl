# ## Инициализация проекта
using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots

gr()
default(fmt = :png, size = (800, 450))

# ## 1. Параметры (Вариант 19)
# w0^2 = 11.0, g = 1.0
w0 = sqrt(11.0)
g = 1.0
x0 = [1.1, 1.1]
tspan = (0.0, 50.0)

# Внешняя сила f(t) = 11*sin(11*t)
f_ext(t) = 11.0 * sin(11.0 * t)

# ## 2. Система уравнений
function oscillator_forced!(du, u, p, t)
    x, y = u
    du[1] = y
    du[2] = -w0^2 * x - 2 * g * y + f_ext(t)
end

# ## 3. Решение и графики
prob3 = ODEProblem(oscillator_forced!, x0, tspan)
sol3 = solve(prob3, Tsit5(), saveat=0.05)

p3_ts = plot(sol3, title="Решение: Внешняя сила (Вариант 19)", xlabel="t", label=["x(t)" "y(t)"])
p3_ph = plot(sol3, vars=(1, 2), title="Фазовый портрет", xlabel="x", ylabel="y", label="Траектория")
layout3 = plot(p3_ts, p3_ph, layout=(1,2))

display(layout3)
