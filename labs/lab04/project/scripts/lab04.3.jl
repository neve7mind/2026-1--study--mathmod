# ## Инициализация проекта
using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots

gr()
default(fmt = :png, size = (800, 450))

# ## 1. Параметры (Вариант 19)
# w0^2 = 11.0, g = 0.0
w0 = sqrt(11.0)
g = 0.0
x0 = [1.1, 1.1] # Начальные условия
tspan = (0.0, 50.0)

# ## 2. Система уравнений
function oscillator_simple!(du, u, p, t)
    x, y = u
    du[1] = y
    du[2] = -w0^2 * x - 2 * g * y
end

# ## 3. Решение и графики
prob1 = ODEProblem(oscillator_simple!, x0, tspan)
sol1 = solve(prob1, Tsit5(), saveat=0.05)

p1_ts = plot(sol1, title="Решение: Без затухания (Вариант 19)", xlabel="t", label=["x(t)" "y(t)"])
p1_ph = plot(sol1, vars=(1, 2), title="Фазовый портрет", xlabel="x", ylabel="y", label="Траектория")
layout1 = plot(p1_ts, p1_ph, layout=(1,2))

display(layout1)
