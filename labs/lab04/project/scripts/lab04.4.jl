# ## Инициализация проекта
using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots

gr()
default(fmt = :png, size = (800, 450))

# ## 1. Параметры (Вариант 19)
# w0^2 = 11.0, g = 11.0
w0 = sqrt(11.0)
g = 11.0
x0 = [1.1, 1.1]
tspan = (0.0, 50.0)

# ## 2. Система уравнений
function oscillator_damped!(du, u, p, t)
    x, y = u
    du[1] = y
    du[2] = -w0^2 * x - 2 * g * y
end

# ## 3. Решение и графики
prob2 = ODEProblem(oscillator_damped!, x0, tspan)
sol2 = solve(prob2, Tsit5(), saveat=0.05)

p2_ts = plot(sol2, title="Решение: С затуханием (Вариант 19)", xlabel="t", label=["x(t)" "y(t)"])
p2_ph = plot(sol2, vars=(1, 2), title="Фазовый портрет", xlabel="x", ylabel="y", label="Траектория")
layout2 = plot(p2_ts, p2_ph, layout=(1,2))

display(layout2)
