using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots
default(fmt = :png)

gr()
default(fmt = :png, size = (800, 450))

w0 = 2.0        # Собственная частота
x0 = [0.0, 1.0] # Начальные условия: x(0)=0, y(0)=1
tspan = (0.0, 15.0)

function oscillator_simple!(du, u, p, t)
    x, y = u
    du[1] = y
    du[2] = -w0^2 * x
end

prob1 = ODEProblem(oscillator_simple!, x0, tspan)
sol1 = solve(prob1, Tsit5(), saveat=0.05)

p1_ts = plot(sol1, title="Решение: Без затухания", xlabel="t", label=["x(t)" "y(t)"])
p1_ph = plot(sol1, vars=(1, 2), title="Фазовый портрет", xlabel="x", ylabel="y", label="Траектория")
layout1 = plot(p1_ts, p1_ph, layout=(1,2))

display(layout1)

mkpath(plotsdir("lab04.0"))
savefig(layout1, plotsdir("lab04.0", "oscillator_simple.png"))
