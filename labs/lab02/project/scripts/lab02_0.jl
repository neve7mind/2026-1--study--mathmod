
# ## Инициализация проекта
using DrWatson
@quickactivate "project"
using DifferentialEquations
using Plots
default(fmt = :png)

## Определение путей для сохранения
script_name = "lab02"
mkpath(plotsdir(script_name))

# ## Постановка задачи
# На море в тумане катер береговой охраны преследует лодку браконьеров.
# Лодка обнаруживается на расстоянии $k = 11.4$ км
# Скорость катера в $n = 4.1$ раза больше скорости лодки
# Необходимо построить траекторию движения для двух случаев

# ## Параметры модели
k = 11.4
n = 4.1
fi = 3*pi/4 # направление движения лодки

# Расстояния, после которых катер начнет двигаться по спирали:
r0_1 = k / (n + 1) # Случай 1
r0_2 = k / (n - 1) # Случай 2

# ## Определение модели
# Дифференциальное уравнение траектории катера в полярных координатах:
# $dr/d\theta = r / \sqrt{n^2 - 1}$
function f(r, p, theta)
    return r / sqrt(n^2 - 1)
end

# ## Решение задачи

# ### Случай 1
# Начальные условия: $\theta_0 = 0, r_0 = k/(n+1)$
tspan1 = (0.0, 2*pi)
prob1 = ODEProblem(f, r0_1, tspan1)
sol1 = solve(prob1, Tsit5(), saveat=0.01)

# ### Случай 2
# Начальные условия: $\theta_0 = -\pi, r_0 = k/(n-1)$
tspan2 = (-pi, pi)
prob2 = ODEProblem(f, r0_2, tspan2)
sol2 = solve(prob2, Tsit5(), saveat=0.01)

# ## Визуализация

# Подготовка данных для траектории лодки (луч под углом fi):
theta_boat = [fi, fi]
r_boat = [0, 15]

# ### Построение графиков
# График для первого случая:
p1 = plot(sol1.t, sol1.u, proj=:polar, lims=(0,15),
          title="Случай 1", label="Катер", lw=2)
plot!(p1, theta_boat, r_boat, label="Лодка", linestyle=:dash, color=:red)

# График для второго случая:
p2 = plot(sol2.t, sol2.u, proj=:polar, lims=(0,15),
          title="Случай 2", label="Катер", lw=2)
plot!(p2, theta_boat, r_boat, label="Лодка", linestyle=:dash, color=:red)

# Объединение графиков:
final_plot = plot(p1, p2, layout=(1,2), size=(1000, 500))

# ## Сохранение результатов
savefig(final_plot, plotsdir(script_name, "pursuit_trajectories.png"))
