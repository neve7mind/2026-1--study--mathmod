# ## Инициализация проекта
using Pkg
Pkg.activate("../project")
using DrWatson
using DifferentialEquations
using Plots

# Настройки графики для корректного рендеринга в отчет
gr()
default(fmt = :png, size = (800, 450), titlefont = font(10))

# ## 1. Параметры Варианта 19
N = 12300            # Общая численность
I0 = 140             # Инфицированные
R0 = 54              # С иммунитетом
S0 = N - I0 - R0     # Восприимчивые (здоровые)
u0 = [S0, I0, R0]

alpha = 0.01         # Коэффициент заболеваемости
beta = 0.02          # Коэффициент выздоровления
tspan = (0.0, 60.0)

# ## 2. Определение функций системы

# Случай 1: I(t) <= I* (Инфекция не распространяется, больные изолированы)
function epidemic_case1!(du, u, p, t)
    S, I, R = u
    du[1] = 0
    du[2] = -beta * I
    du[3] = beta * I
end

# Случай 2: I(t) > I* (Инфекция свободно распространяется)
function epidemic_case2!(du, u, p, t)
    S, I, R = u
    du[1] = -alpha * S
    du[2] = alpha * S - beta * I
    du[3] = beta * I
end

# ## 3. Решение

# Моделирование первого случая
prob1 = ODEProblem(epidemic_case1!, u0, tspan)
sol1 = solve(prob1, Tsit5(), saveat=0.1)

# Моделирование второго случая
prob2 = ODEProblem(epidemic_case2!, u0, tspan)
sol2 = solve(prob2, Tsit5(), saveat=0.1)

# ## 4. Визуализация

# График для первого случая
p1 = plot(sol1,
          title="Эпидемия при I(t) ≤ I* (Вариант 19)",
          xlabel="Время t", ylabel="Численность",
          label=["S(t) - здоровые" "I(t) - больные" "R(t) - иммунитет"],
          lw=2, color=[:blue :red :green])

# График для второго случая
p2 = plot(sol2,
          title="Эпидемия при I(t) > I* (Вариант 19)",
          xlabel="Время t", ylabel="Численность",
          label=["S(t) - здоровые" "I(t) - больные" "R(t) - иммунитет"],
          lw=2, color=[:blue :red :green])

layout = plot(p1, p2, layout=(2,1), size=(800, 800))
display(layout)

# ## 5. Сохранение результатов
mkpath(plotsdir("lab06"))
savefig(layout, plotsdir("lab06", "variant19_epidemic.png"))
