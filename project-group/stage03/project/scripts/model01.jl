using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using LinearAlgebra

# --- 1. Параметры системы ---
# Мы инкапсулируем физические константы в словарь для удобства логирования и варьирования
params = Dict(
    :beta  => 10.0,    # Число Зельдовича (чувствительность скорости к температуре)
    :alpha => 0.8,     # Безразмерный температурный скачок
    :n     => 1.0,     # Порядок реакции
    :L     => 100.0,   # Размер расчетной области
    :dx    => 0.5,     # Шаг пространственной сетки
    :tspan => (0.0, 50.0) # Временной интервал моделирования
    )

# --- 2. Реакционный член w(θ) ---
# Функция описывает локальное выделение тепла согласно закону Аррениуса
function reaction_term(theta, p)
    beta, alpha, n = p[:beta], p[:alpha], p[:n]

    if theta <= 0 || theta >= 1 # Защита от отрицательных температур и значений выше температуры горения
        return 0.0
    end

    w = (beta^2 / 2) * (theta^n) * (1 - theta) * # Формула ЗФК из теоретической базы
        exp(-beta * (1 - theta) / (1 - alpha * (1 - theta)))
    return w
end

# --- 3. Дискретизация уравнения в частных производных ---
# Применяем метод прямых (Method of Lines) для преобразования УЧП в систему ОДУ
function zfk_dynamics!(dtheta, theta, p, t)
    dx = p[:dx]
    N = length(theta) # Внутренние узлы: диффузия (теплопроводность) + источник (реакция)

    for i in 1:N # Граничные условия Дирихле: θ(0, t) = 1 (горячая зона), θ(L, t) = 0 (холодная зона)

        left  = (i == 1) ? 1.0 : theta[i-1] # Аппроксимация второй производной (центральная разность)
        right = (i == N) ? 0.0 : theta[i+1]

        diffusion = (left - 2*theta[i] + right) / dx^2
        reaction  = reaction_term(theta[i], p)

        dtheta[i] = diffusion + reaction
    end
end

# --- 4. Инициализация и решение ---
function run_simulation(p)
    x = 0:p[:dx]:p[:L]
    N = length(x)

    theta_init = [xi < 5.0 ? 1.0 : 0.0 for xi in x] # Начальное условие: резкий перепад температуры (ступенька)

        prob = ODEProblem(zfk_dynamics!, theta_init, p[:tspan], p) # Формируем задачу ОДУ. Система становится жесткой при больших beta, поэтому используем неявный решатель Rosenbrock23

        println("Начало численного интегрирования для beta = $(p[:beta])...")
        sol = solve(prob, Rosenbrock23(), reltol=1e-6, abstol=1e-6)

        return x, sol
    end
    x, sol = run_simulation(params) # --- 5. Визуализация и сохранение ---

    plt = plot(title="Эволюция детерминированного фронта ЗФК", # Построение эволюции фронта во времени
               xlabel="Координата x", ylabel="Температура θ",
               legend=:topright,
               fmt = :png)

    for t in 0:10:params[:tspan][2]
        plot!(plt, x, sol(t), label="t = $t")
    end

    path = plotsdir("zfk_deterministic_front.png")
    safesave(path, plt)
    println("График сохранен в: $path")

    display(plt)
