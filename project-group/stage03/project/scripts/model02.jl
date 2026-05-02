using DrWatson
@quickactivate "project"

using DifferentialEquations
using Plots
using LinearAlgebra
using Trapz

# --- 1. Параметры системы ---
# Инкапсулируем физические константы в словарь для удобства логирования и варьирования
params = Dict(
    :D    => 1.0,      # Коэффициент диффузии
    :L    => 100.0,    # Размер расчётной области
    :dx   => 0.5,      # Шаг пространственной сетки
    :tspan => (0.0, 80.0)  # Временной интервал моделирования
    )

# --- 2. Реакционный член f(u) для модели КПП ---
# Логистический рост: f(u) = u * (1 - u)
function reaction_term(u, p)

    return u * (1 - u) # В модели КПП реакционный член не зависит от дополнительных параметров, но оставляем интерфейс единым для совместимости
end

# --- 3. Дискретизация уравнения в частных производных ---
# Применяем метод прямых (Method of Lines) для преобразования УЧП в систему ОДУ
# Уравнение: ∂u/∂t = D * ∇²u + f(u)
function kpp_dynamics!(du, u, p, t)
    D = p[:D]
    dx = p[:dx]
    N = length(u)

    for i in 1:N # Граничные условия Дирихле: u(0,t) = 1 (прореагировавшая зона), u(L,t) = 0 (несмесь)
        left  = (i == 1) ? 1.0 : u[i-1] # Аппроксимация второй производной (центральная разность)
        right = (i == N) ? 0.0 : u[i+1]

        diffusion = D * (left - 2*u[i] + right) / dx^2
        reaction  = reaction_term(u[i], p)

        du[i] = diffusion + reaction
    end
end

# --- 4. Функция для вычисления положения фронта X(t) = ∫(1 - u(x,t)) dx ---
function front_position(u_vec, x_vec)
    integrand = 1.0 .- u_vec # Интеграл по правилу трапеций
    return trapz(x_vec, integrand)
end

# --- 5. Функция для оценки скорости фронта по линейной аппроксимации ---
function estimate_front_speed(t_vals, X_vals, frac=0.6)
    n = length(t_vals) # Используем последние frac доли времени, где фронт движется равномерно
    idx_start = Int(floor(frac * n))
    t_linear = t_vals[idx_start:end]
    X_linear = X_vals[idx_start:end]

    A = hcat(t_linear, ones(length(t_linear))) # Линейная регрессия (метод наименьших квадратов)
    coeffs = A \ X_linear # Решаем систему: X = v * t + b
    return coeffs[1]  # коэффициент наклона (скорость)
end

# --- 6. Инициализация, решение и анализ ---
function run_simulation(p)
    x = 0:p[:dx]:p[:L]
    N = length(x)

    u_init = [xi < 20.0 ? 1.0 : 0.0 for xi in x] # Начальное условие: плавная ступенька (прореагировавшая область слева)

        for i in 1:N # Дополнительное сглаживание (необязательно, но улучшает устойчивость)
            if x[i] > 18.0 && x[i] < 22.0
                u_init[i] = 1.0 - (x[i] - 18.0) / 4.0
            end
        end


        prob = ODEProblem(kpp_dynamics!, u_init, p[:tspan], p)  # Формируем задачу ОДУ. Система КПП нежёсткая, но для единообразия используем Rosenbrock23

        println("Начало численного интегрирования модели КПП...")
        sol = solve(prob, Rosenbrock23(), reltol=1e-6, abstol=1e-6, saveat=range(p[:tspan][1], p[:tspan][2], length=200))

        X_vals = [front_position(sol.u[i], x) for i in 1:length(sol.t)] # --- Вычисление положения фронта и скорости ---
            v_num = estimate_front_speed(sol.t, X_vals)

            c_KPP = 2 * sqrt(p[:D] * 1.0) # Теоретическая скорость КПП: c_KPP = 2 * sqrt(D * f'(0))

            println("Теоретическая скорость КПП: $c_KPP") # Для f(u)=u(1-u) производная в нуле f'(0)=1
            println("Численная скорость фронта: $v_num")
            println("Относительная ошибка: ", abs(v_num - c_KPP) / c_KPP * 100, "%")

            return x, sol, X_vals, v_num, c_KPP
        end

        x, sol, X_vals, v_num, c_KPP = run_simulation(params) # --- 7. Визуализация и сохранение ---

        plt1 = plot(title="Эволюция детерминированного фронта КПП", # График 1: Эволюция профиля концентрации во времени
                    xlabel="Координата x", ylabel="Концентрация u",
                    legend=:topright,
                    fmt = :png)

        times_to_plot = range(params[:tspan][1], params[:tspan][2], length=6) # Выбираем несколько моментов времени для отображения
        for t in times_to_plot
            plot!(plt1, x, sol(t), label="t = $(round(t, digits=1))", lw=2)
        end

        n = length(sol.t)  # График 2: Положение фронта X(t) и линейная аппроксимация
        idx_start = Int(floor(0.6 * n)) # Вычисляем линейную аппроксимацию для последних 60% времени
        t_linear = sol.t[idx_start:end]
        X_linear = X_vals[idx_start:end]
        A = hcat(t_linear, ones(length(t_linear)))  # Регрессия для прямой на графике
        coeffs = A \ X_linear
        v_reg, b_reg = coeffs

        plt2 = plot(title="Положение фронта КПП",
                    xlabel="Время t", ylabel="Положение фронта X(t)",
                    legend=:topleft,
                    fmt = :png)
        plot!(plt2, sol.t, X_vals, label="Численный X(t)", lw=2, color=:blue)
        plot!(plt2, t_linear, X_linear, label="Линейная аппроксимация (v = $(round(v_num, digits=3)))",
              linestyle=:dash, lw=2, color=:red)
        plot!(plt2, sol.t, c_KPP .* sol.t .- (c_KPP * sol.t[1] - X_vals[1]),
              label="Теоретическая скорость c_KPP = $c_KPP", linestyle=:dot, lw=2, color=:green)

        plt = plot(plt1, plt2, layout=(2,1), size=(800, 600), fmt = :png) # Объединяем графики в один figure

        path = plotsdir("kpp_deterministic_front.png") # Сохранение результатов с использованием DrWatson
        safesave(path, plt)
        println("График сохранён в: $path")

        display(plt)
