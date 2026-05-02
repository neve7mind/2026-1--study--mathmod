using DrWatson
@quickactivate "project"

using DifferentialEquations
using StochasticDiffEq
using Plots
using Statistics
using LinearAlgebra
using Random

params = Dict(
    :γ => 1.0,             # коэффициент трения
    :D => 0.5,             # коэффициент диффузии
    :x0 => 2.0,            # начальное положение частицы
    :tspan => (0.0, 10.0), # временной интервал моделирования
    :trajectories => 200   # количество случайных траекторий
    )

function drift!(du, u, p, t)
    γ = p[:γ]
    du[1] = -γ * u[1]
end

function diffusion!(du, u, p, t)
    D = p[:D]
    du[1] = sqrt(2D)
end

function theoretical_mean(t, p)
    return p[:x0] * exp(-p[:γ] * t)
end

function theoretical_variance(t, p)
    γ = p[:γ]
    D = p[:D]
    return (D / γ) * (1 - exp(-2γ * t))
end

function run_single_trajectory(p)
    u0 = [p[:x0]]
    prob = SDEProblem(drift!, diffusion!, u0, p[:tspan], p)

    println("Начало моделирования одной траектории уравнения Ланжевена...")

    sol = solve(
        prob,
        EM(),
        dt = 0.001,
        saveat = 0.01
        )

    return sol
end

function run_ensemble(p)
    u0 = [p[:x0]]
    prob = SDEProblem(drift!, diffusion!, u0, p[:tspan], p)

    ensemble_prob = EnsembleProblem(prob)

    println("Начало моделирования ансамбля траекторий...")

    ensemble_sol = solve(
        ensemble_prob,
        EM(),
        dt = 0.001,
        saveat = 0.01,
        trajectories = p[:trajectories]
        )

    return ensemble_sol
end

function analyze_ensemble(ensemble_sol, p)
    t_vals = ensemble_sol[1].t
    n_t = length(t_vals)
    n_traj = length(ensemble_sol)

    x_matrix = zeros(n_traj, n_t)

    for i in 1:n_traj
        for j in 1:n_t
            x_matrix[i, j] = ensemble_sol[i].u[j][1]
        end
    end

    mean_vals = [mean(x_matrix[:, j]) for j in 1:n_t]
        var_vals = [var(x_matrix[:, j]) for j in 1:n_t]

            mean_theory = [theoretical_mean(t, p) for t in t_vals]
                var_theory = [theoretical_variance(t, p) for t in t_vals]

                    return t_vals, x_matrix, mean_vals, var_vals, mean_theory, var_theory
                end

                Random.seed!(123)

                single_sol = run_single_trajectory(params)
                ensemble_sol = run_ensemble(params)

                t_vals, x_matrix, mean_vals, var_vals, mean_theory, var_theory =
                    analyze_ensemble(ensemble_sol, params)

                plt1 = plot(
                    single_sol.t,
                    [u[1] for u in single_sol.u],
                        title = "Одна траектория уравнения Ланжевена",
                        xlabel = "Время t",
                        ylabel = "Координата x(t)",
                        label = "Случайная траектория",
                        lw = 2,
                        fmt = :png
                        )

                plt2 = plot(
                    title = "Ансамбль стохастических траекторий",
                    xlabel = "Время t",
                    ylabel = "Координата x(t)",
                    legend = false,
                    fmt = :png
                    )

                for i in 1:min(20, params[:trajectories])
                    plot!(plt2, t_vals, x_matrix[i, :], alpha = 0.5)
                end

                plt3 = plot(
                    t_vals,
                    mean_vals,
                    title = "Среднее значение координаты",
                    xlabel = "Время t",
                    ylabel = "E[x(t)]",
                    label = "Численное среднее",
                    lw = 2,
                    fmt = :png
                    )

                plot!(
                    plt3,
                    t_vals,
                    mean_theory,
                    label = "Теория: x₀ exp(-γt)",
                    linestyle = :dash,
                    lw = 2,
                    fmt = :png
                    )

                plt4 = plot(
                    t_vals,
                    var_vals,
                    title = "Дисперсия координаты",
                    xlabel = "Время t",
                    ylabel = "Var[x(t)]",
                    label = "Численная дисперсия",
                    lw = 2,
                    fmt = :png
                    )

                plot!(
                    plt4,
                    t_vals,
                    var_theory,
                    label = "Теория",
                    linestyle = :dash,
                    lw = 2,
                    fmt = :png
                    )

                plt = plot(
                    plt1,
                    plt2,
                    plt3,
                    plt4,
                    layout = (2, 2),
                    size = (1000, 700),
                    fmt = :png
                    )

                path = plotsdir("langevin_equation_simulation.png")
                safesave(path, plt)

                println("График сохранён в: $path")

                display(plt)
