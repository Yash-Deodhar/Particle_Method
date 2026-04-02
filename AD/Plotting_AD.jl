#using GLMakie #comment if you GLMakie is not compatible with your system
using CairoMakie
using FileIO
include("supporting_functions_AD.jl")

function reconstruction(N, w_p, x, x_p, eps)
    f = zeros(N)
    for i in 1:N
        f[i] = sum(w_p.*phi(x[i] .- x_p, eps))
    end
    return f
end

function initialze_plot(Eqn, FP_Solver, true_sol, t_0, f, x, m)
    time_obs = Observable(t_0)
    f_obs = Observable(f)
    true_sol_obs = Observable(true_sol(t_0, x, m))
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], title = lift(t -> "$Eqn using $(nameof(FP_Solver)) - Time: $(round(t, digits=2))", time_obs), xlabel = "x", ylabel = "f(x)")
    ylims!(ax, 0.0, maximum(true_sol(t_0, x, m))) 
    lines!(ax, x, true_sol_obs, label = "True", color = :black, linestyle = :dash)
    lines!(ax, x, f_obs, label = "Numerical", color = :blue)
    axislegend(ax)
    display(GLMakie.Screen(), fig)
    return time_obs, f_obs, true_sol_obs
end

function update_plot(true_sol, time_obs, f_obs, true_sol_obs, t, f, x, m)
    f_obs[] = f
    true_sol_obs[] = true_sol(t, x, m)
    time_obs[] = t
    yield()
end

function save_plots_AD(eqn_name, AA_iter, AA_error, AA_energy, FPI_iter, FPI_error, FPI_energy, N, L, t_0, t_f, dt, win_size, beta)
    
    folder = joinpath("Figures and Data", eqn_name, "Plots")
    base_AA  = @sprintf("N%d_L%d_t0%g_tf%g_dt%g_w%d_b%.2e", N, L, t_0, t_f, dt, win_size, beta)
    t_axis = collect(range(t_0, t_f, length=round(Int,(t_f-t_0)/dt)))
    quantities = [("iter_hist",   AA_iter,   FPI_iter,   "Iteration Count"), ("error_hist",  AA_error,  FPI_error,  "Error"), ("energy_hist", AA_energy, FPI_energy, "Energy")]

    for (fname, AA_data, FPI_data, ylabel_text) in quantities
        fig = Figure(size = (900, 600))
        ax = Axis(fig[1, 1], title = "$(eqn_name): AA vs FPI $(replace(fname, "_" => " "))", xlabel = "Time", ylabel = ylabel_text)
        lines!(ax, t_axis, AA_data, label = "AA", linewidth = 2)
        lines!(ax, t_axis, FPI_data, label = "FPI", linewidth = 2)
        axislegend(ax, position = :rb)
        save(joinpath(folder, "$(fname)_$(base_AA).png"), fig)
    end
end