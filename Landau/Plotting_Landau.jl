#using GLMakie #comment if you GLMakie is not compatible with your system
include("supporting_functions_Landau.jl")


function reconstruction(N, w_p, vx, Vx_p, vy, Vy_p, eps)
    f = zeros(N,N)
    for i in 1:N
        for j in 1:N
            f[i,j] = sum(w_p.*phi_2d(vx[i,j] .- Vx_p, vy[i,j] .- Vy_p, eps))
        end
    end
    return f
end

function initialize_plot(Eqn, FP_Solver, true_sol, t_0, f, vx, vy)
    time_obs     = Observable(t_0)
    f_obs        = Observable(f)
    true_sol_obs = Observable(true_sol(t_0, vx, vy))
    xs = vx[:, 1]
    ys = vy[1, :]
    fig = Figure(size = (1200, 600))
    title_str = lift(t -> "$Eqn using $(nameof(FP_Solver)) — t = $(round(t, digits=2))", time_obs)
    Label(fig[0, :], title_str, fontsize = 16)
    ax1 = Axis3(fig[1, 1], title = "True Solution",    xlabel = "vx", ylabel = "vy", zlabel = "f")
    ax2 = Axis3(fig[1, 2], title = "Numerical Solution", xlabel = "vx", ylabel = "vy", zlabel = "f")
    surface!(ax1, xs, ys, true_sol_obs, colormap = :viridis)
    surface!(ax2, xs, ys, f_obs,        colormap = :viridis)
    display(fig)
    return time_obs, f_obs, true_sol_obs
end

function update_plot(true_sol, time_obs, f_obs, true_sol_obs, t, f, vx, vy)
    f_obs[] = f
    true_sol_obs[] = true_sol(t, vx, vy)
    time_obs[] = t
    yield()
end

function save_plots_Landau(eqn_name, AA_iter, AA_error, AA_mom_x, AA_mom_y, AA_KE, AA_energy, AA_fisher, AA_diss, FPI_iter, FPI_error, FPI_mom_x, FPI_mom_y, FPI_KE, FPI_energy, FPI_fisher, FPI_diss, N, L, t_0, t_f, dt, win_size, beta)
    
    folder = joinpath("Figures and Data", eqn_name, "Plots")
    base_AA = @sprintf("N%d_L%d_t0%g_tf%g_dt%g_w%d_b%.2e", N, L, t_0, t_f, dt, win_size, beta)
    t_axis = collect(range(t_0, t_f, length=round(Int,(t_f-t_0)/dt)))

    quantities = [
        ("iter_hist",   AA_iter,   FPI_iter,   "Iteration Count"),
        ("error_hist",  AA_error,  FPI_error,  "Error"),
        ("mom_x_hist",  AA_mom_x,  FPI_mom_x,  "Momentum x"),
        ("mom_y_hist",  AA_mom_y,  FPI_mom_y,  "Momentum y"),
        ("KE_hist",     AA_KE,     FPI_KE,     "Kinetic Energy"),
        ("energy_hist", AA_energy, FPI_energy, "Energy"),
        ("fisher_hist", AA_fisher, FPI_fisher, "Fisher Information"),
        ("diss_hist",   AA_diss,   FPI_diss,   "Dissipation"),
    ]

    for (fname, AA_data, FPI_data, ylabel_text) in quantities
        fig = Figure(size = (900, 600))
        ax = Axis(fig[1, 1], title = "$(eqn_name): AA vs FPI $(replace(fname, "_" => " "))", xlabel = "Time", ylabel = ylabel_text)
        lines!(ax, t_axis, AA_data,  label = "AA",  linewidth = 2)
        lines!(ax, t_axis, FPI_data, label = "FPI", linewidth = 2)
        axislegend(ax, position = :rb)
        save(joinpath(folder, "$(fname)_$(base_AA).png"), fig)
    end
end