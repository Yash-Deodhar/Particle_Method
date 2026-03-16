#using GLMakie #comment if you GLMakie is not compatible with your system
include("supporting_functions.jl")


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
