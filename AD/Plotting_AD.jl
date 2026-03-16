#using GLMakie #comment if you GLMakie is not compatible with your system
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

