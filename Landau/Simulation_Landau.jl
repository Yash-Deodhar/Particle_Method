using LinearAlgebra
using MAT
using Printf
#using GLMakie #comment if GLMakie is not compatible with your system
include("supporting_functions_Landau.jl")
include("FP_Solvers_Landau.jl")
include("Plotting_Landau.jl")

function Simulation_Landau(Eqn, FP_Solver; N = 40, L = 4, t_0 = 0.0, t_f = 5.0, dt = 0.01/8, tol = 1e-15, win_size = 1, beta = 1.0, plots = false)
    
    
    gamma, true_sol, Fish_Calc = Equation_lookup[Eqn]
    

    # define constants and mesh

    N_total = N^2
    dv = 2*L/N
    eps = 4*(0.4*(dv)^0.99)^2
    v = collect((-L + dv/2):dv:(L - dv/2))
    vx = [x for x in v, y in v]   
    vy = [y for x in v, y in v]
    Vx = vec(vx)
    Vy = vec(vy)


    # Calculate intial particle positions and weights

    v_p = collect((-L + dv/2):dv:(L - dv/2))
    vx_p = [x for x in v_p, y in v_p]
    vy_p = [y for x in v_p, y in v_p]
    Vx_p = vec(vx_p)
    Vy_p = vec(vy_p)
    f_initial = true_sol(t_0, vx_p, vy_p)
    w_p = vec(dv^2 .* f_initial)


    # Set up plot and plot initial data

    f = reconstruction(N, w_p, vx, Vx_p, vy, Vy_p, eps)
    if plots
        time_obs, f_obs, true_sol_obs = initialize_plot(Eqn, FP_Solver, true_sol, t_0, f, vx, vy)
    end


    # Define Loop parameters and history variables

    max_iter = 100
    n_steps = round(Int,(t_f - t_0)/dt) # number of time steps
    Iter_history = zeros(n_steps)
    error_history = zeros(n_steps)
    Mom_x_history = zeros(n_steps)
    Mom_y_history = zeros(n_steps)
    KE_history = zeros(n_steps)
    Energy_history = zeros(n_steps)
    Fisher_history = zeros(n_steps)
    Dissipation_history = zeros(n_steps)

    Mom_x_initial = sum(w_p .* Vx_p)
    Mom_y_initial = sum(w_p .* Vy_p)
    KE_initial = sum(w_p .* (Vx_p .^ 2 + Vy_p .^ 2))

    # Main loop

    for n in 1:n_steps

        t = t_0 + dt*n
        Ux, Uy, s_int_x, s_int_y, dissipation = RHS_Landau(w_p, Vx_p, Vy_p, Vx_p, Vy_p, Vx, Vy, dv, eps, gamma, N_total)
        Vx_p_new = Vx_p .+ dt*Ux
        Vy_p_new = Vy_p .+ dt*Uy

        Vx_p_new, Vy_p_new, iter_n = FP_Solver(RHS_Landau, max_iter, Vx_p_new, Vy_p_new, w_p, Vx_p, Vy_p, Vx, Vy, dv, eps, N_total, dt, tol, gamma, win_size, beta)
        Iter_history[n] = iter_n
        Vx_p .= Vx_p_new
        Vy_p .= Vy_p_new    

        f = reconstruction(N, w_p, vx, Vx_p, vy, Vy_p, eps)
        error_history[n] = norm(true_sol(t, vx, vy) .- f, 2)
        Mom_x_history[n] = abs(sum(w_p .* Vx_p) - Mom_x_initial)
        Mom_y_history[n] = abs(sum(w_p .* Vy_p) - Mom_y_initial)
        KE_history[n] = abs(sum(w_p .* (Vx_p .^ 2 + Vy_p .^ 2)) - KE_initial)
        Energy_history[n] = dv^2*sum(f .* log.(f))
        Fisher_history[n] = Fish_Calc(f, w_p, s_int_x, s_int_y, Vx_p, Vy_p, vx, vy, dv, eps)
        Dissipation_history[n] = dissipation

        if plots
            update_plot(true_sol, time_obs, f_obs, true_sol_obs, t, f, vx, vy)
        end

        if (100*n)%n_steps == 0
           println(100*n/n_steps, "% done solving $Eqn using $(nameof(FP_Solver))")
           flush(stdout)
        end

    end

    suffix = nameof(FP_Solver) == :FPI ? "" : @sprintf("_w%d_b%.2e", win_size, beta)
    filename = @sprintf("%s_N%d_L%d_t0%g_tf%g_dt%g", nameof(FP_Solver), N, L, t_0, t_f, dt) * suffix * ".mat"
    filepath = joinpath("Figures and Data", Eqn, "Data", filename)
    matwrite(filepath, Dict("Iter_history" => Iter_history, "error_history" => error_history, "Mom_x_history" => Mom_x_history, "Mom_y_history" => Mom_y_history, "KE_history" => KE_history, "Energy_history" => Energy_history, "Fisher_history" => Fisher_history, "Dissipation_history" => Dissipation_history))
    params = (N, L, t_0, t_f, dt, win_size, beta)
    
    return Iter_history, error_history, Mom_x_history, Mom_y_history, KE_history, Energy_history, Fisher_history, Dissipation_history, params
end