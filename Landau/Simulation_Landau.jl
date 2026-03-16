using LinearAlgebra
#using GLMakie #comment if you GLMakie is not compatible with your system
include("supporting_functions.jl")
include("FP_Solvers.jl")
include("Plotting.jl")

function Simulation_Landau(Eqn, FP_Solver; N = 20, L = 4, t_0 = 0.0, t_f = 5.0, dt = 0.02, tol = 1e-15, win_size = 1, beta = 1.0, plots = false)
    
    
    gamma, true_sol, RHS = Equation_lookup[Eqn]
    

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
    error_history = zeros(n_steps+1)
    error_history[1] = norm(f_initial .- f, 2)


    # Main loop

    for n in 1:n_steps

        if (5*n)%n_steps == 0
            println(100*n/n_steps, "% done solving $Eqn using $(nameof(FP_Solver))")
        end

        t = t_0 + dt*n
        Ux, Uy = RHS(w_p, Vx_p, Vy_p, Vx_p, Vy_p, Vx, Vy, dv, eps, gamma, N_total)
        Vx_p_new = Vx_p .+ dt*Ux
        Vy_p_new = Vy_p .+ dt*Uy

        Vx_p_new, Vy_p_new, iter_n = FP_Solver(RHS, max_iter, Vx_p_new, Vy_p_new, w_p, Vx_p, Vy_p, Vx, Vy, dv, eps, N_total, dt, tol, gamma, win_size, beta)
        Iter_history[n] = iter_n
        Vx_p .= Vx_p_new
        Vy_p .= Vy_p_new    

        f = reconstruction(N, w_p, vx, Vx_p, vy, Vy_p, eps)
        error_history[n+1] = norm(true_sol(t, vx, vy) .- f, 2)

        if plots
            update_plot(true_sol, time_obs, f_obs, true_sol_obs, t, f, vx, vy)
        end

    end
    return Iter_history, error_history
end