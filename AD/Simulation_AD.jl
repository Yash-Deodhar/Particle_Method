using LinearAlgebra
using MAT
using Printf
#using GLMakie #comment if you GLMakie is not compatible with your system
include("supporting_functions_AD.jl")
include("FP_Solvers_AD.jl")
include("Plotting_AD.jl")

function Simulation_AD(Eqn, FP_Solver; N = 100, L = 15, t_0 = 2.0, t_f = 3.0, dt = 0.01, tol = 1e-15,  m = 0.0, win_size = 1, beta = 1.0, sim_plot = false)

    true_sol, RHS, energy = Equation_lookup_AD[Eqn]

    # Define important constants and mesh #

    dx = 2*L/N  # Mesh/Cell width
    eps = 4*(0.4*(dx)^0.99)^2   # Regularization parameter
    x = collect((-L + dx/2):dx:(L - dx/2))  # function evaluation points


    # Calculate intial particle positions and weights

    x_p = collect((-L + dx/2):dx:(L - dx/2)) # initial particle positions
    f_initial = true_sol(t_0, x_p,m) # inital condition
    w_p = dx .* f_initial # calculate paricle weights from intial condition

    # Set up plot and plot initial data

    f = reconstruction(N, w_p, x, x_p, eps)
    if sim_plot
        time_obs, f_obs, true_sol_obs = initialze_plot(Eqn, FP_Solver, true_sol, t_0, f, x, m)
    end


    # Define Loop parameters and history variables

    max_iter = 100
    n_steps = round(Int,(t_f - t_0)/dt) # number of time steps
    Iter_history = zeros(n_steps)
    error_history = zeros(n_steps)
    energy_history = zeros(n_steps)
    
    
    # Main loop

    for n in 1:n_steps

        t = t_0 + n*dt
        rhs = RHS(w_p, x_p, x_p, x, dx, eps, N, m) 
        x_p_new = x_p .+ dt*rhs

        x_p_new, iter_n = FP_Solver(RHS, max_iter, x_p_new, w_p, x_p, x, dx, eps, N, dt, tol, m, win_size, beta)
        Iter_history[n] = iter_n
        x_p .= x_p_new

        f = reconstruction(N, w_p, x, x_p, eps)
        error_history[n] = norm(true_sol(t, x, m) .- f, 2)
        energy_history[n] = energy(f, w_p, x, x_p, dx, m)
        
        if sim_plot
            update_plot(true_sol, time_obs, f_obs, true_sol_obs, t, f, x, m)
        end

        # if (5*n)%n_steps == 0
        #     println(100*n/n_steps, "% done solving $Eqn using $(nameof(FP_Solver))")
        # end

    end

    suffix = nameof(FP_Solver) == :FPI ? "" : @sprintf("_w%d_b%.2e", win_size, beta)
    filename = @sprintf("%s_N%d_L%d_t0%g_tf%g_dt%g", nameof(FP_Solver), N, L, t_0, t_f, dt) * suffix * ".mat"
    filepath = joinpath("Figures and Data", Eqn, "Data", filename)
    matwrite(filepath, Dict("Iter_history"   => Iter_history, "error_history"  => error_history, "energy_history" => energy_history))
    params = (N, L, t_0, t_f, dt, win_size, beta)

    return Iter_history, error_history, energy_history, params
end