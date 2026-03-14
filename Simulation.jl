#= Code to solve the 1D heat equation using the paticle method 
with an energy-dissipative and mass conservative scheme =#
# Code written by Yash Deodhar, Plotting code written by Gemini #

using LinearAlgebra
using GLMakie
include("supporting_functions.jl")
include("FP_Solvers.jl")
include("Plotting.jl")

function Simulation(Eqn, FP_Solver; N = 100, L = 15, t_0 = 2.0, t_f = 3.0, dt = 0.01, tol = 1e-15,  m = 0.0, win_size = 1, beta = 1.0, plots = false)

    true_sol, RHS = Equation_lookup[Eqn]

    # Define important constants and mesh #

    dx = 2*L/N  # Mesh/Cell width
    eps = 4*(0.4*(dx)^0.99)^2   # Regularization parameter
    x = collect((-L + dx/2):dx:(L - dx/2))  # function evaluation points
    n_steps = round(Int,(t_f - t_0)/dt) # number of time steps


    # Calculate intial particle positions and weights

    x_p = collect((-L + dx/2):dx:(L - dx/2)) # initial particle positions
    f_initial = true_sol(t_0, x_p,m) # inital condition
    w_p = dx .* f_initial # calculate paricle weights from intial condition

    # Set up plot and plot initial data
    f = reconstruction(N, w_p, x, x_p, eps)
    if plots
        time_obs, f_obs, true_sol_obs = initialze_plot(Eqn, FP_Solver, true_sol, t_0, f, x, m)
    end


    # Define Loop parameters and history variables

    max_iter = 100
    Iter_history = zeros(n_steps)
    error_history = zeros(n_steps+1)
    error_history[1] = norm(true_sol(t_0, x, m) .- f, 2)

    # Call main loop function to run simulation

    for n in 1:n_steps

        if (5*n)%n_steps == 0
            println(100*n/n_steps, "% done solving $Eqn using $(nameof(FP_Solver))")
        end

        t = t_0 + n*dt
        rhs = RHS(w_p, x_p, x_p, x, dx, eps, N, m) 
        x_p_new = x_p .+ dt*rhs

        x_p_new, iter_n = FP_Solver(RHS, max_iter, x_p_new, w_p, x_p, x, dx, eps, N, dt, tol, m, win_size, beta)
        Iter_history[n] = iter_n
        x_p .= x_p_new

        f = reconstruction(N, w_p, x, x_p, eps)
        error_history[n+1] = norm(true_sol(t, x, m) .- f, 2)
        
        if plots
            update_plot(true_sol, time_obs, f_obs, true_sol_obs, t, f, x, m)
        end

    end
    return Iter_history, error_history
end