#using GLMakie #comment if GLMakie is not compatible with your system. 
using CairoMakie
include("AD/Simulation_AD.jl")
include("AD/FP_Solvers_AD.jl")
include("AD/supporting_functions_AD.jl")
include("AD/Plotting_AD.jl")
include("Landau/Simulation_Landau.jl")
include("Landau/FP_Solvers_Landau.jl")
include("Landau/supporting_functions_Landau.jl")
include("Landau/Plotting_Landau.jl")

# AA_iter_hist_heat, AA_error_hist_heat, AA_energy_hist_heat, params_heat = Simulation_AD("Heat", AA, win_size = 3, beta = 0.95); 
# FPI_iter_hist_heat, FPI_error_hist_heat, FPI_energy_hist_heat, params_heat = Simulation_AD("Heat", FPI);
# save_plots_AD("Heat", AA_iter_hist_heat, AA_error_hist_heat, AA_energy_hist_heat, FPI_iter_hist_heat, FPI_error_hist_heat, FPI_energy_hist_heat, params_heat...)
# println("Heat done")
# flush(stdout)

# AA_iter_hist_porous, AA_error_hist_porous, AA_energy_hist_porous, params_porous = Simulation_AD("PorousMedium", AA, L = 8, m = 3/2, win_size = 3, beta = 0.95); 
# FPI_iter_hist_porous, FPI_error_hist_porous, FPI_energy_hist_porous, params_porous = Simulation_AD("PorousMedium", FPI, L = 8, m = 3/2)
# save_plots_AD("PorousMedium", AA_iter_hist_porous, AA_error_hist_porous, AA_energy_hist_porous, FPI_iter_hist_porous, FPI_error_hist_porous, FPI_energy_hist_porous, params_porous...)
# println("PorousMedium done")
# flush(stdout)

# AA_iter_hist_lf, AA_error_hist_lf, AA_energy_hist_lf, params_lf = Simulation_AD("LinearFokker", AA, L = 5, t_0 = 0.5, t_f = 1.0, win_size = 3, beta = 0.95); 
# FPI_iter_hist_lf, FPI_error_hist_lf, FPI_energy_hist_lf, params_lf = Simulation_AD("LinearFokker", FPI, L = 5, t_0 = 0.5, t_f = 1.0)
# save_plots_AD("LinearFokker", AA_iter_hist_lf, AA_error_hist_lf, AA_energy_hist_lf, FPI_iter_hist_lf, FPI_error_hist_lf, FPI_energy_hist_lf, params_lf...)
# println("Linear Fokker done")
# flush(stdout)

# AA_iter_hist_nlf, AA_error_hist_nlf, AA_energy_hist_nlf, params_nlf = Simulation_AD("NLFokker", AA, L = 5, t_0 = 0.5, t_f = 1.0, dt = 0.001, win_size = 3, beta = 0.95); 
# FPI_iter_hist_nlf, FPI_error_hist_nlf, FPI_energy_hist_nlf, params_nlf = Simulation_AD("NLFokker", FPI, L = 5, t_0 = 0.5, t_f = 1.0, dt = 0.001)
# save_plots_AD("NLFokker", AA_iter_hist_nlf, AA_error_hist_nlf, AA_energy_hist_nlf, FPI_iter_hist_nlf, FPI_error_hist_nlf, FPI_energy_hist_nlf, params_nlf...)
# println("NLFokker done")
# flush(stdout)

# FPI_iter_hist_Maxwell, FPI_error_hist_Maxwell, FPI_mom_x_Maxwell, FPI_mom_y_Maxwell, FPI_KE_Maxwell, FPI_Energy_Maxwell, FPI_Fish_Maxwell, FPI_Diss_Maxwell, params_Maxwell = Simulation_Landau("Landau_Maxwell", FPI, N = 40)
# println("Max_FPI done")
# flush(stdout)
# AA_iter_hist_Maxwell, AA_error_hist_Maxwell, AA_mom_x_Maxwell, AA_mom_y_Maxwell, AA_KE_Maxwell, AA_Energy_Maxwell, AA_Fish_Maxwell, AA_Diss_Maxwell, params_Maxwell = Simulation_Landau("Landau_Maxwell", AA, N = 40, win_size = 3, beta = 0.95)
# println("Max_AA done")
# flush(stdout)
# save_plots_Landau("Landau_Maxwell", AA_iter_hist_Maxwell, AA_error_hist_Maxwell, AA_mom_x_Maxwell, AA_mom_y_Maxwell, AA_KE_Maxwell, AA_Energy_Maxwell, AA_Fish_Maxwell, AA_Diss_Maxwell, FPI_iter_hist_Maxwell, FPI_error_hist_Maxwell, FPI_mom_x_Maxwell, FPI_mom_y_Maxwell, FPI_KE_Maxwell, FPI_Energy_Maxwell, FPI_Fish_Maxwell, FPI_Diss_Maxwell, params_Maxwell...)

FPI_iter_hist_Coulomb, FPI_error_hist_Coulomb, FPI_mom_x_Coulomb, FPI_mom_y_Coulomb, FPI_KE_Coulomb, FPI_Energy_Coulomb, FPI_Fish_Coulomb, FPI_Diss_Coulomb, params_Coulomb = Simulation_Landau("Landau_Coulomb", FPI, N = 40, L = 10, t_f = 20, dt = 0.05)
println("Coul_FPI done")
flush(stdout)
AA_iter_hist_Coulomb, AA_error_hist_Coulomb, AA_mom_x_Coulomb, AA_mom_y_Coulomb, AA_KE_Coulomb, AA_Energy_Coulomb, AA_Fish_Coulomb, AA_Diss_Coulomb, params_Coulomb = Simulation_Landau("Landau_Coulomb", AA, N = 40, L = 10, t_f = 20, dt = 0.05, win_size = 3, beta = 0.95)
println("Coul_AA done")
flush(stdout)
save_plots_Landau("Landau_Coulomb", AA_iter_hist_Coulomb, AA_error_hist_Coulomb, AA_mom_x_Coulomb, AA_mom_y_Coulomb, AA_KE_Coulomb, AA_Energy_Coulomb, AA_Fish_Coulomb, AA_Diss_Coulomb, FPI_iter_hist_Coulomb, FPI_error_hist_Coulomb, FPI_mom_x_Coulomb, FPI_mom_y_Coulomb, FPI_KE_Coulomb, FPI_Energy_Coulomb, FPI_Fish_Coulomb, FPI_Diss_Coulomb, params_Coulomb...)
