#using GLMakie #comment if you GLMakie is not compatible with your system. 
include("AD/Simulation_AD.jl")
include("AD/FP_Solvers.jl")
include("AD/supporting_functions.jl")
include("Landau/Simulation_Landau.jl")
include("Landau/FP_Solvers.jl")
include("Landau/supporting_functions.jl")

# AA_iter_hist_heat, AA_error_hist_heat = Simulation_AD("Heat", AA, win_size = 3, beta = 0.95); 
# FPI_iter_hist_heat, FPI_error_hist_heat = Simulation_AD("Heat", FPI);

# AA_iter_hist_porous, AA_error_hist_porous = Simulation_AD("PorousMedium", AA, L = 8, m = 3/2, win_size = 3, beta = 0.95); 
# FPI_iter_hist_porous, FPI_error_hist_porous = Simulation_AD("PorousMedium", FPI, L = 8, m = 3/2)

# AA_iter_hist_lf, AA_error_hist_lf = Simulation_AD("LinearFokker", AA, L = 5, t_0 = 0.5, t_f = 1.0, dt = 0.001, win_size = 3, beta = 0.95); 
# FPI_iter_hist_lf, FPI_error_hist_lf = Simulation_AD("LinearFokker", FPI, L = 5, t_0 = 0.5, t_f = 1.0, dt = 0.001)

# AA_iter_hist_nlf, AA_error_hist_nlf = Simulation_AD("NLFokker", AA, L = 5, t_0 = 0.5, t_f = 1.0, dt = 0.001, win_size = 3, beta = 0.95); 
# FPI_iter_hist_nlf, FPI_error_hist_nlf = Simulation_AD("NLFokker", FPI, L = 5, t_0 = 0.5, t_f = 1.0, dt = 0.001)

FPI_iter_hist_maxwell, FPI_error_hist_maxwell = Simulation_Landau("Maxwell", FPI)
AA_iter_hist_maxwell, AA_error_hist_maxwell = Simulation_Landau("Maxwell", AA, win_size = 3, beta = 0.95)
println(sum(FPI_iter_hist_maxwell), sum(AA_iter_hist_maxwell))


# fig1 = Figure()
# ax1 = Axis(fig1[1,1], xlabel = "Time step", ylabel = "Iterations", title = "Heat iteration history")
# lines!(ax1, 1:length(FPI_iter_hist_heat), FPI_iter_hist_heat, label = "FPI")
# lines!(ax1, 1:length(AA_iter_hist_heat), AA_iter_hist_heat, label = "AA")
# axislegend(ax1)
# display(GLMakie.Screen(), fig1)

# fig2 = Figure()
# ax2 = Axis(fig2[1,1], xlabel = "Time step", ylabel = "Iterations", title = "PorousMedium iteration history")
# lines!(ax2, 1:length(FPI_iter_hist_porous), FPI_iter_hist_porous, label = "FPI")
# lines!(ax2, 1:length(AA_iter_hist_porous), AA_iter_hist_porous, label = "AA")
# axislegend(ax2)
# display(GLMakie.Screen(), fig2)

# fig3 = Figure()
# ax3 = Axis(fig3[1,1], xlabel = "Time step", ylabel = "Iterations", title = "LinearFokker iteration history")
# lines!(ax3, 1:length(FPI_iter_hist_lf), FPI_iter_hist_lf, label = "FPI")
# lines!(ax3, 1:length(AA_iter_hist_lf), AA_iter_hist_lf, label = "AA")
# axislegend(ax3)
# display(GLMakie.Screen(), fig3)

# fig4 = Figure()
# ax4 = Axis(fig4[1,1], xlabel = "Time step", ylabel = "Iterations", title = "NLFokker iteration history")
# lines!(ax4, 1:length(FPI_iter_hist_nlf), FPI_iter_hist_nlf, label = "FPI")
# lines!(ax4, 1:length(AA_iter_hist_nlf), AA_iter_hist_nlf, label = "AA")
# axislegend(ax4)
# display(GLMakie.Screen(), fig4)




