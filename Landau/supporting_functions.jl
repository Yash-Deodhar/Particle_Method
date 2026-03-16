using FastGaussQuadrature

 # Define regularized delta and its derivative

phi_2d(x,y,eps) = 1/(2*pi*eps) .* exp.(- (x .^ 2 .+ y .^ 2) ./ (2*eps))
phi_2d_grad(x, y, eps) = (phi_2d(x, y, eps) .* (-x ./ eps), phi_2d(x, y, eps) .* (-y ./ eps))


# Define true solutions

function true_sol_maxwell(t,x,y) 
    R = 1 - 0.5*exp(-t/8)
    abs_xy = x.^2 .+ y.^2
    return 1/(2*pi*R) .* exp.(-(abs_xy) ./ (2*R)) .* ((2*R-1)/R .+ (1 - R)/(2*R^2) .* abs_xy)
end


# Define functions for calculating right hand side

function RHS_Maxwell(w, Vx_prev_t, Vy_prev_t, Vx_curr_itr, Vy_curr_itr, Vx_cc, Vy_cc, dv, eps, gamma, N_total)
    
    Vmid_x = 0.5*(Vx_curr_itr + Vx_prev_t);
    Vmid_y = 0.5*(Vy_curr_itr + Vy_prev_t);
    Vdif_x = Vx_curr_itr - Vx_prev_t;
    Vdif_y = Vy_curr_itr - Vy_prev_t;

    gn, gw = gausslegendre(4) # 4 point Gauss-Legendre quadrature
    gn = (gn .+ 1) ./ 2 # rescale nodes to [0,1]
    gw = gw ./ 2 # rescale weights to [0,1
    log_sum = zeros(N_total, length(gn))
    s_int_x = zeros(N_total)
    s_int_y = zeros(N_total)

    for i in 1:N_total
        for j in 1:length(gn)
            log_sum[i,j] = sum(w .* phi_2d(Vx_cc[i] .- (Vx_prev_t .+ gn[j] .* Vdif_x), Vy_cc[i] .- (Vy_prev_t .+ gn[j] .* Vdif_y), eps))
        end
    end
    log_sum = log.(log_sum)

    term1_x = zeros(N_total, length(gn))
    term1_y = zeros(N_total, length(gn))
    for i in 1:N_total
        for j in 1:length(gn)
            A, B = phi_2d_grad(Vx_cc .- (Vx_prev_t[i] .+ gn[j] .* Vdif_x[i]), Vy_cc .- (Vy_prev_t[i] .+ gn[j] .* Vdif_y[i]), eps)
            term1_x[i,j] = dv^2*sum(A .* log_sum[:,j])
            term1_y[i,j] = dv^2*sum(B .* log_sum[:,j])
        end
    end

    for i in 1:N_total
        s_int_x[i] = sum(gw .* term1_x[i,:])
        s_int_y[i] = sum(gw .* term1_y[i,:])
    end

    U_x = zeros(N_total)
    U_y = zeros(N_total)
    for i in 1:N_total
        len = sqrt.((Vmid_x[i] .- Vmid_x) .^ 2 + (Vmid_y[i] .- Vmid_y) .^ 2)

        A11 = 1/16 * len .^ gamma .* (len .^ 2 - (Vmid_x[i] .- Vmid_x) .^ 2)
        A12 = -1/16 * len .^ gamma .* (Vmid_x[i] .- Vmid_x) .* (Vmid_y[i] .- Vmid_y)
        A21 = A12
        A22 = 1/16 * len .^ gamma .* (len .^ 2 - (Vmid_y[i] .- Vmid_y) .^ 2)

        termA = s_int_x[i] .- s_int_x
        termB = s_int_y[i] .- s_int_y
        U_x[i] = sum(w .* (A11 .* termA .+ A12 .* termB))
        U_y[i] = sum(w .* (A21 .* termA .+ A22 .* termB))

    end

    return U_x, U_y
end

Equation_lookup = Dict(
    "Maxwell"   => (0, true_sol_maxwell, RHS_Maxwell)
)