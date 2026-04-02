using FastGaussQuadrature


# Define regularized delta and its derivative

phi(x, eps) = 1/sqrt(2*pi*eps) .* exp.(- (x .^ 2) ./ (2*eps))
phi_prime(x, eps) = - x ./ (eps*sqrt(2*pi*eps)) .* exp.(- (x .^ 2) ./ (2*eps))


# Define true solutions (parameter m required in heat and fokker for compatibility with porous)

true_sol_heat(t,x,m) = 1/sqrt(4*pi*t) .* exp.(- (x .^ 2) ./ (4*t))
true_sol_fokker(t,x,m) = (2*pi*(1-exp(-2*t)))^(-0.5) .* exp.( - (x .^ 2) ./ (2*(1-exp(-2*t))))
true_sol_porous(t, x, m) = (1 ./ t .^ (1/(m+1))) .* (max.(0, 1 .- (m-1)/(2*m*(m+1)) .* (x.^2 ./ t .^ (2/(m+1)))) .^ (1/(m-1)))

# Define energy functions

energy_heat(f, w, x, x_p, dx, m) = dx*sum(f .* log.(f))
energy_porous(f, w, x, x_p, dx, m) = dx*sum(f .* log.(f))
energy_linear_fokker(f, w, x, x_p, dx, m) = dx*sum(f .* log.(f)) + 0.5*sum(w .* x_p.^2)
function energy_nl_fokker(f, w, x, x_p, dx, m)
    W_term = zeros(length(x_p))
    for i in 1:length(x_p)
        W_term[i] = sum(0.5*w .* (x_p[i] .- x_p).^2)
    end
    return dx*(sum(f .* log.(f)) + 0.5*sum(w .* W_term))
end


# Define functions for calculating right hand side

function RHS_Heat(w, x_prev_t, x_curr_itr, x_cc, dx, eps, n, m=0)

    # w = particle weights (n x 1)
    # x_prev_t = x_p at previous time  (n x 1)
    # x_curr_itr = x_p at previous iteration at current time step (n x 1)
    # x_cc = cell centers (n x 1)
    # dx = cell width (scalar)
    # eps = parameter for regularized delta function (scalar)
    # n = total number of particles (scalar)
    # m is not used here, required for compatibility with RHS_PorousMedium

    gn, gw = gausslegendre(4) # 4 point Gauss-Legendre quadrature
    gn = (gn .+ 1) ./ 2 # rescale nodes to [0,1]
    gw = gw ./ 2 # rescale weights to [0,1]
    total_dist = x_curr_itr .- x_prev_t
    num_nodes = length(gn)
    log_sum = zeros(num_nodes, n)
    x_int = zeros(num_nodes, n)
    out = zeros(n)

    for i in 1:n
        for j in 1:num_nodes
            log_sum[j,i] = log(sum(w .* phi(x_cc[i] .- (x_prev_t .+ gn[j]*total_dist), eps))) # Calculate log term           
        end
    end

    for i in 1:n
        for j in 1:num_nodes
            prime_term = phi_prime(-(x_prev_t[i] .+ gn[j]*total_dist[i]) .+ x_cc, eps) # calculate phi_prime term
            x_int[j,i] = dx*sum(prime_term .* log_sum[j,:])  # integrate over x
        end
    end

    for i in 1:n
        out[i] = sum(gw .* x_int[:,i]) # integrate over s
    end

    return out
end

function RHS_LinearFokker(w, x_prev_t, x_curr_itr, x_cc, dx, eps, n, m=0)

    # w = particle weights (n x 1)
    # x_prev_t = x_p at previous time  (n x 1)
    # x_curr_itr = x_p at previous iteration at current time step (n x 1)
    # x_cc = cell centers (n x 1)
    # dx = cell width (scalar)
    # eps = parameter for regularized delta function (scalar)
    # n = total number of particles (scalar)
    # m is not used here, required for compatibility with RHS_PorousMedium

    gn, gw = gausslegendre(4) # 4 point Gauss-Legendre quadrature
    gn = (gn .+ 1) ./ 2 # rescale nodes to [0,1]
    gw = gw ./ 2 # rescale weights to [0,1]
    total_dist = x_curr_itr .- x_prev_t
    num_nodes = length(gn)
    log_sum = zeros(num_nodes, n)
    x_int = zeros(num_nodes, n)
    out = zeros(n)

    for i in 1:n
        for j in 1:num_nodes
            log_sum[j,i] = log(sum(w .* phi(x_cc[i] .- (x_prev_t .+ gn[j]*total_dist), eps))) # Calculate log term           
        end
    end

    for i in 1:n
        for j in 1:num_nodes
            prime_term = phi_prime(-(x_prev_t[i] .+ gn[j]*total_dist[i]) .+ x_cc, eps) # calculate phi_prime term
            x_int[j,i] = dx*sum(prime_term .* log_sum[j,:])  # integrate over x
        end
    end

    for i in 1:n
        out[i] = sum(gw .* x_int[:,i]) - 0.5*(x_curr_itr[i] + x_prev_t[i]) # integrate over s and add last term
    end

    return out
end

function RHS_NLFokker(w, x_prev_t, x_curr_itr, x_cc, dx, eps, n, m=0)

    # w = particle weights (n x 1)
    # x_prev_t = x_p at previous time  (n x 1)
    # x_curr_itr = x_p at previous iteration at current time step (n x 1)
    # x_cc = cell centers (n x 1)
    # dx = cell width (scalar)
    # eps = parameter for regularized delta function (scalar)
    # n = total number of particles (scalar)
    # m is not used here, required for compatibility with RHS_PorousMedium

    gn, gw = gausslegendre(4) # 4 point Gauss-Legendre quadrature
    gn = (gn .+ 1) ./ 2 # rescale nodes to [0,1]
    gw = gw ./ 2 # rescale weights to [0,1]
    total_dist = x_curr_itr .- x_prev_t
    num_nodes = length(gn)
    log_sum = zeros(num_nodes, n)
    x_int = zeros(num_nodes, n)
    out = zeros(n)

    for i in 1:n
        for j in 1:num_nodes
            log_sum[j,i] = log(sum(w .* phi(x_cc[i] .- (x_prev_t .+ gn[j]*total_dist), eps))) # Calculate log term           
        end
    end

    for i in 1:n
        for j in 1:num_nodes
            prime_term = phi_prime(-(x_prev_t[i] .+ gn[j]*total_dist[i]) .+ x_cc, eps) # calculate phi_prime term
            x_int[j,i] = dx*sum(prime_term .* log_sum[j,:])  # integrate over x
        end
    end

    for i in 1:n
        out[i] = sum(gw .* x_int[:,i]) - 0.5*sum(w .* (x_prev_t[i] .- x_cc .+ x_curr_itr[i] .- x_curr_itr)) # integrate over s and add last term
    end

    return out
end

function RHS_PorousMedium(w, x_prev_t, x_curr_itr, x_cc, dx, eps, n, m)

    # w = particle weights (n x 1)
    # x_prev_t = x_p at previous time  (n x 1)
    # x_curr_itr = x_p at previous iteration at current time step (n x 1)
    # x_cc = cell centers (n x 1)
    # dx = cell width (scalar)
    # eps = parameter for regularized delta function (scalar)
    # n = total number of particles (scalar)
    # m = porous medium equation parameter (scalar)

    gn, gw = gausslegendre(4) # 4 point Gauss-Legendre quadrature
    gn = (gn .+ 1) ./ 2 # rescale nodes to [0,1]
    gw = gw ./ 2 # rescale weights to [0,1]
    total_dist = x_curr_itr .- x_prev_t
    num_nodes = length(gn)
    f_m = zeros(num_nodes, n)
    x_int = zeros(num_nodes, n)
    out = zeros(n)

    for i in 1:n
        for j in 1:num_nodes
            f_m[j,i] = (sum(w .* phi(x_cc[i] .- (x_prev_t .+ gn[j]*total_dist), eps))) .^ (m-1) # Calculate m-1 power term           
        end
    end

    for i in 1:n
        for j in 1:num_nodes
            prime_term = phi_prime(-(x_prev_t[i] .+ gn[j]*total_dist[i]) .+ x_cc, eps) # calculate phi_prime term
            x_int[j,i] = (m/(m-1))*dx*sum(prime_term .* f_m[j,:])  # integrate over x
        end
    end

    for i in 1:n
        out[i] = sum(gw .* x_int[:,i])  # integrate over s
    end

    return out
end


Equation_lookup_AD = Dict(
    "Heat"     => (true_sol_heat, RHS_Heat, energy_heat),
    "LinearFokker" => (true_sol_fokker, RHS_LinearFokker, energy_linear_fokker),
    "NLFokker" => (true_sol_fokker, RHS_NLFokker, energy_nl_fokker),
    "PorousMedium" => (true_sol_porous, RHS_PorousMedium, energy_porous)
)