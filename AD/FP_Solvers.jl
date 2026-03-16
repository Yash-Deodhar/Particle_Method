include("supporting_functions.jl")  


function FPI(RHS, max_iter, x_p_new, w_p, x_p, x, dx, eps, N, dt, tol, m, win_size, beta)

    # m, beta are not used anywhere but are necessary for compact code in simulations

    i = 0
    # Start Acceleration loop
    while i <= max_iter
        i += 1
        x_p_old = x_p_new
        rhs = RHS(w_p, x_p, x_p_new, x, dx, eps, N, m) 
        x_p_new = x_p .+ dt*rhs
        abs_res = norm(x_p_new .- x_p_old, 2)
        rel_res = abs_res/norm(x_p_new, 2)

        if rel_res < tol
            break
        end

    end
    return x_p_new, i
end

function AA(RHS, max_iter, x_p_new, w_p, x_p, x_cc, dx, eps, N, dt, tol, m, win_size, beta)
    
    # Define arrays for residual (w), iteration history (x) and difference between updates (e) and function (g) for AA
    w = Vector{Vector{Float64}}()
    e = Vector{Vector{Float64}}()
    x = Vector{Vector{Float64}}()
    g(z) = x_p .+ dt .* RHS(w_p, x_p, z, x_cc, dx, eps, N, m)

    # fixed-point uodate step (k = 0)
    push!(x, copy(x_p_new))
    push!(w, g(x[end]) .- x[end])
    push!(x, x[end] .+ beta .* w[end])
    push!(e, x[end] .- x[end-1]) 

    i = 0
    # Start Acceleration loop
    while i <= max_iter
        i += 1
        push!(w, g(x[end]) .- x[end]); # compute w_{k+1}
        m_k = min(length(x)-1, win_size) # k = length(x) - 1
        F_k = zeros(length(x[end]), m_k)
        E_k = zeros(length(x[end]), m_k)
        for j in 1:m_k
            F_k[:,j] = w[end-m_k+j] .- w[end-m_k+j-1] # compute F_k
            E_k[:,j] = e[end-m_k+j] # compute E_k
        end
        gamma_k = F_k\w[end]; # solution to argmin ||F_k*gamma_k - w_{k+1}||
        push!(x, x[end] .+ beta .* w[end] .- (E_k .+ beta .* F_k) * gamma_k) # new iteration update
        push!(e, x[end] .- x[end-1]) # compute e_k
        if norm(e[end], 2)/norm(x[end], 2) < tol
            break
        end
    end
    return x[end], i
end