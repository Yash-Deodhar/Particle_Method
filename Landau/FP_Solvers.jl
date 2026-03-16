include("supporting_functions.jl")  


function FPI(RHS, max_iter, Vx_new, Vy_new, w_p, Vx_prev, Vy_prev, Vx_cc, Vy_cc, dv, eps, N_total, dt, tol, gamma, win_size, beta)

    # win_size, beta are not used anywhere but are necessary for compatibility with AA

    i = 0
    while i <= max_iter
        i += 1

        Vx_old = Vx_new
        Vy_old = Vy_new

        Ux, Uy = RHS(w_p, Vx_prev, Vy_prev, Vx_new, Vy_new, Vx_cc, Vy_cc, dv, eps, gamma, N_total)

        Vx_new = Vx_prev .+ dt .* Ux
        Vy_new = Vy_prev .+ dt .* Uy

        abs_res = sqrt(sum((Vx_new .- Vx_old).^2 .+ (Vy_new .- Vy_old).^2))
        rel_res = abs_res / sqrt(sum(Vx_new.^2 .+ Vy_new.^2))

        if rel_res < tol
            break
        end
    end

    return Vx_new, Vy_new, i
end

function AA(RHS, max_iter, Vx_new, Vy_new, w_p, Vx_prev, Vy_prev, Vx_cc, Vy_cc, dv, eps, N_total, dt, tol, gamma, win_size, beta)

    # Combine the x- and y-particle locations into one state vector
    z0 = vcat(Vx_new, Vy_new)

    # Define arrays for residual (w), iteration history (x),
    # difference between iterates (e), and fixed-point map g
    w = Vector{Vector{Float64}}()
    e = Vector{Vector{Float64}}()
    x = Vector{Vector{Float64}}()

    function g(z)
        Vx_itr = z[1:N_total]
        Vy_itr = z[N_total+1:end]

        Ux, Uy = RHS(w_p, Vx_prev, Vy_prev, Vx_itr, Vy_itr, Vx_cc, Vy_cc, dv, eps, gamma, N_total)

        return vcat(Vx_prev .+ dt .* Ux,
                    Vy_prev .+ dt .* Uy)
    end

    # Fixed-point update step (k = 0)
    push!(x, copy(z0))
    push!(w, g(x[end]) .- x[end])
    push!(x, x[end] .+ beta .* w[end])
    push!(e, x[end] .- x[end-1])

    i = 0
    while i <= max_iter
        i += 1

        push!(w, g(x[end]) .- x[end])   # compute w_{k+1}

        m_k = min(length(x) - 1, win_size)

        F_k = zeros(length(x[end]), m_k)
        E_k = zeros(length(x[end]), m_k)

        for j in 1:m_k
            F_k[:, j] = w[end-m_k+j] .- w[end-m_k+j-1]
            E_k[:, j] = e[end-m_k+j]
        end

        gamma_k = F_k \ w[end]

        push!(x, x[end] .+ beta .* w[end] .- (E_k .+ beta .* F_k) * gamma_k)
        push!(e, x[end] .- x[end-1])

        if norm(e[end], 2)/norm(x[end], 2) < tol
            break
        end
    end

    Vx_final = x[end][1:N_total]
    Vy_final = x[end][N_total+1:end]

    return Vx_final, Vy_final, i
end