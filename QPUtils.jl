using LinearAlgebra, Random, SparseArrays

export genQf

function genQf(n; seed::Integer=0, rank::Real=1.1, conv::Real=1, ecc::Real=0.99, dom::Real=1, box::Real=1)
    #=

        Parameters: 
        - seed: seed for the random number generator
        
        - rank (real, scalar, default: 1.1): if rank < 1 then the matrix is rank-deficient, if rank > 1 then the matrix is full-rank
        
        - conv (real, scalar, default: 1): if conv = 1 then the matrix positive (semi) definite, else if conv = 0 then the matrix is negative (semi) definite. Also the values in between will lead to a matrix that will tend to either positive or negative (semi) definiteness

        - ecc (real, scalar, default: 0.99): eccentricity of Q, i.e. the ratio (λ_1 - λ_n)/(λ_1 + λ_n), where λ_1 and λ_n are the largest and smallest eigenvalues of Q, respectively. Note that with ecc = 0, all the eigenvalues are equal, and as ecc -> 1 all the eigenvalues gets larger w.r.t. the smallest one. (This is computed only if the matrix is positive definite, i.e. if the smallest eigenvalue is positive, since otherwise that would lead to a fixed eccentricity of 1)

        -dom (real, scalar, default: 1): if dom = 1 then the matrix is unchanged, else if dom > 1 then the matrix is scaled up, and if dom < 1 then the matrix is scaled down w.r.t. the diagonal elements. This will lead to a diagonally dominant matrix if dom > 1, and to a diagonally non-dominant matrix if dom < 1

        - box (real, scalar, default: 1): if box > 0 then the problem is feasible, else if box < 0 then the problem is infeasible. This will lead to a feasible problem if box > 0, and to an infeasible problem if box < 0

        Outputs: 
        - Q: a symmetric matrix of size n x n
        - q: a vector of size n
        - v: a scalar, solution to the unconstrained minimum
        - a: a vector of size n, constraint coefficients
        - b: a scalar, constraint constant
        - l: a vector of size n, lower bounds for the box constraints on x
        - u: a vector of size n, upper bounds for the box constraints on x
    =#
    
    if !isa(n, Integer) || n <= 0
        error("n must be a positive integer")
    end

    Random.seed!(seed)

    # --------------------------------------------------------------------------
    # Generate Q
    r = round(Int, rank * n)
    p = round(Int, r * conv)
    Q = zeros(n, n)

    if p > 0
        G = rand(Float64, (p, n))
        Q += G' * G
    end
    if r > p
        G = rand(Float64, (r - p, n))
        Q -= G' * G
    end
    
    # Modify the diagonal if dom != 1
    # Increase or decrease randomly each element in the diagonal by [-1/3, 1/3] of its initial value, to then multiply it by dom
    if dom != 1
        D = diagm(0 => vec(diag(Q)))
        D = D .* (1 .+ (2 .* rand(n) .- 1) ./ 3)
        D = dom .* D
        Q = Q - diagm(0 => vec(diag(Q))) + D
    end
    
    Q = Symmetric(Q)
    vals, V = eigen(Q)
    
    # Compute eigenvalue decomposition and modify eccentricity
    if vals[1] > 1e-14
        D = diagm(0 => vals)
        λ_n = vals[1]
        λ_1 = vals[end]
        l = λ_n * ones(n) + ((λ_n / (λ_1 - λ_n))) * (2 * ecc / (1 - ecc)) * (vals .- λ_n)        
        Q = V * diagm(0 => l) * V'
    end
    # --------------------------------------------------------------------------
    
    # --------------------------------------------------------------------------
    # Generate q if required
    q = nothing
    x = 2 * abs(box) * rand(n) .- abs(box)
    q = -Q * x
    
    # Randomly modify q if box < 0 and Q is rank-deficient
    if box < 0 && vals[1] <= 1e-14
        q .*= (4 / 3) .* rand(n) .- 2 / 3
    end
    # --------------------------------------------------------------------------
    
    
    # --------------------------------------------------------------------------
    # Compute v if required
    v = nothing
    if vals[1] > 1e-14 || (vals[1] > -1e-14 && box > 0)
        v = dot(q, x) / 2
    else
        v = -Inf
    end
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Generate a, b, l, u  
    a = rand(n)
    while norm(a) == 0
        a = rand(n)
    end
    
    b = rand()

    l = rand(n)
    u = rand(n)

    u = l + u
    # --------------------------------------------------------------------------

    return Q, q, v, a, b, l, u
end