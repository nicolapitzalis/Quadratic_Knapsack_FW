function x_bar = linear_subproblem(BCQP, c, x, t)
    n = length(BCQP.a);
    x_bar = zeros(n, 1);
    
    % Compute the trust region bounds
    if nargin > 3
        l_t = max(BCQP.l, x - t .* BCQP.u);
        u_t = min(BCQP.u, x + t .* BCQP.u);
    else
        l_t = BCQP.l;
        u_t = BCQP.u;
    end
    
    % Initialize x_bar based on the sign of c
    for i = 1:n
        if c(i) >= 0
            x_bar(i) = l_t(i);
        else
            x_bar(i) = u_t(i);
        end
    end
    
    % Check if current x_bar satisfies the linear inequality constraint
    if BCQP.a' * x_bar >= BCQP.b
        return;
    end
    
    % Compute the permutation P that orders indices
    [~, P] = sort(c ./ BCQP.a);
    
    % Adjust x_bar based on the permutation
    for p = P'
        x_bar(p) = u_t(p);
        if BCQP.a' * x_bar > BCQP.b
            x_bar(p) = (BCQP.b - sum(BCQP.a .* x_bar) + BCQP.a(p) * x_bar(p)) / BCQP.a(p);
            return;
        end
    end
end
