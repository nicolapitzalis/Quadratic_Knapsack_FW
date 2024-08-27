function x_bar = linear_subproblem(BCQP, c, x, t)
    % Input:
    % a, c - vectors of coefficients
    % b - scalar
    % l, u - lower and upper bounds (vectors)
    
    n = length(BCQP.a);
    x_bar = zeros(n, 1);
    
    if nargin > 2
        BCQP.l = max(BCQP.l, x - t);
        BCQP.u = min(BCQP.u, x + t);
    end

    % Initialize x_bar based on the sign of c
    for i = 1:n
        if c(i) >= 0
            x_bar(i) = BCQP.l(i);
        else
            x_bar(i) = BCQP.u(i);
        end
    end
    
    % Check if current x_bar satisfies the condition
    if BCQP.a' * x_bar >= BCQP.b
        return;
    end
    
    % Compute the permutation P that orders indices
    [~, P] = sort(c ./ BCQP.a);
    
    % Adjust x_bar based on the permutation
    for p = P'
        x_bar(p) = BCQP.u(p);
        if BCQP.a' * x_bar > BCQP.b
            x_bar(p) = (BCQP.b - sum(BCQP.a .* x_bar) + BCQP.a(p) * x_bar(p)) / BCQP.a(p);
            return;
        end
    end
end
