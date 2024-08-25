function x_bar = linear_subproblem(BCQP, c)
    % Input:
    % a, c - vectors of coefficients
    % b - scalar
    % l, u - lower and upper bounds (vectors)
    
    n = length(BCQP.a);
    x_bar = zeros(n, 1);
    
    % Initialize x_bar based on the sign of c
    for i = 1:n
        if c(i) >= 0
            x_bar(i) = BCQP.l(i);
        else
            x_bar(i) = BCQP.u(i);
        end
    end
    
    % Check if current x_bar satisfies the condition
    if BCQP.a' * x_bar >= b
        return;
    end
    
    % Compute the permutation P that orders indices
    [~, P] = sort(c ./ BCQP.a);
    
    % Adjust x_bar based on the permutation
    for p = P'
        x_bar(p) = BCQP.u(p);
        if a' * x_bar > BCQP.b
            x_bar(p) = (BCQP.b - sum(BCQP.a .* x_bar) + BCQP.a(p) * x_bar(p)) / BCQP.a(p);
            return;
        end
    end
end
