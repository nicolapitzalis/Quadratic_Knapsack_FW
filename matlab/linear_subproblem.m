function x_bar = linear_subproblem(BCQP, c)
    % Input:
    % a, c - vectors of coefficients
    % b - scalar
    % l, u - lower and upper bounds (vectors)
    
    n = length(BCQP.a);
    x_bar = zeros(n, 1);
    
    % Step 1-7: Initialize x_bar based on the sign of c
    for i = 1:n
        if c(i) >= 0
            x_bar(i) = BCQP.l(i);
        else
            x_bar(i) = BCQP.u(i);
        end
    end
    
    % Step 8-10: Check if current x_bar satisfies the condition
    if BCQP.a' * x_bar >= 0
        return;
    end
    
    % Step 11: Compute the permutation P that orders indices
    [~, P] = sort(c ./ BCQP.a);
    
    % Step 12-18: Adjust x_bar based on the permutation
    for p = P'
        x_bar(p) = BCQP.u(p);
        if a' * x_bar > BCQP.b
            x_bar(p) = (BCQP.b - sum(BCQP.a .* x_bar) + BCQP.a(p) * x_bar(p)) / BCQP.a(p);
            return;
        end
    end
end
