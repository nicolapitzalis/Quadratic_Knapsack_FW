function x_bar = LP(problem, c)

    % The LP function solves a linear subproblem used in the Frank-Wolfe algorithm, where the goal is to minimize a 
    % linear objective subject to box constraints and a linear inequality constraint in the form: 
    % min_x  c'x
    % subject to:
    % l <= x <= u   and   a'x >= b
    % This function forms part of the inner loop of the Frank-Wolfe algorithm and is used to determine a feasible solution that minimizes the 
    % linearized objective function.
    %
    % Syntax:
    % x_bar = LP(problem, c)
    %
    % Inputs:
    % - problem: (struct) Contains the following fields representing the quadratic knapsack problem:
    %   - l: (vector) Lower bounds for each variable.
    %   - u: (vector) Upper bounds for each variable.
    %   - a: (vector) Coefficients of the linear constraint a'x >= b.
    %   - b: (scalar) Right-hand side of the linear constraint.
    % - c: (vector) Coefficients of the linear objective function c'x that needs to be minimized in this subproblem.
    
    % Outputs:
    % - x_bar: (vector) The solution to the linear subproblem, a vector that satisfies the box constraints and 
    % the linear inequality.
    
    % Example Usage:
    % Define problem parameters
    % problem.l = [0; 0; 0];
    % problem.u = [1; 1; 1];
    % problem.a = [0.5; 0.7; 0.9];
    % problem.b = 1.5;
    % c = [-1; 2; -3];
    % x_bar = LP(problem, c);
    
    % See Also:
    % - Any Frank-Wolfe Algorithm: This function is typically used within the Frank-Wolfe algorithm as part of 
    % the iterative optimization process.
    
    n = length(problem.a);
    x_bar = zeros(n, 1);

    % Initialize x_bar based on the sign of c
    for i = 1:n
        if c(i) >= 0
            x_bar(i) = problem.l(i);
        else
            x_bar(i) = problem.u(i);
        end
    end
    
    % Check if current x_bar satisfies the linear inequality constraint
    if problem.a' * x_bar >= problem.b
        return;
    end
    
    % Compute the permutation P that orders indices
    [~, P] = sort(c ./ problem.a);
    
    % Adjust x_bar based on the permutation
    for p = P'
        x_bar(p) = problem.u(p);
        if problem.a' * x_bar > problem.b
            x_bar(p) = (problem.b - sum(problem.a .* x_bar) + problem.a(p) * x_bar(p)) / problem.a(p);
            return;
        end
    end
end
