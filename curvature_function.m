function cf = curvature_function(problem)

    % The curvature_function computes an approximation of the curvature of the quadratic objective function associated with a convex optimization problem.
    % The curvature is a measure of how much the objective function deviates from its linear approximation and plays a critical role in determining the convergence rates of optimization algorithms.
    % This function uses optimization to maximize the curvature over a sampled region of the feasible space.
    %
    % Syntax:
    % cf = curvature_function(problem)
    %
    % Inputs:
    % - problem: (struct) A structure representing the optimization problem. It must contain the following fields:
    %   - Q: (matrix) The quadratic matrix in the objective function.
    %   - q: (vector) The linear term in the objective function.
    %   - l: (vector) Lower bounds for the variables.
    %   - u: (vector) Upper bounds for the variables.
    %   - a: (vector) Coefficients of the inequality constraint a'x >= b.
    %   - b: (scalar) The right-hand side of the inequality constraint a'x >= b.
    %
    % Outputs:
    % - cf: (scalar) The approximated curvature of the objective function based on the optimization of the curvature function.
    %
    % Example Usage:
    % problem.Q = rand(5,5); problem.Q = problem.Q'*problem.Q; 
    % problem.q = rand(5,1);
    % problem.l = zeros(5,1);
    % problem.u = ones(5,1);
    % problem.a = rand(5,1);
    % problem.b = 1;
    % cf = curvature_function(problem);
    %
    % Notes:
    % - The curvature is computed by solving a non-linear optimization problem using `fmincon`, where the objective is to maximize the curvature function over random initial points and directions.
    % - The function also ensures that the linear inequality constraint a'x >= b is satisfied for both the starting point x and the direction s.
    % - It is used for theoretical bound computation in optimization algorithms.

    n = size(problem.Q, 1);

    % Define the function f(x) = x'Qx + q'x
    f_x = @(x) x' * problem.Q * x + problem.q' * x;

    % Define the gradient of f(x)
    grad_f = @(x) 2 * problem.Q * x + problem.q;

    % Curvature function to maximize
    curv_f = @(x, s, gamma) ...
        2 / gamma^2 * (f_x(x + gamma * (s - x)) - f_x(x) - (x + gamma * (s - x) - x)' * grad_f(x));

    % Optimization options
    options = optimoptions('fmincon', 'Algorithm', 'sqp');

    % Initial guesses for x, s, and gamma
    x0 = problem.l + (problem.u - problem.l) .* rand(size(problem.l));
    s0 = problem.l + (problem.u - problem.l) .* rand(size(problem.l));

    % Ensure the linear constraint a'x >= b is satisfied
    if problem.a' * x0 < problem.b
        lambda = (problem.b - problem.a' * x0) / norm(problem.a)^2;
        x0 = x0 + lambda * problem.a;
    end
    if problem.a' * s0 < problem.b
        lambda = (problem.b - problem.a' * s0) / norm(problem.a)^2;
        s0 = s0 + lambda * problem.a;
    end

    % Initial guess for gamma
    gamma0 = 0.5;

    % Define the full variable vector: [x0; s0; gamma]
    initial_guess = [x0; s0; gamma0];

    % Define the objective function for fmincon (negative because we are maximizing)
    objective = @(var) -curv_f(var(1:n), var(n+1:2*n), var(2*n+1));

    % Define the linear inequality constraint
    A = [-problem.a' zeros(1, n) 0;  % Constraint for x0
        zeros(1, n) -problem.a' 0];  % Constraint for s0
    b_ineq = [-problem.b; -problem.b];

    % Define the bounds: [l; l; 0] <= [x; s; gamma] <= [u; u; 1]
    l = [problem.l; problem.l; 0];
    u = [problem.u; problem.u; 1];

    % Solve the optimization problem using fmincon
    [~, fval] = fmincon(objective, initial_guess, A, b_ineq, [], [], l, u, [], options);
    cf = -fval;
end
