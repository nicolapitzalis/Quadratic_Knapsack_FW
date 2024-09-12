function cf = curvature_function(problem)
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
