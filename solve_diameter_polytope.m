% Define the function and its gradient
syms x1 x2;  % Example for two variables, extend as needed
f_sym = x1^2 + x2^2;  % Example function f(x)
grad_f_sym = gradient(f_sym, [x1, x2]);
f = matlabFunction(f_sym, 'Vars', {[x1, x2]});
grad_f = matlabFunction(grad_f_sym, 'Vars', {[x1, x2]});

% Define the objective function to maximize (convert to minimization)
objective = @(vars) -compute_supremum(vars, f, grad_f);  % Negate for maximization

% Initial guess for x, s, and gamma
initial_guess = [0, 0, 0, 0, 0.5];  % Example initial guess [x1, x2, s1, s2, gamma]

% Define bounds for the optimization (x, s in D and gamma in [0, 1])
lb = [-1, -1, -1, -1, 0];  % Lower bounds (x, s in D, gamma in [0,1])
ub = [1, 1, 1, 1, 1];  % Upper bounds (x, s in D, gamma in [0,1])

% Call fmincon to find the maximum
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
solution = fmincon(objective, initial_guess, [], [], [], [], lb, ub, [], options);

% Display result
disp('Optimized x, s, gamma:');
disp(solution);
disp('Maximum value of the supremum:');
disp(-objective(solution));  % Negate back to get maximum value

% Function to compute the supremum expression for fmincon
function val = compute_supremum(vars, f, grad_f)
    x = vars(1:2);  % Extract x from vars
    s = vars(3:4);  % Extract s from vars
    gamma = vars(5);  % Extract gamma from vars
    
    y = x + gamma * (s - x);  % Compute y
    
    % Compute the terms
    term1 = f(y) - f(x);
    term2 = dot(y - x, grad_f(x));
    
    % Return the value to be maximized (inside the supremum)
    val = (2 / gamma^2) * (term1 - term2);
end
