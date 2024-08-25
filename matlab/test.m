tries = 1000;
failed = 0;

for i = 0:tries-1
    rng(i);  % Set the random seed

    % Set random_bool = true directly as per the Python code
    random_bool = true;

    % Generate data (assuming gen_CQKP is a defined function in MATLAB)
    data = genBCQP(100, 0.4, 1.1, 0.99, i, 8, 12, 1, 5); %, random_bool);

    % Define variables
    Q = data.Q;
    q = data.q;
    a = data.a;
    b = data.b;
    l = data.l;
    u = data.u;
    n = length(q);

    % Convert the problem to a standard form for quadprog
    H = Q; % This is the matrix for the quadratic term
    f = q; % This is the vector for the linear term
    A = -a; % We need to negate the 'a' because quadprog minimizes A*x <= b
    B = -b; % Corresponding negated 'b'
    
    % Define bounds
    lb = l; % Lower bounds
    ub = u; % Upper bounds
    
    % Solve the problem using quadprog
    options = optimoptions('quadprog', 'Display', 'off', 'TolFun', 1e-6);
    [x, fval, exitflag] = quadprog(H, f, A, B, [], [], lb, ub, [], options);

    gold_solution_value = fval;
    gold_solution_x = x;

    % Call the CQKP_FW function (assuming it is defined in MATLAB)
    [y, x_opt, status] = CQKP_FW(data, 'verbose', false, 'eps', 1e-6, 'MaxIter', 1000);
    our_solution_value = y;
    our_solution_x = x_opt;

    % Compare solutions
    if abs(gold_solution_value - our_solution_value) > 1e-6 % || norm(gold_solution_x - our_solution_x, inf) > 1e-6
        fprintf('Test %d failed\n', i);
        fprintf('Gold solution: %f, exitflag: %d\n', gold_solution_value, exitflag);
        fprintf('Our solution: %f, status: %s\n', our_solution_value, status);
        failed = failed + 1;
    end
end

if failed == 0
    fprintf('All %d tests passed\n', tries);
else
    fprintf('%d out of %d tests failed\n', failed, tries);
end