tries = 10;
failed = 0;

for seed = 0:tries-1
    rng(seed); % Set the seed for reproducibility
    n = 5;
    data = genBCQP(n, 0.4, 1.1, 0.99, seed, 8, 12, 1, 5);
    % gen randomly c
    c = -1 + 2 * rand(n, 1);

    % Define the constraints for linprog
    A = -data.a';
    b_ = -data.b;
    lb = data.l;
    ub = data.u;

    % Solve using the custom linear_subproblem function
    solution = linear_subproblem(data, c);

    % Solve using MATLAB's linprog function
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex-highs', 'Display', 'off');
    [golden, ~, exitflag] = linprog(c, A, b_, [], [], lb, ub, options);

    % Check if the solutions are close
    if exitflag == 1
        if ~all(abs(solution - golden) < 1e-6)
            fprintf('Failed for seed %d\n', seed);
            fprintf('Golden solution: %f\n', golden);
            fprintf('Our solution: %f\n', solution);
            failed = failed + 1;
        end
    else
        fprintf('Solver did not converge for seed %d\n', seed);
        fprintf('Failed for seed %d\n', seed);
        fprintf('Golden solution: %f\n', golden);
        fprintf('Our solution: %f\n', solution);
        failed = failed + 1;
    end
end

fprintf('Failed %d out of %d tries\n', failed, tries);
