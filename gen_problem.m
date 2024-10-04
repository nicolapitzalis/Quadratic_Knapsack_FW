function problem = gen_problem(n, actv, density, rc, lmin, lmax, umin, umax, seed)
    
    % The gen_problem function generates a quadratic programming problem with inequality constraints in the form:
    % min_x  x'Qx + q'x
    % subject to:
    % l <= x <= u   and   a'x >= b
    %
    % Syntax: 
    % problem = gen_problem(n, actv, density, rc, lmin, lmax, umin, umax, seed)
    %
    % Inputs:
    % - n: (integer) Size of the problem.
    % - actv: (float) A value between 0 and 1 specifying the fraction of variables that should be initially out of the bounds (either lower or upper bounds).
    % - density: (float) The density of the sparse matrix Q (percentage of non-zero entries).
    % - rc: (float) Reciprocal condition number for the matrix Q, controlling its conditioning.
    % - lmin: (float) Minimum value for the lower bound l on the variables.
    % - lmax: (float) Maximum value for the lower bound l on the variables.
    % - umin: (float) Minimum value for the upper bound u on the variables.
    % - umax: (float) Maximum value for the upper bound u on the variables.
    % - seed: (integer) Seed for the random number generator to ensure reproducibility of results.
    
    % Outputs:
    % - problem: (struct) A struct containing the following fields representing the quadratic knapsack problem:
    %   - l: (vector) Lower bounds for each variable l, with values randomly generated between `lmin` and `lmax`.
    %   - u: (vector) Upper bounds for each variable u, with values randomly generated between `umin` and `umax`.
    %   - Q: (sparse matrix) Symmetric, sparse matrix Q with specified density and conditioning.
    %   - q: (vector) Linear term q in the objective function, generated to create a problem where some variables are initially out of bounds.
    %   - a: (vector) Coefficients of the equality constraint a'x >= b.
    %   - b: (scalar) Right-hand side of the equality constraint b, ensuring a feasible solution within the bounds.
    
    % Example Usage:
    % n = 10;  % Number of variables
    % actv = 0.3;  % 30% of the variables are out of bounds
    % density = 0.2;  % 20% of Q's entries are non-zero
    % rc = 0.01;  % Reciprocal condition number
    % lmin = 0; lmax = 1;  % Lower bounds between 0 and 1
    % umin = 1; umax = 2;  % Upper bounds between 1 and 2
    % seed = 42;  % Random seed
    % problem = gen_problem(n, actv, density, rc, lmin, lmax, umin, umax, seed);
    
    % Notes:
    % - This function is useful for testing quadratic programming solvers, as it ensures that we can generate a variety of problems with different characteristics.
    % - The random seed input ensures the problem is reproducible when running multiple experiments.
    
    % See Also:
    % - `sprandsym`: MATLAB function for generating sparse, symmetric random matrices.

    rng(seed);
    % generate l, u- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    problem.l = lmin * ones( n, 1 ) + ( lmax - lmin ) * rand( n, 1 );
    problem.u = umin * ones( n , 1 ) + ( umax - umin ) * rand( n , 1 );

    % generate Q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    problem.Q = sprandsym(n, density, rc, 1);

    % generate q- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    z = zeros(n, 1);

    % outb( i ) = true if z( i ) will be out of the bounds
    outb = rand(n, 1) <= actv;

    % 50/50 chance of being left of lb or right of ub
    lr = rand(n, 1) <= 0.5;
    l = outb & lr;
    r = outb & ~lr;

    % a random amount left of the lb (0)
    z(l) = problem.l(l) .* (rand(sum(l), 1) - 1);

    % a random amount right of the ub (u)
    z(r) = problem.u(r) .* (1 + rand(sum(r), 1));

    outb = ~outb;  % entries that will be inside the bound
    % pick at random in [0, u]
    z(outb) = rand(sum(outb), 1) .* problem.u(outb);

    problem.q = -2 * problem.Q * z;

    % generate ax=b - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    problem.a = rand(n, 1);
    x0 = problem.l + (problem.u - problem.l) .* rand(size(problem.l));
    problem.b = problem.a' * x0;

end
