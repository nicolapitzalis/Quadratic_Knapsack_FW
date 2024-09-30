function rho = compute_row(problem)
    mu = min(eig(problem.Q));
    L = max(eig(problem.Q));
    delta = computePyramidalWidth(problem.a, problem.b, problem.l, problem.u, 1000);
    D = solve_diameter_polytope(problem);
    rho = (mu / 4*L) * (delta / D)^2;
    return;