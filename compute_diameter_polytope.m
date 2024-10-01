function diameter = solve_diameter_polytope(problem)
diameter = norm(problem.u - problem.l);
return;