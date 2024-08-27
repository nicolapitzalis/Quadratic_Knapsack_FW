difference_norms = zeros(100, 1);
differing_solutions_count = 0;

for i = 1:1000
    n = 1000;
    BCQP = genBCQP(100, 0.2, 1.1, 0.99999999);
    c = randn(100, 1);
    
    x_linprog = linprog(c, -BCQP.a', -BCQP.b, [], [], BCQP.l, BCQP.u);
    my_sol = linear_subproblem(BCQP, c); 
    
    difference_norms(i) = norm(x_linprog - my_sol);

    if norm(x_linprog - my_sol) > 1e-6
        differing_solutions_count = differing_solutions_count + 1;
    end
end

disp('Average difference norm over 1000 experiments:');
disp(mean(difference_norms));
disp('Number of experiments with differing solutions (≥1e-6):');
disp(differing_solutions_count);
