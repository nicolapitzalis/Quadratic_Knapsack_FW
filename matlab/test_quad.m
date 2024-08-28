diff = zeros(100, 1);
differing_solutions_count = 0;

for i = 1:100
    BCQP = genBCQP(100, 0.6, 1.1, 0.99, i);
    
    [X,FVAL,EXITFLAG,OUTPUT] = quadprog(2 * BCQP.Q, BCQP.q, -BCQP.a', -BCQP.b, [], [], BCQP.l, BCQP.u);
    [v, x, s] = FWBCQP ( BCQP );
    
    % Compute the norm of the difference between the two solutions
    diff(i) = abs(FVAL - v);

    % Check if the difference exceeds the threshold (1e-6)
    if diff(i) > 1e-6
        disp(FVAL);
        disp(v);
        differing_solutions_count = differing_solutions_count + 1;
    end
end

% Display the results
disp('Average difference norm over 100 experiments:');
disp(mean(diff));
disp('Number of experiments with differing solutions (≥1e-6):');
disp(differing_solutions_count);
