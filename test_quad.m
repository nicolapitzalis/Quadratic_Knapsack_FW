function test_quad(alg, n, actv, density, max_eig, trials)

    profile on
logs_gaps = [];
times = [];
iterations = [];

epsilon = 1e-16;
MaxIter = 5000;

differing_solutions_count = 0;
lmin = 0;
lmax = 0;

umin = 1;
umax = 1;

if ~ismember(alg, {'quadprog', 'fw_og', 'fw_away_step', 'fw_pairwise'})
    error('alg is not in ["fw_og", "fw_away_step", "fw_pairwise"]');
end

for i = 1:trials
    
    BCQP = gen_problem(n, actv, density, max_eig, lmin, lmax, umin, umax, i);
    
    if strcmp(alg, 'fw_og')
        [v, log, ~] = fw_og(BCQP, epsilon, MaxIter, false);
    elseif strcmp(alg, 'fw_away_step')
        [v, log, ~] = fw_away_step(BCQP, epsilon, MaxIter, false);
    elseif strcmp(alg, 'fw_pairwise')
        [v, log, ~] = fw_pairwise(BCQP, epsilon, MaxIter, false);
    end

    % replace negative values with 0
    log.gaps(log.gaps < 0.0) = 0;
    logs_gaps = [logs_gaps, log.gaps]; % TODO: rename
    times = [times, log.time];
    iterations = [iterations, log.iteration];
    
    % Plot single run log gap
    semilogy(1:length(log.gaps), log.gaps, 'Color', [0.5, 0.5, 0.5]);
    hold on;

    % Confront with quadprog
    [~, FVAL, ~, ~] = quadprog(2 * BCQP.Q, BCQP.q, -BCQP.a', -BCQP.b, [], [], BCQP.l, BCQP.u); 

    % Check if the difference between ours and gold solutions exceeds the threshold (1e-6)
    if abs(FVAL - v) > 1e-4
        disp(FVAL)
        disp(v)
        differing_solutions_count = differing_solutions_count + 1;
    end

end

disp('Differing solutions count:');
disp(differing_solutions_count);

mean_gaps = mean(logs_gaps, 2);
geo_mean_gaps = geomean(logs_gaps, 2);
mean_times = mean(times);
mean_iterations = mean(iterations);

% Plot mean log gap
semilogy(1:length(mean_gaps), mean_gaps);
hold off;

disp('Mean times:');
disp(mean_times);
disp('Mean iterations:');
disp(mean_iterations);

experiment_info = struct('alg', alg, 'n', n, 'actv', actv, 'density', density, 'max_eig', max_eig, 'trials', trials, 'mean_gaps', mean_gaps, 'geo_mean_gaps', geo_mean_gaps, 'mean_times', mean_times, 'mean_iterations', mean_iterations);
experiment_info_json = jsonencode(experiment_info, PrettyPrint=true);
% Save experiment info to a file
file_name = strcat(alg, '_', string(n), '_', string(actv), '_', string(density), '_', string(max_eig), '.json');
fid = fopen(file_name,'w');
% if fid == -1
%     error('Error opening file: %s', file_name);
% else
%     disp(['File opened successfully: ', file_name]);
% end
fprintf(fid,'%s',experiment_info_json);
fclose(fid);

end

