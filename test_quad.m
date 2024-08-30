function test_quad(alg, n, act, rank, ecc, trials)


logs_gaps = [];
times = [];
iterations = [];

epsilon = 1e-15;
MaxIter = 1000;

if ~ismember(alg, {'quadprog', 'fw_og', 'fw_away_step', 'fw_pairwise'})
    error('alg is not in ["fw_og", "fw_away_step", "fw_pairwise"]');
end

for i = 1:trials
    BCQP = genBCQP(n, act, rank, ecc, i);
    
    if strcmp(alg, 'quadprog')
        a = 1; % TODO
    elseif strcmp(alg, 'fw_og')
        [v, log, ~] = fw_og(BCQP, epsilon, MaxIter, false);
    elseif strcmp(alg, 'fw_away_step')
        [v, log, ~] = fw_away_step(BCQP, epsilon, MaxIter, false);
    elseif strcmp(alg, 'fw_pairwise')
        [v, log, ~] = fw_pairwise(BCQP, epsilon, MaxIter, false);
    end

    % [v, x, s] = fw_step_away ( BCQP );
    % [v, x, s] = fw_pairwise ( BCQP );
    logs_gaps = [logs_gaps, log.gaps]; % TODO: rename
    times = [times, log.time];
    iterations = [iterations, log.iteration];

end

mean_gaps = mean(logs_gaps, 2);
geo_mean_gaps = geomean(logs_gaps, 2);
mean_times = mean(times);
mean_iterations = mean(iterations);

disp('Mean times:');
disp(mean_times);
disp('Mean iterations:');
disp(mean_iterations);

experiment_info = struct('alg', alg, 'n', n, 'act', act, 'rank', rank, 'ecc', ecc, 'trials', trials, 'mean_gaps', mean_gaps, 'geo_mean_gaps', geo_mean_gaps, 'mean_times', mean_times, 'mean_iterations', mean_iterations);
experiment_info_json = jsonencode(experiment_info, PrettyPrint=true);
% save
file_name = strcat(alg, '_', string(n), '_', string(act), '_', string(rank), '_', string(ecc), '.json');
fid = fopen(file_name,'w');
fprintf(fid,'%s',experiment_info_json);
fclose(fid);

end

