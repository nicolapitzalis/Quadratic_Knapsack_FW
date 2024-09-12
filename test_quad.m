function test_quad(alg, epsilon, max_iter, n, actv, density, cond_num, trials, path)

    % profile on

    logs_gaps = [];
    logs_primals = [];
    times = [];
    iterations = [];

    differing_solutions_count = 0;
    lmin = 0;
    lmax = 0;

    umin = 1;
    umax = 1;

    fig = figure('Visible', 'off');

    if ~ismember(alg, {'quadprog', 'fw_og', 'fw_away_step', 'fw_pairwise'})
        error('alg is not in ["fw_og", "fw_away_step", "fw_pairwise"]');
    end

    for i = 1:trials
        
        P = gen_problem(n, actv, density, cond_num, lmin, lmax, umin, umax, i);
        
        % Define the options for quadprog with the desired tolerance
        options = optimoptions('quadprog', 'TolFun', 1e-16);
        % Call quadprog with the new options
        
        [~, golden_val, ~, ~] = quadprog(2 * P.Q, P.q, -P.a', -P.b, [], [], P.l, P.u, [], options);
        if strcmp(alg, 'fw_og')
            [val, ~, ~, log] = fw_og(P, epsilon, max_iter, golden_val, false);
        elseif strcmp(alg, 'fw_away_step')
            [val, ~, ~, log] = fw_away_step(P, epsilon, max_iter, golden_val, false);
        elseif strcmp(alg, 'fw_pairwise')
            [val, ~, ~, log] = fw_pairwise(P, epsilon, max_iter, golden_val, false);
        end

        % replace negative values with 0
        log.gaps(log.gaps < 0.0) = 0;
        logs_gaps = [logs_gaps, log.gaps]; % TODO: rename
        logs_primals = [logs_primals, log.primals];
        times = [times, log.time];
        iterations = [iterations, log.iteration];
        
        
        semilogy(1:length(log.gaps), log.gaps, '--', 'LineWidth', 1, 'Color', [0 0 0 0.2]);
        hold on;

        % Confront with quadprog solution
        % Check if the difference between ours and gold solutions exceeds the threshold
        if abs(golden_val - val) > 1e-6
            disp(golden_val)
            disp(val)
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
    semilogy(1:length(mean_gaps), mean_gaps, '-', 'LineWidth', 3, 'Color', 'blue');
    hold off;

    disp('Mean times:');
    disp(mean_times);
    disp('Mean iterations:');
    disp(mean_iterations);

    experiment_info = struct('alg', alg, 'n', n, 'actv', actv, 'density', density, 'cond_num', cond_num, 'trials', trials, 'mean_gaps', mean_gaps, 'geo_mean_gaps', geo_mean_gaps, 'mean_times', mean_times, 'mean_iterations', mean_iterations);
    experiment_info_json = jsonencode(experiment_info, PrettyPrint=true);
    % Save experiment info to a file
    experiment_path = strcat(string(alg), '_', string(n), '_', string(actv), '_', string(density), '_', string(cond_num));
    save_path = strcat(path, '/', experiment_path);

    if ~exist(save_path, 'dir')
        mkdir(save_path); % Create the full path if it doesn't exist
    end
    
    % Save plot
    xlabel('Iterations'); 
    ylabel('Gap');
    saveas(fig, strcat(save_path, '/plot.png'));
    % Save experiment info
    experiment_info_file_name = strcat(save_path, '/log.json');
    fid = fopen(experiment_info_file_name, 'w');
    % if fid == -1
    %     error('Error opening file: %s', file_name);
    % else
    %     disp(['File opened successfully: ', file_name]);
    % end
    fprintf(fid,'%s',experiment_info_json);
    fclose(fid);

end

