function test_quad(alg, epsilon, max_iter, n, actv, density, rc, trials, bound_flag, path)

    % The test_quad function runs multiple trials to test different variants of the Frank-Wolfe algorithm on 
    % instances of the convex quadratic knapsack problem. It collects performance metrics such as dual gaps, 
    % computational times, and iterations, comparing the results with the `quadprog` solver in MATLAB. 
    % The function supports parallel execution to speed up trials and logs results, plots, and information 
    % about the experiment to files.
    %
    % Syntax:
    % test_quad(alg, epsilon, max_iter, n, actv, density, rc, trials, bound_flag, path)
    %
    % Inputs:
    % - alg: (string) The optimization algorithm to test, one of: 'fw_og', 'fw_away_step', 'fw_pairwise', or 'quadprog'.
    % - epsilon: (float) The convergence tolerance for the Frank-Wolfe variants.
    % - max_iter: (integer) Maximum number of iterations allowed for the algorithms.
    % - n: (integer) Problem size.
    % - actv: (float) Fraction of active constraints in the problem, specifying how many variables are initially out of bounds.
    % - density: (float) Density of the quadratic matrix Q, determining the percentage of non-zero elements.
    % - rc: (float) Reciprocal condition number of the matrix Q.
    % - trials: (integer) Number of trials to run for statistical robustness.
    % - bound_flag: (boolean) If true, plots the dual gap along with upper bounds during the experiment.
    % - path: (string) Directory path where the results, plots, and logs will be saved.
    %
    % Outputs:
    % This function does not return any output directly. Instead, it generates several output files:
    % - A JSON file containing experiment information such as algorithm name, mean times, iterations, and dual gaps.
    % - A plot visualizing the dual gaps (and upper bounds if specified) over the iterations for all trials.
    % - A log showing the number of differing solutions between the tested algorithm and `quadprog`.
    
    % Example Usage:
    % alg = 'fw_og';         
    % epsilon = 1e-6;       
    % max_iter = 1500;      
    % n = 100;              
    % actv = 0.5;            
    % density = 0.2;         
    % rc = 0.01;             
    % trials = 10;           
    % bound_flag = true;     
    % path = 'results'; 
    % test_quad(alg, epsilon, max_iter, n, actv, density, rc, trials, bound_flag, path);

    % Notes:
    % - The function leverages MATLAB’s parallel computing capabilities to speed up the trials.
    % - The random seed for each trial is based on the trial number, ensuring reproducibility.
    % - The use of quadpro` as a reference allows for validation of the accuracy of the tested algorithms.
    %
    % See Also:
    % - gen_problem: Generates instances of the convex quadratic knapsack problem.
    % - fw_og, fw_away_step, fw_pairwise: Frank-Wolfe algorithm variants for solving quadratic programming problems.
    % - quadprog: MATLAB's built-in solver for quadratic programming problems.
    
    % profile on
    if isempty(gcp('nocreate'))
        parpool;  % Start the parallel pool if it's not already started
    end

    logs_gaps = zeros(max_iter+1, trials);
    logs_bounds = zeros(max_iter+1, trials);
    times = zeros(1, trials);
    iterations = zeros(1, trials);
    differing_solutions = zeros(1, trials);

    lmin = 0;
    lmax = 0;

    umin = 1;
    umax = 1;

    if ~ismember(alg, {'quadprog', 'fw_og', 'fw_away_step', 'fw_pairwise'})
        error('alg is not in ["fw_og", "fw_away_step", "fw_pairwise"]');
    end
    
    tic % start timer
    parfor i = 1:trials
        
        P = gen_problem(n, actv, density, rc, lmin, lmax, umin, umax, i);
        
        % Define the options for quadprog with the desired tolerance
        options = optimoptions('quadprog', 'TolFun', 1e-16);
        % Call quadprog with the new options
        
        [~, golden_val, ~, ~] = quadprog(2 * P.Q, P.q, -P.a', -P.b, [], [], P.l, P.u, [], options);
        if strcmp(alg, 'fw_og')
            [val, ~, ~, log] = fw_og(P, epsilon, max_iter, false, bound_flag);
        elseif strcmp(alg, 'fw_away_step')
            [val, ~, ~, log] = fw_away_step(P, epsilon, max_iter, false, bound_flag);
        elseif strcmp(alg, 'fw_pairwise')
            [val, ~, ~, log] = fw_pairwise(P, epsilon, max_iter, false, bound_flag);
        end

        % replace negative values with 0
        log.gaps(log.gaps < 0.0) = 0;
        logs_gaps(:, i) = log.gaps;
        logs_bounds(:, i) = log.bounds;
        times(i) = log.time;
        iterations(i) = log.iteration;

        % Confront with quadprog solution
        % Check if the difference between ours and gold solutions exceeds the threshold
        if abs(golden_val - val) > 1e-6
            differing_solutions(i) = 1;
        end

    end
    disp(toc);

    disp('Differing solutions count:');
    disp(sum(differing_solutions));

    mean_gaps = mean(logs_gaps, 2);
    geo_mean_gaps = geomean(logs_gaps, 2);
    mean_times = mean(times);
    mean_iterations = mean(iterations);
    geo_mean_bounds = geomean(logs_bounds, 2);

    disp('Mean times:');
    disp(mean_times);
    disp('Mean iterations:');
    disp(mean_iterations);

    experiment_info = struct('alg', alg, 'n', n, 'actv', actv, 'density', density, 'rc', rc, 'trials', trials, 'mean_gaps', mean_gaps, 'geo_mean_gaps', geo_mean_gaps, 'mean_times', mean_times, 'mean_iterations', mean_iterations);
    experiment_info_json = jsonencode(experiment_info, PrettyPrint=true);
    % Save experiment info to a file
    experiment_path = strcat(string(alg), '_', string(n), '_', string(actv), '_', string(density), '_', string(rc));
    save_path = strcat(path, '/', experiment_path);

    if ~exist(save_path, 'dir')
        mkdir(save_path); % Create the full path if it doesn't exist
    end
    
    % Plotting 
    fig = figure('Visible', 'off');

    if bound_flag
        semilogy(1:length(geo_mean_gaps), geo_mean_gaps, '-', 'LineWidth', 3, 'Color', 'blue', 'DisplayName', 'Dual Gap (geomean)');
        hold on;
        semilogy(1:length(geo_mean_bounds), geo_mean_bounds, '--', 'LineWidth', 3, 'Color', 'blue', 'DisplayName', 'Upper Bound (geomean)');
        hold off;
        legend('Location','southwest');
    else
        semilogy(1:length(logs_gaps), logs_gaps, '--', 'LineWidth', 1, 'Color', [0 0 0 0.2]);
        hold on;
        semilogy(1:length(geo_mean_gaps), geo_mean_gaps, '-', 'LineWidth', 3, 'Color', 'blue', 'DisplayName', 'Dual Gap (geomean)');
        hold off;
    end
    
    % Add a legend to distinguish the lines
    legend('show');
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

