function run_experiments()
   
    base_n = 100; 
    base_actv = 0.7;
    base_density = 1.0;
    base_cond_num = 0.001;
    trials = 2;

    epsilon = 1e-8;
    max_iter = 5000;
    
    test_quad('fw_og', epsilon, max_iter, base_n, base_actv, base_density, base_cond_num, trials, 'results/convergence');
    test_quad('fw_away_step', epsilon, max_iter, base_n, base_actv, base_density, base_cond_num, trials, 'results/convergence');
    test_quad('fw_pairwise', epsilon, max_iter, base_n, base_actv, base_density, base_cond_num, trials, 'results/convergence');

    n_values = [50, 100, 500, 1000, 5000, 10000]; % Varying the size
    for n = n_values
        n_str = strrep(num2str(n), '.', '_');
        experiment_path = strcat('results/n/', n_str);
        test_quad('fw_og', epsilon, max_iter, n, base_actv, base_density, base_cond_num, trials, experiment_path);
        test_quad('fw_away_step', epsilon, max_iter, n, base_actv, base_density, base_cond_num, trials, experiment_path);
        test_quad('fw_pairwise', epsilon, max_iter, n, base_actv, base_density, base_cond_num, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    actv_values = [0.1, 0.3, 0.5, 0.7, 0.9]; % Varying the activation
    for actv = actv_values
        actv_str = strrep(num2str(actv), '.', '_');
        experiment_path = strcat('results/actv/', actv_str);
        test_quad('fw_og',epsilon, max_iter,  base_n, actv, base_density, base_cond_num, trials, experiment_path);
        test_quad('fw_away_step', epsilon, max_iter, base_n, actv, base_density, base_cond_num, trials, experiment_path);
        test_quad('fw_pairwise', epsilon, max_iter, base_n, actv, base_density, base_cond_num, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    density_values = [0.1, 0.3, 0.5, 0.7, 0.9]; % Varying the density
    for density = density_values
        density_str = strrep(num2str(density), '.', '_');
        experiment_path = strcat('results/density/', density_str);
        test_quad('fw_og', epsilon, max_iter, base_n, base_actv, density, base_cond_num, trials, experiment_path);
        test_quad('fw_away_step', epsilon, max_iter, base_n, base_actv, density, base_cond_num, trials, experiment_path);
        test_quad('fw_pairwise', epsilon, max_iter, base_n, base_actv, density, base_cond_num, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    cond_num_values = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]; % Varying the condition number
    for cond_num = cond_num_values
        cond_num_str = strrep(num2str(cond_num), '.', '_');
        experiment_path = strcat('results/cond_num/', cond_num_str);
        test_quad('fw_og', epsilon, max_iter,  base_n, base_actv, base_density, cond_num, trials, experiment_path);
        test_quad('fw_away_step', epsilon, max_iter, base_n, base_actv, base_density, cond_num, trials, experiment_path);
        test_quad('fw_pairwise', epsilon, max_iter, base_n, base_actv, base_density, cond_num, trials, experiment_path);
        plot_comparison(experiment_path);
    end