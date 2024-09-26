function run_experiments()
   
    base_n = 100; 
    base_actv = 0.7;
    base_density = 0.9;
    base_rc = 0.1;
    trials = 100;

    epsilon = 1e-16;
    max_iter = 1500;
    
    test_quad_parallel('fw_og', epsilon, max_iter, base_n, base_actv, base_density, base_rc, trials, 'results/convergence');
    test_quad_parallel('fw_away_step', epsilon, max_iter, base_n, base_actv, base_density, base_rc, trials, 'results/convergence');
    test_quad_parallel('fw_pairwise', epsilon, max_iter, base_n, base_actv, base_density, base_rc, trials, 'results/convergence');
    plot_comparison('results/convergence');

    n_values = [10, 50, 100, 200, 500]; % Varying the size
    for n = n_values
        n_str = num2str(n);
        experiment_path = strcat('results/n/', n_str);
        test_quad_parallel('fw_og', epsilon, max_iter, n, base_actv, base_density, base_rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, n, base_actv, base_density, base_rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, n, base_actv, base_density, base_rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    actv_values = [0.1, 0.3, 0.5, 0.7, 0.9]; % Varying the activation
    for actv = actv_values
        actv_str = num2str(actv);
        experiment_path = strcat('results/actv/', actv_str);
        test_quad_parallel('fw_og',epsilon, max_iter,  base_n, actv, base_density, base_rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, base_n, actv, base_density, base_rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, base_n, actv, base_density, base_rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    density_values = [0.1, 0.3, 0.5, 0.7, 0.9]; % Varying the density
    for density = density_values
        density_str = num2str(density);
        experiment_path = strcat('results/density/', density_str);
        test_quad_parallel('fw_og', epsilon, max_iter, base_n, base_actv, density, base_rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, base_n, base_actv, density, base_rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, base_n, base_actv, density, base_rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    rc_values = [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001]; % Varying the condition number
    for rc = rc_values
        rc_str = num2str(rc);
        experiment_path = strcat('results/rc/', rc_str);
        test_quad_parallel('fw_og', epsilon, max_iter,  base_n, base_actv, base_density, rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, base_n, base_actv, base_density, rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, base_n, base_actv, base_density, rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    big_n_values = [500, 1000];
    density_for_big = 0.1;
    for n = big_n_values
        n_str = num2str(n);
        experiment_path = strcat('results/n_sparse/', n_str);
        test_quad_parallel('fw_og', epsilon, max_iter, n, base_actv, density_for_big, base_rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, n, base_actv, density_for_big, base_rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, n, base_actv, density_for_big, base_rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    huge_n_values = [5000, 10000];
    density_for_huge = 0.01;
    for n = huge_n_values
        n_str = num2str(n);
        experiment_path = strcat('results/n_very_sparse/', n_str);
        test_quad_parallel('fw_og', epsilon, max_iter, n, base_actv, density_for_huge, base_rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, n, base_actv, density_for_huge, base_rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, n, base_actv, density_for_huge, base_rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    huge_n_values = [5000, 10000];
    density_smaller_for_huge = 0.001;
    for n = huge_n_values
        n_str = num2str(n);
        experiment_path = strcat('results/n_extremely_sparse/', n_str);
        test_quad_parallel('fw_og', epsilon, max_iter, n, base_actv, density_smaller_for_huge, base_rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, n, base_actv, density_smaller_for_huge, base_rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, n, base_actv, density_smaller_for_huge, base_rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    big_n_values = [500, 1000];
    density_for_big = 0.1;
    actv = 1.0;
    for n = big_n_values
        n_str = num2str(n);
        experiment_path = strcat('results/n_sparse_all_active/', n_str);
        test_quad_parallel('fw_og', epsilon, max_iter, n, actv, density_for_big, base_rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, n, actv, density_for_big, base_rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, n, actv, density_for_big, base_rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    huge_n_values = [5000, 10000];
    density_for_huge = 0.01;
    actv = 1.0;
    for n = huge_n_values
        n_str = num2str(n);
        experiment_path = strcat('results/n_very_sparse_all_active/', n_str);
        test_quad_parallel('fw_og', epsilon, max_iter, n, actv, density_for_huge, base_rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, n, actv, density_for_huge, base_rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, n, actv, density_for_huge, base_rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    huge_n_values = [5000, 10000];
    density_smaller_for_huge = 0.001;
    actv = 1.0;
    for n = huge_n_values
        n_str = strrep(num2str(n), '.', '_');
        experiment_path = strcat('results/n_extremely_sparse_all_active/', n_str);
        test_quad_parallel('fw_og', epsilon, max_iter, n, actv, density_smaller_for_huge, base_rc, trials, experiment_path);
        test_quad_parallel('fw_away_step', epsilon, max_iter, n, actv, density_smaller_for_huge, base_rc, trials, experiment_path);
        test_quad_parallel('fw_pairwise', epsilon, max_iter, n, actv, density_smaller_for_huge, base_rc, trials, experiment_path);
        plot_comparison(experiment_path);
    end

    winner_density = 0.1;
    winner_n = 1000;
    winner_max_iter = 100000;
    test_quad_parallel('fw_og', epsilon, winner_max_iter, winner_n, base_actv, winner_density, base_rc, trials, 'results/winner');
    test_quad_parallel('fw_away_step', epsilon, winner_max_iter, winner_n, base_actv, winner_density, base_rc, trials, 'results/winner');
    test_quad_parallel('fw_pairwise', epsilon, winner_max_iter, winner_n, base_actv, winner_density, base_rc, trials, 'results/winner');
    plot_comparison('results/winner');


