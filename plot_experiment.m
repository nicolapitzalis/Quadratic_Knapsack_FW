function plot_experiment(paths)
    figure;
    
    for i = 1:length(paths)
        % Load the data
        data = fileread(paths{i});
        data = jsondecode(data);
        % Plot the data
        semilogy(1:length(data.mean_gaps), data.mean_gaps);
        legend(data.alg);
        hold on;
        xlabel('Iterations');
        ylabel('gap');
        grid on;
    end

    hold off;
end