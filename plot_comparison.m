function plot_comparison(experiment_path)

    % The `plot_comparison` function reads and compares the results of three Frank-Wolfe algorithm variants (fw_og, fw_away_step, fw_pairwise)
    % by plotting their geometric mean gaps across iterations. This function is useful for visualizing the performance differences between
    % the algorithms based on logs stored in the experiment_path.
    %
    % Syntax:
    % plot_comparison(experiment_path)
    %
    % Inputs:
    % - experiment_path: (string) The directory path containing the logs of the Frank-Wolfe algorithm variants. The logs must be stored in 
    %   JSON format within subdirectories of the `experiment_path`.
    %
    % Outputs:
    % - This function does not return any values, but it saves a comparison plot (PNG) of the geometric mean gaps for the three algorithms
    %   in the `experiment_path`.
    %
    % Example Usage:
    % experiment_path = 'results/convergence';
    % plot_comparison(experiment_path);

    files = dir(experiment_path);
    
    % Loop through each file in the folder
    for i = 1:length(files)
        % disp(files(i).name);
        % Check if the file is a JSON file and contains the desired substring
        if contains(files(i).name, 'fw_og') && files(i).isdir && ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..')
            % Construct the full file path
            filePath = fullfile(experiment_path, files(i).name);
            jsonFilePath = fullfile(filePath, 'log.json');
            % Read the JSON file content
            fileContent = fileread(jsonFilePath);
            % Decode the JSON content to a MATLAB structure
            fw_og_log = jsondecode(fileContent);
        end

        if contains(files(i).name, 'fw_away_step') && files(i).isdir && ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..')
            % Construct the full file path
            filePath = fullfile(experiment_path, files(i).name);
            jsonFilePath = fullfile(filePath, 'log.json');
            % Read the JSON file content
            fileContent = fileread(jsonFilePath);
            % Decode the JSON content to a MATLAB structure
            fw_away_step_log = jsondecode(fileContent);
        end

        if contains(files(i).name, 'fw_pairwise') && files(i).isdir && ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..')
            % Construct the full file path
            filePath = fullfile(experiment_path, files(i).name);
            jsonFilePath = fullfile(filePath, 'log.json');
            % Read the JSON file content
            fileContent = fileread(jsonFilePath);
            % Decode the JSON content to a MATLAB structure
            fw_pairwise_log = jsondecode(fileContent);
        end
    end

    % Plot the mean gaps
    fig = figure('Visible', 'off');
    semilogy(1:length(fw_og_log.geo_mean_gaps), fw_og_log.geo_mean_gaps, '-', 'LineWidth', 3, 'Color', 'blue', 'DisplayName', 'fw og');
    hold on;
    semilogy(1:length(fw_away_step_log.geo_mean_gaps), fw_away_step_log.geo_mean_gaps, '-', 'LineWidth', 3, 'Color', 'red', 'DisplayName', 'fw away step');
    hold on;
    semilogy(1:length(fw_pairwise_log.geo_mean_gaps), fw_pairwise_log.geo_mean_gaps, '-', 'LineWidth', 3, 'Color', 'green', 'DisplayName', 'fw pairwise');

    % Specify the axis names
    xlabel('Iterations'); 
    ylabel('Gap'); 

    % Add a title
    experiment_str = strrep(experiment_path, 'results/', ' ');
    experiment_str = strrep(experiment_str, '/', ': ');
    experiment_str = strrep(experiment_str, '_', '\_');
    plot_title = strcat('Comparison for ', experiment_str);
    title(plot_title);

    % Add a legend to distinguish the lines
    legend('show');

    % Grid for better visualization
    grid on;

    % Save the plot
    experiment_str = strrep(experiment_str, ': ', '_');
    experiment_str = strrep(experiment_str, '\_', '_');
    experiment_str = strrep(experiment_str, ' ', ''); % delete white spaces
    % experiment_str = strrep(experiment_str, ' ', '_');
    saveas(fig, strcat(experiment_path, '/', experiment_str, '_comparison_plot.png'));
