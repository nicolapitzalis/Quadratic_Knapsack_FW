function plot_comparison(experiment_path)

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
    plot_title = strcat('Comparison for ', experiment_str);
    title(plot_title);

    % Add a legend to distinguish the lines
    legend('show');
    % Save the plot
    saveas(fig, strcat(experiment_path, '/comparison_plot.png'));
