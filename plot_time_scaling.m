exp = 'rc';
% Data from the table
% values = [50, 100, 200, 500];
% fw_og = [0.0756, 0.1600, 0.4511, 5.5692];
% fw_away_step = [0.3048, 0.4505, 1.0227, 6.8633];
% fw_pairwise = [0.3840, 0.7664, 1.5065, 7.9847];

values = [10^-1, 10^-2, 10^-3, 10^-4, 10^-5];
fw_og = [0.2007, 0.1937, 0.2002, 0.2012, 0.2198];
fw_away_step = [0.5771, 0.5789, 0.6231, 0.6751, 0.8159];
fw_pairwise = [0.9470, 1.0260, 1.2220, 1.5697, 2.2624];

% Create the plot
figure;
hold on;
plot(values, fw_og, '-o', 'Color', 'blue', 'DisplayName', 'fw og', 'LineWidth', 2);
plot(values, fw_away_step, '-o', 'Color', 'red', 'DisplayName', 'fw away step', 'LineWidth', 2);
plot(values, fw_pairwise, '-o', 'Color', 'green', 'DisplayName', 'fw pairwise', 'LineWidth', 2);

% Set the y-axis to log scale
set(gca, 'XScale', 'log');

% Labels and title
xlabel(exp);
ylabel('Time (seconds)');
title(strcat('Times varying: ', exp));

% Add a legend
legend('Location', 'southeast'); % Change position of legend bottom right
legend('show');

% Grid for better visualization
grid on;

hold off;
% Save the plot
saveas(gcf, strcat('plot_time_', exp, '.png'));
