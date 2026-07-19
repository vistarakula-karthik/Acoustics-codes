% MATLAB Code to read a 10x11 CSV of SPL data, calculate Directivity Index, and plot.

clear; clc; close all;

%% 1. Define Your Setup
% REPLACE these 11 numbers with the exact 11 angles you extracted from the x-axis
angles = [26, 36, 46, 56, 66, 76, 86, 96, 106, 126, 146]; 

% The 10 frequencies from the legend (in Hz)
frequencies = [10, 20, 50, 100, 200, 501, 1000, 1995, 5012, 10000]; 

num_freqs = length(frequencies);
num_angles = length(angles);

%% 2. Read the CSV File
filename = 'NASAsp8072.csv';

% Check if the file exists before trying to read it
if ~isfile(filename)
    error('File "%s" not found. Make sure it is in the same folder as this script.', filename);
end

% Load the 10x11 data matrix
SPL = readmatrix(filename);

% Verify dimensions just to be safe
[num_rows, num_cols] = size(SPL);
if num_rows ~= 10 || num_cols ~= 11
    error('The CSV file is %dx%d. It must be exactly 10 rows by 11 columns.', num_rows, num_cols);
end

%% 3. Calculate Directivity Index (DI)
% Step A: Convert SPL (dB) to linear acoustic energy
linear_energy = 10.^(SPL / 10);

% Step B: Calculate the spatial average energy across all 11 angles 
% (mean() across dimension 2 averages the columns for each row)
mean_energy = mean(linear_energy, 2); 

% Step C: Convert the averaged linear energy back to dB
average_SPL = 10 * log10(mean_energy);

% Step D: Subtract the average SPL from the original SPL at each angle
DI = SPL - repmat(average_SPL, 1, num_angles);

%% 4. Plot the Directivity Index (DI)
figure('Name', 'Directivity Index Plot', 'Position', [200, 200, 600, 500]);
hold on;

% Define line markers to mimic the provided image
markers = {'-d', '-s', '-^', '-x', '-*', '-o', '-+', '-_', '-p', '-h'};

% Plot each frequency curve
for i = 1:num_freqs
    plot(angles, DI(i, :), markers{i}, 'LineWidth', 1.5, 'MarkerSize', 5);
end
hold off;

% Apply exact formatting to match the reference plot
xlabel('Angle from Exhaust Axis (degrees)', 'FontWeight', 'bold');
ylabel('Directivity Index (dB)', 'FontWeight', 'bold');
xlim([20 160]);
ylim([-20 10]);

% Generate legend labels dynamically and format the legend
leg_labels = arrayfun(@(x) sprintf('%d Hz', x), frequencies, 'UniformOutput', false);
legend(leg_labels, 'NumColumns', 2, 'Location', 'northeast', 'Box', 'off');

% Polish the axes text size
set(gca, 'FontSize', 11);