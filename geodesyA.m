%% Script: Analyze Potential Dependency on Height
%% Author: Adapted to User's Scientific Question
%% Description:
%% Computes and compares gravitational potential values for different heights above the ellipsoid.
%% Outputs plots showing potential variations and differences with height.

%% Workspace Reset
clc;
clear;
close all;

%% Add Required Paths
addpath('./synthesis_m/cptcmap');
addpath('./synthesis_m');
addpath('./variants');
addpath('./ICGEM');
addpath('./results');

%% User Input
institution = 'GFZ';
ddk = 'DDK5';

% Input file paths
pathVARIANTS = './variants/';
pathICGEM = './ICGEM/';
pathRESULTS = './results/';

filenameIN = ['input' institution ddk '.txt'];
list = strsplit(fileread([pathVARIANTS institution '/' filenameIN]));

% Constants
deg2rad = pi/180; % Conversion factor
nmin = 2; % Minimum degree
h_values = [0, 500e3, 1000e3, 1500e3]; % Heights in meters
R = 6378136.3; % Earth's radius in meters
GM = 0.3986004415e+15; % Earth's gravitational constant
a = 6378137; % Semi-major axis
e2 = 0.00669438002290; % Square of first eccentricity

% Grid setup
DlatitudeDEG = 1.0;
DlongitudeDEG = 1.0;
latitude_V_max = 90 - DlatitudeDEG / 2;
latitude_V_min = -90 + DlatitudeDEG / 2;
longitude_V_min = -180 + DlongitudeDEG / 2;
longitude_V_max = 180 - DlongitudeDEG / 2;

latitudeDEG_V = (latitude_V_max:-DlatitudeDEG:latitude_V_min)';
longitudeDEG_V = (longitude_V_min:DlongitudeDEG:longitude_V_max)';

% Convert to radians
latitude_V = latitudeDEG_V * deg2rad;
longitude_V = longitudeDEG_V * deg2rad;

% Initialize storage for potentials
potentials = cell(1, length(h_values));

%% Load Coefficients and Calculate Potentials for Each Height
disp('Loading models and computing potentials...');

modelname1 = list{1}; % First model
filename1 = [pathICGEM institution '/' ddk '/' modelname1];

% Read spherical harmonic coefficients for the first model
[Cnm, Snm, ~, ~, GM1, R1, nmax1, ~, ~, ~] = ReadCoefficientsICGEM(filename1, 0, 0);
[Cnm_rescaled, Snm_rescaled] = rescaleCnm(Cnm, Snm, GM, R, GM1, R1, nmax1);

% Compute potential for each height
for idx = 1:length(h_values)
    h_V = h_values(idx);
    disp(['Calculating potential for height: ', num2str(h_V), ' meters']);
    potentials{idx} = SHS_grid_ell(latitude_V, longitude_V, Cnm_rescaled, Snm_rescaled, ...
                                   GM, h_V, a, e2, nmin, nmax1);
end

disp('Potential computation complete.');

%% Calculate Differences Between Heights
potential_differences = cell(1, length(h_values) - 1);
for i = 1:length(h_values) - 1
    potential_differences{i} = potentials{i+1} - potentials{i};
end

%% Plot Potentials for Different Heights
figure;
hold on;
legendEntries = {}; % For custom legend entries
for idx = 1:length(h_values)
    plot(latitudeDEG_V, potentials{idx}, 'DisplayName', ...
         ['Height = ' num2str(h_values(idx)) ' m']);
    legendEntries{end+1} = ['Height = ' num2str(h_values(idx)) ' m'];
end
title('Potential vs. Height');
xlabel('Latitude (degrees)');
ylabel('Potential (m^2/s^2)');
legend(legendEntries, 'Location', 'Best'); % Use custom legend entries
grid on;
hold off;

%% Plot Differences Between Heights
figure;
hold on;
legendEntries = {}; % Reset legend entries
for idx = 1:length(potential_differences)
    plot(latitudeDEG_V, potential_differences{idx}, 'DisplayName', ...
         ['Difference: ' num2str(h_values(idx)) 'm to ' num2str(h_values(idx+1)) 'm']);
    legendEntries{end+1} = ['Difference: ' num2str(h_values(idx)) 'm to ' num2str(h_values(idx+1)) 'm'];
end
title('Potential Differences vs. Height');
xlabel('Latitude (degrees)');
ylabel('Potential Difference (m^2/s^2)');
legend(legendEntries, 'Location', 'Best'); % Use custom legend entries
grid on;
hold off;

%% Save Results
outputFolder = './results/';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

for idx = 1:length(h_values)
    temp_potential = potentials{idx}; % Store data temporarily
    save([outputFolder 'potential_height_' num2str(h_values(idx)) '.mat'], 'temp_potential');
end

disp('Results saved.');
disp('Script execution completed.');
