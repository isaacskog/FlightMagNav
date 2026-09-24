close all; clc; 

% Dataset folders.
Root = fullfile('..','..');
parsedRoot = fullfile(Root,'parsed');

% Support functions
addpath(fullfile(Root,'matlab/utils'));
addpath(fullfile(Root,'matlab/plottings'));

% Load data
filename = fullfile(parsedRoot,'FlightMagNav');
load(filename);

% Get settings
settings=get_settings();

% Fit model
model = fit_mag_map_model(data,settings);

% Plot flight paths
plot_flight_paths(data,settings)


% Plot estimated map
plot_mag_map(model)

fprintf('Map negative log evidence: %.3f\n',model.negative_log_evidence);
for ii = 1:numel(model.validInfo)
    fprintf('Flight %d: validation noise %.3f nT\n', ...
        model.validInfo(ii).flight_index,model.validInfo(ii).sigma_validation);
end
