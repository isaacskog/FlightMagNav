close all; clc; 

% Dataset folders.
Root = fullfile('..','..');
rawRoot = fullfile(Root,'raw');
metadataFile = fullfile(Root,'metadata','flights.txt');
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


% Plot filter results
plot_filter_parameters(model.obs,model.filterInfo)

% Plot estimated map, including uncertinaty 
plot_mag_map(model,'model')

% Plot validation results
plot_filter_parameters(model.val_obs,model.validInfo)

for ii=1:numel(model.validInfo)
disp(model.validInfo(ii).sigma_validation)
end

% % Save the results
% timestamp = string(datetime('now','Format','yyyyMMdd_HHmmss'));
% save("GeneratedModel_" + timestamp + ".mat", "model");