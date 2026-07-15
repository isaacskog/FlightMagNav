close all; clc; 

addpath(fullfile('..','..','Matlab','utils/'))
addpath(fullfile('..','..','Matlab','plottings/'))
% Load data
filename = fullfile('..','..','Matlab','SkovdeFlightData');
load(filename);

% Get settings
settings=get_settings();

% Fit model
model = fit_mag_map_model(data,settings);


% Plot heights
plot_trajectory_heights(model);

% Plot filter results
plot_filter_parameters(model.obs,model.filterInfo)

% Plot estimated map, including uncertinaty 
plot_mag_map(model,'model')

% Plot map according to SGU
plot_mag_map(model,'sgu','sgudata.mat')


% Plot validation results
plot_filter_parameters(model.val_obs,model.validInfo)


disp(model.validInfo.sigma_validation)

% Save the results
save("GeneratedModel","model");