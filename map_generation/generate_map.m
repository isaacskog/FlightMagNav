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




% Plot filter results
plot_filter_parameters(model)

% Plot estimated map, including uncertinaty 
plot_mag_map(model,'model')

% Plot map according to SGU
plot_mag_map(model,'sgu','sgudata.mat')
