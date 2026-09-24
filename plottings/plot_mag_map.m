function plot_mag_map(model,saveFigures)
%PLOT_MAG_MAP Plot the learned magnetic anomaly map and its uncertainty.
%
% Usage:
%   plot_mag_map(model)
%   plot_mag_map(model,true)
%
% Inputs:
%   model
%       Output from fit_mag_map_model().
%
%   saveFigures
%       Logical flag. If true, only the learned magnetic anomaly map is
%       exported as a tightly cropped PDF for publication.
%       Default: false.

if nargin < 2 || isempty(saveFigures)
    saveFigures = false;
end

fontSize = 13;
gridRes = 5;
markerSize = 10;
mapAlpha = 0.85;
anomalyClim = [-70 70];
uncertaintyClim = [0.5 2];
maxMapStd = 1.2; % Maximum posterior standard deviation shown in map [nT]
baseMap = 'topographic';
pathFigureSize = [17.8 12];

outputFolder = fullfile( ...
    fileparts(fileparts(mfilename('fullpath'))),'publication_figures');

% Create a local evaluation grid.
centers = model.basis.map.centers;

north = min(centers(:,1)):gridRes:max(centers(:,1));
east  = min(centers(:,2)):gridRes:max(centers(:,2));
[N,E] = meshgrid(north,east);

% Predict map mean and variance.
rNed = [N(:),E(:),zeros(numel(N),1)];
[B,sigma2] = model.predict_map(rNed);
sigma = sqrt(sigma2);

% Convert the evaluation grid to latitude and longitude.
[latitude,longitude] = local_ned_to_lla( ...
    rNed,model.settings.reference_lla);

%% Learned magnetic anomaly map

% Only plot map values where the posterior standard deviation is below
% the selected threshold.
idxModel = ...
    ~isnan(B) & ...
    ~isnan(sigma) & ...
    sigma < maxMapStd;

Bmodel = B(idxModel);

figModel = create_map_figure(pathFigureSize);
gxModel = geoaxes(figModel);
hold(gxModel,'on');

geobasemap(gxModel,baseMap);

geoscatter( ...
    gxModel, ...
    latitude(idxModel), ...
    longitude(idxModel), ...
    markerSize, ...
    Bmodel, ...
    'filled', ...
    'MarkerEdgeColor','none', ...
    'MarkerFaceAlpha',mapAlpha);

plot_flight_tracks( ...
    gxModel,model.obs,model.settings.reference_lla);

format_geoaxes( ...
    gxModel, ...
    'Learned magnetic anomaly map', ...
    'Magnetic anomaly [nT]', ...
    anomalyClim, ...
    fontSize);

geotickformat -dd
%% Posterior standard deviation

% The uncertainty figure is shown for analysis, but is not exported.
idxUncertainty = ~isnan(sigma);

figUncertainty = create_map_figure(pathFigureSize);
gxUncertainty = geoaxes(figUncertainty);
hold(gxUncertainty,'on');

geobasemap(gxUncertainty,baseMap);

geoscatter( ...
    gxUncertainty, ...
    latitude(idxUncertainty), ...
    longitude(idxUncertainty), ...
    markerSize, ...
    sigma(idxUncertainty), ...
    'filled', ...
    'MarkerEdgeColor','none', ...
    'MarkerFaceAlpha',mapAlpha);

plot_flight_tracks( ...
    gxUncertainty,model.obs,model.settings.reference_lla);

format_geoaxes( ...
    gxUncertainty, ...
    'Uncertainty in the learned magnetic map', ...
    'Posterior standard deviation [nT]', ...
    uncertaintyClim, ...
    fontSize);

geotickformat -dd

%% Export publication figure

if saveFigures
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    export_figure_pdf( ...
        gxModel, ...
        fullfile(outputFolder,'learned_magnetic_anomaly_map.pdf'));
end

end


function fig = create_map_figure(pathFigureSize)
%CREATE_MAP_FIGURE Create a consistently sized publication figure.

fig = figure( ...
        'Units','centimeters', ...
        'Position',[2 2 pathFigureSize], ...
        'Color','white');
end


function format_geoaxes( ...
    gx,plotTitle,colorbarLabel,colorLimits,fontSize)
%FORMAT_GEOAXES Apply common formatting to a geographic map figure.

title( ...
    gx, ...
    plotTitle, ...
    'FontSize',fontSize, ...
    'FontWeight','normal');

gx.FontSize = fontSize;

colormap(gx,turbo);
clim(gx,colorLimits);

cb = colorbar(gx);
cb.Label.String = colorbarLabel;
cb.FontSize = fontSize;
cb.Label.FontSize = fontSize;
cb.FontWeight = 'normal';

drawnow;
end


function export_figure_pdf(gx,fileName)
%EXPORT_FIGURE_PDF Export a tightly cropped publication-quality PDF.

exportgraphics( ...
    gx, ...
    fileName, ...
    'ContentType','vector', ...
    'BackgroundColor','white', ...
    'Resolution',300);
end


function plot_flight_tracks(gx,obs,referenceLLA)
%PLOT_FLIGHT_TRACKS Overlay flight tracks in latitude and longitude.

for ii = 1:numel(obs)
    [latitude,longitude] = local_ned_to_lla( ...
        obs(ii).r_ned,referenceLLA);

    valid = ~isnan(latitude) & ~isnan(longitude);

    geoplot( ...
        gx, ...
        latitude(valid), ...
        longitude(valid), ...
        'Color',[0.4 0.4 0.4], ...
        'LineWidth',0.6);
end
end


function [latitude,longitude,altitude] = ...
    local_ned_to_lla(rNed,referenceLLA)
%LOCAL_NED_TO_LLA Convert local NED coordinates to WGS84 coordinates.
%
% This is the inverse of the small-area conversion used when learning
% the map.

R = 6378137;

lat0 = referenceLLA(1);
lon0 = referenceLLA(2);
alt0 = referenceLLA(3);

north = rNed(:,1);
east = rNed(:,2);

if size(rNed,2) >= 3
    down = rNed(:,3);
else
    down = zeros(size(north));
end

latitude = lat0 + rad2deg(north/R);
longitude = lon0 + rad2deg(east/(R*cosd(lat0)));
altitude = alt0 - down;
end
