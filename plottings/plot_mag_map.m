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
%       exported as a tightly cropped vector PDF for publication.
%       Default: false.

if nargin < 2 || isempty(saveFigures)
    saveFigures = false;
end

fontSize = 14;
gridRes = 5;
anomalyClim = [-70 70];
uncertaintyClim = [0.5 2];
maxMapStd = 1.2; % Maximum posterior standard deviation shown in map [nT]

outputFolder = fullfile( ...
    '..','..','Matlab','publication_figures');

% Collect flight-track points in local NED coordinates.
rTrack = [];
for ii = 1:numel(model.obs)
    rTrack = [rTrack; model.obs(ii).r_ned(:,1:2)]; %#ok<AGROW>
end

% Convex hull of the area covered by the mapping flights.
kHull = convhull(rTrack(:,1),rTrack(:,2));

% Create a local evaluation grid.
centers = model.basis.map.centers;

nMin = min(centers(:,1));
nMax = max(centers(:,1));
eMin = min(centers(:,2));
eMax = max(centers(:,2));

north = nMin:gridRes:nMax;
east  = eMin:gridRes:eMax;
[N,E] = meshgrid(north,east);

inside = inpolygon( ...
    N(:), ...
    E(:), ...
    rTrack(kHull,1), ...
    rTrack(kHull,2));

% Predict map mean and variance.
r = [N(:),E(:),zeros(numel(N),1)];
[B,sigma2] = model.predict_map(r);
sigma = sqrt(sigma2);

%% Learned magnetic anomaly map

% Only plot map values where the posterior standard deviation is below
% the selected threshold.
idxModel = ...
    inside & ...
    ~isnan(B) & ...
    ~isnan(sigma) & ...
    sigma < maxMapStd;

Bmodel = B(idxModel);
Bmodel = Bmodel - mean(Bmodel,'omitnan');

figModel = create_map_figure();
axModel = axes(figModel);
hold(axModel,'on');

scatter( ...
    axModel, ...
    E(idxModel), ...
    N(idxModel), ...
    10, ...
    Bmodel, ...
    'filled');

plot_flight_tracks(axModel,model.obs);

format_map_axes( ...
    axModel, ...
    'Learned magnetic anomaly map', ...
    'Magnetic anomaly [nT]', ...
    anomalyClim, ...
    fontSize, ...
    [eMin eMax], ...
    [nMin nMax]);

%% Posterior standard deviation

% The uncertainty figure is shown for analysis, but is not exported.
idxUncertainty = inside & ~isnan(sigma);

figUncertainty = create_map_figure();
axUncertainty = axes(figUncertainty);
hold(axUncertainty,'on');

scatter( ...
    axUncertainty, ...
    E(idxUncertainty), ...
    N(idxUncertainty), ...
    10, ...
    sigma(idxUncertainty), ...
    'filled');

plot_flight_tracks(axUncertainty,model.obs);

format_map_axes( ...
    axUncertainty, ...
    'Uncertainty in the learned magnetic map', ...
    'Posterior standard deviation [nT]', ...
    uncertaintyClim, ...
    fontSize, ...
    [eMin eMax], ...
    [nMin nMax]);

%% Export publication figure

if saveFigures
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    export_figure_pdf( ...
        figModel, ...
        fullfile(outputFolder,'learned_magnetic_anomaly_map.pdf'));
end

end


function fig = create_map_figure()
%CREATE_MAP_FIGURE Create a consistently sized publication figure.

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 16 14]);
end


function format_map_axes( ...
    ax,plotTitle,colorbarLabel,colorLimits,fontSize,xLimits,yLimits)
%FORMAT_MAP_AXES Apply common formatting to a magnetic-map figure.

xlabel(ax,'East [m]');
ylabel(ax,'North [m]');
title(ax,plotTitle);

axis(ax,'equal');
xlim(ax,xLimits);
ylim(ax,yLimits);
grid(ax,'on');
box(ax,'on');

ax.FontSize = fontSize;
ax.Layer = 'top';

colormap(ax,turbo);

if ~isempty(colorLimits)
    clim(ax,colorLimits);
end

cb = colorbar(ax);
cb.Label.String = colorbarLabel;
cb.FontSize = fontSize;
cb.Label.FontSize = fontSize;

drawnow;
end


function export_figure_pdf(fig,fileName)
%EXPORT_FIGURE_PDF Export a tightly cropped vector PDF.

exportgraphics( ...
    fig, ...
    fileName, ...
    'ContentType','vector', ...
    'BackgroundColor','white');
end


function plot_flight_tracks(ax,obs)
%PLOT_FLIGHT_TRACKS Overlay flight tracks.

for ii = 1:numel(obs)
    east  = obs(ii).r_ned(:,2);
    north = obs(ii).r_ned(:,1);

    plot( ...
        ax, ...
        east, ...
        north, ...
        'Color',[0.4 0.4 0.4], ...
        'LineWidth',0.6);
end
end
