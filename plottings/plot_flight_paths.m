function plot_flight_paths(data,settings,saveAsPdf)
%PLOT_FLIGHT_PATHS Plot map-building and validation flight paths.
%
% Usage:
%   plot_flight_paths(data,settings)
%   plot_flight_paths(data,settings,true)
%
% Inputs:
%   data
%       Struct array containing data(ii).GPSaidedINS.
%
%   settings
%       Structure containing
%           settings.idx
%           settings.idx_validation_data_set
%           settings.altitude_min
%
%   saveAsPdf
%       Logical flag. If true, the four figures are saved as publication-
%       ready PDF files in the folder specified below. Default is false.

    if nargin < 3 || isempty(saveAsPdf)
        saveAsPdf = false;
    end

    % ---------------------------------------------------------------------
    % Figure settings.
    % ---------------------------------------------------------------------
    fontSize = 13;
    lineWidth = 1.5;

    % Full-width IEEE figure dimensions [width height] in centimetres.
    pathFigureSize = [17.8 12];
    heightFigureSize = [pathFigureSize(1) pathFigureSize(2)/2];


    outputFolder =fullfile('..','..','Matlab','publication_figures');

    mapPathFile = fullfile(outputFolder,'map_building_flight_paths.pdf');
    mapHeightFile = fullfile(outputFolder,'map_building_flight_heights.pdf');
    validationPathFile = fullfile(outputFolder,'validation_flight_paths.pdf');
    validationHeightFile = fullfile(outputFolder,'validation_flight_heights.pdf');

    if saveAsPdf && ~isfolder(outputFolder)
        mkdir(outputFolder);
    end

    % ---------------------------------------------------------------------
    % Flight paths used for map building.
    % ---------------------------------------------------------------------
    figMapPath = figure( ...
        'Units','centimeters', ...
        'Position',[2 2 pathFigureSize], ...
        'Color','white');

    gx1 = geoaxes(figMapPath);
    hold(gx1,'on');

    for k = settings.idx
        lat = data(k).GPSaidedINS.latitude;
        lon = data(k).GPSaidedINS.longitude;
        valid = ~isnan(lat) & ~isnan(lon);

        geoplot(gx1,lat(valid),lon(valid),'LineWidth',lineWidth);
    end

    geobasemap(gx1,'topographic');
    title(gx1,'Flight trajectories used for map learning', ...
        'FontSize',fontSize,'FontWeight','normal');

    lgd = legend(gx1,"Flight " + string(settings.idx), ...
        'Location','northeast');
    lgd.FontSize = fontSize;

    gx1.FontSize = fontSize;

    drawnow;
    [latLim,lonLim] = geolimits(gx1);
    geotickformat -dd

    if saveAsPdf
        export_publication_pdf(gx1,mapPathFile);
    end

    % ---------------------------------------------------------------------
    % Heights of flights used for map building.
    % ---------------------------------------------------------------------
    figMapHeight = figure( ...
        'Units','centimeters', ...
        'Position',[2 2 heightFigureSize], ...
        'Color','white');

    ax1 = axes(figMapHeight);
    hold(ax1,'on');

    for k = settings.idx
        time = seconds( ...
            data(k).GPSaidedINS.time - data(k).GPSaidedINS.time(1));
        height = data(k).GPSaidedINS.altitude - ...
            data(k).GPSaidedINS.altitude(1);
        valid = ~isnan(height);

        plot(ax1,time(valid),height(valid),'LineWidth',lineWidth);
    end
    plot(ax1,[0 1100],settings.altitude_min*ones(1,2), ...
        'k','LineWidth',lineWidth);
    grid(ax1,'minor');
    xlabel(ax1,'Time [s]','FontSize',fontSize);
    ylabel(ax1,'Height [m]','FontSize',fontSize);
    title(ax1,'Height above ground versus time', ...
        'FontSize',fontSize,'FontWeight','normal');

    xlim(ax1,[0 1100]);
    ylim(ax1,[0 110]);
    ax1.FontSize = fontSize;
    box on;
    text(550,settings.altitude_min-10,'Height threshold map learning','FontSize',fontSize,'HorizontalAlignment','center','Color','k')

    drawnow;

    if saveAsPdf
        export_publication_pdf(ax1,mapHeightFile);
    end

    % ---------------------------------------------------------------------
    % Flight paths used for validation.
    % ---------------------------------------------------------------------
    figValidationPath = figure( ...
        'Units','centimeters', ...
        'Position',[2 2 pathFigureSize], ...
        'Color','white');

    gx2 = geoaxes(figValidationPath);
    hold(gx2,'on');

    for k = settings.idx_validation_data_set
        lat = data(k).GPSaidedINS.latitude;
        lon = data(k).GPSaidedINS.longitude;
        valid = ~isnan(lat) & ~isnan(lon);

        geoplot(gx2,lat(valid),lon(valid),'LineWidth',lineWidth);
    end

    geobasemap(gx2,'topographic');
    geolimits(gx2,latLim,lonLim);

    title(gx2,'Flight trajectories used for map validation', ...
        'FontSize',fontSize,'FontWeight','normal');

    lgd = legend(gx2, ...
        "Flight " + string(settings.idx_validation_data_set), ...
        'Location','northeast');
    lgd.FontSize = fontSize;

    gx2.FontSize = fontSize;

    geotickformat -dd
    drawnow;


    if saveAsPdf
        export_publication_pdf(gx2,validationPathFile);
    end

    % ---------------------------------------------------------------------
    % Heights of validation flights.
    % ---------------------------------------------------------------------
    figValidationHeight = figure( ...
        'Units','centimeters', ...
        'Position',[2 2 heightFigureSize], ...
        'Color','white');

    ax2 = axes(figValidationHeight);
    hold(ax2,'on');

    for k = settings.idx_validation_data_set
        time = seconds( ...
            data(k).GPSaidedINS.time - data(k).GPSaidedINS.time(1));
        height = data(k).GPSaidedINS.altitude - ...
            data(k).GPSaidedINS.altitude(1);
        valid = ~isnan(height);

        plot(ax2,time(valid),height(valid),'LineWidth',lineWidth);
    end

    grid(ax2,'minor');
    xlabel(ax2,'Time [s]','FontSize',fontSize);
    ylabel(ax2,'Height [m]','FontSize',fontSize);
    title(ax2,'Height above ground versus time', ...
        'FontSize',fontSize,'FontWeight','normal');

    xlim(ax2,[0 700]);
    ylim(ax2,[0 110]);
    ax2.FontSize = fontSize;
    box on;

    drawnow;

    if saveAsPdf
        export_publication_pdf(ax2,validationHeightFile);
    end
end


function export_publication_pdf(graphicsHandle,fileName)
%EXPORT_PUBLICATION_PDF Export tightly cropped publication-quality PDF.
%
% Vector output is used when possible. The resolution setting controls any
% raster content, such as geographic basemap tiles.

    exportgraphics(graphicsHandle,fileName, ...
        'ContentType','vector', ...
        'BackgroundColor','white', ...
        'Resolution',600);
end
