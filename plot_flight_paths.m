function plot_flight_paths(data,idx)
%PLOT_FLIGHT_PATHS Plot flight paths for selected data sets.
%
% Usage:
%   plot_flight_paths(data,1:4)
%
% Input:
%   data - struct array with field data(ii).GPSaidedINS
%   idx  - indices of flights to plot, e.g. 1:4
%
% Assumes:
%   data(ii).GPSaidedINS.latitude
%   data(ii).GPSaidedINS.longitude

figure;
gx = geoaxes;
hold(gx,'on');

for k = idx
    lat = data(k).GPSaidedINS.latitude;
    lon = data(k).GPSaidedINS.longitude;

    valid = ~isnan(lat) & ~isnan(lon);

    geoplot(gx,lat(valid),lon(valid),'LineWidth',1.5);

    % Label close to start point
    if any(valid)
        ii0 = find(valid,1,'first');
        text(gx,lat(ii0),lon(ii0),sprintf('%d',k), ...
            'FontSize',12, ...
            'FontWeight','bold');
    end
end

geobasemap(gx,'topographic');
title(gx,'Flight paths');
legend(gx,"Flight " + string(idx),'Location','best');
end