function plot_mag_map(model,mapSource,sguMatFile)
%PLOT_MAG_MAP Plot estimated magnetic map or SGU magnetic data.
%
% Usage:
%   plot_mag_map(model)
%   plot_mag_map(model,'model')
%   plot_mag_map(model,'sgu','SGU_data.mat')
%
% Inputs:
%   model
%       Output from fit_mag_map_model().
%
%   mapSource
%       'model' : plot the estimated map from model.predict_map().
%       'sgu'   : plot SGU data loaded from sguMatFile.
%
%   sguMatFile
%       MAT-file containing variable y with columns:
%           y(:,1) = projected easting  [m]
%           y(:,2) = projected northing [m]
%           y(:,3) = magnetic field     [nT]

if nargin < 2 || isempty(mapSource)
    mapSource = 'model';
end

mapSource = lower(string(mapSource));

% Collect flight-track points in local NED.
rTrack = [];

for ii = 1:numel(model.obs)
    rTrack = [rTrack; model.obs(ii).r_ned(:,1:2)]; %#ok<AGROW>
end

% Convex hull of mapped area in local NED coordinates.
kHull = convhull(rTrack(:,1),rTrack(:,2));

% Use same local grid for model and SGU map.
centers = model.basis.map.centers;

nMin = min(centers(:,1));
nMax = max(centers(:,1));
eMin = min(centers(:,2));
eMax = max(centers(:,2));

gridRes = 5;

north = nMin:gridRes:nMax;
east  = eMin:gridRes:eMax;

[N,E] = meshgrid(north,east);

% Grid points inside flight-covered area.
inside = inpolygon( ...
    N(:), ...
    E(:), ...
    rTrack(kHull,1), ...
    rTrack(kHull,2));

figure;
clf;
hold on;

switch mapSource

    case "model"

        r = [N(:),E(:),zeros(numel(N),1)];

        [B,sigma2] = model.predict_map(r);

        idx = inside & ~isnan(B);

        Bplot = B(idx);
        Bplot = Bplot - mean(Bplot,'omitnan');


        scatter( ...
            E(idx), ...
            N(idx), ...
            10, ...
            Bplot, ...
            'filled');

        cb = colorbar;
        cb.Label.String = 'Predicted map [nT]';

        title('Estimated magnetic map');
        clim([-70 70])


        colormap turbo

        % Overlay flight tracks.
        plot_flight_tracks(model.obs)

        xlabel('East [m]');
        ylabel('North [m]');

        axis equal
        grid on



        figure
        clf;
        hold on;
        sigma=sqrt(sigma2(idx));

        scatter( ...
            E(idx), ...
            N(idx), ...
            10, ...
            sigma, ...
            'filled');

        cb = colorbar;
        cb.Label.String = 'Uncertinaty predicted map, mean removed [nT]';


        colormap turbo

        % Overlay flight tracks.
        plot_flight_tracks(model.obs)

        xlabel('East [m]');
        ylabel('North [m]');

        axis equal
        grid on
         title('Uncertinaty in estimated magnetic map');

    case "sgu"

        if nargin < 3 || isempty(sguMatFile)
            error('For SGU plotting, provide sguMatFile.');
        end

        S = load(sguMatFile);
        y = S.y;

        % Convert SGU projected coordinates to local east/north.
        refXY = get_reference_sgu_xy(model);

        eastSGU  = y(:,1) - refXY(1);
        northSGU = y(:,2) - refXY(2);
        Bsgu     = y(:,3);

        % Keep only SGU points in a bounding box around the mapped area.
        % This avoids building scatteredInterpolant from the full SGU
        % data set, which can be slow.
        bboxMargin = 2*model.settings.map.center_spacing;

        nBoxMin = min(rTrack(:,1)) - bboxMargin;
        nBoxMax = max(rTrack(:,1)) + bboxMargin;
        eBoxMin = min(rTrack(:,2)) - bboxMargin;
        eBoxMax = max(rTrack(:,2)) + bboxMargin;

        idxBox = ...
            northSGU >= nBoxMin & northSGU <= nBoxMax & ...
            eastSGU  >= eBoxMin & eastSGU  <= eBoxMax;

        eastSGU  = eastSGU(idxBox);
        northSGU = northSGU(idxBox);
        Bsgu     = Bsgu(idxBox);

        % Build SGU interpolant using only nearby SGU data.
        F = scatteredInterpolant( ...
            eastSGU, ...
            northSGU, ...
            Bsgu, ...
            'linear', ...
            'none');

        % Interpolate SGU data only at grid points inside the mapped
        % area. This avoids evaluating the interpolant on unnecessary
        % grid points.
        Einside = E(inside);
        Ninside = N(inside);

        Binside = F(Einside,Ninside);

        idx = ~isnan(Binside);

        Bplot = Binside(idx);
        Bplot = Bplot - mean(Bplot,'omitnan');

        scatter( ...
            Einside(idx), ...
            Ninside(idx), ...
            10, ...
            Bplot, ...
            'filled');

        cb = colorbar;
        cb.Label.String = 'SGU magnetic field, mean removed [nT]';

        title('Low resolution SGU magnetic map');
        clim([-70 70])


        colormap turbo


        % Overlay flight tracks.
        plot_flight_tracks(model.obs)

        xlabel('East [m]');
        ylabel('North [m]');

        axis equal
        grid on



    otherwise

        error('Unknown mapSource. Use ''model'' or ''sgu''.');
end




end


function refXY = get_reference_sgu_xy(model)
%GET_REFERENCE_SGU_XY Project reference_lla to SGU coordinate system.
%
% Uses the user's sweref() function.
%
% Input:
%   model.settings.reference_lla = [lat0 lon0 alt0]
%
% Output:
%   refXY = [E0 N0]

lat0 = model.settings.reference_lla(1);
lon0 = model.settings.reference_lla(2);

[NArray,EArray,~] = sweref(lat0,lon0,0);

refXY = [EArray NArray];
end


function plot_flight_tracks(obs)
%PLOT_FLIGHT_TRACKS Plot flight tracks with distinct dark colors and labels.

nFlight = numel(obs);



for ii = 1:nFlight

    east  = obs(ii).r_ned(:,2);
    north = obs(ii).r_ned(:,1);

    plot( ...
        east, ...
        north, ...
        'k', ...
        'LineWidth',1.5);

    % Put label near the middle of the track.
    kk = round(numel(east)/2);

    text( ...
        east(kk), ...
        north(kk), ...
        sprintf('%d',ii), ...
        'FontWeight','bold', ...
        'FontSize',10, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'BackgroundColor','w', ...
        'Margin',1);

end
end