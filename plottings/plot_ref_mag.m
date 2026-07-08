
function plot_ref_mag(data,model,settings,doFiltering)
%PLOT_REF_MAG Plot reference magnetometer data in tiled subplots.
%
% Usage:
%   plot_ref_mag(data,1:4)
%   plot_ref_mag(data,1:4,true)
%
% Inputs:
%   data
%       Struct array created by the main parsing script.
%
%   idx
%       Indices of the flights to plot, e.g. 1:4.
%
%   doFiltering
%       Optional logical flag.
%
%       If false, the raw parsed ref_mag data is plotted.
%
%       If true, the signal is:
%         1) low-pass filtered at 0.05 Hz,
%         2) downsampled/resampled to 1 Hz.
%
%       This mode is intended for extracting/visualizing only the slow
%       temporal variation in the reference magnetometer. This is useful
%       when the reference sensor is only used to compensate for slow
%       background-field variations when creating magnetic maps.
%
% Assumes:
%   data(ii).ref_mag.time
%   data(ii).ref_mag.tot_field
%
% Notes:
%   - No notch/comb filtering is applied.
%   - The function first tries to use butter/filtfilt if available.
%   - If Signal Processing Toolbox is unavailable, it falls back to a
%     zero-phase moving-average low-pass approximation.

if nargin < 3
    doFiltering = false;
end

nPlot = numel(settings.idx)+1;

nCols = ceil(sqrt(nPlot));
nRows = ceil(nPlot/nCols);

figure;
tl = tiledlayout(nRows,nCols, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for ii = 1:nPlot

    if ii<nPlot
        k = settings.idx(ii);
        obs=model.obs(ii);
    else
        k=settings.idx_validation_data_set;
        obs=model.val_obs;
    end

    nexttile;
    hold on;
    grid on;

    t = data(k).ref_mag.time;
    B = data(k).ref_mag.tot_field;
    t_ups=data(k).UPS_ref_mag.time;
    B_ups=data(k).UPS_ref_mag.tot_field;
    


    valid = ~isnat(t) & ~isnan(B);
    t = t(valid);
    B = B(valid);

    if doFiltering
        [tPlot,BPlot] = filter_ref_mag_for_plotting(t,B);

        idx=tPlot>obs.time(1) & tPlot<obs.time(end);
        tPlot=tPlot(idx);
        BPlot=BPlot(idx);


        idx=t_ups>obs.time(1) & t_ups<obs.time(end);
        t_ups=t_ups(idx);
        B_ups=B_ups(idx);

        % Plot relative variation to make slow drift easier to compare.
        plot(tPlot,BPlot - BPlot(1),'k','LineWidth',1.2);
         plot(t_ups,B_ups - B_ups(1),'r','LineWidth',1.2);
        title(sprintf('Flight %d, LP 0.05 Hz, 1 Hz',k));
        ylabel('\Delta field [nT]');


    else
        plot(t,B,'LineWidth',1.2);
        title(sprintf('Flight %d',k));
        ylabel('Field [nT]');
    end

    if ii > (nRows-1)*nCols
        xlabel('Time (UTC)');
    end
    ylim([-5 5])
    legend('Reference magnetometers','Uppsala reference magnetometer')
end

if doFiltering
    title(tl,'Reference magnetometer, slow variation');
else
    title(tl,'Reference magnetometer');
end
end


function [tOut,BOut] = filter_ref_mag_for_plotting(t,B)
%FILTER_REF_MAG_FOR_PLOTTING Extract slow reference-field variation.
%
% Processing:
%   1) Low-pass filter at 0.05 Hz.
%   2) Resample/downsample to 1 Hz.
%
% Output:
%   tOut - datetime vector at 1 Hz
%   BOut - low-pass filtered magnetic field [nT]

t = t(:);
B = B(:);

% Work in seconds relative to first sample.
ts = seconds(t - t(1));

% Remove duplicate timestamps, if any.
[ts,idxUnique] = unique(ts,'stable');
t = t(idxUnique);
B = B(idxUnique);

% Estimate input sample rate.
fsIn = 1/median(diff(ts),'omitnan');

if ~isfinite(fsIn) || fsIn <= 0
    error('Could not estimate input sample rate.');
end

% Fill short NaN gaps before filtering.
valid = ~isnan(B);
if ~all(valid)
    B = interp1(ts(valid),B(valid),ts,'linear','extrap');
end

% ---------------------------------------------------------------------
% 1) Low-pass filter at 0.05 Hz.
% ---------------------------------------------------------------------
fc = 0.05;          % Hz
filterOrder = 3;

if fsIn <= 2*fc
    warning('Input sample rate %.2f Hz is too low for %.3f Hz low-pass filtering. Skipping low-pass.',fsIn,fc);
    BLp = B;
else
    try
        Wn = fc/(fsIn/2);
        [b,a] = butter(filterOrder,Wn,'low');
        BLp = filtfilt(b,a,B);
    catch
        % Toolbox-free fallback: zero-phase moving average.
        %
        % For fc = 0.05 Hz, one period is 20 s. A moving-average window
        % on this order strongly suppresses the observed ~7.5 s
        % disturbance while keeping only slow background variation.
        winLength = max(3,round(fsIn/fc));

        if mod(winLength,2) == 0
            winLength = winLength + 1;
        end

        kernel = ones(winLength,1)/winLength;

        BLp = conv(B,kernel,'same');
        BLp = flipud(BLp);
        BLp = conv(BLp,kernel,'same');
        BLp = flipud(BLp);
    end
end

% ---------------------------------------------------------------------
% 2) Resample to 1 Hz on a uniform datetime grid.
% ---------------------------------------------------------------------
fsOut = 1;                  % Hz
dtOut = seconds(1/fsOut);

tOut = (t(1):dtOut:t(end)).';
tsOut = seconds(tOut - t(1));

BOut = interp1(ts,BLp,tsOut,'linear');
end
