function [tCorrected,delayInfo] = estimate_quspin_time( ...
    tIn,angularRateQuSpin,angularRateRef)
%ESTIMATE_QUSPIN_TIME Estimate and correct the QuSpin time scale.
%
% Usage:
%   [tCorrected,delayInfo] = estimate_quspin_time( ...
%       tIn,angularRateQuSpin,angularRateRef)
%
% Inputs:
%   tIn
%       Nx1 assumed common time scale [s] for both angular-rate signals.
%
%   angularRateQuSpin
%       Nx3 QuSpin angular-rate data [rad/s].
%
%   angularRateRef
%       Nx3 reference angular-rate data [rad/s].
%
% Outputs:
%   tCorrected
%       Nx1 corrected time scale [s] for the QuSpin data.
%
%   delayInfo.time
%       Center time of each delay-estimation window [s].
%
%   delayInfo.delay_raw
%       Raw delay estimate in each window [s].
%
%   delayInfo.delay
%       Filtered delay estimate in each window [s].
%
%   delayInfo.correlation
%       Peak normalized cross-correlation in each window.
%
%   delayInfo.valid
%       Logical vector indicating accepted delay estimates.
%
%   delayInfo.delay_at_input_time
%       Estimated delay interpolated to every sample in tIn [s].
%
% Sign convention:
%   The delay is obtained from xcorr(reference,QuSpin), and the corrected
%   QuSpin time scale is calculated as
%
%       tCorrected = tIn + delay.
%
%   A negative delay therefore moves the QuSpin samples earlier in time.
%
% Notes:
%   - Both signals are assumed to be sampled at 100 Hz on the same assumed
%     input time scale.
%   - Each estimate is based on a window centered at delayInfo.time.
%   - All settings are defined inside this function.

    % ---------------------------------------------------------------------
    % Settings.
    % ---------------------------------------------------------------------
    fs = 100;                 % Sample rate [Hz]
    windowLength = 60;        % Centered window length [s]
    stepLength = 1;           % Time between delay estimates [s]
    maxDelay = 3;             % Maximum searched delay [s]
    lowpassCutoff = 1;        % Gyroscope-norm low-pass cutoff [Hz]
    filterOrder = 5;          % Butterworth filter order
    minVariance = 0.02;       % Minimum variance in reference window
    minCorrelation = 0.95;    % Minimum accepted peak correlation
    medianLength = 21;        % Median-filter length in estimates

    % ---------------------------------------------------------------------
    % Form and low-pass filter the angular-rate norms.
    % ---------------------------------------------------------------------
    uRef = vecnorm(angularRateRef,2,2);
    uQuSpin = vecnorm(angularRateQuSpin,2,2);

    [b,a] = butter(filterOrder,lowpassCutoff/(fs/2),'low');

    uRef = filtfilt(b,a,uRef);
    uQuSpin = filtfilt(b,a,uQuSpin);

    % ---------------------------------------------------------------------
    % Define centered delay-estimation windows.
    % ---------------------------------------------------------------------
    nWindowSamples = round(windowLength*fs);

    if mod(nWindowSamples,2) == 0
        nWindowSamples = nWindowSamples + 1;
    end

    halfWindow = (nWindowSamples - 1)/2;
    stepSamples = round(stepLength*fs);
    maxLagSamples = round(maxDelay*fs);

    centerIdx = (1 + halfWindow):stepSamples: ...
        (numel(tIn) - halfWindow);

    nEstimates = numel(centerIdx);

    delayRaw = NaN(nEstimates,1);
    peakCorrelation = zeros(nEstimates,1);
    valid = false(nEstimates,1);

    % ---------------------------------------------------------------------
    % Estimate the delay in each centered window.
    % ---------------------------------------------------------------------
    for ii = 1:nEstimates
        idxWindow = centerIdx(ii) + (-halfWindow:halfWindow);

        [delayRaw(ii),peakCorrelation(ii)] = estimate_delay( ...
            uRef(idxWindow),uQuSpin(idxWindow), ...
            fs,maxLagSamples,minVariance,minCorrelation);

        valid(ii) = ~isnan(delayRaw(ii));
    end

    delayTime = tIn(centerIdx);

    % ---------------------------------------------------------------------
    % Fill rejected estimates and remove isolated outliers.
    % ---------------------------------------------------------------------
    delay = fillmissing(delayRaw,'linear','EndValues','nearest');
    delay = medfilt1(delay,medianLength);

    % ---------------------------------------------------------------------
    % Interpolate the delay to every sample on the input time scale.
    % ---------------------------------------------------------------------
    delayAtInputTime = interp1( ...
        delayTime,delay,tIn,'linear','extrap');

    tCorrected = tIn + delayAtInputTime;

    % ---------------------------------------------------------------------
    % Return diagnostic information.
    % ---------------------------------------------------------------------
    delayInfo = struct();
    delayInfo.time = delayTime;
    delayInfo.delay_raw = delayRaw;
    delayInfo.delay = delay;
    delayInfo.correlation = peakCorrelation;
    delayInfo.valid = valid;
    delayInfo.delay_at_input_time = delayAtInputTime;
end


function [delay,peakCorrelation] = estimate_delay( ...
    uRef,uQuSpin,fs,maxLagSamples,minVariance,minCorrelation)
%ESTIMATE_DELAY Estimate the QuSpin delay in one centered window.

    x = detrend(uRef);
    y = detrend(uQuSpin);

    [c,lags] = xcorr(x,y,maxLagSamples,'coeff');
    [peakCorrelation,idxMax] = max(c);

    if var(x) > minVariance && peakCorrelation > minCorrelation
        delay = lags(idxMax)/fs;
    else
        delay = NaN;
    end
end
