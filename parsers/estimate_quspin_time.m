function [tCorrected,delayInfo] = estimate_quspin_time( ...
    tIn,angularRateQuSpin,angularRateRef)
%ESTIMATE_QUSPIN_TIME Estimate an affine QuSpin clock correction.
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
%       Raw delay estimate in each accepted or rejected window [s].
%
%   delayInfo.delay
%       Delay predicted by the fitted affine clock model at each window
%       center [s].
%
%   delayInfo.correlation
%       Peak normalized cross-correlation in each window.
%
%   delayInfo.valid
%       Logical vector indicating which raw delay estimates were accepted.
%
%   delayInfo.delay_at_input_time
%       Affine delay correction evaluated at every sample in tIn [s].
%
%   delayInfo.alpha
%       Estimated QuSpin clock scale.
%
%   delayInfo.beta
%       Estimated QuSpin clock offset at delayInfo.t0 [s].
%
%   delayInfo.t0
%       Time origin used for the centered affine clock model [s].
%
%   delayInfo.fit_residual
%       Difference between each raw delay estimate and the fitted affine
%       delay model [s].
%
% Clock model:
%
%   tCorrected = t0 + alpha*(tIn - t0) + beta
%
% Equivalently, the estimated delay is
%
%   delay(t) = beta + (alpha - 1)*(t - t0).
%
% Sign convention:
%   The local delay is obtained from xcorr(reference,QuSpin). A negative
%   delay therefore moves the QuSpin timestamps earlier in time.
%
% Notes:
%   - Both signals are assumed to be sampled at 100 Hz on the same assumed
%     input time scale.
%   - Delay estimates use non-overlapping 60-second centered windows.
%   - A robust bisquare least-squares fit is used to estimate alpha and
%     beta from the accepted local delay estimates.
%   - All settings are defined inside this function.

    % ---------------------------------------------------------------------
    % Settings.
    % ---------------------------------------------------------------------
    fs = 100;                 % Sample rate [Hz]
    windowLength = 60;        % Centered window length [s]
    stepLength = 60;          % Time between delay estimates [s]
    maxDelay = 3;             % Maximum searched delay [s]
    lowpassCutoff = 1;        % Gyroscope-norm low-pass cutoff [Hz]
    filterOrder = 5;          % Butterworth filter order
    minVariance = 0.02;       % Minimum variance in reference window
    minCorrelation = 0.95;    % Minimum accepted peak correlation

    % ---------------------------------------------------------------------
    % Form and low-pass filter the angular-rate norms.
    % ---------------------------------------------------------------------
    uRef = vecnorm(angularRateRef,2,2);
    uQuSpin = vecnorm(angularRateQuSpin,2,2);

    [b,a] = butter(filterOrder,lowpassCutoff/(fs/2),'low');

    uRef = filtfilt(b,a,uRef);
    uQuSpin = filtfilt(b,a,uQuSpin);

    % ---------------------------------------------------------------------
    % Define non-overlapping centered delay-estimation windows.
    % ---------------------------------------------------------------------
    nWindowSamples = round(windowLength*fs);

    % Use an odd number of samples to obtain one unique center sample.
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

    if sum(valid) < 2
        error('Need at least two valid delay estimates to fit the clock model.');
    end

    % ---------------------------------------------------------------------
    % Estimate the affine clock model using robust least squares.
    %
    % Centering the time avoids numerical problems when tIn contains large
    % absolute values, such as POSIX time.
    % ---------------------------------------------------------------------
    t0 = delayTime(1);
    tau = delayTime - t0;

    [clockFit,fitStats] = robustfit( ...
        tau(valid),delayRaw(valid),'cauchy');

    beta = clockFit(1);
    clockDrift = clockFit(2);
    alpha = 1 + clockDrift;

    % Delay predicted by the affine model at the estimation times.
    delay = beta + clockDrift*tau;

    % Delay predicted by the affine model at every input timestamp.
    delayAtInputTime = beta + clockDrift*(tIn - t0);

    % Apply the affine clock correction.
    tCorrected = tIn + delayAtInputTime;

    % Residuals are only defined for accepted raw delay estimates.
    fitResidual = NaN(size(delayRaw));
    fitResidual(valid) = delayRaw(valid) - delay(valid);

    robustWeight = NaN(size(delayRaw));
    robustWeight(valid) = fitStats.w;

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
    delayInfo.alpha = alpha;
    delayInfo.beta = beta;
    delayInfo.t0 = t0;
    delayInfo.fit_residual = fitResidual;
    delayInfo.robust_weight = robustWeight;
end


function [delay,peakCorrelation] = estimate_delay( ...
    uRef,uQuSpin,fs,maxLagSamples,minVariance,minCorrelation)
%ESTIMATE_DELAY Estimate the QuSpin delay in one centered window.
%
% The estimate is accepted only when the reference signal contains enough
% variation and the normalized cross-correlation peak is sufficiently high.

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
