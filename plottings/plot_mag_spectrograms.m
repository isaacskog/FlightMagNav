
function plot_mag_spectrograms(data,idx)
%PLOT_MAG_SPECTROGRAMS Plot spectrograms for reference and front magnetometers.
%
% Usage:
%   plot_mag_spectrograms(data,1:4)
%
% Processing:
%   1. Low-pass filter total field to approximately 4 Hz.
%   2. Resample to a uniform 10 Hz grid.
%   3. Remove mean.
%   4. Plot toolbox-free STFT spectrograms.
%
% The function does not use spectrogram, butter, filtfilt or pwelch.

    if nargin < 2 || isempty(idx)
        idx = 1:numel(data);
    end

    fsOut = 10;
    fc = 4;

    nPlot = numel(idx);

    figure;
    tl = tiledlayout(nPlot,2,'TileSpacing','compact','Padding','compact');

    for ii = 1:nPlot
        k = idx(ii);

        [tRef,BRef] = prepare_mag_signal_for_spectrogram( ...
            data(k).ref_mag.time, data(k).ref_mag.tot_field, fsOut, fc);

        nexttile;
        local_spectrogram_plot(tRef,BRef,fsOut);
        title(sprintf('Flight %d: ref\\_mag',k));

        [tFront,BFront] = prepare_mag_signal_for_spectrogram( ...
            data(k).front_mag.time, data(k).front_mag.tot_field, fsOut, fc);

        nexttile;
        local_spectrogram_plot(tFront,BFront,fsOut);
        title(sprintf('Flight %d: front\\_mag',k));
    end

    figure;
    plot(tRef(10:end-10),BRef(10:end-10))

    title(tl,'Magnetometer total-field spectrograms');
end


function [tOut,BOut] = prepare_mag_signal_for_spectrogram(t,B,fsOut,fc)

    t = t(:);
    B = B(:);

    valid = ~isnat(t) & ~isnan(B);
    t = t(valid);
    B = B(valid);

    if numel(t) < 10
        error('Too few valid samples for spectrogram.');
    end

    [t,idxSort] = sort(t);
    B = B(idxSort);

    ts = seconds(t - t(1));

    [ts,idxUnique] = unique(ts,'stable');
    t = t(idxUnique);
    B = B(idxUnique);

    fsIn = 1/median(diff(ts),'omitnan');

    validB = ~isnan(B);
    if ~all(validB)
        B = interp1(ts(validB),B(validB),ts,'linear','extrap');
    end

    Bf = lowpass_moving_average_zero_phase(B,fsIn,fc);

    dtOut = seconds(1/fsOut);
    tOut = (t(1):dtOut:t(end)).';
    tsOut = seconds(tOut - t(1));

    BOut = interp1(ts,Bf,tsOut,'linear');
    BOut = BOut - mean(BOut,'omitnan');
end


function y = lowpass_moving_average_zero_phase(x,fs,fc)

    x = x(:);

    if fs <= 2*fc
        warning('Input sample rate %.2f Hz is too low for %.1f Hz low-pass filtering. Filtering skipped.',fs,fc);
        y = x;
        return;
    end

    winLength = max(3,round(fs/(2*fc)));

    if mod(winLength,2) == 0
        winLength = winLength + 1;
    end

    h = ones(winLength,1)/winLength;

    y = conv(x,h,'same');
    y = flipud(y);
    y = conv(y,h,'same');
    y = flipud(y);
end


function local_spectrogram_plot(t,x,fs)

    x = x(:);

    winSec = 30;
    winLength = round(winSec*fs);
    winLength = min(winLength,numel(x));

    if winLength < 32
        error('Signal is too short for spectrogram.');
    end

    hopLength = max(1,round(winLength/4));
    nfft = 2^nextpow2(winLength);

    n = (0:winLength-1).';
    w = 0.5 - 0.5*cos(2*pi*n/(winLength-1));

    startIdx = 1:hopLength:(numel(x)-winLength+1);
    nFrames = numel(startIdx);

    S = zeros(nfft/2+1,nFrames);
    tFrame = NaT(1,nFrames,'TimeZone',t.TimeZone);

    for jj = 1:nFrames
        ind = startIdx(jj):(startIdx(jj)+winLength-1);

        xw = x(ind).*w;
        X = fft(xw,nfft);

        S(:,jj) = abs(X(1:nfft/2+1)).^2;
        tFrame(jj) = t(ind(round(end/2)));
    end

    f = (0:nfft/2)'*fs/nfft;
    SdB = 10*log10(S + eps);

    imagesc(tFrame,f,SdB);
    axis xy;
    ylim([0 5]);
    grid on;

    xlabel('Time (UTC)');
    ylabel('Frequency [Hz]');
    cb = colorbar;
    cb.Label.String = 'Power [dB]';

    hold on;
    yline(1,'--','1 Hz');
    yline(2,'--','2 Hz');
    yline(3,'--','3 Hz');
    yline(1/8.5,'--','8.5 s');
end
