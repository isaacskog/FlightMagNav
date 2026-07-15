function data = parse_front_mag_data(filename,GPSaidedINS)
%PARSE_FRONT_MAG_DATA Parse front magnetometer data and sample to INS clock.
%
% Usage:
%   data = parse_front_mag_data(filename,timeGPSaidedINS)
%
% Inputs:
%   filename
%       Path to the front magnetometer log file.
%
%   timeGPSaidedINS
%       Nx1 datetime vector from the GPS-aided INS parser. This function
%       samples the front magnetometer data at exactly these time instants.
%       The time vector should have TimeZone = 'UTC'.
%
% Output:
%   data is a struct with the following fields:
%
%   data.time
%       Nx1 datetime vector. This is the subset of timeGPSaidedINS that is
%       inside the magnetometer time span.
%
%   data.tot_field
%       Nx1 total magnetic field [nT], sampled at data.time.
%
%   data.vec_field
%       Nx3 vector magnetic field matrix [nT], ordered as [X, Y, Z].
%       The log reports vector components cyclically. Therefore, for each
%       output sample only the component that is closest in time is filled,
%       while the other two components are set to NaN.
%
%   data.angular_rate
%       Nx3 angular-rate matrix, ordered as [i, j, k], linearly
%       interpolated to data.time. Values are returned in the units used
%       in the log file.
%
%   data.acceleration
%       Nx3 acceleration matrix, ordered as [a, b, c], linearly
%       interpolated to data.time. Values are returned in the units used
%       in the log file.
%
%   data.quality
%       Nx1 quality indicator sampled to data.time using nearest-neighbour
%       interpolation.
%
% What the function does:
%   1. Reads front magnetometer data lines of the form
%          internalTimestamp!totField_XvecX=>...sQuality...
%          internalTimestamp!totField_YvecY=>...sQuality...
%          internalTimestamp!totField_ZvecZ=>...sQuality...
%
%      Example:
%          118743331!51318.312_X17130.334=>141949s054v076i0.43j-0.49k0.49
%
%      Here 51318.312 is the total magnetic field and X17130.334 is one
%      vector component. The vector component cycles between X, Y and Z.
%
%   2. Reads GNSS timing lines of the form
%          internalTimestamp,GNSSFIX,...,year,month,day,hour,min,sec,...
%
%   3. Uses the GNSS timing lines to map the magnetometer internal clock to
%      absolute UTC time.
%
%   4. Samples the total field and the quality indicator to the GPS-aided
%      INS time vector.
%
%   5. Samples the cyclic vector-component readings to the GPS-aided INS
%      time vector using nearest-neighbour association. Only the reported
%      component is filled at each output sample; the other components are
%      NaN. No low-pass filtering is applied to the vector components.
%
%   6. Parses the IMU angular-rate samples identified by i, j and k, and
%      acceleration samples identified by a, b and c. These are linearly
%      interpolated to the 100 Hz GPS-aided INS time vector.
%
% Notes:
%   - The total-field magnetometer data is low-pass filtered at 45 Hz before
%     interpolation to the 100 Hz INS clock. This reduces possible aliasing
%     from 50 Hz power-line disturbances.
%   - No low-pass filtering is applied to the vector components, since they
%     are reported cyclically at approximately 125/3 Hz.
%   - Samples in timeGPSaidedINS outside the magnetometer time span are
%     removed from the output.


    % Make sure the INS time vector is interpreted as UTC.
    timeGPSaidedINS=GPSaidedINS.time;
    timeGPSaidedINS.TimeZone = 'UTC';

    % ---------------------------------------------------------------------
    % Read file as raw bytes. Remove possible null bytes before parsing.
    % ---------------------------------------------------------------------
    fid = fopen(filename,'r');
    if fid == -1
        error('Could not open file: %s',filename);
    end

    raw = fread(fid,Inf,'*uint8');
    fclose(fid);

    raw(raw == 0) = [];

    txt = char(raw.');
    lines = regexp(txt,'\r\n|\n|\r','split');

    nLines = numel(lines);

    % ---------------------------------------------------------------------
    % Preallocate arrays.
    % ---------------------------------------------------------------------
    tMagInternal = zeros(nLines,1);
    totField     = zeros(nLines,1);
    vecFieldRaw  = NaN(nLines,3);
    quality      = zeros(nLines,1);
    nMag = 0;

    tGyroInternal = zeros(nLines,1);
    angularRateRaw = zeros(nLines,3);
    nGyro = 0;

    tAccInternal = zeros(nLines,1);
    accelerationRaw = zeros(nLines,3);
    nAcc = 0;

    tGpsInternal = zeros(nLines,1);
    tGpsDatetime = NaT(nLines,1,'TimeZone','UTC');
    nGps = 0;

    % Data line pattern:
    %   timestamp!total_componentValue=>...sQuality...
    %
    % Example:
    %   118743331!51318.312_X17130.334=>141949s054v076i...
    %
    % Tokens:
    %   1 timestamp
    %   2 total field
    %   3 component label X/Y/Z
    %   4 component value
    %   5 quality after 's'
    dataPattern = ['^(\d+)!([+-]?\d*\.?\d+)_([XYZ])', ...
                   '([+-]?\d*\.?\d+)=>.*?s(\d+)'];

    % IMU data appended to the magnetic-field line:
    %   angular rate: i<value>j<value>k<value>
    %   acceleration: a<value>b<value>c<value>
    gyroPattern = ['i([+-]?\d*\.?\d+)', ...
                   'j([+-]?\d*\.?\d+)', ...
                   'k([+-]?\d*\.?\d+)'];

    accPattern = ['a([+-]?\d*\.?\d+)', ...
                  'b([+-]?\d*\.?\d+)', ...
                  'c([+-]?\d*\.?\d+)'];

    % ---------------------------------------------------------------------
    % Parse magnetic field samples and GNSS timing messages.
    % ---------------------------------------------------------------------
    for i = 1:nLines
        line = strtrim(lines{i});

        if isempty(line)
            continue;
        end

        % Magnetic field data line.
        tok = regexp(line,dataPattern,'tokens','once');

        if ~isempty(tok)
            nMag = nMag + 1;

            tMagInternal(nMag) = str2double(tok{1});
            totField(nMag)     = str2double(tok{2});
            compLabel          = tok{3};
            compValue          = str2double(tok{4});
            quality(nMag)      = str2double(tok{5});

            switch compLabel
                case 'X'
                    vecFieldRaw(nMag,1) = compValue;
                case 'Y'
                    vecFieldRaw(nMag,2) = compValue;
                case 'Z'
                    vecFieldRaw(nMag,3) = compValue;
            end

            % Parse angular-rate data, if present on this line.
            gyroTok = regexp(line,gyroPattern,'tokens','once');
            if ~isempty(gyroTok)
                nGyro = nGyro + 1;
                tGyroInternal(nGyro) = tMagInternal(nMag);
                angularRateRaw(nGyro,:) = [ ...
                    str2double(gyroTok{1}), ...
                    str2double(gyroTok{2}), ...
                    str2double(gyroTok{3})];
            end

            % Parse acceleration data, if present on this line.
            accTok = regexp(line,accPattern,'tokens','once');
            if ~isempty(accTok)
                nAcc = nAcc + 1;
                tAccInternal(nAcc) = tMagInternal(nMag);
                accelerationRaw(nAcc,:) = [ ...
                    str2double(accTok{1}), ...
                    str2double(accTok{2}), ...
                    str2double(accTok{3})];
            end

            continue;
        end

        % GNSS timing line.
        if contains(line,'GNSSFIX')
            parts = strsplit(line,',');

            if numel(parts) >= 10
                nGps = nGps + 1;

                tGpsInternal(nGps) = str2double(parts{1});

                yearVal  = str2double(parts{5});
                monthVal = str2double(parts{6});
                dayVal   = str2double(parts{7});
                hourVal  = str2double(parts{8});
                minVal   = str2double(parts{9});
                secVal   = str2double(parts{10});

                tGpsDatetime(nGps) = datetime( ...
                    yearVal,monthVal,dayVal, ...
                    hourVal,minVal,secVal, ...
                    'TimeZone','UTC');
            end
        end
    end


    % Trim unused preallocated values.
    tMagInternal = tMagInternal(1:nMag);
    totField     = totField(1:nMag);
    vecFieldRaw  = vecFieldRaw(1:nMag,:);
    quality      = quality(1:nMag);

    tGyroInternal = tGyroInternal(1:nGyro);
    angularRateRaw = angularRateRaw(1:nGyro,:);

    tAccInternal = tAccInternal(1:nAcc);
    accelerationRaw = accelerationRaw(1:nAcc,:);

    tGpsInternal = tGpsInternal(1:nGps);
    tGpsDatetime = tGpsDatetime(1:nGps);

    if nMag == 0
        error('No front magnetometer samples found.');
    end

    if nGyro < 2
        error('Need at least two angular-rate samples.');
    end

    if nAcc < 2
        error('Need at least two acceleration samples.');
    end

    if nGps < 2
        error('Need at least two GNSSFIX records for GPS time interpolation.');
    end

    % ---------------------------------------------------------------------
    % Clean and sort GNSS timing data.
    % ---------------------------------------------------------------------
    validGps = ~isnan(tGpsInternal) & ~isnat(tGpsDatetime);

    tGpsInternal = tGpsInternal(validGps);
    tGpsDatetime = tGpsDatetime(validGps);

    [tGpsInternal,idxGps] = sort(tGpsInternal);
    tGpsDatetime = tGpsDatetime(idxGps);

    [tGpsInternal,idxUniqueGps] = unique(tGpsInternal,'stable');
    tGpsDatetime = tGpsDatetime(idxUniqueGps);

    if numel(tGpsInternal) < 2
        error('Need at least two unique GNSSFIX timestamps.');
    end

    % ---------------------------------------------------------------------
    % Sort magnetometer samples by internal timestamp.
    % ---------------------------------------------------------------------
    [tMagInternal,idxMag] = sort(tMagInternal);
    totField    = totField(idxMag);
    vecFieldRaw = vecFieldRaw(idxMag,:);
    quality     = quality(idxMag);

    % ---------------------------------------------------------------------
    % Map internal magnetometer timestamps to absolute UTC time.
    % ---------------------------------------------------------------------
    tGpsSeconds = posixtime(tGpsDatetime);

    tMagSeconds = interp1( ...
        tGpsInternal, ...
        tGpsSeconds, ...
        tMagInternal, ...
        'linear', ...
        'extrap');


    tGyroSeconds = interp1( ...
        tGpsInternal, ...
        tGpsSeconds, ...
        tGyroInternal, ...
        'linear', ...
        'extrap');

    tAccSeconds = interp1( ...
        tGpsInternal, ...
        tGpsSeconds, ...
        tAccInternal, ...
        'linear', ...
        'extrap');

    % Sort by absolute time.
    [tMagSeconds,idxTime] = sort(tMagSeconds);
    totField    = totField(idxTime);
    vecFieldRaw = vecFieldRaw(idxTime,:);
    quality     = quality(idxTime);

    [tGyroSeconds,idxGyro] = sort(tGyroSeconds);
    angularRateRaw = angularRateRaw(idxGyro,:);

    [tAccSeconds,idxAcc] = sort(tAccSeconds);
    accelerationRaw = accelerationRaw(idxAcc,:);

    % Remove possible duplicate IMU timestamps by averaging each axis.
    [tGyroSeconds,angularRateRaw] = merge_duplicate_vector_samples( ...
        tGyroSeconds,angularRateRaw);

    [tAccSeconds,accelerationRaw] = merge_duplicate_vector_samples( ...
        tAccSeconds,accelerationRaw);


    % ---------------------------------------------------------------------
    % Remove duplicate timestamps.
    % Total field and quality are merged. For vector components, each
    % component is merged separately, ignoring NaN.
    % ---------------------------------------------------------------------
    [tMagSecondsUnique,~,groupIdx] = unique(tMagSeconds,'stable');

    if numel(tMagSecondsUnique) < numel(tMagSeconds)
        totFieldMerged = accumarray(groupIdx,totField,[],@mean);
        qualityMerged  = accumarray(groupIdx,quality,[],@(x) round(mean(x)));

        vecFieldMerged = NaN(numel(tMagSecondsUnique),3);
        for jj = 1:3
            for kk = 1:numel(tMagSecondsUnique)
                values = vecFieldRaw(groupIdx == kk,jj);
                values = values(~isnan(values));
                if ~isempty(values)
                    vecFieldMerged(kk,jj) = mean(values);
                end
            end
        end

        tMagSeconds = tMagSecondsUnique;
        totField    = totFieldMerged;
        vecFieldRaw = vecFieldMerged;
        quality     = qualityMerged;
    else
        tMagSeconds = tMagSecondsUnique;
    end

    % ---------------------------------------------------------------------
    % Estimate total-field magnetometer sample rate from the GPS-timestamped
    % scalar samples.
    % ---------------------------------------------------------------------
    dtMag = median(diff(tMagSeconds),'omitnan');
    fsMag = 1/dtMag;

    if ~isfinite(fsMag) || fsMag <= 0
        error('Could not estimate a valid magnetometer sample rate.');
    end

    % ---------------------------------------------------------------------
    % Low-pass filter the total field before resampling to the 100 Hz INS
    % clock. The cutoff is 45 Hz to reduce possible aliasing of the 50 Hz
    % power-line disturbance.
    %
    % First try butter/filtfilt from Signal Processing Toolbox. If that is
    % not available, use fixed coefficients for a 3rd-order Butterworth
    % low-pass filter designed for approximately fs = 125 Hz and fc = 45 Hz,
    % together with a simple forward-backward filter.
    % ---------------------------------------------------------------------
    fc = 45;       % Low-pass cutoff frequency [Hz]
    filterOrder = 3;

    if fsMag <= 2*fc
        warning(['Estimated magnetometer sample rate is %.2f Hz, which is ', ...
                 'too low for a %.1f Hz low-pass cutoff. Filtering skipped.'], ...
                 fsMag,fc);
        totFieldFilt = totField;
    else
        try
            Wn = fc/(fsMag/2);
            [b,a] = butter(filterOrder,Wn,'low');
            totFieldFilt = filtfilt(b,a,totField);
        catch
            b = [0.401885920660598, 1.205657761981794, 1.205657761981794, 0.401885920660598];
            a = [1.000000000000000, 1.280264479365736, 0.775064545154888, 0.159758340764158];

            totFieldFilt = zero_phase_filter_simple(b,a,totField);
        end
    end

    
    % ---------------------------------------------------------------------
    % Sample to the GPS-aided INS clock.
    % ---------------------------------------------------------------------
    tInsSeconds = posixtime(timeGPSaidedINS);

    tStart = max([tMagSeconds(1),tGyroSeconds(1),tAccSeconds(1)]);
    tEnd = min([tMagSeconds(end),tGyroSeconds(end),tAccSeconds(end)]);

    idxInside = tInsSeconds >= tStart & tInsSeconds <= tEnd;

    tOut = timeGPSaidedINS(idxInside);
    tOutSeconds = tInsSeconds(idxInside);

    % Total field is interpolated linearly to the INS clock after
    % low-pass filtering.
    totFieldOut = interp1(tMagSeconds,totFieldFilt,tOutSeconds,'linear');

    % Quality is a discrete indicator and is sampled using nearest-neighbour.
    qualityOut = interp1(tMagSeconds,quality,tOutSeconds,'nearest');

    % Vector components are reported cyclically. Associate each INS time
    % sample with the nearest magnetometer vector reading, and keep only the
    % reported component. The other components remain NaN.
    idxNearest = nearest_sample_indices(tMagSeconds,tOutSeconds);
    vecFieldOut = vecFieldRaw(idxNearest,:);

    % Angular rate and acceleration are linearly interpolated to the same
    % 100 Hz time vector as the other outputs.
    angularRateOut = interp1( ...
        tGyroSeconds,deg2rad(angularRateRaw),tOutSeconds,'linear');

    accelerationOut = interp1( ...
        tAccSeconds,accelerationRaw,tOutSeconds,'linear'); % Unclear unit


    %---------------------------------------------------------------------%
    % Find error/offset in clock of QuSpin sensor and re-interpolate
    %---------------------------------------------------------------------%
    gyrRef = GPSaidedINS.gyr(idxInside,:);
    [tCorrected,delayInfo] = estimate_quspin_time(tOutSeconds,...
        angularRateOut,gyrRef); 

    % Total field
    totFieldOut = interp1( ...
        tCorrected,totFieldOut,tOutSeconds,'linear','extrap');

    qualityOut = interp1( ...
        tCorrected,qualityOut,tOutSeconds,'nearest','extrap');

    idxNearest = nearest_sample_indices(tCorrected,tOutSeconds);
    vecFieldOut = vecFieldOut(idxNearest,:);

    angularRateOut = interp1( ...
        tCorrected,angularRateOut,tOutSeconds,'linear','extrap');

    accelerationOut = interp1( ...
        tCorrected,accelerationOut,tOutSeconds,'linear','extrap');


    % ---------------------------------------------------------------------
    % Return output struct.
    % ---------------------------------------------------------------------
    data = struct();
    data.time = tOut;
    data.tot_field = totFieldOut;
    data.vec_field = vecFieldOut;
    data.angular_rate = angularRateOut;
    data.acceleration = accelerationOut;
    data.quality = qualityOut;
    data.delay_info = delayInfo;
end



function y = zero_phase_filter_simple(b,a,x)
%ZERO_PHASE_FILTER_SIMPLE Simple forward-backward filtering without filtfilt.
%
% This approximates filtfilt by filtering forward, reversing the signal,
% filtering again, and reversing back.
%
% Note:
%   This does not use the same initial-condition handling as filtfilt, so
%   edge transients may be larger near the beginning and end of the signal.

    x = x(:);

    y = filter(b,a,x);
    y = flipud(y);
    y = filter(b,a,y);
    y = flipud(y);
end


function [tUnique,xUnique] = merge_duplicate_vector_samples(t,x)
%MERGE_DUPLICATE_VECTOR_SAMPLES Average vector samples with equal times.

    [tUnique,~,groupIdx] = unique(t,'stable');

    if numel(tUnique) == numel(t)
        xUnique = x;
        return;
    end

    xUnique = zeros(numel(tUnique),size(x,2));

    for jj = 1:size(x,2)
        xUnique(:,jj) = accumarray(groupIdx,x(:,jj),[],@mean);
    end
end


function idxNearest = nearest_sample_indices(tSample,tQuery)
%NEAREST_SAMPLE_INDICES Return index of nearest sample for each query time.
%
% Inputs:
%   tSample - Mx1 sorted sample times [seconds]
%   tQuery  - Nx1 query times [seconds]
%
% Output:
%   idxNearest - Nx1 indices into tSample.

    idxNearest = interp1( ...
        tSample, ...
        (1:numel(tSample)).', ...
        tQuery, ...
        'nearest', ...
        'extrap');

    idxNearest = round(idxNearest);
    idxNearest = max(1,min(numel(tSample),idxNearest));
end
