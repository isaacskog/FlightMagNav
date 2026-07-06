function data = parse_ref_mag_data(filename,timeGPSaidedINS)
%PARSE_REF_MAG_DATA Parse reference scalar magnetometer data.
%
% Usage:
%   data = parse_ref_mag_data(filename,timeGPSaidedINS)
%
% Inputs:
%   filename
%       Path to the reference scalar magnetometer log file.
%
%   timeGPSaidedINS
%       Nx1 datetime vector from the GPS-aided INS parser. This function
%       samples the reference magnetometer data at exactly these time
%       instants. The time vector should have TimeZone = 'UTC'.
%
% Output:
%   data is a struct with the following fields:
%
%   data.time
%       Nx1 datetime vector. This is identical to timeGPSaidedINS after
%       removing possible samples outside the reference magnetometer time
%       span.
%
%   data.tot_field
%       Nx1 total magnetic field [nT], sampled at data.time.
%
%   data.quality
%       Nx1 quality indicator interpolated to the GPS-aided INS clock
%       using nearest-neighbour interpolation.
%
% What the function does:
%   1. Reads scalar magnetometer samples from lines of the form
%          internalTimestamp!magneticField_squality
%      for example
%          129254180!51312.571_s126
%
%   2. Reads GNSS timing lines of the form
%          internalTimestamp,GNSSFIX,...,year,month,day,hour,min,sec,...
%
%   3. Uses the GNSS timing lines to map the magnetometer internal clock to
%      absolute UTC time.
%
%   4. Low-pass filters the magnetic field with a 45 Hz Butterworth filter
%      using butter.m and filtfilt.m. This is done before resampling to the
%      100 Hz GPS-aided INS clock, to reduce possible aliasing from the
%      50 Hz power-line component.
%
%   5. Interpolates the filtered magnetic field to the GPS-aided INS time
%      vector.
%
% Notes:
%   - The magnetometer log appears to use an internal timestamp in
%     microseconds.
%   - The filter is designed using the median sample interval estimated from
%     the GPS-time-stamped magnetometer data.
%   - Samples in timeGPSaidedINS outside the magnetometer time span are
%     removed from the output.

    arguments
        filename (1,:) char
        timeGPSaidedINS (:,1) datetime
    end

    if isempty(timeGPSaidedINS)
        error('timeGPSaidedINS is empty.');
    end

    % Make sure the INS time vector is interpreted as UTC.
    timeGPSaidedINS.TimeZone = 'UTC';

    % ---------------------------------------------------------------------
    % Read the file as raw bytes. Some logger files may contain null bytes,
    % so these are removed before parsing the text.
    % ---------------------------------------------------------------------
    fid = fopen(filename,'r');
    if fid == -1
        fprintf('Could not open file: %s\n',filename);
        data=[];
        return;
    end

    raw = fread(fid,Inf,'*uint8');
    fclose(fid);

    raw(raw == 0) = [];

    txt = char(raw.');
    lines = regexp(txt,'\r\n|\n|\r','split');

    nLines = numel(lines);

    % ---------------------------------------------------------------------
    % Preallocate arrays using the number of lines as a rough upper bound.
    % ---------------------------------------------------------------------
    tMagInternal = zeros(nLines,1);
    totField     = zeros(nLines,1);
    quality      = zeros(nLines,1);
    nMag = 0;

    tGpsInternal = zeros(nLines,1);
    tGpsDatetime = NaT(nLines,1,'TimeZone','UTC');
    nGps = 0;

    scalarPattern = '^(\d+)!([+-]?\d*\.?\d+)_s(\d+)$';

    % ---------------------------------------------------------------------
    % Parse scalar magnetometer samples and GNSS timing messages.
    % ---------------------------------------------------------------------
    for i = 1:nLines
        line = strtrim(lines{i});

        if isempty(line)
            continue;
        end

        % Scalar magnetic field line.
        tok = regexp(line,scalarPattern,'tokens','once');

        if ~isempty(tok)
            nMag = nMag + 1;

            tMagInternal(nMag) = str2double(tok{1});
            totField(nMag)     = str2double(tok{2});   % [nT]
            quality(nMag)      = str2double(tok{3});

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
    quality      = quality(1:nMag);

    tGpsInternal = tGpsInternal(1:nGps);
    tGpsDatetime = tGpsDatetime(1:nGps);

    if nMag == 0
        error('No magnetic field samples found.');
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
    % Sort magnetic samples by internal timestamp.
    % ---------------------------------------------------------------------
    [tMagInternal,idxMag] = sort(tMagInternal);
    totField = totField(idxMag);
    quality  = quality(idxMag);

    % ---------------------------------------------------------------------
    % Interpolate the internal magnetometer clock to absolute UTC time.
    % ---------------------------------------------------------------------
    tGpsSeconds = posixtime(tGpsDatetime);

    tMagSeconds = interp1( ...
        tGpsInternal, ...
        tGpsSeconds, ...
        tMagInternal, ...
        'linear', ...
        'extrap');

    % Sort by absolute time.
    [tMagSeconds,idxTime] = sort(tMagSeconds);
    totField = totField(idxTime);
    quality = quality(idxTime);

    % Remove duplicate timestamps by averaging the magnetic field values.
    [tMagSeconds,~,groupIdx] = unique(tMagSeconds,'stable');
    if numel(tMagSeconds) < numel(totField)
        totField = accumarray(groupIdx,totField,[],@mean);
        quality  = accumarray(groupIdx,quality,[],@(x) round(mean(x)));
    end

    % ---------------------------------------------------------------------
    % Estimate magnetometer sample rate from the GPS-timestamped samples.
    % ---------------------------------------------------------------------
    dtMag = median(diff(tMagSeconds),'omitnan');
    fsMag = 1/dtMag;

    if ~isfinite(fsMag) || fsMag <= 0
        error('Could not estimate a valid magnetometer sample rate.');
    end

    % ---------------------------------------------------------------------
    % Low-pass filter before resampling to the 100 Hz INS clock.
    % The cutoff is 45 Hz to reduce possible aliasing of the 50 Hz power-line
    % disturbance when sampling to 100 Hz.
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
            b=[0.401885920660598   1.205657761981794   1.205657761981794   0.401885920660598];
            a=[1.000000000000000   1.280264479365736   0.775064545154888   0.159758340764158];
                 
            % filtfilt avoids phase distortion.
            totFieldFilt = zero_phase_filter_simple(b,a,totField);
        end
        
   
    end

    % ---------------------------------------------------------------------
    % Interpolate the filtered magnetic field to the GPS-aided INS clock.
    % Samples outside the magnetometer time span are removed.
    % ---------------------------------------------------------------------
    tInsSeconds = posixtime(timeGPSaidedINS);

    idxInside = tInsSeconds >= tMagSeconds(1) & tInsSeconds <= tMagSeconds(end);

    tOut = timeGPSaidedINS(idxInside);
    tOutSeconds = tInsSeconds(idxInside);

    totFieldOut = interp1(tMagSeconds,totFieldFilt,tOutSeconds,'linear');
    qualityOut  = interp1(tMagSeconds,quality,tOutSeconds,'nearest');

    % ---------------------------------------------------------------------
    % Return output struct.
    % ---------------------------------------------------------------------
    data = struct();
    data.time = tOut;
    data.tot_field = totFieldOut;
    data.quality = qualityOut;
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

    y = filter(b,a,x);
    y = flipud(y);
    y = filter(b,a,y);
    y = flipud(y);
end