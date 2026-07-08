function data = parse_ups_ref_mag_data(filename,timeGPSaidedINS)
%PARSE_UPS_REF_MAG_DATA Parse Uppsala observatory IAGA-2002 .sec data.
%
% Usage:
%   data = parse_ups_ref_mag_data(filename,timeGPSaidedINS)
%
% Inputs:
%   filename
%       Path to Uppsala observatory .sec file in IAGA-2002 format.
%
%   timeGPSaidedINS
%       Nx1 datetime vector from the GPS-aided INS parser.
%
% Output:
%   data.time
%       Datetime vector, restricted to the overlap with the observatory data.
%
%   data.tot_field
%       Total magnetic field UPSF [nT], interpolated to data.time.
%
%   data.X, data.Y, data.Z
%       Vector magnetic field components [nT], interpolated to data.time.

    arguments
        filename (1,:) char
        timeGPSaidedINS (:,1) datetime
    end

    if isempty(timeGPSaidedINS)
        error('timeGPSaidedINS is empty.');
    end

    timeGPSaidedINS.TimeZone = 'UTC';

    fid = fopen(filename,'r');
    if fid == -1
        error('Could not open file: %s',filename);
    end

    raw = fread(fid,Inf,'*char')';
    fclose(fid);

    lines = regexp(raw,'\r\n|\n|\r','split');

    isDataLine = false(numel(lines),1);
    for ii = 1:numel(lines)
        isDataLine(ii) = ~isempty(regexp(lines{ii}, ...
            '^\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2}\.\d{3}', ...
            'once'));
    end

    dataLines = lines(isDataLine);

    if isempty(dataLines)
        error('No IAGA-2002 data lines found in file.');
    end

    txt = strjoin(dataLines,newline);

    C = textscan(txt,'%s %s %f %f %f %f %f', ...
        'Delimiter',' ', ...
        'MultipleDelimsAsOne',true);

    dateStr = C{1};
    timeStr = C{2};

    X = C{4};
    Y = C{5};
    Z = C{6};
    F = C{7};

    tObs = datetime(strcat(dateStr,{' '},timeStr), ...
        'InputFormat','yyyy-MM-dd HH:mm:ss.SSS', ...
        'TimeZone','UTC');

    % Remove missing/sentinel values.
    valid = ~isnat(tObs) & ...
            isfinite(X) & isfinite(Y) & isfinite(Z) & isfinite(F) & ...
            abs(X) < 90000 & abs(Y) < 90000 & ...
            abs(Z) < 90000 & abs(F) < 90000;

    tObs = tObs(valid);
    X = X(valid);
    Y = Y(valid);
    Z = Z(valid);
    F = F(valid);

    if numel(tObs) < 2
        error('Need at least two valid observatory samples.');
    end

    % Sort and remove duplicate timestamps.
    [tObs,idxSort] = sort(tObs);
    X = X(idxSort);
    Y = Y(idxSort);
    Z = Z(idxSort);
    F = F(idxSort);

    tObsSec = posixtime(tObs);

    [tObsSec,~,groupIdx] = unique(tObsSec,'stable');
    if numel(tObsSec) < numel(F)
        X = accumarray(groupIdx,X,[],@mean);
        Y = accumarray(groupIdx,Y,[],@mean);
        Z = accumarray(groupIdx,Z,[],@mean);
        F = accumarray(groupIdx,F,[],@mean);
    end

    tInsSec = posixtime(timeGPSaidedINS);

    idxInside = tInsSec >= tObsSec(1) & tInsSec <= tObsSec(end);

    tOut = timeGPSaidedINS(idxInside);
    tOutSec = tInsSec(idxInside);

    data = struct();
    data.time = tOut;
    data.tot_field = interp1(tObsSec,F,tOutSec,'linear');
    data.X = interp1(tObsSec,X,tOutSec,'linear');
    data.Y = interp1(tObsSec,Y,tOutSec,'linear');
    data.Z = interp1(tObsSec,Z,tOutSec,'linear');
    data.source = 'UPS';
end