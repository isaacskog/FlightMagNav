function data = parse_gps_aided_ins_data(filename,startTime,stopTime)
%PARSE_GPS_AIDED_INS_DATA Parse GPS-aided INS data from Xsens/MTi text export.
%
% Usage:
%   data = parse_gps_aided_ins_data(filename,startTime,stopTime)
%
% Output:
%   data.time is a UTC datetime vector uniformly sampled at 100 Hz. The time
%   stamps are forced to lie on exact hundredths of a second.
%
% Notes:
%   - The original Xsens/MTi file is nominally 100 Hz, but can contain small
%     timestamp jitter.
%   - All signals are resampled to a clean 100 Hz clock using nearest
%     neighbour. This is deliberate: for quaternions, nearest neighbour is
%     safer than naive component-wise linear interpolation.
%   - Velocity is converted from ENU in the file to NED in the output.

    arguments
        filename (1,:) char
        startTime (1,1) datetime
        stopTime (1,1) datetime
    end

    if stopTime <= startTime
        error('stopTime must be later than startTime.');
    end

    startTime.TimeZone = 'UTC';
    stopTime.TimeZone  = 'UTC';

    % Find header line.
    fid = fopen(filename,'r');
    if fid == -1
        error('Could not open file: %s',filename);
    end

    headerLineNumber = [];
    lineNumber = 0;

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end

        lineNumber = lineNumber + 1;

        if startsWith(strtrim(line),'PacketCounter')
            headerLineNumber = lineNumber;
            break;
        end
    end

    fclose(fid);

    if isempty(headerLineNumber)
        error('Could not find CSV header line starting with PacketCounter.');
    end

    opts = detectImportOptions(filename, ...
        'FileType','text', ...
        'Delimiter',',', ...
        'NumHeaderLines',headerLineNumber-1);

    opts.VariableNamingRule = 'preserve';
    T = readtable(filename,opts);
    T.Properties.VariableNames = matlab.lang.makeValidName(T.Properties.VariableNames);

    requiredVars = { ...
        'UTC_Year','UTC_Month','UTC_Day', ...
        'UTC_Hour','UTC_Minute','UTC_Second','UTC_Nano', ...
        'Acc_X','Acc_Y','Acc_Z', ...
        'Gyr_X','Gyr_Y','Gyr_Z', ...
        'Quat_q0','Quat_q1','Quat_q2','Quat_q3', ...
        'Roll','Pitch','Yaw', ...
        'Latitude','Longitude','Altitude', ...
        'Vel_E','Vel_N','Vel_U'};

    for i = 1:numel(requiredVars)
        if ~ismember(requiredVars{i},T.Properties.VariableNames)
            error('Required variable "%s" was not found in the file.',requiredVars{i});
        end
    end

    % Build UTC datetime vector.
    t = datetime( ...
        T.UTC_Year, ...
        T.UTC_Month, ...
        T.UTC_Day, ...
        T.UTC_Hour, ...
        T.UTC_Minute, ...
        T.UTC_Second, ...
        'TimeZone','UTC');

    t = t + seconds(T.UTC_Nano*1e-9);

    % Crop.
    validTime = ~isnat(t);
    idx = validTime & t >= startTime & t <= stopTime;

    T = T(idx,:);
    t = t(idx);

    % Sort and remove duplicate timestamps.
    [t,sortIdx] = sort(t);
    T = T(sortIdx,:);

    [t,uniqueIdx] = unique(t,'stable');
    T = T(uniqueIdx,:);

    % Clean 100 Hz clock on exact hundredths of a second.
    fsOut = 100;
    dtOut = seconds(1/fsOut);

    tStart = round_datetime_to_sample_grid(t(1),fsOut,'ceil');
    tStop  = round_datetime_to_sample_grid(t(end),fsOut,'floor');

    tOut = (tStart:dtOut:tStop).';

    % Use seconds relative to first original timestamp.
    tSec = seconds(t - t(1));
    tOutSec = seconds(tOut - t(1));

    % Nearest-neighbour resampling.
    latitude  = interp1(tSec,T.Latitude, tOutSec,'nearest','extrap');
    longitude = interp1(tSec,T.Longitude,tOutSec,'nearest','extrap');
    altitude  = interp1(tSec,T.Altitude, tOutSec,'nearest','extrap');

    velE = interp1(tSec,T.Vel_E,tOutSec,'nearest','extrap');
    velN = interp1(tSec,T.Vel_N,tOutSec,'nearest','extrap');
    velU = interp1(tSec,T.Vel_U,tOutSec,'nearest','extrap');

    quat = zeros(numel(tOut),4);
    quat(:,1) = interp1(tSec,T.Quat_q0,tOutSec,'nearest','extrap');
    quat(:,2) = interp1(tSec,T.Quat_q1,tOutSec,'nearest','extrap');
    quat(:,3) = interp1(tSec,T.Quat_q2,tOutSec,'nearest','extrap');
    quat(:,4) = interp1(tSec,T.Quat_q3,tOutSec,'nearest','extrap');

    euler = zeros(numel(tOut),3);
    euler(:,1) = interp1(tSec,T.Roll, tOutSec,'nearest','extrap');
    euler(:,2) = interp1(tSec,T.Pitch,tOutSec,'nearest','extrap');
    euler(:,3) = interp1(tSec,T.Yaw,  tOutSec,'nearest','extrap');

    acc = zeros(numel(tOut),3);
    acc(:,1) = interp1(tSec,T.Acc_X,tOutSec,'nearest','extrap');
    acc(:,2) = interp1(tSec,T.Acc_Y,tOutSec,'nearest','extrap');
    acc(:,3) = interp1(tSec,T.Acc_Z,tOutSec,'nearest','extrap');

    gyr = zeros(numel(tOut),3);
    gyr(:,1) = interp1(tSec,T.Gyr_X,tOutSec,'nearest','extrap');
    gyr(:,2) = interp1(tSec,T.Gyr_Y,tOutSec,'nearest','extrap');
    gyr(:,3) = interp1(tSec,T.Gyr_Z,tOutSec,'nearest','extrap');

    % Normalize quaternion rows after resampling.
    quatNorm = sqrt(sum(quat.^2,2));
    quat = quat./quatNorm;

    data = struct();
    data.time = tOut;
    data.latitude = latitude;
    data.longitude = longitude;
    data.altitude = altitude;
    data.velocity_ned = [velN, velE, -velU];
    data.quaternion = quat;
    data.euler = euler;
    data.acc = acc;
    data.gyr = gyr;
end


function tRound = round_datetime_to_sample_grid(t,fs,mode)
%ROUND_DATETIME_TO_SAMPLE_GRID Round datetime to multiples of 1/fs seconds.
%
% For fs = 100 this gives exact hundredths of a second. Rounding is done
% relative to the start of the UTC day.

    dayStart = dateshift(t,'start','day');
    sec = seconds(t - dayStart);

    sample = sec*fs;

    switch mode
        case 'ceil'
            sample = ceil(sample);
        case 'floor'
            sample = floor(sample);
        case 'round'
            sample = round(sample);
        otherwise
            error('Unknown rounding mode.');
    end

    tRound = dayStart + seconds(sample/fs);
end
