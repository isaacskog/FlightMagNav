% Main script for parsing all data from the data collection in Skovde.
% All data are synchronized to the GPS-aided INS and converted to a 100 Hz
% sample rate.
%
% Notes:
%
% 1) Scalar magnetometer data are low-pass filtered at 45 Hz to avoid
% aliasing of the 50 Hz components when downsampling to 100 Hz.
%
% 2) Vector readings from the magnetometers are kept at 125/3 Hz.
%
% Isaac Skog and ChatGPT, 2026
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preamble

% Dataset folders.
datasetRoot = fullfile('..','..');

rawRoot = fullfile(datasetRoot,'raw');
metadataFile = fullfile(datasetRoot,'metadata','flights.txt');
parsedRoot = fullfile(datasetRoot,'parsed');

% Read flight metadata.
flightInfo = read_flight_metadata(metadataFile);

nFlights = numel(flightInfo);

%% Parse all data

template = struct( ...
    'flight_id',[], ...
    'start_time',[], ...
    'stop_time',[], ...
    'raw_folder','', ...
    'purpose','', ...
    'ref_mag',[], ...
    'front_mag',[], ...
    'back_mag',[], ...
    'GPSaidedINS',[], ...
    'UPS_ref_mag',[]);

data = repmat(template,1,nFlights);

parfor ii = 1:nFlights

    disp(['Flight ' num2str(flightInfo(ii).flight) ...
        ' out of ' num2str(nFlights)]);

    flightData = template;

    % ---------------------------------------------------------------------
    % Flight metadata.
    % ---------------------------------------------------------------------
    flightData.flight_id = flightInfo(ii).flight;
    flightData.start_time = flightInfo(ii).start_time;
    flightData.stop_time = flightInfo(ii).stop_time;
    flightData.raw_folder = flightInfo(ii).raw_folder;
    flightData.purpose = flightInfo(ii).purpose;

    % ---------------------------------------------------------------------
    % Determine raw-data folder.
    %
    % Example:
    %   raw/2026-06-10/acquisition_01/
    % ---------------------------------------------------------------------
    dateFolder = datestr(flightInfo(ii).start_time,'yyyy-mm-dd');

    acquisitionFolder = fullfile( ...
        rawRoot, ...
        dateFolder, ...
        flightInfo(ii).raw_folder);

    % ---------------------------------------------------------------------
    % GPS-aided INS.
    %
    % The time base of this system is used as the common time base for all
    % other sensors.
    % ---------------------------------------------------------------------
    filename = fullfile(acquisitionFolder,'GPSINS.txt');

    flightData.GPSaidedINS = parse_gps_aided_ins_data( ...
        filename, ...
        flightInfo(ii).start_time, ...
        flightInfo(ii).stop_time);

    % ---------------------------------------------------------------------
    % Local reference magnetometer.
    % ---------------------------------------------------------------------
    if flightInfo(ii).local_ref_ok

        filename = fullfile(acquisitionFolder,'RefMag.txt');

        flightData.ref_mag = parse_ref_mag_data( ...
            filename, ...
            flightData.GPSaidedINS.time);
    end

    % ---------------------------------------------------------------------
    % Uppsala geomagnetic reference station.
    %
    % Example:
    %   raw/reference/ups_ref_20260610.sec
    % ---------------------------------------------------------------------
    dateString = datestr(flightInfo(ii).start_time,'yyyymmdd');

    filename = fullfile( ...
        rawRoot, ...
        'reference', ...
        ['ups_ref_' dateString '.sec']);

    flightData.UPS_ref_mag = parse_ups_ref_mag_data( ...
        filename, ...
        flightData.GPSaidedINS.time);

    % ---------------------------------------------------------------------
    % Front magnetometer.
    % ---------------------------------------------------------------------
    filename = fullfile(acquisitionFolder,'FrontMag.txt');

    flightData.front_mag = parse_front_mag_data( ...
        filename, ...
        flightData.GPSaidedINS);

    % ---------------------------------------------------------------------
    % Back magnetometer.
    % ---------------------------------------------------------------------
    filename = fullfile(acquisitionFolder,'BackMag.txt');

    flightData.back_mag = parse_front_mag_data( ...
        filename, ...
        flightData.GPSaidedINS);

    % Store parsed flight.
    data(ii) = flightData;
end

%% Save parsed data

if ~exist(parsedRoot,'dir')
    mkdir(parsedRoot);
end

filename = fullfile(parsedRoot,'FlightData.mat');

save(filename,'data','-v7.3');


%% Local functions

function flightInfo = read_flight_metadata(filename)
%READ_FLIGHT_METADATA Read metadata describing the individual flights.
%
% Expected format:
%
% flight raw_folder start_utc stop_utc local_ref_ok purpose

    lines = readlines(filename);

    % Remove header and empty lines.
    lines = lines(2:end);
    lines = lines(strlength(strtrim(lines)) > 0);

    nFlights = numel(lines);

    flightInfo = repmat(struct( ...
        'flight',[], ...
        'raw_folder','', ...
        'start_time',[], ...
        'stop_time',[], ...
        'local_ref_ok',false, ...
        'purpose',''), ...
        1,nFlights);

    for ii = 1:nFlights

        parts = regexp(strtrim(lines(ii)),'\s+','split');

        flightInfo(ii).flight = str2double(parts{1});
        flightInfo(ii).raw_folder = parts{2};

        flightInfo(ii).start_time = datetime( ...
            [parts{3} ' ' parts{4}], ...
            'InputFormat','yyyy-MM-dd HH:mm:ss', ...
            'TimeZone','UTC');

        flightInfo(ii).stop_time = datetime( ...
            [parts{5} ' ' parts{6}], ...
            'InputFormat','yyyy-MM-dd HH:mm:ss', ...
            'TimeZone','UTC');

        flightInfo(ii).local_ref_ok = strcmpi(parts{7},'true');
        flightInfo(ii).purpose = parts{8};
    end
end

