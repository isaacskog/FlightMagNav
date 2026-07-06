% Main script for parsing all the data from the data collection in Skovde. 
% All data is synchronized to the GPS aided INS and converted into a 100 Hz
% sample rate.
%
% Note: 
% 
% 1) Scalar magnetometer data is low-pass filtred at 45 Hz to avoid
% aliasing of the 50 Hz components when downsampling to 100 Hz. 
%
% 2) Vector readings from the magnetometers are keeped at 125/3 Hz. 
%
% Isaac Skog and ChatGPT, 2026-06-15
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Preamble 

% Folders
folders={'Flygning1','Flygning2','Flygning3','Flygning4','Flygning5','Flygning6','Flygning6','Flygning7'}; % Map 6 containes data from 2 flights

% Start and stop times according to the protocols
startandstoptimes=...
    [datetime(2026,6,10,10,47,0,0,"TimeZone","UTC") datetime(2026,6,10,11,2,0,0,"TimeZone","UTC");... % 1
    datetime(2026,6,10,13,5,0,0,"TimeZone","UTC") datetime(2026,6,10,13,23,0,0,"TimeZone","UTC");... % 2
    datetime(2026,6,10,14,4,0,0,"TimeZone","UTC") datetime(2026,6,10,14,22,0,0,"TimeZone","UTC"); ... % 3
    datetime(2026,6,11,8,5,0,0,"TimeZone","UTC") datetime(2026,6,11,8,22,0,0,"TimeZone","UTC"); ... % 4
    datetime(2026,6,11,9,40,0,0,"TimeZone","UTC") datetime(2026,6,11,9,59,0,0,"TimeZone","UTC");... % 5
    datetime(2026,6,11,12,15,0,0,"TimeZone","UTC") datetime(2026,6,11,12,25,0,0,"TimeZone","UTC");... % 6.1
        datetime(2026,6,11,12,31,0,0,"TimeZone","UTC") datetime(2026,6,11,12,42,0,0,"TimeZone","UTC");... % 6.2
                datetime(2026,6,11,12,52,0,0,"TimeZone","UTC") datetime(2026,6,11,13,03,0,0,"TimeZone","UTC");... % 7
    ];


%% Parse all data
data=repmat(struct('ref_mag',[],'front_mag',[],'back_mag',[],'GPSaidedINS',[]),1,numel(folders));
for ii=1:numel(folders)
    
    disp(['Data set ' num2str(ii) ' out of ' num2str(numel(folders))]);

    % GPS aided INS. The time base of this system will be the time base for
    % all other sensors
    filename = fullfile('..','..',folders{ii},'GPSINS.txt');
    data(ii).GPSaidedINS=parse_gps_aided_ins_data(filename,startandstoptimes(ii,1),startandstoptimes(ii,2));
 
    % Reference magnetometers
    filename = fullfile('..','..',folders{ii},'RefMag.txt');
    data(ii).ref_mag = parse_ref_mag_data(filename,data(ii).GPSaidedINS.time);

    % Front magnetometer
    filename = fullfile('..','..',folders{ii},'FrontMag.txt');
    data(ii).front_mag = parse_front_mag_data(filename,data(ii).GPSaidedINS.time);

    % Back magnetometer
    filename = fullfile('..','..',folders{ii},'BackMag.txt');
    data(ii).back_mag = parse_front_mag_data(filename,data(ii).GPSaidedINS.time);
end

%% Save the data
filename = fullfile('..','..','Matlab','SkovdeFlightData');
save(filename,"data");


