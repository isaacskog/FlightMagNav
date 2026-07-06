function obs = build_observation_data(data,settings)
%BUILD_OBSERVATION_DATA Build low-rate front/back magnetometer observations.
%
% obs(ii) contains one flight. obs(ii).y has two columns:
%   column 1: front magnetometer total field
%   column 2: back magnetometer total field

obs = repmat(struct( ...
    'flight_idx',[], ...
    'time',[], ...
    't_sec',[], ...
    'y',[], ...
    'r_ned',[], ...
    'q',[], ...
    'yaw',[]), ...
    1,numel(settings.idx));

% Design low pass filter    
fsIn = 100;
fc = 0.45*settings.fs_map;
filterOrder = 3;
[b,a] = butter(filterOrder,fc/(fsIn/2),'low');

% Extract data
for ii = 1:numel(settings.idx)

    kk = settings.idx(ii);
    ins = data(kk).GPSaidedINS;

    tIn = seconds(ins.time - ins.time(1));
    tOut = (0:1/settings.fs_map:tIn(end)).';

    % Low pass filter using zero-phase filtering 
    yFront = filtfilt(b,a,data(kk).front_mag.tot_field);
    yBack  = filtfilt(b,a,data(kk).back_mag.tot_field);

    yFront = interp1(tIn,yFront,tOut,'nearest');
    yBack  = interp1(tIn,yBack, tOut,'nearest');

    rBody = lla_to_local_ned( ...
        interp1(tIn,ins.latitude, tOut,'nearest'), ...
        interp1(tIn,ins.longitude,tOut,'nearest'), ...
        interp1(tIn,ins.altitude, tOut,'nearest'), ...
        settings.reference_lla);

    q = interp1(tIn,ins.quaternion,tOut,'nearest');

    % Get the index of the data points above the specified minimum
    % height
    idx = -rBody(:,3) >= settings.altitude_min;

    % STore the data
    obs(ii).flight_idx = kk;
    obs(ii).time = ins.time(1) + seconds(tOut(idx));
    obs(ii).t_sec = tOut(idx);
    obs(ii).y = [yFront(idx) yBack(idx)]-mean([yFront(idx); yBack(idx)]);
    obs(ii).r_ned = rBody(idx,:);
    obs(ii).q = q(idx,:);

end
end