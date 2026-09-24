function settings = get_settings()

settings.tune_hyper_par=false;                                              % Turn on optimization of hyperparameters (in this case the measurement noise variance)
settings.idx=2:5;                                                           % Flights to use for the map creation. Flights 1-5 follow a lawnmower pattern 
settings.idx_validation_data_set=6:8;                                       % Data set used for validation of the quality of the learned map   
settings.reference_lla=[58.2952919010000;...                                % Reference location (lat,lon, altitude)
    13.9020338060000;...
    166.931992000000];
settings.fs_map = 1;                                                        % Map and validation sampling rate [Hz].
settings.map_cutoff_hz = 0.45;                                              % Low-pass cutoff before downsampling [Hz].
settings.altitude_min = 50;                                                 % Min altitude [m] relative the starting altitude for the data used in the map creation     
settings.pos_front_mag = -1e-3*[958/2+52/2-30 0 0];                         % Position of front sensor in the platform coordinate [m]
settings.pos_back_mag = 1e-3*[958/2-52/2-30 0 0];                           % Position of back sensor in the platform coordinate [m]
settings.B_e_ned = wrldmagm( ...                                            % Fixed part of the magnetic field
    settings.reference_lla(3)/1000, ...
    settings.reference_lla(1), ...
    settings.reference_lla(2), ...
    decyear(datetime(2026,6,10)), 'Custom', ...
    fullfile(fileparts(mfilename('fullpath')),'WMM.COF'));
settings.map.center_spacing = 50;                                           % Spacing between the basis functions of the map [m]
settings.map.margin = settings.map.center_spacing;                          % How much should the grid of basis function extend beyond the area covered by the flight paths
settings.map.sigma = 25;                                                    % Prior on the basis function weights [nT]
settings.map.length_scale =-2*pi*50/log(0.01);                              % Length scale used in the basis functions [m]
settings.time.sigma_g0 = 22.5;                                              % Prior uncertainty (std) on the temporal variation bias [nT]. This also captures total changes in the field level between the flights. 
settings.noise.sigma = [0.74;0.88;0.83;1.41]';                              % Measurement noise/model error standard deviation [nT] for the different flights 
settings.calibration.sigma_ori = 1.34e-4*norm(settings.B_e_ned);
settings.noise.sigma_validation=3;

end
