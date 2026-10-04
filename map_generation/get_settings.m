function settings = get_settings()

settings.tune_hyper_par=true;                                              % Estimate noise and prior scales by model evidence/marginal likelihood.
settings.idx=2:5;                                                           % Flights to use for the map creation. Flights 1-5 follow a lawnmower pattern 
settings.idx_validation_data_set=6:8;                                       % Data set used for validation of the quality of the learned map   
settings.reference_lla=[58.2952919010000;...                                % Reference location (lat,lon, altitude)
    13.9020338060000;...
    166.931992000000];
settings.fs_map = 1;                                                        % Map and validation sampling rate [Hz].
settings.altitude_min = 50;                                                 % Min altitude [m] relative the starting altitude for the data used in the map creation     
settings.pos_front_mag = -1e-3*[958/2+52/2-30 0 0];                         % Position of front sensor in the platform coordinate [m]
settings.pos_back_mag = 1e-3*[958/2-52/2-30 0 0];                           % Position of back sensor in the platform coordinate [m]
settings.m0 =1.0e+04*[  1.574159830418450; ...                              % Nominal Earth-field vector in NED [nT] according to WMM-2025, WMM.COF, 2026-06-10, at reference_lla  
   0.165375165966756;...
   4.890952638323985]';                               
settings.map.center_spacing = 50;                                           % Spacing between the basis functions of the map [m]
settings.map.margin = settings.map.center_spacing;                          % How much should the grid of basis function extend beyond the area covered by the flight paths
settings.map.sigma = 25;                                                    % Prior on the basis function weights [nT]
settings.map.length_scale =60;                                              % Length scale used in the basis functions [m]
settings.noise.sigma = [0.6 0.8 0.8 1.3]';                                  % Measurement noise/model error standard deviation [nT] for the different flights 
settings.calibration.sigma_xi = 10;                                         % Initial prior standard deviation for the calibration coefficients [nT].
settings.noise.sigma_validation=3;                                          % Initial value for the noise standard deviation in optimization done in the validation 

end
