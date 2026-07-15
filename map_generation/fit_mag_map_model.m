function model = fit_mag_map_model(data,settings)
%FIT_MAG_MAP_MODEL Estimate a magnetic-field map using a Kalman filter.
%
% The state contains:
%   - one common static map,
%   - one calibration block per flight,
%   - one temporal random-walk state per flight.


% Extract the observations need to build the map
[obs,validation_obs] = build_observation_data(data,settings);

% Build the basis
basis = build_map_basis(obs,settings);

% Construct information about parameter structure of the state vector
paramInfo = build_parameter_info(settings,basis);

% Check if hyperparamters should be tuned
if(settings.tune_hyper_par)

    % Set up optimization algorithm
    options = optimoptions('fminunc', ...
        'Algorithm','quasi-newton', ...
        'Display','iter', ...
        'MaxIterations',10, ...
        'MaxFunctionEvaluations',500, ...
        'UseParallel',true, ...
        'PlotFcn','optimplotfval');

    % Get initial conditions
    theta0 = get_init(settings);

    % Run optimization
    theta_hat = fminunc( ...
        @(theta) run_kalman_filter(obs,basis,paramInfo,complete_settings(settings,theta)), ...
        theta0, ...
        options);

    % Store the optimized parameters in the settings
    settings = complete_settings(settings,theta_hat);
end

% Run the map learning 
[NlogL,theta,P,filterInfo] = run_kalman_filter(obs,basis,paramInfo,settings); %#ok<ASGLU>


% Save the result into a struct
model = struct();
model.theta = theta;
model.P = P;
model.paramInfo = paramInfo;
model.basis = basis;
model.settings = settings;
model.obs = obs;
model.val_obs = validation_obs;
model.filterInfo = filterInfo;
model.predict_map = @(r_ned) predict_map(r_ned,theta,P,basis,paramInfo);
model.predict_map_lla = @(lat,lon,alt) ...
    predict_map(lla_to_local_ned(lat,lon,alt,settings.reference_lla), ...
    theta,P,basis,paramInfo);


% Run the validation 
model.validInfo=run_validation(validation_obs,model);



end


function theta0 = get_init(settings)
%Initial log-hyperparameter vector theta.
    theta0 = [log(settings.noise.sigma); log(settings.map.sigma); log(settings.time.sigma_g0); log(settings.calibration.sigma_ori)];
end

function settings = complete_settings(settings,theta)
% Update settings with the values in theta 
        settings.noise.sigma=exp(theta(1:4));
        settings.map.sigma=exp(theta(5));
        settings.time.sigma_g0=exp(theta(6));
        settings.calibration.sigma_ori=exp(theta(7));
end


function paramInfo = build_parameter_info(settings,basis)
%BUILD_PARAMETER_INFO Define full parameter/state layout.
%
% State:
%   x = [
%       theta_map
%       theta_ori_front_1
%       theta_ori_back_1
%       ...
%       theta_ori_front_M
%       theta_ori_back_M
%       g_1
%       ...
%       g_M
%   ]
%
% The map and calibration states are static. The g_i states are random
% walks, but only the active flight's g_i receives process noise.

    nMap = size(basis.map.centers,1);
    nFlights = numel(settings.idx);

    paramInfo = struct();
    paramInfo.idx_map = 1:nMap;

    next = nMap + 1;

    for ii = 1:nFlights
        paramInfo.flight(ii).flight_idx = settings.idx(ii);

        paramInfo.flight(ii).idx_ori_front = next:(next+2);
        next = next + 3;

        paramInfo.flight(ii).idx_ori_back = next:(next+2);
        next = next + 3;
    end

    for ii = 1:nFlights
        paramInfo.flight(ii).idx_g = next;
        next = next + 1;
    end

    paramInfo.n_static = next - 1;
    paramInfo.n_state = next - 1;
end





















