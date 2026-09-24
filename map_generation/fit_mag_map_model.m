function model = fit_mag_map_model(data,settings)
%FIT_MAG_MAP_MODEL Estimate a magnetic-field map using a Kalman filter.
%
% The state contains a common static map and an eight-element calibration
% vector xi for each flight: three orientation terms and one constant bias
% for each of the front and back sensors.


% Extract the observations need to build the map
[obs,validation_obs] = build_observation_data(data,settings);

% Build the basis
basis = build_map_basis(obs,settings);

% Construct information about parameter structure of the state vector
paramInfo = build_parameter_info(settings,basis);

% Check if hyperparamters should be tuned
if(settings.tune_hyper_par)

       fprintf('\nTunes hyperparameters \n');

    % Set up optimization algorithm
    options = optimoptions('fminunc', ...
        'Algorithm','quasi-newton', ...
        'Display','iter', ...
        'MaxIterations',50, ...
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
 fprintf('\n Learns map \n');
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

% Run the validation 
 fprintf('\n Runs validation \n');
model.validInfo=run_validation(model);

end


function theta0 = get_init(settings)
%Initial log-hyperparameter vector theta.
    theta0 = [log(settings.noise.sigma(:)); ...
        log(settings.map.sigma); log(settings.calibration.sigma_xi)];
end

function settings = complete_settings(settings,theta)
% Update settings with the values in theta 
        nFlights = numel(settings.idx);
        settings.noise.sigma = exp(theta(1:nFlights));
        settings.map.sigma = exp(theta(nFlights+1));
        settings.calibration.sigma_xi = exp(theta(nFlights+2));
end


function paramInfo = build_parameter_info(settings,basis)
%BUILD_PARAMETER_INFO Define full parameter/state layout.
%
% State: theta = [w; xi^(1); ...; xi^(I)]. Each xi contains
% [orientation_front(3); bias_front; orientation_back(3); bias_back].

    nMap = size(basis.map.centers,1);
    nFlights = numel(settings.idx);

    paramInfo = struct();
    paramInfo.idx_map = 1:nMap;

    next = nMap + 1;

    for ii = 1:nFlights
        paramInfo.flight(ii).flight_idx = settings.idx(ii);

        paramInfo.flight(ii).idx_xi = next:(next+7);
        paramInfo.flight(ii).idx_xi_front = next:(next+3);
        paramInfo.flight(ii).idx_xi_back = (next+4):(next+7);
        next = next + 8;
    end

    paramInfo.n_state = next - 1;
end



















