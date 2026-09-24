function model = fit_mag_map_model(data,settings)
%FIT_MAG_MAP_MODEL Learn the map and flight calibrations in one batch.
%
% The parameter vector is theta = [w; xi^(1); ...; xi^(I)], where each
% eight-element xi contains separate orientation terms and offsets for the
% front and back sensors.

[obs,validationObs] = build_observation_data(data,settings);
basis = build_map_basis(obs,settings);
paramInfo = build_parameter_info(settings,basis);
[H,y,flightId] = build_batch_map_data(obs,basis,paramInfo,settings);

if settings.tune_hyper_par
    options = optimoptions('fminunc', ...
        'Algorithm','quasi-newton', ...
        'Display','iter', ...
        'MaxIterations',50, ...
        'MaxFunctionEvaluations',500,'UseParallel',true);

    eta0 = log([settings.noise.sigma(:); ...
        settings.map.sigma; settings.calibration.sigma_xi]);
    etaHat = fminunc(@(eta) infer_map( ...
        H,y,flightId,paramInfo,complete_settings(settings,eta)), ...
        eta0,options);
    settings = complete_settings(settings,etaHat);
end

[negativeLogMargLikelihood,theta,P] = infer_map( ...
    H,y,flightId,paramInfo,settings);

model = struct();
model.theta = theta;
model.P = P;
model.negative_log_evidence = negativeLogMargLikelihood;
model.paramInfo = paramInfo;
model.basis = basis;
model.settings = settings;
model.obs = obs;
model.val_obs = validationObs;
model.predict_map = @(r_ned) predict_map(r_ned,theta,P,basis,paramInfo);
model.validInfo = run_validation(model);
end


function [negativeLogEvidence,theta,P] = infer_map( ...
    H,y,flightId,paramInfo,settings)
%MAP_EVIDENCE Evaluate (19) and, when requested, the posterior in (16).

p = paramInfo.n_state;
priorStd = zeros(p,1);
priorStd(paramInfo.idx_map) = settings.map.sigma;
for ii = 1:numel(paramInfo.flight)
    priorStd(paramInfo.flight(ii).idx_xi) = settings.calibration.sigma_xi;
end

noiseStd = settings.noise.sigma(flightId);
if nargout == 1
    negativeLogEvidence = linear_gaussian_batch( ...
        H,y,zeros(p,1),diag(priorStd),noiseStd);
else
    [negativeLogEvidence,theta,P] = linear_gaussian_batch( ...
        H,y,zeros(p,1),diag(priorStd),noiseStd);
end
end


function settings = complete_settings(settings,eta)
%COMPLETE_SETTINGS Convert log standard deviations back to model settings.

nFlights = numel(settings.idx);
settings.noise.sigma = exp(eta(1:nFlights));
settings.map.sigma = exp(eta(nFlights+1));
settings.calibration.sigma_xi = exp(eta(nFlights+2));
end


function paramInfo = build_parameter_info(settings,basis)
%BUILD_PARAMETER_INFO Define indices in theta = [w; xi^(1); ...; xi^(I)].

nMap = size(basis.map.centers,1);
paramInfo.idx_map = 1:nMap;
next = nMap + 1;

for ii = 1:numel(settings.idx)
    paramInfo.flight(ii).flight_idx = settings.idx(ii);
    paramInfo.flight(ii).idx_xi = next:(next+7);
    paramInfo.flight(ii).idx_xi_front = next:(next+3);
    paramInfo.flight(ii).idx_xi_back = (next+4):(next+7);
    next = next + 8;
end
paramInfo.n_state = next - 1;
end
