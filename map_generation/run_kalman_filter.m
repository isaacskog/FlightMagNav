function [NlogL_tot,x,P,filter_par] = run_kalman_filter(obs,basis,paramInfo,settings)
%RUN_KALMAN_FILTER Estimate map and flight-specific calibration states.

x = zeros(paramInfo.n_state,1);
P = zeros(paramInfo.n_state,paramInfo.n_state);
P(paramInfo.idx_map,paramInfo.idx_map) = ...
    settings.map.sigma^2 *eye(numel(paramInfo.idx_map));

% Flight-specific calibration priors.
for ii = 1:numel(obs)
    P(paramInfo.flight(ii).idx_ori_front,paramInfo.flight(ii).idx_ori_front) = ...
        settings.calibration.sigma_ori^2*eye(3);

    P(paramInfo.flight(ii).idx_ori_back,paramInfo.flight(ii).idx_ori_back) = ...
        settings.calibration.sigma_ori^2*eye(3);

    P(paramInfo.flight(ii).idx_g,paramInfo.flight(ii).idx_g) = settings.time.sigma_g0^2;

    if ii == 1
        % First flight defines the reference level.
        P(paramInfo.flight(ii).idx_g,paramInfo.flight(ii).idx_g) = 0;
    end
end

% Total negative log likelihood      
NlogL_tot = 0;

% Structure of output data beyond the map parameters
filter_par = repmat(struct('NIS',[],'g',[],'g_var',[],'ori_bias_front',[],'ori_bias_front_cov_diag',[],'ori_bias_back',[],'ori_bias_back_cov_diag',[]),1,numel(obs));

for ii = 1:numel(obs)

    % Allocate memory and store initial values
    filter_par(ii).NIS = zeros(size(obs(ii).y,1),1);
    filter_par(ii).g = zeros(size(obs(ii).y,1),1);
    filter_par(ii).g_var = zeros(size(obs(ii).y,1),1);
    filter_par(ii).ori_bias = zeros(size(obs(ii).y,1),3);
    filter_par(ii).ori_bias_cov_diag = zeros(size(obs(ii).y,1),3);
    filter_par(ii).g(1) = x(paramInfo.flight(ii).idx_g);
    filter_par(ii).g_var(1) = P(paramInfo.flight(ii).idx_g,paramInfo.flight(ii).idx_g);
    filter_par(ii).ori_bias_front(1,:) = x(paramInfo.flight(ii).idx_ori_front).';
    filter_par(ii).ori_bias_front_cov_diag(1,:) = diag(P(paramInfo.flight(ii).idx_ori_front,paramInfo.flight(ii).idx_ori_front)).';
    filter_par(ii).ori_bias_back(1,:) = x(paramInfo.flight(ii).idx_ori_back).';
    filter_par(ii).ori_bias_back_cov_diag(1,:) = diag(P(paramInfo.flight(ii).idx_ori_back,paramInfo.flight(ii).idx_ori_back)).';


    % Set the measurement covariance for the specific flight
    ssm.R = settings.noise.sigma(ii)^2*eye(2);

    % Run the Kalman filter 
    for kk = 2:numel(obs(ii).t_sec)

        % Construct the measurement matrix    
        ssm.H = get_mesaurement_matrix( ...
            obs(ii).r_ned(kk,:), ...
            obs(ii).q(kk,:), ...
            ii, ...
            basis, ...
            paramInfo, ...
            settings);

        % Do one step of the Kalman filter algorithm
        [x,P,logL,NIS] = step_kf(obs(ii).y(kk,:).',x,P,ssm);

        % Accumulate the negative log likelihood
        NlogL_tot = NlogL_tot - logL;

        % Store data
        filter_par(ii).NIS(kk) = NIS;
        filter_par(ii).g(kk) = x(paramInfo.flight(ii).idx_g);
        filter_par(ii).g_var(kk) = P(paramInfo.flight(ii).idx_g,paramInfo.flight(ii).idx_g);
        filter_par(ii).ori_bias_front(kk,:) = x(paramInfo.flight(ii).idx_ori_front).';
        filter_par(ii).ori_bias_front_cov_diag(kk,:) = diag(P(paramInfo.flight(ii).idx_ori_front,paramInfo.flight(ii).idx_ori_front)).';
        filter_par(ii).ori_bias_back(kk,:) = x(paramInfo.flight(ii).idx_ori_back).';
        filter_par(ii).ori_bias_back_cov_diag(kk,:) = diag(P(paramInfo.flight(ii).idx_ori_back,paramInfo.flight(ii).idx_ori_back)).';
    end
end    
end