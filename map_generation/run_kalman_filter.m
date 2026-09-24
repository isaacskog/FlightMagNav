function [NlogL_tot,theta,P,filter_par] = run_kalman_filter( ...
    obs,basis,paramInfo,settings)
%RUN_KALMAN_FILTER Infer the static map and flight-specific calibration vectors.

theta = zeros(paramInfo.n_state,1);
P = zeros(paramInfo.n_state,paramInfo.n_state);
P(paramInfo.idx_map,paramInfo.idx_map) = ...
    settings.map.sigma^2*eye(numel(paramInfo.idx_map));

for ii = 1:numel(obs)
    idx_xi = paramInfo.flight(ii).idx_xi;
    P(idx_xi,idx_xi) = settings.calibration.sigma_xi^2*eye(8);
end

NlogL_tot = 0;
filter_par = repmat(struct( ...
    'flight_index',[], ...
    'NIS',[], ...
    'xi_front',[], ...
    'xi_front_cov_diag',[], ...
    'xi_back',[], ...
    'xi_back_cov_diag',[]),1,numel(obs));

for ii = 1:numel(obs)
    n = size(obs(ii).y,1);
    idx_front = paramInfo.flight(ii).idx_xi_front;
    idx_back = paramInfo.flight(ii).idx_xi_back;

    filter_par(ii).flight_index = settings.idx(ii);
    filter_par(ii).NIS = zeros(n,1);
    filter_par(ii).xi_front = zeros(n,4);
    filter_par(ii).xi_front_cov_diag = zeros(n,4);
    filter_par(ii).xi_back = zeros(n,4);
    filter_par(ii).xi_back_cov_diag = zeros(n,4);

    ssm.R = settings.noise.sigma(ii)^2*eye(2);

    for kk = 1:n
        ssm.H = get_mesaurement_matrix( ...
            obs(ii).r_ned(kk,:),obs(ii).q(kk,:), ...
            ii,basis,paramInfo,settings);

        [theta,P,logL,NIS] = step_kf(obs(ii).y(kk,:).',theta,P,ssm);
        NlogL_tot = NlogL_tot - logL;

        filter_par(ii).NIS(kk) = NIS;
        filter_par(ii).xi_front(kk,:) = theta(idx_front).';
        filter_par(ii).xi_front_cov_diag(kk,:) = diag(P(idx_front,idx_front)).';
        filter_par(ii).xi_back(kk,:) = theta(idx_back).';
        filter_par(ii).xi_back_cov_diag(kk,:) = diag(P(idx_back,idx_back)).';
    end
end
end
