function filter_par = run_validation(obs,model)
%RUN_VALIDATION Validate the learned map and tune the effective noise level.
%
% The validation noise standard deviation settings.noise.sigma_validation
% is estimated by maximizing the validation-data log likelihood.
%
% The optimization is performed over log(sigma_validation), which ensures
% that the estimated standard deviation remains positive.
%
% Output:
%   filter_par
%       Kalman-filter results and validation diagnostics. In addition to
%       the state estimates and NIS values, the structure contains
%
%       filter_par.sigma_validation
%           Maximum-likelihood estimate of the effective validation-noise
%           standard deviation [nT].
%
%       filter_par.log_likelihood
%           Maximized validation-data log likelihood.
%
%       filter_par.optimization
%           Information returned by fminsearch.

    % ---------------------------------------------------------------------
    % Optimization settings.
    % ---------------------------------------------------------------------
    optimOptions = optimset( ...
        'Display','on', ...
        'TolX',1e-4, ...
        'TolFun',1e-3, ...
        'MaxIter',60, ...
        'MaxFunEvals',120);

    % ---------------------------------------------------------------------
    % Extract important variables.
    % ---------------------------------------------------------------------
    theta = model.theta;
    Pmodel = model.P;
    basis = model.basis;
    settings = model.settings;

    paramInfo = build_parameter_info_validation();

    % Number of two-sensor observations.
    N = size(obs.y,1);

    % Stack the front measurements followed by the back measurements.
    yStacked = obs.y(:);

    % Build the map and calibration regression matrices.
    [Hmap,Hcal] = get_mesaurement_matrix_validation( ...
        obs.r_ned,obs.q,basis,paramInfo,settings);

    % ---------------------------------------------------------------------
    % Quantities that do not depend on sigma_validation.
    % ---------------------------------------------------------------------
    thetaMap = theta(model.paramInfo.idx_map);
    Pmap = Pmodel( ...
        model.paramInfo.idx_map,model.paramInfo.idx_map);

    mapResidual = yStacked - Hmap*thetaMap;
    mapCovariance = Hmap*Pmap*Hmap.';

    % ---------------------------------------------------------------------
    % Tune sigma_validation by maximizing the log likelihood.
    %
    % fminsearch performs minimization, so the objective is the negative
    % log likelihood. Optimizing log(sigma) guarantees sigma > 0.
    % ---------------------------------------------------------------------
    logSigma0 = log(settings.noise.sigma_validation);

    objective = @(logSigma) negative_validation_log_likelihood( ...
        logSigma,mapResidual,mapCovariance,Hcal, ...
        paramInfo,settings,N);

    [logSigmaOpt,negativeLogLikelihood,exitflag,output] = ...
        fminsearch(objective,logSigma0,optimOptions);

    sigmaValidation = exp(logSigmaOpt);
    settings.noise.sigma_validation = sigmaValidation;

    % ---------------------------------------------------------------------
    % Run the validation filter once more using the optimized noise level
    % and store the complete filter output.
    % ---------------------------------------------------------------------
    [filter_par,logLikelihood] = run_validation_filter( ...
        sigmaValidation,mapResidual,mapCovariance,Hcal, ...
        paramInfo,settings,N,true);

    filter_par.sigma_validation = sigmaValidation;
    filter_par.log_likelihood = logLikelihood;

    filter_par.optimization = struct();
    filter_par.optimization.exitflag = exitflag;
    filter_par.optimization.output = output;
    filter_par.optimization.negative_log_likelihood = ...
        negativeLogLikelihood;
end


function negativeLogLikelihood = negative_validation_log_likelihood( ...
    logSigma,mapResidual,mapCovariance,Hcal,paramInfo,settings,N)
%NEGATIVE_VALIDATION_LOG_LIKELIHOOD Objective used by fminsearch.

    sigmaValidation = exp(logSigma);

    [~,logLikelihood] = run_validation_filter( ...
        sigmaValidation,mapResidual,mapCovariance,Hcal, ...
        paramInfo,settings,N,false);

    negativeLogLikelihood = -logLikelihood;
end


function [filter_par,logLikelihood] = run_validation_filter( ...
    sigmaValidation,mapResidual,mapCovariance,Hcal, ...
    paramInfo,settings,N,storeOutput)
%RUN_VALIDATION_FILTER Run the prewhitened validation Kalman filter.
%
% The original stacked measurement model is
%
%   y - Hmap*thetaMap = Hcal*x + e,
%
% where
%
%   cov(e) = sigmaValidation^2*I + Hmap*Pmap*Hmap'.
%
% If Sigma = L*L', prewhitening gives
%
%   L\(y - Hmap*thetaMap) = (L\Hcal)*x + w,
%
% with w ~ N(0,I). Both the measurements and Hcal must therefore be
% transformed by L.

    % ---------------------------------------------------------------------
    % Prewhiten the validation measurements and calibration regressors.
    % ---------------------------------------------------------------------
    Sigma = sigmaValidation^2*eye(2*N) + mapCovariance;

    [L,cholFlag] = chol(Sigma,'lower');

    if cholFlag ~= 0
        logLikelihood = -Inf;
        filter_par = [];
        return;
    end

    yTilde = L\mapResidual;
    HcalTilde = L\Hcal;

    % The Kalman filter below evaluates the likelihood in the whitened
    % domain. Subtracting log(det(L)) converts it back to the likelihood of
    % the original measurements.
    logJacobian = sum(log(diag(L)));

    % ---------------------------------------------------------------------
    % Initialize the validation states.
    % ---------------------------------------------------------------------
    [x,P] = get_initial_states(paramInfo,settings);

    ssm.F = eye(paramInfo.nstate);
    ssm.Q = zeros(paramInfo.nstate);
    ssm.R = eye(2);

    logLikelihoodWhitened = 0;

    if storeOutput
        filter_par = struct();

        filter_par.NIS = zeros(N,1);
        filter_par.normalizedInnovation=zeros(N,2);
        filter_par.g = zeros(N,1);
        filter_par.g_var = zeros(N,1);

        filter_par.ori_bias_front = zeros(N,3);
        filter_par.ori_bias_front_cov_diag = zeros(N,3);

        filter_par.ori_bias_back = zeros(N,3);
        filter_par.ori_bias_back_cov_diag = zeros(N,3);
    else
        filter_par = [];
    end

    % ---------------------------------------------------------------------
    % Run the Kalman filter.
    % ---------------------------------------------------------------------
    for nn = 1:N

        % Front and back measurements at the current sample.
        y = [yTilde(nn); yTilde(nn+N)];

        % The calibration matrix must be prewhitened using the same
        % Cholesky factor as the measurements.
        ssm.H = [HcalTilde(nn,:); HcalTilde(nn+N,:)];

        [x,P,logL,NIS,normalizedInnovation] = step_kf(y,x,P,ssm);

        logLikelihoodWhitened = logLikelihoodWhitened + logL;

        if storeOutput
            filter_par.NIS(nn) = NIS;
            filter_par.normalizedInnovation(nn,:) = normalizedInnovation';

            filter_par.g(nn) = x(paramInfo.idx_g);
            filter_par.g_var(nn) = ...
                P(paramInfo.idx_g,paramInfo.idx_g);

            filter_par.ori_bias_front(nn,:) = ...
                x(paramInfo.idx_ori_front).';

            filter_par.ori_bias_front_cov_diag(nn,:) = ...
                diag(P( ...
                paramInfo.idx_ori_front, ...
                paramInfo.idx_ori_front)).';

            filter_par.ori_bias_back(nn,:) = ...
                x(paramInfo.idx_ori_back).';

            filter_par.ori_bias_back_cov_diag(nn,:) = ...
                diag(P( ...
                paramInfo.idx_ori_back, ...
                paramInfo.idx_ori_back)).';
        end
    end

    % Convert the likelihood from the whitened measurements back to the
    % original measurement domain.
    logLikelihood = logLikelihoodWhitened - logJacobian;
end


function paramInfo = build_parameter_info_validation()
%BUILD_PARAMETER_INFO_VALIDATION Define validation-state layout.
%
% State:
%   x = [
%       theta_ori_front
%       theta_ori_back
%       g
%   ]

    paramInfo = struct();
    paramInfo.idx_ori_front = 1:3;
    paramInfo.idx_ori_back = 4:6;
    paramInfo.idx_g = 7;
    paramInfo.nstate = 7;
end


function [x,P] = get_initial_states(paramInfo,settings)
%GET_INITIAL_STATES Construct validation-state prior.

    x = zeros(paramInfo.nstate,1);
    P = zeros(paramInfo.nstate,paramInfo.nstate);

    P(paramInfo.idx_ori_front,paramInfo.idx_ori_front) = ...
        settings.calibration.sigma_ori^2*eye(3);

    P(paramInfo.idx_ori_back,paramInfo.idx_ori_back) = ...
        settings.calibration.sigma_ori^2*eye(3);

    P(paramInfo.idx_g,paramInfo.idx_g) = ...
        settings.time.sigma_g0^2;
end


function [Hmap,Hcal] = get_mesaurement_matrix_validation( ...
    r_ned,q,basis,paramInfo,settings)
%GET_MESAUREMENT_MATRIX_VALIDATION Build validation regression matrices.
%
% Measurements are stacked as
%
%   y_stacked = [y_front; y_back].

    N = size(r_ned,1);
    M = size(basis.map.centers,1);

    Hmap = zeros(2*N,M);
    Hcal = zeros(2*N,paramInfo.nstate);

    for nn = 1:N

        [PhiFront,PhiBack,z] = get_measurement_regressors( ...
            r_ned(nn,:),q(nn,:),basis,settings);

        Hmap(nn,:) = PhiFront;
        Hmap(N+nn,:) = PhiBack;

        Hcal(nn,paramInfo.idx_ori_front) = z;
        Hcal(nn,paramInfo.idx_g) = 1;

        Hcal(N+nn,paramInfo.idx_ori_back) = z;
        Hcal(N+nn,paramInfo.idx_g) = 1;
    end
end
