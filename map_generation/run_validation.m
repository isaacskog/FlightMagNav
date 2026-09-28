function validInfo = run_validation(model)
%RUN_VALIDATION Estimate effective noise on independent flights using (14).
%
% The predictive mean and covariance follow (15b)-(15c). Each flight is
% evaluated separately, conditional on the posterior map from learning.

obs = model.val_obs;
validInfo = repmat(struct( ...
    'flight_index',[], ...
    'sigma_validation',[], ...
    'log_likelihood',[], ...
    'xi_front',[], ...
    'xi_back',[], ...
    'xi_cov',[], ...
    'optimization',[]),size(obs));

idxMap = model.paramInfo.idx_map;
muMap = model.theta(idxMap);
Pmap = model.P(idxMap,idxMap);
Lmap = chol((Pmap + Pmap.')/2,'lower');
sigmaXi = model.settings.calibration.sigma_xi;

options = optimset( ...
    'Display','off', ...
    'TolX',1e-4, ...
    'TolFun',1e-3, ...
    'MaxIter',60, ...
    'MaxFunEvals',120);

for ii = 1:numel(obs)
    [A,B,y] = build_validation_data(obs(ii),model.basis,model.settings);

    % Equation (15): y ~ N(meanY, F*F' + sigma^2*I).
    meanY = A*muMap;
    F = [A*Lmap, sigmaXi*B];
    objective = @(logSigma) -log_gaussian_marginal( ...
        y,meanY,F,exp(logSigma));

    [logSigma,~,exitflag,output] = fminsearch( ...
        objective,log(model.settings.noise.sigma_validation),options);
    sigma = exp(logSigma);
    [logDensity,muU,U] = log_gaussian_marginal(y,meanY,F,sigma);
    [xi,Pxi] = validation_calibration_posterior( ...
        muU,U,sigmaXi,numel(idxMap));

    validInfo(ii).flight_index = model.settings.idx_validation_data_set(ii);
    validInfo(ii).sigma_validation = sigma;
    validInfo(ii).log_likelihood = logDensity;
    validInfo(ii).xi_front = xi(1:4);
    validInfo(ii).xi_back = xi(5:8);
    validInfo(ii).xi_cov = Pxi;
    validInfo(ii).optimization.exitflag = exitflag;
    validInfo(ii).optimization.output = output;
end
end


function [xi,Pxi] = validation_calibration_posterior(muU,U,sigmaXi,nMap)
%VALIDATION_CALIBRATION_POSTERIOR Recover the new flight's calibration.

idxXi = nMap+(1:8);
xi = sigmaXi*muU(idxXi);
C = U\eye(size(U));
Cxi = sigmaXi*C(idxXi,:);
Pxi = Cxi*Cxi.';
end


function [A,B,y] = build_validation_data(obs,basis,settings)
%BUILD_VALIDATION_DATA Construct (15b)-(15c) for one validation flight.

n = size(obs.y,1);
nMap = size(basis.map.centers,1);
A = zeros(2*n,nMap);
B = zeros(2*n,8);
y = obs.y(:);

for kk = 1:n
    [phiFront,phiBack,z] = get_measurement_regressors( ...
        obs.r_ned(kk,:),obs.q(kk,:),basis,settings);

    A(kk,:) = phiFront;
    A(n+kk,:) = phiBack;
    B(kk,1:4) = [z 1];
    B(n+kk,5:8) = [z 1];
end
end
