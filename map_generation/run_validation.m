function validInfo = run_validation(model)
%RUN_VALIDATION Estimate effective noise on independent flights using (22).
%
% Each flight is evaluated separately, conditional on the learned map.
% The map weights and the new flight's calibration xi are marginalized
% when estimating its effective measurement-noise standard deviation.

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
L0 = blkdiag(Lmap,sigmaXi*eye(8));
mu0 = [muMap;zeros(8,1)];

options = optimset( ...
    'Display','off', ...
    'TolX',1e-4, ...
    'TolFun',1e-3, ...
    'MaxIter',60, ...
    'MaxFunEvals',120);

for ii = 1:numel(obs)
    [H,y] = build_validation_data(obs(ii),model.basis,model.settings);
    F = H*L0;
    residual = y - H*mu0;
    gram = F.'*F;
    rhs = F.'*residual;
    residualNorm2 = residual.'*residual;
    objective = @(logSigma) validation_evidence( ...
        logSigma,gram,rhs,residualNorm2,numel(y));

    [logSigma,~,exitflag,output] = fminsearch( ...
        objective,log(model.settings.noise.sigma_validation),options);
    sigma = exp(logSigma);
    [negativeLogEvidence,muU,U] = validation_evidence( ...
        logSigma,gram,rhs,residualNorm2,numel(y));
    mu = mu0 + L0*muU;
    C = L0/U;
    P = C*C.';

    idxXi = numel(idxMap)+(1:8);
    validInfo(ii).flight_index = model.settings.idx_validation_data_set(ii);
    validInfo(ii).sigma_validation = sigma;
    validInfo(ii).log_likelihood = -negativeLogEvidence;
    validInfo(ii).xi_front = mu(idxXi(1:4));
    validInfo(ii).xi_back = mu(idxXi(5:8));
    validInfo(ii).xi_cov = P(idxXi,idxXi);
    validInfo(ii).optimization.exitflag = exitflag;
    validInfo(ii).optimization.output = output;
end
end


function [negativeLogEvidence,muU,U] = validation_evidence( ...
    logSigma,gram,rhs,residualNorm2,n)
%VALIDATION_EVIDENCE Apply (23) through a parameter-space Cholesky factor.

sigma2 = exp(2*logSigma);
p = size(gram,1);
precision = eye(p) + gram/sigma2;
precision = (precision + precision.')/2;
U = chol(precision,'upper');
h = rhs/sigma2;
muU = U\(U.'\h);

logDetPyy = n*log(sigma2) + 2*sum(log(diag(U)));
quadratic = residualNorm2/sigma2 - h.'*muU;
negativeLogEvidence = 0.5*(n*log(2*pi) + logDetPyy + quadratic);
end


function [H,y] = build_validation_data(obs,basis,settings)
%BUILD_VALIDATION_DATA Construct [A B] and stack [front; back] observations.

n = size(obs.y,1);
nMap = size(basis.map.centers,1);
H = zeros(2*n,nMap+8);
y = obs.y(:);

for kk = 1:n
    [phiFront,phiBack,z] = get_measurement_regressors( ...
        obs.r_ned(kk,:),obs.q(kk,:),basis,settings);

    H(kk,1:nMap) = phiFront;
    H(n+kk,1:nMap) = phiBack;
    H(kk,nMap+(1:4)) = [z 1];
    H(n+kk,nMap+(5:8)) = [z 1];
end
end
