function [logDensity,muU,U] = log_gaussian_marginal(y,meanY,F,noiseStd)
%LOG_GAUSSIAN_MARGINAL Log N(y; meanY, F*F' + diag(noiseStd.^2)).
%
% F is a factor of the parameter-induced covariance. A scalar noiseStd
% gives isotropic noise; a vector allows different noise for each sample.
% The optional outputs describe the posterior of u when
%   y = meanY + F*u + noise,  u ~ N(0,I):
%   u | y ~ N(muU, inv(U'*U)).

residual = y(:) - meanY(:);
n = numel(residual);
noiseStd = noiseStd(:);
if isscalar(noiseStd)
    noiseStd = repmat(noiseStd,n,1);
end
if numel(y) ~= numel(meanY) || size(F,1) ~= n || ...
        numel(noiseStd) ~= n || any(~isfinite(noiseStd)) || ...
        any(noiseStd <= 0)
    error('Incompatible dimensions or nonpositive measurement noise.');
end

G = F./noiseStd;
r = residual./noiseStd;
gram = G.'*G;
rhs = G.'*r;

precision = eye(size(F,2)) + gram;
precision = (precision + precision.')/2;
U = chol(precision,'upper');
muU = U\(U.'\rhs);

logDetCov = 2*sum(log(noiseStd)) + 2*sum(log(diag(U)));
quadratic = r.'*r - rhs.'*muU;
logDensity = -0.5*(n*log(2*pi) + logDetCov + quadratic);
end
