function [negativeLogEvidence,muU,U] = linear_gaussian_batch_factors( ...
    F,residual,noiseStd)
%LINEAR_GAUSSIAN_BATCH_FACTORS Evidence with theta = mu0 + L0*u.
%
% The supplied F = H*L0 and residual = y - H*mu0 can be reused when only
% the measurement-noise standard deviation changes.

residual = residual(:);
noiseStd = noiseStd(:);
n = numel(residual);
if isscalar(noiseStd)
    noiseStd = repmat(noiseStd,n,1);
end
if size(F,1) ~= n || numel(noiseStd) ~= n || ...
        any(~isfinite(noiseStd)) || any(noiseStd <= 0)
    error('Incompatible dimensions or nonpositive measurement noise.');
end

G = F./noiseStd;
r = residual./noiseStd;
[logDensity,muU,U] = log_gaussian_marginal( ...
    G.'*G,G.'*r,r.'*r,2*sum(log(noiseStd)),n);
negativeLogEvidence = -logDensity;
end
