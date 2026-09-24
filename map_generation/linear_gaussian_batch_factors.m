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
p = size(F,2);
precision = eye(p) + G.'*G;
precision = (precision + precision.')/2;
U = chol(precision,'upper');
h = G.'*r;
muU = U\(U.'\h);

% Determinant lemma and Woodbury identity for P_yy = R + H*P0*H'.
logDetPyy = 2*sum(log(noiseStd)) + 2*sum(log(diag(U)));
quadratic = r.'*r - h.'*muU;
negativeLogEvidence = 0.5*(n*log(2*pi) + logDetPyy + quadratic);
end
