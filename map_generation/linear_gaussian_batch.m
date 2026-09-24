function [negativeLogEvidence,mu,P] = linear_gaussian_batch( ...
    H,y,mu0,L0,noiseStd)
%LINEAR_GAUSSIAN_BATCH Posterior and evidence for a linear Gaussian model.
%
%   y | theta ~ N(H*theta,diag(noiseStd.^2))
%   theta     ~ N(mu0,L0*L0')
%
% The result is equivalent to (16) and (19) in the paper. The factorized
% calculation avoids forming the large measurement covariance P_yy.

y = y(:);
mu0 = mu0(:);
p = numel(mu0);
if size(H,1) ~= numel(y) || size(H,2) ~= p || ~isequal(size(L0),[p p])
    error('Incompatible sizes of H, y, mu0 and L0.');
end

F = H*L0;
residual = y - H*mu0;
[negativeLogEvidence,muU,U] = linear_gaussian_batch_factors( ...
    F,residual,noiseStd);

if nargout > 1
    mu = mu0 + L0*muU;
    if nargout > 2
        C = L0/U;
        P = C*C.';
        P = (P + P.')/2;
    end
end
end
