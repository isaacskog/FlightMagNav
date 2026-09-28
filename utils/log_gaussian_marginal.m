function [logDensity,muU,U] = log_gaussian_marginal( ...
    gram,rhs,residualNorm2,logDetNoise,n)
%LOG_GAUSSIAN_MARGINAL Log N(y; H*mu0, R + F*F') without forming its covariance.
%
% For residual = y - H*mu0 and F = H*L0, supply the sufficient statistics
%   gram          = F'*(R\F)
%   rhs           = F'*(R\residual)
%   residualNorm2 = residual'*(R\residual)
%   logDetNoise   = log(det(R)).
% The optional outputs give the posterior of u in theta = mu0 + L0*u:
%   u | y ~ N(muU, inv(U'*U)).

p = size(gram,1);
precision = eye(p) + gram;
precision = (precision + precision.')/2;
U = chol(precision,'upper');
muU = U\(U.'\rhs);

logDetCov = logDetNoise + 2*sum(log(diag(U)));
quadratic = residualNorm2 - rhs.'*muU;
logDensity = -0.5*(n*log(2*pi) + logDetCov + quadratic);
end
