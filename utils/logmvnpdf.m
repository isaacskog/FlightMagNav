function logp = logmvnpdf(x,mu,Sigma)
% Log density of N(mu,Sigma) at x.

D = size(x,1);
xc = x-mu;

logp = -0.5*(D*log(2*pi) + logdet(Sigma) + xc'*(Sigma\xc));
end


function y = logdet(A)
% Log determinant for symmetric positive definite matrix.

U = chol(A);
y = 2*sum(log(diag(U)));
end