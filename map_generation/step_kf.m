function [x,P,logL,NIS] = step_kf(y,x,P,ssm)


% Innovation covariance
S = ssm.H*P*ssm.H' + ssm.R;

% Innovation 
z = y - ssm.H*x;

% Normalized innovation squared 
NIS = z'*(S\z);

% Kalman gain
K = (P*ssm.H')/S;

% Filter update
x = x + K*z;
P = (eye(size(P))-K*ssm.H)*P;

% Log of marginalized likelihood
logL =logmvnpdf(z,0,S);
end