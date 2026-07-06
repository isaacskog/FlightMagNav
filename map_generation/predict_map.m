function [m,var_m] = predict_map(r_ned,theta,P,basis,paramInfo)
% Predict static map at local NED positions.

Phi = spatial_rbf_2d( ...
    r_ned(:,1:2), ...
    basis.map.centers, ...
    basis.map.length_scale);

idx = paramInfo.idx_map;
theta_map = theta(idx);
P_map = P(idx,idx);
m = Phi*theta_map;
var_m = sum((Phi*P_map).*Phi,2);
end