function [Phi_front,Phi_back,z] = get_measurement_regressors(r_ned,q,basis,settings)

% Get rotation matrices
Rb2n_enu = quat2rotm(q);                  % Note that Xsense uses ENU. Body to Navigation frame 
Rb2n=[0 1 0; 1 0 0; 0 0 -1]*Rb2n_enu;
Rn2b = Rb2n';

% Equation (7): use the direction of the nominal field, not its magnitude.
m0 = settings.m0(:);
z = (Rn2b*(m0/norm(m0))).';

% Get position of front and back sensor in NED coordinates
r_front = r_ned + (Rb2n*settings.pos_front_mag(:)).';
r_back  = r_ned + (Rb2n*settings.pos_back_mag(:)).';

% Get basis function regressors
Phi_front = spatial_rbf_2d(r_front(1:2),basis.map.centers,basis.map.length_scale);
Phi_back  = spatial_rbf_2d(r_back(1:2), basis.map.centers,basis.map.length_scale);

end
