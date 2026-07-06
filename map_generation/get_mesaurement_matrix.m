function H = get_mesaurement_matrix(r_ned,q,flightNumber,basis,paramInfo,settings)
% Create two-row measurement matrix.
%
% Measurement:
%   y = [y_front; y_back]
%
% Only the parameter block for flightNumber is active.

% Get rotation matrices
Rb2n_enu = quat2rotm(q);                  % Note that Xsense uses ENU. Body to Navigation frame 
Rb2n=[0 1 0; 1 0 0; 0 0 -1]*Rb2n_enu;
Rn2b = Rb2n';

% Create regressor for the orientation dependent bias
z = (Rn2b*settings.B_e_ned).';

% Get position of front and back sensor in NED coordinates
r_front = r_ned + (Rb2n*settings.pos_front_mag(:)).';
r_back  = r_ned + (Rb2n*settings.pos_back_mag(:)).';

% Get basis function regressors
Phi_front = spatial_rbf_2d(r_front(1:2),basis.map.centers,basis.map.length_scale);
Phi_back  = spatial_rbf_2d(r_back(1:2), basis.map.centers,basis.map.length_scale);

% Build the matrix
H = zeros(2,paramInfo.n_state);
H(1,paramInfo.idx_map) = Phi_front;
H(1,paramInfo.flight(flightNumber).idx_ori_front) = z;
H(1,paramInfo.flight(flightNumber).idx_g) = 1;
H(2,paramInfo.idx_map) = Phi_back;
H(2,paramInfo.flight(flightNumber).idx_ori_back) = z;
H(2,paramInfo.flight(flightNumber).idx_back_bias) = 1;
H(2,paramInfo.flight(flightNumber).idx_g) = 1;
end