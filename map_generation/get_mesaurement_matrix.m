function H = get_mesaurement_matrix(r_ned,q,flightNumber,basis,paramInfo,settings)
% Create two-row measurement matrix.
%
% Measurement:
%   y = [y_front; y_back]
%
% Only the parameter block for flightNumber is active.

% Get regressors
[Phi_front,Phi_back,z] = get_measurement_regressors(r_ned,q,basis,settings);

% Build the matrix
H = zeros(2,paramInfo.n_state);
H(1,paramInfo.idx_map) = Phi_front;
H(1,paramInfo.flight(flightNumber).idx_xi_front) = [z 1];
H(2,paramInfo.idx_map) = Phi_back;
H(2,paramInfo.flight(flightNumber).idx_xi_back) = [z 1];
end
