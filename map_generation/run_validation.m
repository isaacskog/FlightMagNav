function filter_par=run_validation(obs,model)

% Extract important variables 
theta     = model.theta;
P         = model.P;
basis     = model.basis;
settings  = model.settings;

paramInfo = build_parameter_info_validation();

% Number of observations
N=size(obs.y,1);

% Stack all observations
y_stacked=obs.y(:);

% Build regression matrices
[H_map,H_cal]=get_mesaurement_matrix_validation(obs.r_ned,obs.q,basis,paramInfo,settings);


% Prewhitning
theta_map = theta(model.paramInfo.idx_map);
Pmap = P(model.paramInfo.idx_map,model.paramInfo.idx_map);
Sigma  =settings.noise.sigma_validation^2*eye(2*N) + H_map*Pmap*H_map';
L  = chol(Sigma,'lower');
y_tilde = L\(y_stacked-H_map*theta_map);

% Get initial states
[x,P]=get_initial_states(paramInfo,settings);


% Structure of output data beyond the map parameters
filter_par = struct('NIS',[],'g',[],'g_var',[],'bias',[],'bias_var',[],'ori_bias_front',[],'ori_bias_front_cov_diag',[],'ori_bias_back',[],'ori_bias_back_cov_diag',[]);

% Allocate memory and store initial values
filter_par.NIS = zeros(N,1);
filter_par.g = zeros(N,1);
filter_par.g_var = zeros(N,1);
filter_par.bias = zeros(N,1);
filter_par.bias_var = zeros(N,1);
filter_par.ori_bias = zeros(N,3);
filter_par.ori_bias_cov_diag = zeros(N,3);
filter_par.g(1) = x(paramInfo.idx_g);
filter_par.g_var(1) = P(paramInfo.idx_g,paramInfo.idx_g);
filter_par.bias(1) = x(paramInfo.idx_back_bias);
filter_par.bias_var(1) = P(paramInfo.idx_back_bias,paramInfo.idx_back_bias);
filter_par.ori_bias_front(1,:) = x(paramInfo.idx_ori_front).';
filter_par.ori_bias_front_cov_diag(1,:) = diag(P(paramInfo.idx_ori_front,paramInfo.idx_ori_front)).';
filter_par.ori_bias_back(1,:) = x(paramInfo.idx_ori_back).';
filter_par.ori_bias_back_cov_diag(1,:) = diag(P(paramInfo.idx_ori_back,paramInfo.idx_ori_back)).';

% Run Kalman filter
for nn=2:N

    % Build measurement and measurement model 
    y=[y_tilde(nn) y_tilde(nn+N)]';
    ssm.H=[H_cal(nn,:); H_cal(nn+N,:)];
    ssm.R=eye(2);
    
    % Get the process noise variance (only temporal variations)
    ssm.Q=zeros(paramInfo.nstate,paramInfo.nstate);
    dt = obs.t_sec(nn)-obs.t_sec(nn-1);
    ssm.q = dt^2*settings.time.sigma_q^2;
    ssm.q_idx=paramInfo.idx_g;

    % Do one step of the Kalman filter algorithm
    [x,P,~,NIS] = step_kf(y,x,P,ssm);  


    % Store data
    filter_par.NIS(nn) = NIS;
    filter_par.g(nn) = x(paramInfo.idx_g);
    filter_par.g_var(nn) = P(paramInfo.idx_g,paramInfo.idx_g);
    filter_par.bias(nn) = x(paramInfo.idx_back_bias);
    filter_par.bias_var(nn) = P(paramInfo.idx_back_bias,paramInfo.idx_back_bias);
    filter_par.ori_bias_front(nn,:) = x(paramInfo.idx_ori_front).';
    filter_par.ori_bias_front_cov_diag(nn,:) = diag(P(paramInfo.idx_ori_front,paramInfo.idx_ori_front)).';
    filter_par.ori_bias_back(nn,:) = x(paramInfo.idx_ori_back).';
    filter_par.ori_bias_back_cov_diag(nn,:) = diag(P(paramInfo.idx_ori_back,paramInfo.idx_ori_back)).';
end

end



function paramInfo = build_parameter_info_validation()
%BUILD_PARAMETER_INFO Define full parameter/state layout.
%
% State:
%   x = [
%       theta_ori_front_1
%       theta_ori_back_1
%       b_back
%       g]
%

paramInfo = struct();
paramInfo.idx_ori_front = 1:3;
paramInfo.idx_ori_back = 4:6;
paramInfo.idx_back_bias =7;
paramInfo.idx_g = 8;
paramInfo.nstate=8;

end


% Get initial states
function [x,P]=get_initial_states(paramInfo,settings)

x=zeros(paramInfo.nstate,1);
P=zeros(paramInfo.nstate,paramInfo.nstate);
P(paramInfo.idx_ori_front,paramInfo.idx_ori_front)=settings.calibration.sigma_ori^2*eye(3);
P(paramInfo.idx_ori_back,paramInfo.idx_ori_back)=settings.calibration.sigma_ori^2*eye(3);
P(paramInfo.idx_back_bias,paramInfo.idx_back_bias)=settings.calibration.sigma_back_bias^2;
P(paramInfo.idx_g,paramInfo.idx_g)=settings.time.sigma_g0^2;
end



function [Hmap,Hcal]=get_mesaurement_matrix_validation(r_ned,q,basis,paramInfo,settings)
% Create two-row measurement matrix.
%
% Measurement:
%   y_stacked =obs.y(:)
%

% Allocate memory
N=size(r_ned,1);
M=size(basis.map.centers,1);
L=8;
Hmap=zeros(2*N,M);
Hcal=zeros(2*N,L);


for n=1:N

    % Get regressors
    [Phi_front,Phi_back,z] = get_measurement_regressors(r_ned(n,:),q(n,:),basis,settings);

    % Build the matrix
    Hmap(n,:)=Phi_front;
    Hmap(N+n,:)=Phi_back;
    Hcal(n,paramInfo.idx_ori_front)=z;      % Front orientation bias
    Hcal(n,paramInfo.idx_g)=1;              % Time random walk
    Hcal(N+n,paramInfo.idx_ori_back)=z;     % Back orientation bias
    Hcal(N+n,paramInfo.idx_back_bias)=1;    % Back bias    
    Hcal(N+n,paramInfo.idx_g)=1;            % Time random walk
end

end