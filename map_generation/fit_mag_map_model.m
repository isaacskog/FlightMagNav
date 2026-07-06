function model = fit_mag_map_model(data,settings)
%FIT_MAG_MAP_MODEL Estimate a magnetic-field map using a Kalman filter.
%
% The state contains:
%   - one common static map,
%   - one calibration block per flight,
%   - one temporal random-walk state per flight.


% Extract the observations need to build the map
obs = build_observation_data(data,settings);

% Build the basis
basis = build_map_basis(obs,settings);

% Construct information about parameter structure of the state vector
paramInfo = build_parameter_info(settings,basis);

% Check if hyperparamters should be tuned
if(settings.tune_hyper_par)

    % Set up optimization algorithm
    options = optimoptions('fminunc', ...
        'Algorithm','quasi-newton', ...
        'Display','iter', ...
        'MaxIterations',5, ...
        'MaxFunctionEvaluations',500, ...
        'UseParallel',true, ...
        'PlotFcn','optimplotfval');

    % Get initial conditions
    theta0 = get_init(settings);

    % Run optimization
    theta_hat = fminunc( ...
        @(theta) run_kalman_filter(obs,basis,paramInfo,complete_settings(settings,theta)), ...
        theta0, ...
        options);

    % Store the optimized parameters in the settings
    settings = complete_settings(settings,theta_hat);
end

% Run the map learning 
[NlogL,theta,P,filterInfo] = run_kalman_filter(obs,basis,paramInfo,settings); %#ok<ASGLU>

% Save the result into a struct
model = struct();
model.theta = theta;
model.P = P;
model.paramInfo = paramInfo;
model.basis = basis;
model.settings = settings;
model.obs = obs;
model.filterInfo = filterInfo;
model.predict_map = @(r_ned) predict_map(r_ned,theta,P,basis,paramInfo);
model.predict_map_lla = @(lat,lon,alt) ...
    predict_map(lla_to_local_ned(lat,lon,alt,settings.reference_lla), ...
    theta,P,basis,paramInfo);
end


function theta0 = get_init(settings)
%Initial log-hyperparameter vector theta.
    theta0 = log(settings.noise.sigma);
end

function settings = complete_settings(settings,theta)
% Update settings with the values in theta 
        settings.noise.sigma=exp(theta);
end


function obs = build_observation_data(data,settings)
%BUILD_OBSERVATION_DATA Build low-rate front/back magnetometer observations.
%
% obs(ii) contains one flight. obs(ii).y has two columns:
%   column 1: front magnetometer total field
%   column 2: back magnetometer total field

    obs = repmat(struct( ...
        'flight_idx',[], ...
        'time',[], ...
        't_sec',[], ...
        'y',[], ...
        'r_ned',[], ...
        'q',[], ...
        'yaw',[]), ...
        1,numel(settings.idx));

    % Design low pass filter    
    fsIn = 100;
    fc = 0.45*settings.fs_map;
    filterOrder = 3;
    [b,a] = butter(filterOrder,fc/(fsIn/2),'low');

    % Extract data
    for ii = 1:numel(settings.idx)

        kk = settings.idx(ii);
        ins = data(kk).GPSaidedINS;

        tIn = seconds(ins.time - ins.time(1));
        tOut = (0:1/settings.fs_map:tIn(end)).';

        % Low pass filter using zero-phase filtering 
        yFront = filtfilt(b,a,data(kk).front_mag.tot_field);
        yBack  = filtfilt(b,a,data(kk).back_mag.tot_field);

        yFront = interp1(tIn,yFront,tOut,'nearest');
        yBack  = interp1(tIn,yBack, tOut,'nearest');

        rBody = lla_to_local_ned( ...
            interp1(tIn,ins.latitude, tOut,'nearest'), ...
            interp1(tIn,ins.longitude,tOut,'nearest'), ...
            interp1(tIn,ins.altitude, tOut,'nearest'), ...
            settings.reference_lla);

        q = interp1(tIn,ins.quaternion,tOut,'nearest');

        % Get the index of the data points above the specified minimum
        % height
        idx = -rBody(:,3) >= settings.altitude_min;
 
        % STore the data
        obs(ii).flight_idx = kk;
        obs(ii).time = ins.time(1) + seconds(tOut(idx));
        obs(ii).t_sec = tOut(idx);
        obs(ii).y = [yFront(idx) yBack(idx)]-mean([yFront(idx); yBack(idx)]);
        obs(ii).r_ned = rBody(idx,:);
        obs(ii).q = q(idx,:);
      
    end
end


function basis = build_map_basis(obs,settings)
%BUILD_MAP_BASIS Build rectangular RBF grid enclosing all measured points.

    allR = [];
    for ii = 1:numel(obs)
        allR = [allR; obs(ii).r_ned]; %#ok<AGROW>
    end

    north = allR(:,1);
    east  = allR(:,2);

    h = settings.map.center_spacing;
    margin = settings.map.margin;

    nMin = floor((min(north)-margin)/h)*h;
    nMax = ceil( (max(north)+margin)/h)*h;
    eMin = floor((min(east)-margin)/h)*h;
    eMax = ceil( (max(east)+margin)/h)*h;

    [N,E] = meshgrid(nMin:h:nMax,eMin:h:eMax);

    basis = struct();
    basis.map.centers = [N(:),E(:)];
    basis.map.length_scale = settings.map.length_scale;
end


function paramInfo = build_parameter_info(settings,basis)
%BUILD_PARAMETER_INFO Define full parameter/state layout.
%
% State:
%   x = [
%       theta_map
%       theta_ori_front_1
%       theta_ori_back_1
%       b_back_1
%       ...
%       theta_ori_front_M
%       theta_ori_back_M
%       b_back_M
%       g_1
%       ...
%       g_M
%   ]
%
% The map and calibration states are static. The g_i states are random
% walks, but only the active flight's g_i receives process noise.

    nMap = size(basis.map.centers,1);
    nFlights = numel(settings.idx);

    paramInfo = struct();
    paramInfo.idx_map = 1:nMap;

    next = nMap + 1;

    for ii = 1:nFlights
        paramInfo.flight(ii).flight_idx = settings.idx(ii);

        paramInfo.flight(ii).idx_ori_front = next:(next+2);
        next = next + 3;

        paramInfo.flight(ii).idx_ori_back = next:(next+2);
        next = next + 3;

        paramInfo.flight(ii).idx_back_bias = next;
        next = next + 1;
    end

    for ii = 1:nFlights
        paramInfo.flight(ii).idx_g = next;
        next = next + 1;
    end

    paramInfo.n_static = next - 1;
    paramInfo.n_state = next - 1;
end


function [NlogL_tot,x,P,filter_par] = run_kalman_filter(obs,basis,paramInfo,settings)
%RUN_KALMAN_FILTER Estimate map and flight-specific calibration states.

    x = zeros(paramInfo.n_state,1);
    P = zeros(paramInfo.n_state,paramInfo.n_state);
    P(paramInfo.idx_map,paramInfo.idx_map) = ...
    settings.map.sigma^2 *eye(numel(paramInfo.idx_map));

    % Flight-specific calibration priors.
    for ii = 1:numel(obs)
        P(paramInfo.flight(ii).idx_ori_front,paramInfo.flight(ii).idx_ori_front) = ...
            settings.calibration.sigma_ori^2*eye(3);

        P(paramInfo.flight(ii).idx_ori_back,paramInfo.flight(ii).idx_ori_back) = ...
            settings.calibration.sigma_ori^2*eye(3);

        P(paramInfo.flight(ii).idx_back_bias,paramInfo.flight(ii).idx_back_bias) = ...
            settings.calibration.sigma_back_bias^2;
        P(paramInfo.flight(ii).idx_g,paramInfo.flight(ii).idx_g) = settings.time.sigma_g0^2;

        if ii == 1
            % First flight defines the reference level.
            P(paramInfo.flight(ii).idx_g,paramInfo.flight(ii).idx_g) = 0;
        end
    end

    % Total negative log likelihood      
    NlogL_tot = 0;

    % Structure of output data beyond the map parameters
    filter_par = repmat(struct('NIS',[],'g',[],'g_var',[],'bias',[],'bias_var',[],'ori_bias_front',[],'ori_bias_front_cov_diag',[],'ori_bias_back',[],'ori_bias_back_cov_diag',[]),1,numel(obs));

    for ii = 1:numel(obs)

        % Allocate memory and store initial values
        filter_par(ii).NIS = zeros(size(obs(ii).y,1),1);
        filter_par(ii).g = zeros(size(obs(ii).y,1),1);
        filter_par(ii).g_var = zeros(size(obs(ii).y,1),1);
        filter_par(ii).bias = zeros(size(obs(ii).y,1),1);
        filter_par(ii).bias_var = zeros(size(obs(ii).y,1),1);
        filter_par(ii).ori_bias = zeros(size(obs(ii).y,1),3);
        filter_par(ii).ori_bias_cov_diag = zeros(size(obs(ii).y,1),3);
        filter_par(ii).g(1) = x(paramInfo.flight(ii).idx_g);
        filter_par(ii).g_var(1) = P(paramInfo.flight(ii).idx_g,paramInfo.flight(ii).idx_g);
        filter_par(ii).bias(1) = x(paramInfo.flight(ii).idx_back_bias);
        filter_par(ii).bias_var(1) = P(paramInfo.flight(ii).idx_back_bias,paramInfo.flight(ii).idx_back_bias);
        filter_par(ii).ori_bias_front(1,:) = x(paramInfo.flight(ii).idx_ori_front).';
        filter_par(ii).ori_bias_front_cov_diag(1,:) = diag(P(paramInfo.flight(ii).idx_ori_front,paramInfo.flight(ii).idx_ori_front)).';
        filter_par(ii).ori_bias_back(1,:) = x(paramInfo.flight(ii).idx_ori_back).';
        filter_par(ii).ori_bias_back_cov_diag(1,:) = diag(P(paramInfo.flight(ii).idx_ori_back,paramInfo.flight(ii).idx_ori_back)).';


        % Set the measurement covariance for the specific flight
        ssm.R = settings.noise.sigma(ii)^2*eye(2);

        % Run the Kalman filter 
        for kk = 2:numel(obs(ii).t_sec)

            % Construct the measurement matrix    
            ssm.H = get_mesaurement_matrix( ...
                obs(ii).r_ned(kk,:), ...
                obs(ii).q(kk,:), ...
                ii, ...
                basis, ...
                paramInfo, ...
                settings);

            % Get the process noise variance (only temporal variations)
            dt = obs(ii).t_sec(kk)-obs(ii).t_sec(kk-1);
            ssm.q = dt^2*settings.time.sigma_q^2;
            ssm.q_idx=paramInfo.flight(ii).idx_g;
            
            % Do one step of the Kalman filter algorithm
            [x,P,logL,NIS] = step_kf(obs(ii).y(kk,:).',x,P,ssm);

            % Accumulate the negative log likelihood
            NlogL_tot = NlogL_tot - logL;

            % Store data
            filter_par(ii).NIS(kk) = NIS;
            filter_par(ii).g(kk) = x(paramInfo.flight(ii).idx_g);
            filter_par(ii).g_var(kk) = P(paramInfo.flight(ii).idx_g,paramInfo.flight(ii).idx_g);
            filter_par(ii).bias(kk) = x(paramInfo.flight(ii).idx_back_bias);
            filter_par(ii).bias_var(kk) = P(paramInfo.flight(ii).idx_back_bias,paramInfo.flight(ii).idx_back_bias);
            filter_par(ii).ori_bias_front(kk,:) = x(paramInfo.flight(ii).idx_ori_front).';
            filter_par(ii).ori_bias_front_cov_diag(kk,:) = diag(P(paramInfo.flight(ii).idx_ori_front,paramInfo.flight(ii).idx_ori_front)).';
            filter_par(ii).ori_bias_back(kk,:) = x(paramInfo.flight(ii).idx_ori_back).';
            filter_par(ii).ori_bias_back_cov_diag(kk,:) = diag(P(paramInfo.flight(ii).idx_ori_back,paramInfo.flight(ii).idx_ori_back)).';
        end
    end    
end


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


function r = lla_to_local_ned(lat,lon,alt,refLLA)
% Small-area WGS84 to local NED conversion.

    R = 6378137;
    lat0 = refLLA(1);
    lon0 = refLLA(2);
    alt0 = refLLA(3);
    north = deg2rad(lat - lat0).*R;
    east  = deg2rad(lon - lon0).*R.*cosd(lat0);
    down  = -(alt - alt0);
    r = [north(:),east(:),down(:)];
end


function Phi = spatial_rbf_2d(r,centers,lengthScale)
% Isotropic Gaussian RBFs in the horizontal NED plane.

    D2 = pdist2_local(r,centers).^2;
    Phi = exp(-0.5*D2/lengthScale^2);
end


function D = pdist2_local(A,B)
% Pairwise Euclidean distances.

    AA = sum(A.^2,2);
    BB = sum(B.^2,2).';

    D2 = AA + BB - 2*(A*B.');
    D2 = max(D2,0);

    D = sqrt(D2);
end


function [x,P,logL,NIS] = step_kf(y,x,P,ssm)

    % Time update of state covariance
    P(ssm.q_idx,ssm.q_idx) = P(ssm.q_idx,ssm.q_idx)+ ssm.q;
    
    
    
    % Innovation covariance
    S = ssm.H*P*ssm.H' + ssm.R;
    
    % Innovation 
    z = y - ssm.H*x;

    % Normalized innovation squared 
    NIS = z'*(S\z);

    % Kalman gain
    K = (P*ssm.H')/S;

    % Filter upadte
    x = x + K*z;
    P = (eye(size(P))-K*ssm.H)*P;

    % Log of marganalized likelihood
    logL =logmvnpdf(z,0,S);
end


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
