function [H,y,flightId] = build_batch_map_data(obs,basis,paramInfo,settings)
%BUILD_BATCH_MAP_DATA Stack all map-learning measurements and regressors.
%
% Each consecutive pair of rows corresponds to the front and back sensors
% at one time instant, as in the per-flight matrices in (9)-(13).

samples = arrayfun(@(flight) size(flight.y,1),obs);
nRows = 2*sum(samples);
H = zeros(nRows,paramInfo.n_state);
y = zeros(nRows,1);
flightId = zeros(nRows,1);

row = 1;
for ii = 1:numel(obs)
    for kk = 1:samples(ii)
        rows = row:(row+1);
        H(rows,:) = get_mesaurement_matrix( ...
            obs(ii).r_ned(kk,:),obs(ii).q(kk,:), ...
            ii,basis,paramInfo,settings);
        y(rows) = obs(ii).y(kk,:).';
        flightId(rows) = ii;
        row = row + 2;
    end
end
end
