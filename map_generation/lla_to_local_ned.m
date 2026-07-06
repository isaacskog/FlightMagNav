function r = lla_to_local_ned(lat,lon,alt,refLLA)
%LLA_TO_LOCAL_NED Small-area WGS84 to local NED conversion.
%
% Input:
%   lat, lon [deg], alt [m]
%   refLLA = [lat0 lon0 alt0]
%
% Output:
%   r = [north east down] [m]
%
% This approximation is sufficient for small survey areas and does not
% require Mapping Toolbox.

    R = 6378137;

    lat0 = refLLA(1);
    lon0 = refLLA(2);
    alt0 = refLLA(3);

    north = deg2rad(lat - lat0).*R;
    east  = deg2rad(lon - lon0).*R.*cosd(lat0);
    down  = -(alt - alt0);

    r = [north(:),east(:),down(:)];
end
