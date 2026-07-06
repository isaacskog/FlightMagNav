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