function D_08_station_map(res_root, depth_km)
% D_08_STATION_MAP
% Map of stations on MATLAB's built-in geobasemap (no shapefiles
% required), coloured by posterior mean log10(rho) at depth_km. Each
% station's marker shape encodes its tectonic_group from
% station_metadata.csv; station names are labelled next to the marker.

% ---- defaults ----------------------------------------------------------
if nargin < 1 || isempty(res_root)
    base = fileparts(fileparts(mfilename('fullpath')));
    res_root = fullfile(base, 'output', 'results_custom');
end
if nargin < 2 || isempty(depth_km); depth_km = 100; end
fprintf('  D_08: scanning %s (depth %g km)\n', res_root, depth_km);

ZOOM_EXTENT = [10 40 -35 -12]; % [lon_min lon_max lat_min lat_max]

% ---- read metadata + posterior values --------------------------------
meta = cpc_read_metadata();
if isempty(meta)
    fprintf('  D_08: no station_metadata.csv found.\n');
    return;
end

n = height(meta);
lon=nan(n,1); lat=nan(n,1); val=nan(n,1);
names=strings(n,1); group=strings(n,1); keep=false(n,1);

for i = 1:n
    st     = char(meta.station(i));
    folder = fullfile(res_root, st);
    fc = dir(fullfile(folder, '*_Stat_info.mat'));
    if isempty(fc); continue; end
    Lf = load(fullfile(folder, fc(1).name), 'S');
    S  = Lf.S;
    zkm = S.zPlot(:) / 1000;
    [~, ix] = min(abs(zkm - depth_km));
    if ~isnan(S.pmean(ix))
        val(i)   = S.pmean(ix);
        lon(i)   = double(meta.lon(i));
        lat(i)   = double(meta.lat(i));
        names(i) = meta.station(i);
        group(i) = meta.tectonic_group(i);
        keep(i)  = true;
    end
end
lon=lon(keep); lat=lat(keep); val=val(keep);
names=names(keep); group=group(keep);
if isempty(lon)
    fprintf('  D_08: no valid stations at %g km.\n', depth_km);
    return;
end

% ---- figure with built-in basemap ------------------------------------
fig = figure('Color','w','Position',[60 60 1150 1000]);
gx  = geoaxes(fig);
hold(gx,'on');
geobasemap(gx, 'grayland');
geolimits(gx, [ZOOM_EXTENT(3) ZOOM_EXTENT(4)], ...
              [ZOOM_EXTENT(1) ZOOM_EXTENT(2)]);

% ---- per-group geoscatter (shape encodes tectonic class) -------------
% Use the same global colour scale across groups so the colorbar makes
% sense; clip extremes to the actual val range.
clim_lo = min(val); clim_hi = max(val);
if abs(clim_hi - clim_lo) < 1e-6
    clim_lo = clim_lo - 0.5; clim_hi = clim_hi + 0.5;
end

unique_groups = unique(group);
markers = {'o','s','^','d','v','p','h','>','<','*'};
hg     = gobjects(0);
labels = {};
for g = 1:numel(unique_groups)
    grp = unique_groups(g);
    gi  = (group == grp);
    mk  = markers{mod(g-1, numel(markers)) + 1};
    h = geoscatter(gx, lat(gi), lon(gi), 160, val(gi), 'filled', ...
                   'Marker', mk, ...
                   'MarkerEdgeColor', 'k', 'LineWidth', 0.9);
    hg(end+1)     = h;                                      %#ok<AGROW>
    labels{end+1} = char(grp);                              %#ok<AGROW>
end

% ---- station-name labels (dark, bold, with white halo) ---------------
for i = 1:numel(lon)
    text(gx, lat(i)+0.18, lon(i)+0.18, char(names(i)), ...
        'FontSize', 9, 'Color', 'k', 'FontWeight','bold', ...
        'BackgroundColor','w', 'Margin', 1, ...
        'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 0.3, ...
        'Interpreter','none', ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle');
end

% ---- colormap + colorbar (high rho = bright -- cratons "glow") -------
colormap(gx, cpc_viridis(256));
clim(gx, [clim_lo clim_hi]);
cb = colorbar(gx);
cb.Label.String   = sprintf('log_{10} \\rho at %g km (\\Omega·m)', depth_km);
cb.Label.FontSize = 12;
cb.FontSize       = 10;
cb.Color          = 'k';

% ---- legend (marker shape per tectonic_group) ------------------------
lgd = legend(hg, labels, ...
             'Location','southwest', ...
             'FontSize', 9, 'Interpreter','none', ...
             'Color','w', 'EdgeColor', [0.5 0.5 0.5], ...
             'TextColor', 'k');
title(lgd, 'Tectonic unit', 'FontSize', 10, 'Color','k');

% ---- title -----------------------------------------------------------
title(gx, sprintf(['Posterior mean log_{10}\\rho at %g km' ...
                   '\nMcCourt (2013) tectonic units'], depth_km), ...
      'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');

% Force readable axis text on geoaxes (latitude/longitude tick labels)
gx.FontSize         = 11;
gx.LatitudeAxis.Color  = 'k';
gx.LongitudeAxis.Color = 'k';
gx.GridColor           = [0.4 0.4 0.4];
gx.GridAlpha           = 0.6;

% ---- export ----------------------------------------------------------
out = fullfile(res_root, sprintf('D_08_station_map_%dkm', round(depth_km)));
exportgraphics(fig, [out '.jpg'], 'Resolution', 350);
fprintf('  D_08: saved %s.jpg\n', out);
end
