function D_07_multi_station_summary(res_root)
% D_07_MULTI_STATION_SUMMARY  Small-multiples grid of posterior PDFs.
%
% Reads Custom_Process_Chains/station_metadata.csv to know which stations
% to include and how to group them by tectonic_group, then loads each
% station's *_Stat_info.mat from res_root/<station>/ and lays out one
% posterior pcolor per station, sorted by group.
%
% Usage (called from run_custom_processing):
%   D_07_multi_station_summary(fullfile('output','results_custom'))

if nargin < 1 || isempty(res_root)
    base = fileparts(fileparts(mfilename('fullpath')));   % project root
    res_root = fullfile(base, 'output', 'results_custom');
end
if ~isfolder(res_root)
    % fall back to the project-root location if a relative path was given
    base = fileparts(fileparts(mfilename('fullpath')));
    alt = fullfile(base, res_root);
    if isfolder(alt); res_root = alt; end
end
fprintf('  D_07: scanning %s\n', res_root);

meta = cpc_read_metadata();
if isempty(meta)
    fprintf('  D_07: skipping (no station_metadata.csv).\n');
    return;
end

% sort by tectonic_group then station
[~, ord] = sortrows([string(meta.tectonic_group), string(meta.station)]);
meta = meta(ord, :);

% load each station's Stat_info
n = height(meta);
panels = struct('station',{}, 'group',{}, 'S',{}, 'ref',{});
for i = 1:n
    st = char(meta.station(i));
    folder = fullfile(res_root, st);
    fcandidates = dir(fullfile(folder, '*_Stat_info.mat'));
    if isempty(fcandidates)
        fprintf('  D_07: %s has no *_Stat_info.mat -- skipping.\n', st);
        continue;
    end
    Lf = load(fullfile(folder, fcandidates(1).name), 'S');
    here = pwd; cd(folder);
    [~, ~, ~, ~, ref] = cpc_load_data();   % uses station folder name
    cd(here);
    panels(end+1).station = st;                   %#ok<AGROW>
    panels(end).group     = char(meta.tectonic_group(i));
    panels(end).S         = Lf.S;
    panels(end).ref       = ref;
end
np = numel(panels);
if np == 0
    fprintf('  D_07: no usable stations -- aborting.\n');
    return;
end

% layout: columns capped at 4
ncols = min(4, np);
nrows = ceil(np / ncols);

fig = figure('Color','w','Position',[40 40 350*ncols 360*nrows]);
for i = 1:np
    p = panels(i);
    S = p.S;
    PDF = log10(S.posteriorPDF + 1e-10);
    sampled = PDF(S.posteriorPDF > 0);
    if isempty(sampled); vmax = 0; else; vmax = double(max(sampled)); end
    vmin = vmax - 2.5;

    ax = subplot(nrows, ncols, i); hold on; box on;
    h = pcolor(S.rhoPlot, S.zPlot/1000, PDF);
    set(h,'EdgeColor','none');
    colormap(ax, cpc_viridis(256));
    caxis(ax, [vmin vmax]);
    set(ax, 'YDir','reverse', 'Layer','top', ...
            'FontSize', 9, 'TickDir','out');
    xlim([S.rhoMin S.rhoMax]);
    ylim([0 double(S.zMax)/1000]);

    % overlay true model if synthetic
    if p.ref.has_ref
        [zz_m, rr] = cpc_true_stairs(p.ref, double(S.zMax));
        plot(rr, zz_m/1000, '--w', 'LineWidth', 2.0);
        plot(rr, zz_m/1000, '--k', 'LineWidth', 1.0);
    end
    title(sprintf('%s\n[%s]', p.station, p.group), ...
          'Interpreter','none', 'FontSize', 10);
    if mod(i-1, ncols) == 0; ylabel('Depth (km)'); end
    if i > (nrows-1)*ncols;  xlabel('log_{10}\rho (\Omega·m)'); end
end

sgtitle('Multi-station posterior summary (grouped by tectonic class)', ...
        'FontSize', 13);
cpc_force_light_theme(fig);

out = fullfile(res_root, 'D_07_multi_station_summary');
print(fig, out, '-djpeg', '-r350');
fprintf('  D_07: saved %s.jpg\n', out);
end
