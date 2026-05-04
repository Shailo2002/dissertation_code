function D_06_layers_trace()
% D_06_LAYERS_TRACE  Trans-dimensional trace of the number of subsurface
% layers per chain across the full MCMC, plus the post-burn-in posterior
% histogram on the same figure. The trace is the standard convergence
% diagnostic for RJ-MCMC dimension jumps.

candidates = { ...
    'MT_TD_Chain_',     'MT'; ...
    'DC_TD_Chain_',     'DC'; ...
    'MT_DC_TD_Chain_',  'MT+DC'};
suffix = ''; variant = '';
for i = 1:size(candidates,1)
    if exist([candidates{i,1} '001.mat'], 'file')
        suffix = candidates{i,1}; variant = candidates{i,2}; break;
    end
end
if isempty(suffix)
    error('D_06_layers_trace: no chain .mat files visible in path.');
end

S0 = load([suffix '001.mat']);
nC = S0.CData.nChains;
nS = S0.CData.nsteps;
ns = S0.CData.nsamples;
total = nS * ns;
burnin_x = floor(nS * 0.5) * ns;

ng_trace = nan(nC, total);
keep_chain = false(nC, 1);

for ic = 1:nC
    fname = sprintf('%s%03d.mat', suffix, ic);
    if ~exist(fname, 'file'); continue; end
    Lf = load(fname, 'Samples_Chain', 'CData');
    if Lf.CData.temperature(ic) ~= 1; continue; end
    keep_chain(ic) = true;
    SC = Lf.Samples_Chain;
    for is = 1:numel(SC)
        nc = SC(is).ncells(:)';
        i1 = (is-1)*numel(nc) + 1;
        i2 = i1 + numel(nc) - 1;
        if i2 > total; i2 = total; end
        ng_trace(ic, i1:i2) = nc(1:i2-i1+1);
    end
end
chains_used = find(keep_chain);
if isempty(chains_used)
    error('D_06_layers_trace: no T=1 chains found.');
end

% post-burnin processed sample (already loaded by cpc_load_data, but we
% reuse the file lookup rather than the full helper to avoid re-loading
% the heavy *_Stat_info struct just for ngrid).
proc_files = dir('*_Processed.mat');
P = load(proc_files(1).name, 'ngrid');
ng_post = P.ngrid;

% --- figure -------------------------------------------------------------
fig = figure('Color','w','Position',[80 80 1200 700]);
cmap = lines(numel(chains_used));

% (a) trace
ax1 = subplot(1,2,1); hold on; box on;
for j = 1:numel(chains_used)
    ic = chains_used(j);
    plot(ax1, 1:total, ng_trace(ic,:), '-', ...
         'Color', cmap(j,:), 'LineWidth', 0.5);
end
xline(ax1, burnin_x, ':', 'burn-in', 'Color', [0.55 0.10 0.10], ...
      'LineWidth', 1.3, 'LabelVerticalAlignment','top', 'FontSize', 9);
set(ax1, 'XScale','log', 'FontSize', 11, 'TickDir','out');
xlim(ax1, [1 total]);
xlabel(ax1, '# samples');
ylabel(ax1, '# subsurface layers');
title(ax1, 'Trans-dimensional trace', 'FontSize', 12);
leg_t = arrayfun(@(c) sprintf('chain %d', c), chains_used(:)', ...
                 'UniformOutput', false);
legend(ax1, leg_t, 'Location','best', 'FontSize', 9, 'Box','off');

% (b) post-burnin histogram
ax2 = subplot(1,2,2); box on;
edges = (min(ng_post)-0.5):1:(max(ng_post)+0.5);
histogram(ng_post, edges, 'Normalization','probability', ...
          'FaceColor', [0.275 0.510 0.706], 'EdgeColor','w');
xlabel(ax2, '# subsurface layers');
ylabel(ax2, 'probability');
title(ax2, 'Post-burn-in posterior', 'FontSize', 12);
set(ax2, 'FontSize', 11, 'TickDir','out');

sgtitle(sprintf('Number of layers (%s)', variant), 'FontSize', 13);
cpc_force_light_theme(fig);

[~, station] = fileparts(pwd);
print(fig, sprintf('%s_%s_D06_layers_trace', station, suffix(1:end-1)), ...
      '-djpeg', '-r400');
end
