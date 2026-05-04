function D_03_chain_convergence()
% D_03_CHAIN_CONVERGENCE  White-background MCMC convergence diagnostic.
% Mirrors Process_All_Chains/A_01_see_chain_convergence.m but produces a
% thesis-ready figure: each T=1 chain plotted in a distinct colour, three
% stacked panels (nRMS, log-likelihood, noise sigma) on a log-x axis,
% with the burn-in cutoff marked as a vertical dashed line.

% ---- locate chain files -------------------------------------------------
candidates = { ...
    'MT_TD_Chain_',     'MT'; ...
    'DC_TD_Chain_',     'DC'; ...
    'MT_DC_TD_Chain_',  'MT+DC'};

suffix = ''; variant = ''; inversion_flag = -1;
for i = 1:size(candidates,1)
    if exist([candidates{i,1} '001.mat'], 'file')
        suffix = candidates{i,1};
        variant = candidates{i,2};
        inversion_flag = i - 1;       % 0=MT, 1=DC, 2=MT+DC
        break;
    end
end
if isempty(suffix)
    error('D_03_chain_convergence: no chain .mat files visible in path.');
end

S0 = load([suffix '001.mat']);
CData = S0.CData;
nChains = CData.nChains;
nsteps  = CData.nsteps;
ns      = CData.nsamples;
total   = nsteps * ns;
step_discard = floor(nsteps * 0.5);
burnin_x = step_discard * ns;          % sample index where burn-in ends

% ---- preallocate per-chain traces ---------------------------------------
like_all = nan(nChains, total);
nrms_all = nan(nChains, total);
% sigma is per-component; collapse to one trace by taking the first comp.
sigma_all = nan(nChains, total);
keep_chain = false(nChains, 1);

for ic = 1:nChains
    fname = sprintf('%s%03d.mat', suffix, ic);
    if ~exist(fname, 'file'); continue; end
    Lf = load(fname, 'Samples_Chain', 'CData');
    if Lf.CData.temperature(ic) ~= 1
        continue;             % only T=1 chains contribute to posterior
    end
    keep_chain(ic) = true;
    SC = Lf.Samples_Chain;
    nstep_actual = numel(SC);
    for istep = 1:nstep_actual
        samp = SC(istep);
        npt  = numel(samp.like);
        i1 = (istep-1)*npt + 1;
        i2 = i1 + npt - 1;
        if i2 > total; i2 = total; end
        m = i2 - i1 + 1;
        like_all(ic, i1:i2) = samp.like(1:m);
        if size(samp.misfit,2) >= 1
            nrms_all(ic, i1:i2) = samp.misfit(1:m, 1)';
        end
        if isfield(samp,'sigma') && ~isempty(samp.sigma)
            sigma_all(ic, i1:i2) = samp.sigma(1:m, 1)';
        end
    end
end
chains_used = find(keep_chain);
if isempty(chains_used)
    error('D_03_chain_convergence: no T=1 chains found.');
end

% ---- figure -------------------------------------------------------------
fig = figure('Color','w','Position',[80 80 1000 820]);
cmap = lines(numel(chains_used));
xq = 1:total;

% (1) nRMS
ax1 = subplot(3,1,1); hold(ax1,'on'); box(ax1,'on');
for k = 1:numel(chains_used)
    ic = chains_used(k);
    plot(ax1, xq, nrms_all(ic,:), '-', 'Color', cmap(k,:), 'LineWidth', 0.6);
end
yline(ax1, 1, '--', 'target nRMS = 1', 'Color', [0.30 0.30 0.30], ...
      'LineWidth', 1, 'LabelHorizontalAlignment','left', 'FontSize', 9);
xline(ax1, burnin_x, ':', 'burn-in', 'Color', [0.55 0.10 0.10], ...
      'LineWidth', 1.2, 'LabelVerticalAlignment','top', 'FontSize', 9);
set(ax1, 'XScale','log', 'FontSize',11, 'TickDir','out');
xlim(ax1, [1 total]);
ylim(ax1, [0 max(5, min(10, prctile(nrms_all(:), 99)))]);
xlabel(ax1, '# samples'); ylabel(ax1, 'nRMS');
title(ax1, sprintf('Chain convergence (%s)', variant), 'FontSize', 12);

% (2) log-likelihood
ax2 = subplot(3,1,2); hold(ax2,'on'); box(ax2,'on');
for k = 1:numel(chains_used)
    ic = chains_used(k);
    plot(ax2, xq, like_all(ic,:), '-', 'Color', cmap(k,:), 'LineWidth', 0.6);
end
xline(ax2, burnin_x, ':', '', 'Color', [0.55 0.10 0.10], 'LineWidth', 1.2);
set(ax2, 'XScale','log', 'FontSize', 11, 'TickDir','out');
xlim(ax2, [1 total]);
xlabel(ax2, '# samples'); ylabel(ax2, 'log-likelihood');

% (3) sigma factor
ax3 = subplot(3,1,3); hold(ax3,'on'); box(ax3,'on');
any_sigma = false;
for k = 1:numel(chains_used)
    ic = chains_used(k);
    if any(~isnan(sigma_all(ic,:)))
        plot(ax3, xq, sigma_all(ic,:), '-', 'Color', cmap(k,:), 'LineWidth', 0.6);
        any_sigma = true;
    end
end
if ~any_sigma
    text(ax3, 0.5, 0.5, 'sigma trace not stored in chain output', ...
         'Units','normalized', 'HorizontalAlignment','center');
end
yline(ax3, 1, '--', '\sigma = 1', 'Color', [0.30 0.30 0.30], ...
      'LineWidth', 1, 'FontSize', 9);
xline(ax3, burnin_x, ':', '', 'Color', [0.55 0.10 0.10], 'LineWidth', 1.2);
set(ax3, 'XScale','log', 'FontSize', 11, 'TickDir','out');
xlim(ax3, [1 total]);
xlabel(ax3, '# samples'); ylabel(ax3, 'noise \sigma factor');

% ---- light theme + save -----------------------------------------------
cpc_force_light_theme(fig);

[~, station] = fileparts(pwd);
fname = sprintf('%s_%s_D03_chain_convergence', station, suffix(1:end-1));
print(fig, fname, '-djpeg', '-r400');
end
