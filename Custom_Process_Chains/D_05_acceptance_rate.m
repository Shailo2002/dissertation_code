function D_05_acceptance_rate()
% D_05_ACCEPTANCE_RATE  Per-step acceptance-rate traces for each move
% type recorded by the RJ-MCMC sampler.
%
% Each Samples_Chain(istep).acceptance_all is a 5-vector with
% [Overall, Birth, Death, Move, Change-rho] acceptance percentages.
% Plot one panel per move type, all T=1 chains overlaid in colour, with
% a green band marking the empirically-recommended 23-44% sweet spot.

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
    error('D_05_acceptance_rate: no chain .mat files visible in path.');
end

S0 = load([suffix '001.mat']);
nC = S0.CData.nChains;
nS = S0.CData.nsteps;

% AR(:,:,k) holds rate for move type k=1..5
AR = nan(nC, nS, 5);
keep_chain = false(nC, 1);

for ic = 1:nC
    fname = sprintf('%s%03d.mat', suffix, ic);
    if ~exist(fname, 'file'); continue; end
    Lf = load(fname, 'Samples_Chain', 'CData');
    if Lf.CData.temperature(ic) ~= 1; continue; end
    keep_chain(ic) = true;
    SC = Lf.Samples_Chain;
    for is = 1:numel(SC)
        if isfield(SC(is),'acceptance_all') && numel(SC(is).acceptance_all) >= 5
            AR(ic, is, :) = SC(is).acceptance_all(1:5);
        end
    end
end
chains_used = find(keep_chain);
if isempty(chains_used)
    error('D_05_acceptance_rate: no T=1 chains found.');
end

labels = {'Overall', 'Birth', 'Death', 'Move', 'Change \rho'};
fig = figure('Color','w','Position',[80 80 1100 780]);
cmap = lines(numel(chains_used));

for k = 1:5
    ax = subplot(3,2,k); hold on; box on;
    % target band 23-44%
    xl = [1 nS];
    fill([xl fliplr(xl)], [23 23 44 44], [0.85 0.95 0.85], ...
         'EdgeColor','none','FaceAlpha',0.5);
    for j = 1:numel(chains_used)
        ic = chains_used(j);
        plot(1:nS, squeeze(AR(ic,:,k)), '-', ...
             'Color', cmap(j,:), 'LineWidth', 0.9);
    end
    set(ax, 'FontSize', 10, 'TickDir','out');
    title(labels{k}, 'FontSize', 11);
    xlabel('step'); ylabel('acceptance %');
    ylim([0 100]); xlim(xl);
end

% legend in 6th tile
ax_leg = subplot(3,2,6); axis(ax_leg,'off'); hold(ax_leg,'on');
legh = gobjects(numel(chains_used),1);
for j = 1:numel(chains_used)
    legh(j) = plot(ax_leg, NaN, NaN, '-', ...
                   'Color', cmap(j,:), 'LineWidth', 1.6);
end
ph = patch(ax_leg, NaN(1,4), NaN(1,4), [0.85 0.95 0.85], 'EdgeColor','none');
leg_labels = [arrayfun(@(c) sprintf('chain %d', c), chains_used(:)', ...
                       'UniformOutput', false), {'target 23-44%'}];
legend(ax_leg, [legh; ph], leg_labels, 'Location','west', 'Box','off', ...
       'FontSize', 10);

sgtitle(sprintf('RJ-MCMC acceptance rates (%s)', variant), 'FontSize', 13);
cpc_force_light_theme(fig);

[~, station] = fileparts(pwd);
print(fig, sprintf('%s_%s_D05_acceptance', station, suffix(1:end-1)), ...
      '-djpeg', '-r400');
end
