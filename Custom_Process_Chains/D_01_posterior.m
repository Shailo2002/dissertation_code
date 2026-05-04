function D_01_posterior()
% D_01_POSTERIOR  3-panel posterior summary, MATLAB port of the Python
% reference postprocess/plot_posterior.py used by the Python rewrite.
%
% Panels:
%   1. Posterior PDF in (depth x log10 rho) space; viridis pcolor with
%      colour range fixed to [vmax-2.5, vmax] so unsampled cells fall to
%      the darkest viridis tone instead of leaving white holes;
%      red 5th/95th-percentile step lines; optional true model drawn as
%      a white-under-black double-stroke staircase.
%   2. Interface probability density vs depth; steel-blue posterior,
%      black dashed prior (vertical), red dashed true-interface depths.
%   3. Histogram of posterior on number of layers (steel-blue bars with
%      white edges) plus a red dashed line at the true layer count.

[S, P, prefix, variant, ref] = cpc_load_data();

% Synthetic stations are the only ones with a defensible "truth" -- for
% real-data stations any model_1D file is just a deterministic guess and
% must not be drawn as "True" markers on the interface / layer panels.
[~, station] = fileparts(pwd);
is_synth = startsWith(lower(station), 'synth');

% ---- prep posterior PDF ------------------------------------------------
% Mirror the Python: log10(pdf + 1e-10) so empty cells clip to vmin.
PDF     = S.posteriorPDF;
pdf_log = log10(PDF + 1e-10);
sampled = pdf_log(PDF > 0);
if isempty(sampled)
    vmax = 0; vmin = -2.5;
else
    vmax = double(max(sampled));
    vmin = vmax - 2.5;
end

zm   = S.zPlot(:);              % depth in metres
zmax = double(S.zMax);
ytk  = 0:50000:zmax;

% ---- figure ------------------------------------------------------------
fig = figure('Color','w','Position',[60 60 1500 740]);

% steel-blue (matplotlib's "steelblue" #4682B4)
SB   = [0.2745, 0.5098, 0.7059];
% red used for percentile lines and true interfaces
RED  = [0.85,   0.10,   0.10];

% ---- (1) Posterior PDF -------------------------------------------------
ax1 = subplot('Position',[0.05 0.10 0.27 0.78]); hold on; box on;
h = pcolor(S.rhoPlot, zm, pdf_log);
set(h, 'EdgeColor', 'none');
colormap(ax1, cpc_viridis(256));
caxis(ax1, [vmin vmax]);
set(ax1, 'YDir','reverse', 'Layer','top', 'FontSize',11, 'TickDir','out', ...
         'YTick', ytk);
try; ax1.YAxis.Exponent = 0; catch; end

% 5th / 95th percentile (step ~ MATLAB stairs(...,where='mid'))
hperc = stairs(S.p5,  zm, '-', 'Color', RED, 'LineWidth', 1.5);
        stairs(S.p95, zm, '-', 'Color', RED, 'LineWidth', 1.5);

% True model: white outline (thick) + black dashes on top -- matches python
htrue = [];
if ref.has_ref
    [zz_m, rr] = cpc_true_stairs(ref, zmax);   % zmax > 1000 -> metres
    plot(rr, zz_m, '--w', 'LineWidth', 2.5);
    htrue = plot(rr, zz_m, '--k', 'LineWidth', 1.5);
end

xlim([S.rhoMin S.rhoMax]);
ylim([0 zmax]);
xlabel('log_{10}(\rho) (ohm-m)', 'FontSize', 11);
ylabel('Depth (m)', 'FontSize', 11);
title('Posterior PDF', 'FontSize', 11);

c = colorbar(ax1, 'eastoutside');
c.Label.String   = 'log_{10} (PDF)';
c.Label.FontSize = 10;

leg_h = hperc; leg_t = {'5th/95th %ile'};
if ~isempty(htrue); leg_h(end+1) = htrue; leg_t{end+1} = 'True model'; end
legend(leg_h, leg_t, 'Location','southeast', 'FontSize', 9);

% ---- (2) Interface probability ----------------------------------------
ax2 = subplot('Position',[0.39 0.10 0.27 0.78]); hold on; box on;
kPDF   = sum(S.kSamples, 2, 'omitnan') ./ (sum(S.kSamples(:),'omitnan') * S.dz);
kPrior = mean(kPDF, 'omitnan');
kPDF(1) = 0;       % suppress surface interface (always present)
zk_m = zm(1:numel(kPDF));

hpost  = plot(kPDF, zk_m, '-', 'Color', SB,  'LineWidth', 2);
hprior = xline(kPrior, '--k', 'LineWidth', 1.5);

% True interface depths (horizontal lines) -- synthetic stations only
htrue_if = [];
if is_synth && ref.has_ref && numel(ref.z_km) > 1
    interfaces_m = ref.z_km(2:end) * 1000;   % skip surface
    for ii = 1:numel(interfaces_m)
        h = yline(interfaces_m(ii), '--', '', 'Color', RED, 'LineWidth', 1.5);
        if isempty(htrue_if); htrue_if = h; end
    end
end

set(ax2, 'YDir','reverse', 'FontSize', 11, 'TickDir','out', 'YTick', ytk, ...
         'YGrid','on', 'GridAlpha', 0.3);
try; ax2.YAxis.Exponent = 0; catch; end
ylim([0 zmax]);
xlabel('Probability density', 'FontSize', 11);
ylabel('Depth (m)', 'FontSize', 11);
title('Interface Probability', 'FontSize', 11);

leg_h = [hpost hprior]; leg_t = {'Posterior','Prior'};
if ~isempty(htrue_if); leg_h(end+1) = htrue_if; leg_t{end+1} = 'True interface'; end
legend(leg_h, leg_t, 'Location','southeast', 'FontSize', 9);

% ---- (3) Number of layers ---------------------------------------------
ax3 = subplot('Position',[0.72 0.10 0.25 0.78]); hold on; box on;
ng = P.ngrid;
edges = (min(ng)-0.5):1:(max(ng)+0.5);
histogram(ng, edges, 'Normalization','pdf', ...
          'FaceColor', SB, 'EdgeColor', 'w', 'LineWidth', 0.5);

if is_synth && ref.has_ref
    yl = ylim;
    htrue_n = xline(ref.n_layers, '--', sprintf('True = %d layers', ref.n_layers), ...
        'Color', RED, 'LineWidth', 2, ...
        'LabelHorizontalAlignment','right', ...
        'LabelVerticalAlignment','top', ...
        'FontSize', 9);
    ylim(yl);
    legend(htrue_n, {sprintf('True = %d layers', ref.n_layers)}, ...
           'Location','northeast', 'FontSize', 9);
end
set(ax3, 'FontSize', 11, 'TickDir','out', 'YGrid','on', 'GridAlpha', 0.3);
xlabel('# subsurface layers', 'FontSize', 11);
ylabel('Probability density', 'FontSize', 11);
title('Number of Layers', 'FontSize', 11);

% ---- super-title ------------------------------------------------------
sgt = sprintf('1D %s Transdimensional Bayesian Inversion\nPosterior Distribution', variant);
sgtitle(sgt, 'FontSize', 13);

% ---- light theme + save -----------------------------------------------
cpc_force_light_theme(fig);

fname = sprintf('%s_%s_D01_posterior', station, prefix);
print(fig, fname, '-djpeg', '-r400');
end
