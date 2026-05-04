function D_02_synthetic_compare()
% D_02_SYNTHETIC_COMPARE  True-vs-Inversion comparison plot for synthetic
% stations. Shows the posterior mean resistivity profile, the 5th-95th
% percentile band as a shaded corridor, the true 1D model as blue
% stairs, and the model-domain RMS misfit between mean and truth.
%
% Skipped (with a printed note) when no reference 1D model is present.

[S, ~, prefix, variant, ref] = cpc_load_data();

if ~ref.has_ref
    fprintf('  D_02: no reference 1D model -- skipping synthetic compare.\n');
    return;
end

% ---- per-depth statistics ----------------------------------------------
zm  = S.zPlot(:);            % metres
zkm = zm / 1000;             % km
nZ  = size(S.rhoSamples,1);

p05  = nan(nZ,1);
p95  = nan(nZ,1);
pmu  = nan(nZ,1);
for i = 1:nZ
    row = S.rhoSamples(i,:); row = row(~isnan(row));
    if isempty(row); continue; end
    pmu(i) = mean(row);
    q = prctile(row, [5 95]);
    p05(i) = q(1); p95(i) = q(2);
end

% ---- true model on the same depth grid ---------------------------------
true_rho = nan(nZ,1);
zlay_km  = ref.z_km(:);
rlay     = ref.rho_log10(:);
for i = 1:nZ
    idx = find(zlay_km <= zkm(i), 1, 'last');
    if isempty(idx); idx = 1; end
    if idx > numel(rlay); idx = numel(rlay); end
    true_rho(i) = rlay(idx);
end

% ---- model-domain RMS (log10 ohm-m) ------------------------------------
ok  = ~isnan(pmu) & ~isnan(true_rho);
rms = sqrt(mean((pmu(ok) - true_rho(ok)).^2));

% ---- true-model stair points -------------------------------------------
[zz_true, rr_true] = cpc_true_stairs(ref, zkm(end));   % km because zmax<1000

% ---- figure ------------------------------------------------------------
fig = figure('Color','w','Position',[80 80 720 820]);
ax  = axes('NextPlot','add','YDir','reverse','Box','on', ...
           'FontSize',11,'XGrid','on','YGrid','on', ...
           'GridAlpha',0.3,'TickDir','out');

% Shaded P5-P95 corridor
ok_b = ~isnan(p05) & ~isnan(p95);
zb   = zkm(ok_b);
hband = fill([p05(ok_b); flipud(p95(ok_b))], [zb; flipud(zb)], ...
             [0.95 0.55 0.55], 'EdgeColor','none', 'FaceAlpha', 0.55);

% Inversion mean
hmean = plot(pmu, zkm, '-', 'Color', [0.85 0.10 0.10], 'LineWidth', 1.8);

% True model
htrue = plot(rr_true, zz_true, '-', 'Color', [0.10 0.20 0.85], 'LineWidth', 2.0);

xlim([S.rhoMin S.rhoMax]);
ylim([zkm(1) zkm(end)]);
xlabel('log_{10}(Resistivity) (ohm-m)', 'FontSize', 12);
ylabel('Depth (km)', 'FontSize', 12);

[~, station] = fileparts(pwd);
title(sprintf('Synthetic Test: True vs Inversion (%s)', variant), ...
      'FontSize', 13);

legend([hmean hband htrue], ...
       {'Inversion mean','5th-95th %ile','True model'}, ...
       'Location','northeast', 'FontSize', 10);

% RMS annotation box (bottom-left of the figure)
text(0.03, 0.03, sprintf('RMS = %.3f', rms), ...
     'Units','normalized', ...
     'FontSize', 11, ...
     'BackgroundColor','w', 'EdgeColor', [0.4 0.4 0.4], ...
     'Margin', 4, ...
     'VerticalAlignment','bottom', ...
     'HorizontalAlignment','left');

% ---- force light theme so saved JPEG isn't dark-mode --------------------
cpc_force_light_theme(fig);

fname = sprintf('%s_%s_D02_synthetic_compare', station, prefix);
print(fig, fname, '-djpeg', '-r400');
end
