function D_04_data_fit()
% D_04_DATA_FIT  Posterior-predictive data-fit plot.
%
% Top panel : observed apparent resistivity vs period with error bars,
%             overlaid with the posterior median forward response and a
%             5th-95th percentile predictive band.
% Bot panel : same for phase (degrees).
%
% Posterior models are subsampled to keep the forward-modelling cost
% bounded (~300 calls). Only the MT branch is implemented -- DC stations
% are skipped with a printed note.

[~, P, prefix, variant, ~] = cpc_load_data();
CData = P.CData;

if ~isfield(CData,'MT') || ~isfield(CData.MT,'period') ...
        || ~isfield(CData.MT,'dobs_Z')
    fprintf('  D_04: no MT impedance data -- skipping data-fit plot.\n');
    return;
end

period   = CData.MT.period(:);
Zobs     = CData.MT.dobs_Z(:);
err_Z    = CData.MT.err_Z(:);

mu0   = 4*pi*1e-7;
omega = 2*pi ./ period;

% --- observed app_res / phase + propagated errors -----------------------
app_obs   = abs(Zobs).^2 ./ (mu0 * omega);
phase_obs = atan2(imag(Zobs), real(Zobs)) * 180/pi;
err_app   = 2 * app_obs .* err_Z ./ abs(Zobs);          % linear ohm-m
err_phase = err_Z ./ abs(Zobs) * 180/pi;                % degrees, small-angle

% --- subsample posterior and compute predicted responses ----------------
nKeep = size(P.z_all, 1);
nSub  = min(300, nKeep);
idx   = round(linspace(1, nKeep, nSub));

nT = numel(period);
app_pred = nan(nSub, nT);
phs_pred = nan(nSub, nT);

for ii = 1:nSub
    k = idx(ii);
    n = P.ngrid(k);
    z = P.z_all(k, 1:n);
    rho_lin = 10 .^ P.rho_all(k, 1:n);
    if n > 1
        thick = diff(z);                 % m, top n-1 layers
    else
        thick = 1;                       % halfspace placeholder
    end
    try
        [~, app_pred(ii,:), phs_pred(ii,:)] = MT1DFW(rho_lin, thick, period);
    catch ME
        % Forward call failed for this sample (rare numerical edge case).
        fprintf('  D_04: forward call failed at sample %d: %s\n', k, ME.message);
    end
end

app_p05 = prctile(app_pred, 5,  1);
app_p50 = prctile(app_pred, 50, 1);
app_p95 = prctile(app_pred, 95, 1);
phs_p05 = prctile(phs_pred, 5,  1);
phs_p50 = prctile(phs_pred, 50, 1);
phs_p95 = prctile(phs_pred, 95, 1);

% --- figure -------------------------------------------------------------
fig = figure('Color','w','Position',[80 80 900 820]);

BAND = [0.95 0.55 0.55];
RED  = [0.85 0.10 0.10];
BLUE = [0.10 0.20 0.85];

% (a) apparent resistivity
ax1 = subplot(2,1,1); hold on; box on;
fill([period; flipud(period)], [app_p05'; flipud(app_p95')], ...
     BAND, 'EdgeColor','none','FaceAlpha', 0.5);
hmed = plot(period, app_p50, '-', 'Color', RED, 'LineWidth', 1.8);
hobs = errorbar(period, app_obs, err_app, 'o', ...
                'Color', BLUE, 'MarkerFaceColor', BLUE, ...
                'MarkerSize', 5, 'LineWidth', 0.8);
set(ax1, 'XScale','log', 'YScale','log', 'FontSize', 11, ...
         'TickDir','out', 'XGrid','on','YGrid','on','GridAlpha', 0.3);
xlabel('Period (s)', 'FontSize', 12);
ylabel('Apparent resistivity (\Omega·m)', 'FontSize', 12);
title(sprintf('Data fit (%s)', variant), 'FontSize', 12);
legend([hobs hmed], {'observed \pm err','posterior median'}, ...
       'Location','best', 'FontSize', 10);

% (b) phase
ax2 = subplot(2,1,2); hold on; box on;
fill([period; flipud(period)], [phs_p05'; flipud(phs_p95')], ...
     BAND, 'EdgeColor','none','FaceAlpha', 0.5);
plot(period, phs_p50, '-', 'Color', RED, 'LineWidth', 1.8);
errorbar(period, phase_obs, err_phase, 'o', ...
         'Color', BLUE, 'MarkerFaceColor', BLUE, ...
         'MarkerSize', 5, 'LineWidth', 0.8);
set(ax2, 'XScale','log', 'FontSize', 11, 'TickDir','out', ...
         'XGrid','on','YGrid','on','GridAlpha', 0.3);
xlabel('Period (s)', 'FontSize', 12);
ylabel('Phase (degrees)', 'FontSize', 12);
ylim([0 90]);

cpc_force_light_theme(fig);
[~, station] = fileparts(pwd);
print(fig, sprintf('%s_%s_D04_data_fit', station, prefix), '-djpeg', '-r400');
end
