clc; clear all; close all;

% Generate synthetic 1D MT data of varying difficulty for Bayesian inversion.
% Output mirrors what A2_extract_data.m and A_001_extract_1D_profile.m
% produce for real SAMTEX stations:
%   data/observed_data/<name>.dat        Period Real(Z_SI) Imag(Z_SI) Error
%   data/model_1D/<name>_model_1D.mat    struct model_1D with .z (km) and
%                                        .rho (log10 ohm-m) per layer
% Station names are written to ../sites_required.dat.

base = fileparts(fileparts(mfilename('fullpath')));   % Bayesian_code/
addpath(genpath(fullfile(base,'MT_DC_function')));    % MT1DFW.m

obs_dir = fullfile(base,'data','observed_data');
mod_dir = fullfile(base,'data','model_1D');
if ~exist(obs_dir,'dir'); mkdir(obs_dir); end
if ~exist(mod_dir,'dir'); mkdir(mod_dir); end

mu0 = 4*pi*1e-7;

% Periods: log-spaced 0.1 s -> 30000 s, broad enough to be sensitive
% from sediments (~hundreds of m) down to ~350 km in resistive lithosphere
% (skin depth in 100 Ohm-m at 30000 s ~ 870 km).
T  = logspace(log10(0.1), log10(30000), 30)';
Nf = length(T);

% ------------------------------------------------------------------
%  Define 4 synthetic models of increasing difficulty.
%  Inversion is run with z up to 350 km, so each model places its
%  deepest interface inside the 0-350 km range (halfspace continues
%  below). Thinking in geological terms:
%    upper crust / lower crust / lithospheric mantle / asthenosphere.
%
%   Each model:  layer_top (m, top of each layer; first must be 0)
%                layer_rho (Ohm-m, one per layer; last layer = halfspace)
%                noise_pct (Gaussian noise, % of |Z|)
% ------------------------------------------------------------------
models = struct();

% Model 1 (easy): 2-layer, resistive lithosphere over conductive asthenosphere.
% Single deep interface at 200 km -- well within the 350 km inversion window.
models(1).name      = 'synth01_easy_2layer';
models(1).layer_top = [0; 200000];
models(1).layer_rho = [300; 10];
models(1).noise_pct = 2;

% Model 2 (moderate): 3-layer R-C-R
% resistive crust / mid-crustal conductor / resistive lithospheric mantle,
% with the deepest interface inside the 350 km range.
models(2).name      = 'synth02_mod_3layer_RCR';
models(2).layer_top = [0; 20000; 150000];
models(2).layer_rho = [500; 20; 2000];
models(2).noise_pct = 3;

% Model 3 (hard): 4-layer
% sediments / resistive crust / resistive lithospheric mantle / conductive
% asthenosphere; deepest interface near 300 km.
models(3).name      = 'synth03_hard_4layer';
models(3).layer_top = [0; 5000; 40000; 300000];
models(3).layer_rho = [100; 1000; 5000; 30];
models(3).noise_pct = 5;

% Model 4 (very hard): 5-layer with a thin lower-crustal conductor and
% deep asthenospheric conductor at ~250 km; the thin layer is a known
% resolution challenge for MT inversion.
models(4).name      = 'synth04_vhard_5layer_thin';
models(4).layer_top = [0; 8000; 25000; 30000; 250000];
models(4).layer_rho = [200; 2000; 5; 3000; 50];
models(4).noise_pct = 5;

% Reproducible noise
rng(42);

new_sites = {};
for k = 1:numel(models)
    M = models(k);

    thicknesses = [M.layer_top(2:end) - M.layer_top(1:end-1); 1.5*M.layer_top(end)];

    [Z, appres, phase] = MT1DFW(M.layer_rho, thicknesses, T);

    % Add Gaussian noise (real & imag independently), 5% error floor like A2
    noise_frac = 0.01 * M.noise_pct;
    Z_noisy = Z + noise_frac * abs(Z) .* (randn(Nf,1) + 1i*randn(Nf,1))/sqrt(2);
    err     = 0.05 * abs(Z_noisy);

    % --- Observed-data file (same format as A2_extract_data.m output) ---
    fname = fullfile(obs_dir, [M.name, '.dat']);
    fid   = fopen(fname, 'w');
    for jf = 1:Nf
        fprintf(fid, '%10.5f %+10.5E %+10.5E %+10.5E\n', ...
            T(jf), real(Z_noisy(jf)), imag(Z_noisy(jf)), err(jf));
    end
    fclose(fid);
    fprintf('Wrote %s  (%d periods, %d-layer)\n', fname, Nf, numel(M.layer_rho));

    new_sites{end+1} = M.name; %#ok<SAGROW>

    % --- True 1D model file (same format as A_001_extract_1D_profile.m) ---
    %   model_1D.z   : depths in km, top of each layer
    %   model_1D.rho : log10(resistivity), one value per layer
    model_1D = struct();
    model_1D.z   = M.layer_top / 1000;          % m -> km
    model_1D.rho = log10(M.layer_rho);          % ohm-m -> log10
    save(fullfile(mod_dir, [M.name, '_model_1D.mat']), 'model_1D');
    fprintf('Wrote %s\n', fullfile(mod_dir, [M.name, '_model_1D.mat']));

    % --- Preview plot of forward response ---
    figure('Visible','off');
    subplot(2,1,1);
    loglog(T, appres, '-ko', 'LineWidth', 1.2); grid on;
    xlabel('Period (s)'); ylabel('App. res. (\Omega m)');
    title(sprintf('%s', M.name), 'Interpreter','none');
    subplot(2,1,2);
    semilogx(T, phase, '-ko', 'LineWidth', 1.2); grid on;
    xlabel('Period (s)'); ylabel('Phase (deg)');
    print(gcf, fullfile(obs_dir, [M.name, '_preview']), '-djpeg', '-r150');
    close(gcf);
end

% ------------------------------------------------------------------
%  Overwrite sites_required.dat with the synthetic site names only
% ------------------------------------------------------------------
sites_file = fullfile(base, 'sites_required.dat');
fid = fopen(sites_file,'w');
for k = 1:numel(new_sites)
    fprintf(fid, '%s\n', new_sites{k});
end
fclose(fid);
fprintf('Wrote %d synthetic site names to %s\n', numel(new_sites), sites_file);

fprintf('\nDone. %d synthetic stations ready for inversion.\n', numel(new_sites));
