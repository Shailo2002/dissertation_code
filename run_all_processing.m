clc; close all;

base = fileparts(mfilename('fullpath'));
proc = fullfile(base, 'Process_All_Chains');
addpath(proc);

inv_root = fullfile(base, 'output', 'inversions');
res_root = fullfile(base, 'output', 'results');
if ~exist(res_root, 'dir'); mkdir(res_root); end

if ~exist(inv_root, 'dir')
    error('Inversion folder not found: %s\nRun Bayesian_1D_MTDC_J_looped.m first.', inv_root);
end

station_dirs = dir(inv_root);
station_dirs = station_dirs([station_dirs.isdir]);
station_dirs = station_dirs(~ismember({station_dirs.name}, {'.','..'}));

if isempty(station_dirs)
    error('No station folders found in %s', inv_root);
end

for k = 1:length(station_dirs)
    station    = station_dirs(k).name;
    inv_folder = fullfile(inv_root, station);
    res_folder = fullfile(res_root, station);
    if ~exist(res_folder, 'dir'); mkdir(res_folder); end

    fprintf('\n===== Processing station: %s =====\n', station);

    % Make chain .mat files visible via the MATLAB path so scripts that
    % `load('MT_TD_Chain_001.mat')` resolve them, while cwd stays at the
    % result folder so any `save(...)` lands there.
    addpath(inv_folder);
    cd(res_folder);
    close all;

    % --- A_01 : chain convergence
    try
        cd(res_folder);
        fprintf('  [A_01] cwd = %s\n', pwd);
        run_script('A_01_see_chain_convergence');
        cd(res_folder);
        annotate_and_save(gcf, station, 'convergence', res_folder);
    catch ME
        fprintf('  A_01 failed: %s\n', ME.message);
    end
    close all;

    % --- A_02 : process chains (no figure produced)
    try
        cd(res_folder);
        fprintf('  [A_02] cwd = %s\n', pwd);
        run_script('A_02_process_chains');
        cd(res_folder);
    catch ME
        fprintf('  A_02 failed: %s\n', ME.message);
        rmpath(inv_folder);
        continue;
    end

    % --- A_03 : posterior plots (linear depth)
    try
        cd(res_folder);
        fprintf('  [A_03] cwd = %s\n', pwd);
        run_script('A_03_plot_posteriori_3_things');
        cd(res_folder);
        annotate_and_save(gcf, station, 'posterior', res_folder);
    catch ME
        fprintf('  A_03 failed: %s\n', ME.message);
    end
    close all;

    % --- A_04 : noise hyperparameter diagnostic
    try
        cd(res_folder);
        fprintf('  [A_04] cwd = %s\n', pwd);
        run_script('A_04_plot_noise');
        cd(res_folder);
        annotate_and_save(gcf, station, 'noise', res_folder);
    catch ME
        fprintf('  A_04 failed: %s\n', ME.message);
    end
    close all;

    % --- Z_03 : posterior plots on log-spaced depth axis
    try
        cd(res_folder);
        fprintf('  [Z_03] cwd = %s\n', pwd);
        run_script('Z_03_plot_posteriori_3_things');
        cd(res_folder);
        annotate_and_save(gcf, station, 'logposterior', res_folder);
    catch ME
        fprintf('  Z_03 failed: %s\n', ME.message);
    end
    close all;

    % --- AB_01_cross_plots_3_layer : depth-band resistivity cross-plots
    try
        cd(res_folder);
        fprintf('  [AB_01] cwd = %s\n', pwd);
        run_script('AB_01_cross_plots_3_layer');
        cd(res_folder);
        annotate_and_save(gcf, station, 'crossplot', res_folder);
    catch ME
        fprintf('  AB_01_cross_plots_3_layer failed: %s\n', ME.message);
    end
    close all;

    % Remove duplicate generic-named JPEGs written by the individual
    % processing scripts (we already saved station-prefixed copies).
    generic_jpegs = { ...
        'MT_TD_Chain.jpg', 'DC_TD_Chain.jpg', 'MT_DC_TD_Chain.jpg', ...
        'MT_Post_image.jpg', 'DC_Post_image.jpg', ...
        'MT_LogPost_image.jpg', 'DC_LogPost_image.jpg', 'MT_DC_LogPost_image.jpg', ...
        'MT_noise_image.jpg', 'DC_noise_image.jpg', ...
        'Cross_Plot_3layer.jpg'};
    for gi = 1:numel(generic_jpegs)
        gp = fullfile(res_folder, generic_jpegs{gi});
        if exist(gp, 'file'); delete(gp); end
    end

    rmpath(inv_folder);
end

cd(base);
fprintf('\nAll done.\n');
fprintf('Inversion outputs: %s\n', inv_root);
fprintf('Processed results: %s\n', res_root);

% ----------------------------------------------------------------------
function run_script(scriptname)
    eval(scriptname);
end

function annotate_and_save(figh, station, tag, res_folder)
    if ~isgraphics(figh); return; end
    figure(figh);
    try
        sgtitle(station, 'Interpreter', 'none', 'FontWeight', 'bold');
    catch
    end
    fname = sprintf('%s_%s', station, tag);
    print(figh, fullfile(res_folder, fname), '-djpeg', '-r300');
end
