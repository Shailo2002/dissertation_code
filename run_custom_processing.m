function run_custom_processing()
% run_custom_processing
% ------------------------------------------------------------------
% Driver for the dissertation-specific post-processing pipeline.
% Produces a distinct set of figures (D_01..D_05) under
%   <inv_root_parent>/results_custom/<station>/
%
% NOTE: this is a *function* (not a script) on purpose -- A_02 and Z_03
% start with `clear`, which would otherwise wipe the loop variables.
% ------------------------------------------------------------------
clc; close all;

base       = fileparts(mfilename('fullpath'));
proc_std   = fullfile(base, 'Process_All_Chains');
proc_cust  = fullfile(base, 'Custom_Process_Chains');
fwd_dir    = fullfile(base, 'MT_DC_function');     % needed for MT1DFW (D_04)
addpath(proc_std);
addpath(proc_cust);
if exist(fwd_dir,'dir'); addpath(fwd_dir); end

% ---- which inversion outputs to process -------------------------------
inv_root = fullfile(base, 'output', 'inversions');
res_root = fullfile(base, 'output', 'results_custom');

if ~exist(inv_root,'dir')
    error('Inversion folder not found: %s', inv_root);
end
if ~exist(res_root,'dir'); mkdir(res_root); end

station_dirs = dir(inv_root);
station_dirs = station_dirs([station_dirs.isdir]);
station_dirs = station_dirs(~ismember({station_dirs.name},{'.','..'}));
if isempty(station_dirs)
    error('No station folders in %s', inv_root);
end

custom_steps = { ...
    'D_01_posterior',          ...
    'D_02_synthetic_compare',  ...
    'D_03_chain_convergence',  ...
    'D_04_data_fit',           ...
    'D_05_acceptance_rate',    ...
    'D_06_layers_trace'};

for k = 1:numel(station_dirs)
    station    = station_dirs(k).name;
    inv_folder = fullfile(inv_root, station);
    res_folder = fullfile(res_root, station);
    if ~exist(res_folder,'dir'); mkdir(res_folder); end

    fprintf('\n===== Custom processing: %s =====\n', station);
    addpath(inv_folder);
    cd(res_folder);
    close all;

    % --- ensure A_02 has produced *_Processed.mat / *_Stat_info.mat ----
    have_proc = ~isempty(dir('*_Processed.mat'));
    try
        if ~have_proc
            fprintf('  running A_02_process_chains...\n');
            run_isolated_script('A_02_process_chains', res_folder);
        end
    catch ME
        fprintf('  prerequisite step failed for %s: %s\n', station, ME.message);
        rmpath(inv_folder);
        continue;
    end
    cd(res_folder);
    close all;

    % Defensive cleanup: if Z_03 had been run earlier and dropped a
    % MT_LogPost_image.jpg in this folder, remove it -- we don't want it.
    stale_imgs = {'MT_LogPost_image.jpg','DC_LogPost_image.jpg', ...
                  'MT_DC_LogPost_image.jpg'};
    for ii = 1:numel(stale_imgs)
        f = fullfile(res_folder, stale_imgs{ii});
        if exist(f,'file'); delete(f); end
    end

    % --- run each custom plot
    for s = 1:numel(custom_steps)
        sname = custom_steps{s};
        try
            cd(res_folder);
            fprintf('  [%s]\n', sname);
            feval(sname);
            close all;
        catch ME
            fprintf('  %s failed: %s\n', sname, ME.message);
        end
    end

    rmpath(inv_folder);
end

% ---- cross-station figures (D_07 multi-station summary, D_08 map) -----
cd(base);
try
    fprintf('\n[D_07] multi-station posterior summary...\n');
    D_07_multi_station_summary(res_root);
catch ME
    fprintf('  D_07 failed: %s\n', ME.message);
end
try
    fprintf('\n[D_08] station map at 100 km...\n');
    D_08_station_map(res_root, 100);
catch ME
    fprintf('  D_08 failed: %s\n', ME.message);
end

cd(base);
fprintf('\nCustom processing complete.\n');
fprintf('Figures: %s\n', res_root);
end

% ----------------------------------------------------------------------
function run_isolated_script(scriptname, work_folder)
% Execute a script in this helper's local workspace so any `clear`
% inside the script (A_02 / Z_03 both do `clc; clear;`) cannot reach
% the driver's loop variables.
cd(work_folder);
run(scriptname);
end
