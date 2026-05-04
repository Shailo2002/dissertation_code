clc;    clear all;      close all;

% Add Path of the Folders and all the subfolders having all the subroutines
% restoredefaultpath
subroutinefolder1 = 'MT_DC_function';
addpath(genpath([pwd,'/',subroutinefolder1]));

CData = GetDefaultParameters();

% load('merged_data.mat','dd','cd')
fid = fopen('sites_required.dat', 'r');
indx = textscan(fid,'%s');
for i = 1:length(indx{1})
    sites_names{i} = strtrim(char(indx{1}(i)));
end
fclose(fid);

inv_root = fullfile(pwd, 'output', 'inversions');
if ~exist(inv_root, 'dir'); mkdir(inv_root); end

% Open one combined acceptance-rate summary file for all stations
AR_filename = fullfile(inv_root, 'Acceptance_Rate_Summary_AllStations.txt');
fid_ar = fopen(AR_filename, 'w');
fprintf(fid_ar, '%s\n', repmat('=',1,70));
fprintf(fid_ar, '   ACCEPTANCE RATE SUMMARY (all stations, all chains & steps)\n');
fprintf(fid_ar, '%s\n', repmat('=',1,70));
fprintf(fid_ar, 'Generated: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

for i = 1:length(sites_names)
    filename = sites_names{i};

    destination_folder = fullfile(inv_root, filename);
    if ~exist(destination_folder, 'dir'); mkdir(destination_folder); end

    MT_datafile = fullfile('data', 'observed_data', [filename,'.dat']);
    DC_datafile = ' .dat';
    CData.LogFile = fullfile(destination_folder, [filename,'_LogFile.dat']);

    % which data we nned to invert
    CData.inversion_method = 'MT';   % MT DC MT_DC
    CData.MT.datatype = 'Z';

    [model_m0, model_LF0, data_sigma] = initialise_chain(CData);

    CData = read_data(MT_datafile, DC_datafile, CData);

    fprintf('%s \n','MCMC Procedure Starts ');
    fid = fopen(CData.LogFile,'w');
    swapCount = [];

    for is = 1:CData.nsteps
        for ic = 1:CData.nChains
            % Pick LF0 and model0 after parallel tempering
            model0 = model_m0{ic};  LF0 = model_LF0{ic}; noise = data_sigma{ic};

            [Samples, model0, LF0, noise] =  bayesain_hier(model0, LF0,...
                noise, CData, ic, is);

            % update the models and likehood
            model_m0{ic} = model0;  model_LF0{ic} = LF0; data_sigma{ic} = noise;

            % Let us store the last results only, we will keep updating the
            % results.
            results_all(ic,is) = Samples;
        end

        status = save_chain(results_all, is, CData, swapCount, destination_folder);

        % Begin parallel tempering we will flip the model and likelihood
        if any(CData.temperature > 1)
            [model_LF0, model_m0, data_sigma, swapCount] = swap_temperatures...
                (model_LF0, model_m0, data_sigma, CData, is, swapCount);
        end
    end
    fclose(fid);
    fprintf('%s %s\n','MCMC Procedure Ends for ',filename);

    % =====================================================================
    %  ACCEPTANCE RATE SUMMARY for this station -> append to combined file
    % =====================================================================
    all_AR    = zeros(CData.nChains, CData.nsteps);
    all_AR_BD = zeros(CData.nChains, CData.nsteps);
    all_AR_D  = zeros(CData.nChains, CData.nsteps);
    all_AR_M  = zeros(CData.nChains, CData.nsteps);
    all_AR_C  = zeros(CData.nChains, CData.nsteps);
    for ic = 1:CData.nChains
        for is = 1:CData.nsteps
            ar = results_all(ic, is).acceptance_all;
            all_AR(ic,is)    = ar(1);
            all_AR_BD(ic,is) = ar(2);
            all_AR_D(ic,is)  = ar(3);
            all_AR_M(ic,is)  = ar(4);
            all_AR_C(ic,is)  = ar(5);
        end
    end
    labels  = {'Overall', 'Birth  ', 'Death  ', 'Move   ', 'Change '};
    AR_mats = {all_AR, all_AR_BD, all_AR_D, all_AR_M, all_AR_C};

    fprintf(fid_ar, '%s\n', repmat('-',1,70));
    fprintf(fid_ar, 'Station: %s\n', filename);
    fprintf(fid_ar, '%s\n', repmat('-',1,70));
    fprintf(fid_ar, '%-10s  %8s  %8s  %8s  %8s\n', ...
        'Type','Min(%)','Max(%)','Mean(%)','Std(%)');
    fprintf(fid_ar, '%s\n', repmat('-',1,50));
    for k = 1:length(labels)
        v = AR_mats{k}(:);
        fprintf(fid_ar, '%-10s  %8.2f  %8.2f  %8.2f  %8.2f\n', ...
            labels{k}, min(v), max(v), mean(v), std(v));
    end
    fprintf(fid_ar, '\n%-10s  %8s  %8s  %8s\n', 'Chain','Min(%)','Max(%)','Mean(%)');
    fprintf(fid_ar, '%s\n', repmat('-',1,38));
    for ic = 1:CData.nChains
        v = all_AR(ic,:);
        fprintf(fid_ar, 'Chain %-4d  %8.2f  %8.2f  %8.2f\n', ...
            ic, min(v), max(v), mean(v));
    end
    fprintf(fid_ar, '\n');

    pause(2)
end
fclose(fid_ar);
fprintf('Combined acceptance rate summary saved to: %s\n', AR_filename);
function CData = GetDefaultParameters( )

% CData.temperature = [ones(1,5) logspace(0, log10(50), 5)];
% CData.temperature = ones(10,1);
CData.temperature = [1 1 1 1 1];

% number of chains we want to apply parallel tempering in one go.
CData.nchain_for_PT = min(length(find(CData.temperature==1)),...
    length(find(CData.temperature~=1)));

CData.nChains = length(CData.temperature);
%  The jump type between different temperatures:
%  0, randomly choose a temperature;
%  1, randomly choose a temperature between the two neighbours
%  2, choose the nearest temperature
CData.jumptype = 0;

% log domain for z or not
CData.logdomain = true;
CData.log_normal_noise = false;

CData.scale = 1;   % scale to convert model to km, keep it as 1

% Minimum and maximum value of nodes for creation of new layer
if CData.logdomain == true
    CData.min_z = log10(1000);           CData.max_z = log10(350000);
else
    % linear domain
    CData.min_z = 0;                   CData.max_z = 200;
end
CData.minumum_layer_thickness = 100;     % specify in mts;

CData.proposal = [0.05 0.4 0.7 1];

% No of steps and number of samples in each step
CData.nsteps = 100;       CData.nsamples = 100;

CData.resmin = 1.0E-01;    CData.resmax = 1.0E+05;

% updating birth rates for multiple chains
CData.nChains_atT1 = length(find(CData.temperature==1));  % chains at T == 1
CData.nChains_atoT = CData.nChains - CData.nChains_atT1;  % chains at T /= 1

CData.sigma_rho = [0.1*ones(1,CData.nChains_atT1) linspace(0.1, 0.2, CData.nChains_atoT)];
CData.sigma_rho_birth = [0.2*ones(1,CData.nChains_atT1) linspace(0.2, 0.3, CData.nChains_atoT)];
CData.sigma_rho_delayed = CData.sigma_rho/4;

% Position values range and deivations
CData.sigma_loc_z = [1*ones(1,CData.nChains_atT1) linspace(1, 2, CData.nChains_atoT)];
CData.sigma_loc_z_delayed = CData.sigma_loc_z/4;

% Noise values range and deivations
CData.sigma_Z = ([1 5]);
CData.sigma_app_res = ([0.8 5]);
CData.sigma_phase = ([0.8 5]);
CData.sigma_noise = 0.01;


% Kernel type; 0 for gaussian and 1 for prior;
CData.kernel = 0;            CData.eps = 1.0E-09;

% Minimum and the maximum number of nodes
CData.minnodes = 2;         CData.maxnodes = 30;

% ==================== No need to look further down ===================== %
CData.min_res_log = log10(CData.resmin);
CData.max_res_log = log10(CData.resmax);  % values on the log10 scale

CData.rho = 'log10';       % logrho
end