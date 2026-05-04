function [Slin, P, prefix, variant, ref] = cpc_load_data()
% CPC_LOAD_DATA  Locate processed/stat .mat files in cwd and load them.
%   Returns:
%     Slin    - linear-depth statistics struct (from *_Stat_info.mat)
%     P       - processed samples struct with fields:
%                 z_all, rho_all, nrms_all, like_all, ngrid, step,
%                 sigma, CData
%     prefix  - one of 'MT_TD_Chain' / 'DC_TD_Chain' / 'MT_DC_TD_Chain'
%     variant - 'MT', 'DC', or 'MT+DC' (for plot titles)
%     ref     - struct with fields .has_ref (logical) and, if true,
%                 .z_km, .rho_log10, .n_layers (true reference 1D model)

candidates = { ...
    'MT_TD_Chain',     'MT'; ...
    'DC_TD_Chain',     'DC'; ...
    'MT_DC_TD_Chain',  'MT+DC'};

prefix = ''; variant = '';
for i = 1:size(candidates,1)
    if exist([candidates{i,1} '_Processed.mat'], 'file')
        prefix  = candidates{i,1};
        variant = candidates{i,2};
        break;
    end
end
if isempty(prefix)
    error('cpc_load_data: no *_Processed.mat in %s. Run A_02 first.', pwd);
end

statfile = [prefix, '_Stat_info.mat'];
if ~exist(statfile,'file')
    error('cpc_load_data: %s not found. Run A_02 first.', statfile);
end

P = load([prefix, '_Processed.mat']);
L = load(statfile);
Slin = L.S;

% --- reference 1D model lookup (synthetic stations only) -----------------
ref = struct('has_ref', false);
[~, station] = fileparts(pwd);
proj_root = fileparts(fileparts(mfilename('fullpath')));
candidates_m1d = { ...
    fullfile(proj_root, 'data', [station, '_model_1D.mat']), ...
    fullfile(proj_root, 'data', 'model_1D', [station, '_model_1D.mat'])};
for ii = 1:numel(candidates_m1d)
    if exist(candidates_m1d{ii},'file')
        m = load(candidates_m1d{ii},'model_1D');
        ref.has_ref   = true;
        ref.z_km      = m.model_1D.z(:);     % layer-top depths, km
        ref.rho_log10 = m.model_1D.rho(:);   % log10 ohm-m per layer
        ref.n_layers  = numel(ref.rho_log10);
        break;
    end
end
end
