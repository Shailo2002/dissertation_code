function [zz, rr] = cpc_true_stairs(ref, zmax)
% CPC_TRUE_STAIRS  Build (depth, log10-rho) stair points for plotting a
% reference 1D model on a posterior depth-resistivity plot.
%
%   ref   : struct from cpc_load_data with fields z_km, rho_log10
%   zmax  : maximum plot depth in the SAME units as the returned zz
%   zz    : depth values (input zmax units; convert ref.z_km*1000 first
%           if you want metres)
%
% Returns column vectors. The caller is expected to have already chosen
% the unit for zmax; this function multiplies ref.z_km by 1000 if zmax
% is greater than 1000 (heuristic for metres) -- otherwise keeps km.

if zmax > 1000
    zlayer = ref.z_km(:) * 1000;   % metres
else
    zlayer = ref.z_km(:);          % km
end
r = ref.rho_log10(:);

zz = []; rr = [];
for i = 1:numel(r)
    z_top = zlayer(i);
    if i < numel(r)
        z_bot = zlayer(i+1);
    else
        z_bot = zmax;
    end
    zz = [zz; z_top; z_bot];   %#ok<AGROW>
    rr = [rr; r(i); r(i)];     %#ok<AGROW>
end
end
