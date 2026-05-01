clc; clear all; close all;

% extract 1D profile

% CData = [];
% [CData, MT] = Read_Data_MT('USA_ZVTF_All_sites_108Bath.dat', CData);
% 
% [rho, x, y, z] = ReadModel_iso('Merged_iso_model.dat');

load('model_data.mat',"z","y","x","rho","MT","CData");

out_dir = fullfile('data', 'model_1D');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

save(fullfile(out_dir, 'USA_Model_Data.mat'),"z","y","x","rho","MT","CData");

% load('merged_data.mat','dd','cd')
fid = fopen('sites_required.dat', 'r');
indx = textscan(fid,'%s');
for i = 1:length(indx{1})
   sites_names{i} = strtrim(char(indx{1}(i)));
end
fclose(fid);

code_list = char(CData.MT.Stat.code);
for j = 1:length(CData.MT.Stat.lon)
    trimmed = strtrim(code_list(j,:));
    for i = 1:length(sites_names)
        site = sites_names{i};
        nmin = min(length(trimmed), length(site));
        if nmin >= 4 && strncmpi(trimmed, site, nmin)
            fprintf('%s %s\n','Site selected',char(CData.MT.Stat.code(j,:)));
            site_x = CData.MT.cord.x(j);
            site_y = CData.MT.cord.y(j);

            % find the cell position where the site is sitting
            ix = find(x <= site_x);  ix = ix(end);
            iy = find(y <= site_y);  iy = iy(end);

            model_1D.z = z(1:end-1)/1000;
            model_1D.rho = log10(squeeze(rho(ix, iy, :)));

            save(fullfile(out_dir, [sites_names{i}, '_model_1D.mat']),'model_1D');
        end
    end
end