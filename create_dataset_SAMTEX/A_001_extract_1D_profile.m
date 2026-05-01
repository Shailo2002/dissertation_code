clc; clear all; close all;

% extract 1D profile

CData = [];
[CData, MT] = Read_Data_MT('Samtex_ObservedData_Max_out.dat', CData);

[rho, x, y, z] = ReadModel_iso('ModEM_SA2022_Max.rho');

out_dir = fullfile('data', 'model_1D');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

save(fullfile(out_dir, 'SAMTEX_Max_Data.mat'),"z","y","x","rho","MT","CData");

% load('merged_data.mat','dd','cd')
fid = fopen('sites_required.dat', 'r');
indx = textscan(fid,'%s');
for i = 1:length(indx{1})
   sites_names{i} = strtrim(char(indx{1}(i))); 
end
fclose(fid);

code_list = char(CData.MT.Stat.code);
for j = 1:length(CData.MT.Stat.lon)
    for i = 1:length(sites_names)
        if strcmpi(strtrim(code_list(j,:)), sites_names{i})
            fprintf('%s %s\n','Site selected',char(CData.MT.Stat.code(j,:)));
            site_x = CData.MT.cord.x(j);
            site_y = CData.MT.cord.y(j);

            % find the cell position where the site is sitting
            indx = find(x <= site_x);  indx = indx(end);
            indy = find(y <= site_y);  indy = indy(end);
            
            model_1D.z = z(1:end-1)/1000;
            model_1D.rho = log10(exp(squeeze(rho(indx, indy, :))));

            save(fullfile(out_dir, [sites_names{i}, '_model_1D.mat']),'model_1D');
        end
    end
end