clc; clear all; close all;

% extract 1D profile

% CData = [];
% [CData, MT] = Read_Data_MT('USA_ZVTF_All_sites_108Bath.dat', CData);
% 
% [rho, x, y, z] = ReadModel_iso('Merged_iso_model.dat');

load('model_data.mat',"z","y","x","rho","MT","CData");

% load('merged_data.mat','dd','cd')
fid = fopen('sites_required.dat', 'r');
indx = textscan(fid,'%s');
for i = 1:length(indx{1})
   indxx(i,:) = strtrim(char(indx{1}(i))); 
end
fclose(fid);

code_list = char(CData.MT.Stat.code);
for j = 1:length(CData.MT.Stat.lon)
    for i = 1:size(indxx,1)
        if strcmpi(code_list(j,1:5), indxx(i,:))
            fprintf('%s %s\n','Site selected',char(CData.MT.Stat.code(j,:)));
            site_x = CData.MT.cord.x(j);
            site_y = CData.MT.cord.y(j);

            % find the cell position where the site is sitting
            indx = find(x <= site_x);  indx = indx(end);
            indy = find(y <= site_y);  indy = indy(end);
            
            model_1D.z = z(1:end-1)/1000;
            model_1D.rho = log10(squeeze(rho(indx, indy, :)));

            save([indxx(i,:), '_model_1D.mat'],'model_1D');
        end
    end
end