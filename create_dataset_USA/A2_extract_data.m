clc; clear all;

ref_lat = 39.7134567;      ref_long = -95.9637527;
% CData.rho = [];
% [CData, dd] = Read_Data_MT('SA_array_gr3s_108Bath_NAD83_moreSites.dat', CData);
load('USArray.mat','CData',"dd");

out_dir = fullfile('data', 'observed_data');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

save(fullfile(out_dir, 'USArray.mat'),'CData',"dd");

% load('merged_data.mat','dd','cd')
fid = fopen('sites_required.dat', 'r');
indx = textscan(fid,'%s');
for i = 1:length(indx{1})
   site_names{i} = strtrim(char(indx{1}(i)));
end
fclose(fid);

mu0 = 4*pi*1e-7;   % magnetic permeability
Nf  = length(CData.MT.T);

Zdet     = zeros(1, Nf);
Zdet_SI  = zeros(1, Nf);
rhoa     = zeros(1, Nf);
phase    = zeros(1, Nf);

code_list = char(CData.MT.Stat.code);
for j = 1:length(CData.MT.Stat.lon)
    flag = 0;
    trimmed = strtrim(code_list(j,:));
    for i = 1:length(site_names)
       site = site_names{i};
       nmin = min(length(trimmed), length(site));
       if nmin >= 4 && strncmpi(trimmed, site, nmin)
           fprintf('%s %s\n','Site selected',char(CData.MT.Stat.code(j,:)));

           data = squeeze(dd.data(j,:,:));
           data = conj(data);
           for jf = 1:Nf

               % Extract tensor components
               Zxx = data(1,jf);   Zxy = data(2,jf);
               Zyx = data(3,jf);   Zyy = data(4,jf);

               % ---- Step 1: determinant impedance ----
               Zdet(jf) = sqrt(Zxx*Zyy - Zxy*Zyx);

               % ---- Step 2: convert to SI (E/B) ----
               Zdet_SI(jf) = Zdet(jf)*mu0;   % mV/km/nT -> V/m/T

               % ---- Step 3: apparent resistivity (E/B form) ----
               omega = 2*pi/CData.MT.T(jf);
               rhoa(jf) = abs(Zdet_SI(jf)) * abs(Zdet_SI(jf)) / (mu0*omega);

               % ---- Step 4: phase ----
               phase(jf) = atan2(imag(Zdet_SI(jf)), real(Zdet_SI(jf))) * 180/pi;

           end

           % write the data
           fid = fopen(fullfile(out_dir, [char(site_names{i}), '.dat']),'w');
           for ifreq = 1:Nf
               if abs(Zdet_SI(ifreq)) == 0
                   continue
               end
               fprintf(fid, '%10.5f %+10.5E %+10.5E %+10.5E\n', CData.MT.T(ifreq), ...
                   real(Zdet_SI(ifreq)), imag(Zdet_SI(ifreq)), 0.05*abs(Zdet_SI(ifreq)));
           end
           fclose(fid);
       end
    end
end

subplot(2,1,1)
plot(CData.MT.T, rhoa, 'ok');
set(gca, 'yscale', 'log', 'Xscale', 'log', 'ylim', [0.1 1.0E6])
subplot(2,1,2)
plot(CData.MT.T, phase, 'ok');
set(gca, 'xscale', 'log')