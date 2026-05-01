clc; clear all;

try
    load('MT_TD_Chain_Processed.mat');
    fprintf('Plotting MT chains (log-depth) \n')
    out_prefix = 'MT_TD_Chain';
    image_name = 'MT_LogPost_image';
catch
    try
        load('DC_TD_Chain_Processed.mat');
        fprintf('Plotting DC chains (log-depth) \n')
        out_prefix = 'DC_TD_Chain';
        image_name = 'DC_LogPost_image';
    catch
        load('MT_DC_TD_Chain_Processed.mat');
        fprintf('Plotting MT+DC chains (log-depth) \n')
        out_prefix = 'MT_DC_TD_Chain';
        image_name = 'MT_DC_LogPost_image';
    end
end

filename = [out_prefix, '_LogStat.mat'];

fprintf('%s','Loading .mat data ...');
fprintf('%s \n','Done');

S.zMin = log10(1); S.zMax = log10(350000);
S.dz = 0.02;

S.nZbins = ceil(( S.zMax - S.zMin )/S.dz);  % number of depth bins
S.zPlot = S.zMin+S.dz/2+(0:S.nZbins-1)*S.dz;  % Depth axis of our PDF plots (midpoints of depth bins)
S.zPlot = 10.^S.zPlot;

S.drho = 0.05; % in log10(rho) (ohm-m)
S.rhoMin = -1;  S.rhoMax = 5;   

S.rho = S.rhoMin : S.drho : S.rhoMax;
S.nRho = length(S.rho) - 1;
S.prior = (1/(S.rhoMax-S.rhoMin)) * ones(S.nRho,1);

S.rhoBinEdges = S.rhoMin:S.drho:S.rhoMax;
S.nRhobins = length(S.rhoBinEdges) - 1;
S.rhoPlot = S.rhoMin + S.drho/2 + (0:S.nRhobins-1)*S.drho; % rho axis of our PDF plots (midpoints of rho bins)

S.nSamples = size(ngrid,1);
S.rhoSamples = nan(S.nZbins,S.nSamples);
S.kSamples = nan(S.nZbins,S.nSamples);

iProgress = 1;
for iSample = 1:S.nSamples
    n_grids = ngrid(iSample);   % Number of interfaces
    z = [z_all(iSample,1:n_grids) S.zMax];
    temp = z(2:end,1)-z(1:end-1,1);
    if min(temp) <= 0
        disp('here');
    end
    rho = rho_all(iSample,1:n_grids);
    % set the counting for all the models
    for j = 1:n_grids   % Till the last layer
        iZbin1 = find(S.zPlot >= z(j), 1, 'first');
        iZbin2 = find(S.zPlot < z(j+1), 1, 'last');
        if ~isempty(iZbin1)
            S.rhoSamples(iZbin1:iZbin2, iSample) = rho(j);
            S.kSamples(iZbin1, iSample) = 1;
        end
    end
    if(mod(iSample,floor(iProgress*S.nSamples/10)) == 0 )
        fprintf('%d %% samples processed...\n',floor(iProgress*10))
        iProgress = iProgress + 1;
    end
end
S.kSamples = S.kSamples(1:end-1,:);

% save([cfile(1:end-4),'_rho.mat'],'rhoall','znew');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    making the histograms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S.posteriorPDF = zeros(S.nZbins,S.nRhobins);
S.p5 = zeros(S.nZbins,1);
S.p95 = zeros(S.nZbins,1);
S.KLd = zeros(S.nZbins,1);
S.pmean = zeros(S.nZbins,1);
iProgress = 1;
for iZbin=1:S.nZbins
    a = histcounts(S.rhoSamples(iZbin,:),S.rhoBinEdges,'Normalization','pdf');
%     a = histogram(S.rhoSamples(iZbin,:),S.rhoBinEdges,'Normalization','pdf');
    S.posteriorPDF(iZbin,:) = a;
    S.p5(iZbin) = prctile(S.rhoSamples(iZbin,:),5);
    S.p95(iZbin) = prctile(S.rhoSamples(iZbin,:),95);
    S.pmean(iZbin) = mean(S.rhoSamples(iZbin,:));
    S.KLd(iZbin) = KLdivergence(S.posteriorPDF(iZbin,:),S.prior);

    if(mod(iZbin,floor(iProgress*S.nZbins/10)) == 0 )
        fprintf('%d %% complete...\n',floor(iProgress*10))
        iProgress = iProgress + 1;
    end
end
S.KLd(1) = 0;  % To make the first interface as invalid as its always there

S.p5 = [S.p5(1); S.p5(1:end-1)];
S.p95 = [S.p95(1); S.p95(1:end-1)];

save(filename,'S','nrms_all','ngrid','step','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnt_headg = 14; fnt_ticks = 12;     fnt_label = 12;    marksize = 8;

xwidth = 800; ywidth = 600;

figure1 = figure;
set(gcf,'PaperPositionMode','auto')
set(figure1, 'Renderer','zbuffer', 'Position', [0 0 xwidth ywidth])

S.zPlot = S.zPlot/1000;

subplot(1,3,1)

h = pcolor(S.rhoPlot,S.zPlot,log10(S.posteriorPDF));
set(h,'EdgeColor','none')
hold on
stairs(S.p5,S.zPlot,'-r','linewidth',1)
stairs(S.p95,S.zPlot,'-r','linewidth',1)
[~, station] = fileparts(pwd);
proj_root = fileparts(fileparts(mfilename('fullpath')));
m1d_candidates = { ...
    fullfile(proj_root, 'data', [station, '_model_1D.mat']), ...
    fullfile(proj_root, 'data', 'model_1D', [station, '_model_1D.mat'])};
m1d_file = '';
for ii = 1:numel(m1d_candidates)
    if exist(m1d_candidates{ii}, 'file'); m1d_file = m1d_candidates{ii}; break; end
end
if ~isempty(m1d_file)
   load(m1d_file, 'model_1D')
   plotModel1D(model_1D)
else
   fprintf('  [Z_03] no reference 1D model found for %s\n', station);
end
xlabel('log(\rho) (ohm-m)')
ylabel('log10 Depth (m)')
set(gca,'YDir','reverse','Xtick',-1:1:6)
ylim([S.zPlot(1) S.zPlot(end-1)])
set(gca,'FontSize',fnt_ticks)
c = colorbar;
c.Label.String = 'log10 (PDF)';

subplot(1,3,2)
kPDF = sum(S.kSamples,2,'omitnan')./(sum(S.kSamples(:),'omitnan')*S.dz);
kPrior = mean(kPDF);
kPDF(1) = 0;  % To make the first interface as invalid as its always there
plot(kPDF,S.zPlot(1:end-1),'linewidth',2)
hold on
plot([kPrior kPrior],[S.zPlot(2) S.zPlot(end-1)],'--k')
set(gca,'YDir','reverse')
set(gca,'FontSize',14)
ylim([S.zPlot(1) S.zPlot(end-1)])
ylabel('log10 Depth (m)')
xlabel('probability density')

subplot(1,3,3)
histogram(ngrid,'Normalization','pdf')
xlabel('# of subsurface layers')
ylabel('probability density')
set(gca,'FontSize',fnt_ticks)

% set(gcf,'position',[50  250  1092  477])

f=gcf;
print(f,image_name,'-djpeg',['-r',num2str(600)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    KL divergence function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [KLd] = KLdivergence(P,Q)

if(length(P) ~= length(Q) )
    fprintf('length of input arrays must be the same\n')
    keyboard
    return
end

KLd = 0.0;

for j=1:length(P)
    if( P(j) > 0 && Q(j) > 0 )
        KLd = KLd + P(j)*(log(P(j)) - log(Q(j)));
    end
end

end

function [] = plotModel1D(x)% plotModel1D(b,l,x)

z = x.z;
r = x.rho;

for i=1:length(r)
   k = 2*i-1;
   rr(k) = r(i);
   rr(k+1) = r(i);
end
for i=1:length(z)-1
    k = 2*i - 1;
   zz(k) = z(i);
   zz(k+1) = z(i+1);
end
zz(k+2) = z(end);
zz(k+3) = 2*z(end);

plot(rr,zz,'--k','linewidth',2)
set (gca,'ydir','reverse')
hold on
end