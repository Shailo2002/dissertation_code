clc; clear all;

cfile = 'MT_DC_TD_Chain_Combined.mat';
filename = 'MT_DC_exp3.mat';

fprintf('%s','Loading .mat data ...');
load(cfile);
fprintf('%s \n','Done');

S.layers = zeros(S.nSamples,1);
S.thkt_info = zeros(S.nSamples,30);
S.rho_info = zeros(S.nSamples,30);

iProgress = 1;
for iSample = 1:S.nSamples
    n_grids = ngrid(iSample);   % Number of interfaces
    z = [z_all(iSample,1:n_grids) S.zMax];
    temp = z(2:end,1)-z(1:end-1,1);
    
    S.layers = length(temp);
    if min(temp) <= 0
        disp('here');
    end
    S.layers
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

% figure
% h = pcolor(S.nRhobins,S.zPlot,log10(S.posteriorPDF));
% set(h,'EdgeColor','none')
% hold on
% stairs(S.p5,S.zPlot,'-r','linewidth',1)
% stairs(S.p95,S.zPlot,'-r','linewidth',1)
% xlabel('log(\rho) (ohm-m)')
% ylabel('depth (m)')
% set(gca,'YDir','reverse','Xtick',-1:1:6)
% ylim([S.zPlot(1) S.zPlot(end-1)])
% set(gca,'FontSize',14)
% c = colorbar;

xwidth = 800; ywidth = 600;

figure1 = figure;
set(gcf,'PaperPositionMode','auto')
set(figure1, 'Renderer','zbuffer', 'Position', [0 0 xwidth ywidth])

% S.zPlot = 10.^S.zPlot/1000;

subplot(1,3,1)

h = pcolor(S.rhoPlot,S.zPlot,log10(S.posteriorPDF));
set(h,'EdgeColor','none')
hold on
stairs(S.p5,S.zPlot,'-r','linewidth',1)
stairs(S.p95,S.zPlot,'-r','linewidth',1)
if exist('MT_Example_3.mat','file')
   load('MT_Example_3.mat')
   y.z = [0; cumsum(MT.thicknesses)]; y.z = y.z(1:end-1);
   y.rho = log10(MT.resistivities);
   y.z(1) = -1.;
   y.z(2:end) = log10(y.z(2:end));
   plotModel1D(y)
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
print(f,cfile(1:end-4),'-djpeg',['-r',num2str(600)]);

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