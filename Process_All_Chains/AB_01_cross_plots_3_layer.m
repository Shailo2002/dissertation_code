clc; clear; close all;

% Resolve the log-depth stats file produced by Z_03 for whichever
% inversion variant is present in the current station folder.
if exist('MT_TD_Chain_LogStat.mat','file')
    statfile = 'MT_TD_Chain_LogStat.mat';
elseif exist('DC_TD_Chain_LogStat.mat','file')
    statfile = 'DC_TD_Chain_LogStat.mat';
elseif exist('MT_DC_TD_Chain_LogStat.mat','file')
    statfile = 'MT_DC_TD_Chain_LogStat.mat';
else
    error('AB_01_cross_plots_3_layer: no *_LogStat.mat in current folder. Run Z_03 first.');
end
fprintf('Loading %s\n', statfile);
load(statfile);
zgrid = S.zPlot/1000;

% Find the indices of the desired location of the plots
z = [5 15; 120 200; 250 340];

for j = 1:3
    ind1 = find(zgrid >= z(j,1));
    ind2 = find(zgrid >= z(j,2));
    index(j,:) = [ind1(1) ind2(1)];
end

load(statfile)
for j = 1:3
    temp = S.rhoSamples(index(j,1):index(j,2),:);
    rho_rmt{j} = temp(:);
end

% Model limits for
yl = -1;   yr = 5;   zu = -1;  zd = 5;
% Ticks labels
ytickks = -1:6;     ztickks = -1:6;
% ylabel to right hand side, working as headings.
fnt_ticks = 10;  fnt_headg = 10;
%(Nh,Nw,gap, marg_h, marg_w)
[ha, pos] = tight_subplot(3,3,[.07 .07],[.15 .02],[.06 .02]);
close all

rhoBinEdges = S.rhoMin:S.drho:S.rhoMax;
rhoBinEdgesplot = rhoBinEdges(1:end-1) + 0.5*(rhoBinEdges(2:end)-rhoBinEdges(1:end-1));
nRhobins = length(rhoBinEdges) - 1;

figure1 = figure;
set(gcf,'PaperPositionMode','auto')
set(figure1, 'Renderer','zbuffer', 'Position', [0 0 800 800])

k = 0; skip = 50;
z_value = [1 3 1 3];
y_value = [3 1 3 1];

for j = 1:3
    for i = 1:3
        k = k + 1;
               
        pv = pos{k};
        if j > i; continue; end
        if i ~= j
            subplot = axes('Position',pv,'Parent',figure1,...
                'XTick',ytickks,'XMinorTick','off',...
                'YTick',ztickks,'YMinorTick','off',...
                'TickDir','in','MinorGridLineStyle','-',...
                'FontSize',fnt_ticks);
            xlim(subplot,[yl yr]);    ylim(subplot,[zu zd]);
            box(subplot,'on');        hold(subplot,'all');
        else
            subplot = axes('Position',pv,'Parent',figure1,...
                'XTick',ytickks,'XMinorTick','off',...
                'YMinorTick','off',...
                'TickDir','in','MinorGridLineStyle','-',...
                'FontSize',fnt_ticks);
            xlim(subplot,[yl yr]);    
            box(subplot,'on');        hold(subplot,'all');
        end
        
        if i == j
            temp1 = rho_rmt{i};
            post1 = histcounts(temp1,rhoBinEdges,'Normalization','pdf');            
            plot(rhoBinEdgesplot,post1,'-r');
            xlabel(['log_1_0 \rho',num2str(j)],'FontWeight','bold','FontSize',fnt_headg)
            ylabel('p(m)','FontWeight','bold','FontSize',fnt_headg)
            
            % if i == 1
            %      plot([log10(250) log10(250)],[0 6],'--k');
            % elseif i == 2
            %     plot([log10(15) log10(15)],[0 3],'--k');
            % elseif i == 3
            %     plot([log10(5000) log10(5000)],[0 1.5],'--k');
            % end
        else
            i1 = i;  i2 = j;
            temp1 = rho_rmt{i1}; 
            temp4 = rho_rmt{i2}; 
            % disp([j i i1 i2])
            [n1,~] = size(rho_rmt{i1});    [n2,~] = size(rho_rmt{i2});      
            n = min(n1, n2);
            plot(temp1(1:skip:n),temp4(1:skip:n),'r.','markersize',0.5);
            % if j == 1
            %     if i == 2
            %         plot(log10(15),log10(250),'ko','markersize',5,'linewidth',1.5);
            %     elseif i == 3
            %         plot(log10(5000),log10(250),'ko','markersize',5,'linewidth',1.5);
            %     end
            % elseif j == 2
            %     if i == 3
            %         plot(log10(5000),log10(15),'ko','markersize',5,'linewidth',1.5);
            %     end
            % elseif j == 3
            % 
            % end
            xlabel(['log_1_0 \rho_',num2str(i)],'FontWeight','bold','FontSize',fnt_headg)
            ylabel(['log_1_0 \rho_',num2str(j)],'FontWeight','bold','FontSize',fnt_headg)
        end        
    end
end

subplot = axes('Position',[0.4 0.2 0.15 0.1],'Parent',figure1,...
    'XTick',[],'XMinorTick','off',...
    'YTick',[],'YMinorTick','off',...
    'TickDir','in','MinorGridLineStyle','none');
box(subplot,'on');        hold(subplot,'all');
xlim(subplot,[0 1]);      ylim(subplot,[0 1]);

% h = rectangle; h.LineWidth = 1;
% h.FaceColor = 'r'; h.EdgeColor = 'r'; h.Position = [0.1 0.8 0.25 0.05];
% h = rectangle; h.LineWidth = 1;
% h.FaceColor = 'b'; h.EdgeColor = 'b'; h.Position = [0.1 0.5 0.25 0.05];
% h = rectangle; h.LineWidth = 1;
% h.FaceColor = 'g'; h.EdgeColor = 'g'; h.Position = [0.1 0.2 0.25 0.05];
% text(0.4,0.85,'MT','FontWeight','bold','FontSize',fnt_headg)
% text(0.4,0.55,'DCR','FontWeight','bold','FontSize',fnt_headg)
% text(0.4,0.25,'Joint','FontWeight','bold','FontSize',fnt_headg)


% MT.thicknesses = [800 250 950];
% MT.resistivities = [250 15 5000];
% y.z = [0 cumsum(MT.thicknesses)]; y.z = y.z(1:end-1);
% y.rho = log10(MT.resistivities);

% plot the posteriori plots
load(statfile)
% Model limits for
yl = -1;    yr = 5;     zu = 0;  zd = 350;
% Ticks labels
ytickks = [0 2 4];     ztickks = 0:50:2000;

pv = [0.08 0.1 0.18 0.55];
subplot = axes('Position',pv,'Parent',figure1,...
    'XTick',ytickks,'XMinorTick','off',...
    'YTick',ztickks,'YMinorTick','off',...
    'TickDir','in','MinorGridLineStyle','-',...
    'ydir','reverse','Clim',[-5 0],'FontSize',fnt_ticks);
xlim(subplot,[yl yr]);    ylim(subplot,[zu zd]);
box(subplot,'on');        hold(subplot,'all');

S.zPlot = S.zPlot/1000;      S.posteriorPDF = log10(S.posteriorPDF);
h = pcolor(S.rhoPlot, S.zPlot,S.posteriorPDF);
set(h,'EdgeColor','none')
set(gca,'layer','top');
for j = 1:3
    plot([-1 6], [z(j,1) z(j,1)], '-w', 'linewidth', 1.5);
    plot([-1 6], [z(j,2) z(j,2)], '-w', 'linewidth', 1.5);
end
xlabel('log_1_0 (\Omega-m)','FontWeight','bold','FontSize',fnt_headg)
ylabel('Depth (km)','FontWeight','bold','FontSize',fnt_headg)

% plotModel1D(y);

f = gcf; %f is the handle of the figure you want to export
print(f,fullfile(pwd,'Cross_Plot_3layer'),'-djpeg',['-r',num2str(600)]) %save file

function [] = plotModel1D(x)% plotModel1D(b,l,x)
z = x.z;   r = x.rho;
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
plot(rr,zz,'--r','linewidth',2)
% set (gca,'ydir','reverse')
% hold on
end
