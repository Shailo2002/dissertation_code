clc; clear all; close all

fnt_headg = 10; fnt_ticks = 10;     fnt_label = 10;    marksize = 8;
                        % Nh, Nw, gap,    marg_h,   marg_w
[ha, pos] = tight_subplot(3,3,[.06 .07],[.06 .04],[.1 .11]);
close all

figure1 = figure;
set(gcf,'PaperPositionMode','auto')
set(figure1, 'Position', [0 0 500 1000])

xmodel.z = [0 800 1050]/1000;
xmodel.rho = log10([250 15 5000]);

loglims = [-4 0.5];
strng = ['(a)';'(b)';'(c)'];
for i = 1:3
    if i == 1
        load('MT_DC_exp3.mat','S','nrms_all','ngrid','step')     
    elseif i== 2
        load('MT_DC_exp3.mat','S','nrms_all','ngrid','step')     
    elseif i == 3
        load('MT_DC_exp3.mat','S','nrms_all','ngrid','step')     
    end
    
    pv = pos{i};
    subplot = axes('Position',pv,'Parent',figure1,...
        'XTick',0:6,'XMinorTick','on',...
        'YTick',0:0.5:2,'YMinorTick','on','YDir','reverse',...
        'TickDir','in','MinorGridLineStyle','none',...
        'FontSize',fnt_ticks, 'CLim', loglims);
    xlim(subplot,[0 6]);     ylim(subplot,[0 2]);
    box(subplot,'on');        hold(subplot,'all');
    
    S.zPlot = 10.^S.zPlot/1000;
    h = pcolor(S.rhoPlot,S.zPlot,log10(S.posteriorPDF));
    set(h,'EdgeColor','none')
    
    stairs(S.p5,S.zPlot,'--r','linewidth',1)
    stairs(S.p95,S.zPlot,'--r','linewidth',1)
    stairs(S.pmean,S.zPlot,'-r','linewidth',1)
    if i == 1
        ylabel('Depth (km)','FontWeight','bold','FontSize',fnt_headg)
    end
    xlabel('log_1_0 \rho (\Omega-m)','FontWeight','bold','FontSize',fnt_headg)
        
    if i == 3
        c = colorbar;
        set(c, 'Position', [0.90 pv(2) 0.025 pv(4)])
        c.Label.String = 'log10 (PDF)';
        c.Ticks = -4:0;
    end
    
    xp = -1.5; yp = -0.15;
    text(xp, yp,strng(i,:),'FontWeight','bold','FontSize',fnt_headg)
    
    plotModel1D(xmodel);
end

strng = ['(d)';'(e)';'(f)'];
for i = 1:3
    if i == 1
        load('MT_exp3.mat','S','nrms_all','ngrid','step')     
    elseif i== 2
        load('DC_exp3.mat','S','nrms_all','ngrid','step')     
    elseif i == 3
        load('MT_DC_exp3.mat','S','nrms_all','ngrid','step')     
    end
    
    pv = pos{i+3};
    subplot = axes('Position',pv,'Parent',figure1,...
        'XMinorTick','on',...
        'YTick',0:0.5:2,'YMinorTick','on','YDir','reverse',...
        'TickDir','in','MinorGridLineStyle','none',...
        'FontSize',fnt_ticks);
    xlim(subplot,[0 4]);     
    ylim(subplot,[0 2]);
    box(subplot,'on');        hold(subplot,'all');
    
    xp = -1.20; yp = -0.15;
    text(xp, yp,strng(i,:),'FontWeight','bold','FontSize',fnt_headg)
    
    kPDF = sum(S.kSamples,2,'omitnan')./(sum(S.kSamples(:),'omitnan')*S.dz);
    kPrior = mean(kPDF);
    kPDF(1) = 0;  % To make the first interface as invalid as its always there
    S.zPlot = 10.^S.zPlot/1000;
    plot(kPDF,S.zPlot(1:end-1),'linewidth',2)
    hold on
%     plot([kPrior kPrior],[S.zPlot(2) S.zPlot(end-1)],'--k')
    
    plot([0 4],[0.8 0.8],'--k', 'linewidth', 1);
    plot([0 4],[1.05 1.05],'--k', 'linewidth', 1);
    
    if i == 1
        ylabel('Depth (km)','FontWeight','bold','FontSize',fnt_headg)
    end
    xlabel('Int. prob.','FontWeight','bold','FontSize',fnt_headg)
end

strng = ['(g)';'(h)';'(i)'];

for i = 1:3
    if i == 1
        load('MT_exp3.mat','ngrid')     
    elseif i== 2
        load('DC_exp3.mat','ngrid')     
    elseif i == 3
        load('MT_DC_exp3.mat','ngrid')     
    end
    
    pv = pos{i+6};
    subplot = axes('Position',pv,'Parent',figure1,...
        'YTick',0:.1:.3,'XMinorTick','on',...
        'XTick',0:5:50,'YMinorTick','on',...
        'TickDir','in','MinorGridLineStyle','none',...
        'FontSize',fnt_ticks);
    xlim(subplot,[0 15]);     
    ylim(subplot,[0 0.3]);
    box(subplot,'on');        hold(subplot,'all');
    
    histogram(ngrid,'Normalization','pdf')
    
    xp = -5; yp = .325;
    text(xp, yp,strng(i,:),'FontWeight','bold','FontSize',fnt_headg)
    
    if i == 1
        ylabel('probability density','FontWeight','bold','FontSize',fnt_headg)
    end
    xlabel('# layers','FontWeight','bold','FontSize',fnt_headg)
end

f = gcf; %f is the handle of the figure you want to export
print(f,fullfile(pwd,'Exp_3_Post'),'-djpeg',['-r',num2str(600)]) %save file


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
plot(rr,zz,'--k','linewidth',1)
hold on
end