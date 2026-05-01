clc; clear all; close all;

[ha, pos] = tight_subplot(1,3,[.07 .07],[.10 .06],[.09 .02]);
close all

% % Model limits for
yl = 0;    yr = 8;   zu = 0;  zd = 1.2;  
% % Ticks labels
ytickks = 0:1:8;      ztickks = 0:0.2:1.2;
% % Marker size propoerties
% % ylabel to right hand side, working as headings.
fnt_ticks = 12;  fnt_label = 10;  fnt_headg = 14;

figure1 = figure;
set(gcf,'PaperPositionMode','auto')
set(figure1, 'Renderer','zbuffer', 'Position', [0 0 800 600])

    
for i = 1:3
    pv = pos{i};
    subplot = axes('Position',pv,'Parent',figure1,...
        'XTick',ytickks,'XMinorTick','off',...
        'YTick',ztickks,'YMinorTick','off',...
        'TickDir','in','MinorGridLineStyle','-',...
        'FontSize',fnt_ticks);
    xlim(subplot,[yl yr]);    ylim(subplot,[zu zd]);
    box(subplot,'on');        hold(subplot,'all');
    if i == 1
        load('MT_TD_Chain_Combined.mat')
        [a, edges] = histcounts(sigma(:),'Normalization','pdf');
        plot(edges(1:end-1)-0.5*(edges(2)-edges(1)), a, '-r','linewidth',1.5);
        % histogram(sigma(:),'Normalization','pdf')
        legend ('MT')
    elseif i == 2
        load('DC_TD_Chain_Combined.mat')
        [a, edges] = histcounts(sigma(:),'Normalization','pdf');
        plot(edges(1:end-1)-0.5*(edges(2)-edges(1)), a, '-b','linewidth',1.5);
        legend ('DCR')
    else
        load('MT_DC_TD_Chain_Combined.mat')
        temp = sigma(:,:,1);
        [a, edges] = histcounts(temp(:),'Normalization','pdf');
        plot(edges(1:end-1)-0.5*(edges(2)-edges(1)), a, '-r','linewidth',1.5);
        temp = sigma(:,:,2);
        [a, edges] = histcounts(temp(:),'Normalization','pdf');
        plot(edges(1:end-1)-0.5*(edges(2)-edges(1)), a, '-b','linewidth',1.5);
        legend ('MT','DCR')
    end
    legend boxoff     
    if i == 1
        text(-2,1.25,'(a)','FontWeight','bold','FontSize',fnt_headg)
    elseif i==2
        text(-2,1.25,'(b)','FontWeight','bold','FontSize',fnt_headg)
    elseif i==3
        text(-2,1.25,'(c)','FontWeight','bold','FontSize',fnt_headg)
    end
    xlabel('Lambda factor','FontWeight','bold','FontSize',fnt_headg)
    if i == 1
        ylabel('p(\lambda|d)','FontWeight','bold','FontSize',fnt_headg)
    end
    
    set(gca,'layer','top');
   
end
f = gcf; %f is the handle of the figure you want to export
print(f,fullfile(pwd,'Lambda_histogram'),'-djpeg',['-r',num2str(600)]) %save file
