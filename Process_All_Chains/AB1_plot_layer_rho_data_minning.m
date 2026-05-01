clc; clear all; close all

z_max = log10(350000);

fnt_headg = 10; fnt_ticks = 10;     fnt_label = 10;    marksize = 8;
                        % Nh, Nw, gap,    marg_h,   marg_w
[ha, pos] = tight_subplot(3,3,[.06 .07],[.06 .04],[.1 .11]);
close all

figure1 = figure;
set(gcf,'PaperPositionMode','auto')
set(figure1, 'Position', [0 0 1000 1000])

loglims = [-4 0.5];

load('MT_TD_Chain_Combined.mat');
% Find all the layers having 3 layers
indx = find(ngrid == 3);   length(indx)
layer_thickness = zeros(length(indx),3);
layer_rho = zeros(length(indx),3);

for i = 1:length(indx)
    k = indx(i);   % model number    
    n_grids = ngrid(k);   % Number of interfaces
    z = 10.^[z_all(k,1:n_grids) z_max];  z(1) = 0;
    layer_thickness(i,:) = log10(z(2:end)-z(1:end-1));
    layer_rho(i,:) = rho_all(k,1:n_grids);
end
corr(layer_rho(:,2), layer_thickness(:,2))

strng = ['(a)';'(b)';'(c)'];
for i = 1:3
    pv = pos{i};
    subplot = axes('Position',pv,'Parent',figure1,...
        'XMinorTick','on',...
        'YMinorTick','on',...
        'TickDir','in','MinorGridLineStyle','none',...
        'FontSize',fnt_ticks);
%     xlim(subplot,[1 5]);      ylim(subplot,[0.5 4]);
    box(subplot,'on');        hold(subplot,'all');
    
    h = plot(layer_rho(:,i), layer_thickness(:,i), '.r');
    if i == 1
        ylabel('log (thickness)','FontWeight','bold','FontSize',fnt_headg)
    end
    xlabel('log_1_0 \rho (\Omega-m)','FontWeight','bold','FontSize',fnt_headg)
            
%     xp = -1.5; yp = -0.15;
%     text(xp, yp,strng(i,:),'FontWeight','bold','FontSize',fnt_headg)
end

load('DC_TD_Chain_Combined.mat');
% Find all the layers having 3 layers
indx = find(ngrid == 3);   length(indx)
layer_thickness = zeros(length(indx),3);
layer_rho = zeros(length(indx),3);
for i = 1:length(indx)
    k = indx(i);   % model number    
    n_grids = ngrid(k);   % Number of interfaces
    z = 10.^[z_all(k,1:n_grids) z_max];  z(1) = 0;
    layer_thickness(i,:) = log10(z(2:end)-z(1:end-1));
    layer_rho(i,:) = rho_all(k,1:n_grids);
end

strng = ['(d)';'(e)';'(f)'];
for i = 1:3
    pv = pos{i+3};
    subplot = axes('Position',pv,'Parent',figure1,...
        'XMinorTick','on',...
        'YMinorTick','on',...
        'TickDir','in','MinorGridLineStyle','none',...
        'FontSize',fnt_ticks);
%     xlim(subplot,[1 5]);      ylim(subplot,[0.5 4]);
    box(subplot,'on');        hold(subplot,'all');
    
    h = plot(layer_rho(:,i), layer_thickness(:,i), '.r');
    if i == 1
        ylabel('log (thickness)','FontWeight','bold','FontSize',fnt_headg)
    end
    xlabel('log_1_0 \rho (\Omega-m)','FontWeight','bold','FontSize',fnt_headg)
            
%     xp = -1.5; yp = -0.15;
%     text(xp, yp,strng(i,:),'FontWeight','bold','FontSize',fnt_headg)
end

load('MT_DC_TD_Chain_Combined.mat');
% Find all the layers having 3 layers
indx = find(ngrid == 3);   length(indx)
layer_thickness = zeros(length(indx),3);
layer_rho = zeros(length(indx),3);
for i = 1:length(indx)
    k = indx(i);   % model number    
    n_grids = ngrid(k);   % Number of interfaces
    z = 10.^[z_all(k,1:n_grids) z_max];  z(1) = 0;
    layer_thickness(i,:) = log10(z(2:end)-z(1:end-1));
    layer_rho(i,:) = rho_all(k,1:n_grids);
end
clear subplot
% figure
% subplot(3,1,1);
% histogram(layer_thickness(:,2));
% subplot(3,1,2);
% histogram(layer_rho(:,2));
% subplot(3,1,3);
% histogram(10.^layer_thickness(:,2)./(10.^layer_rho(:,2)));
% corr(layer_rho(:,2), layer_thickness(:,2))
% pause

strng = ['(g)';'(h)';'(i)'];
for i = 1:3
    pv = pos{i+6};
    subplot = axes('Position',pv,'Parent',figure1,...
        'XMinorTick','on',...
        'YMinorTick','on',...
        'TickDir','in','MinorGridLineStyle','none',...
        'FontSize',fnt_ticks);
%     xlim(subplot,[1 5]);      ylim(subplot,[0.5 4]);
    box(subplot,'on');        hold(subplot,'all');
    
    h = plot(layer_rho(:,i), layer_thickness(:,i), '.r');
    if i == 1
        ylabel('log (thickness)','FontWeight','bold','FontSize',fnt_headg)
    end
    xlabel('log_1_0 \rho (\Omega-m)','FontWeight','bold','FontSize',fnt_headg)
            
%     xp = -1.5; yp = -0.15;
%     text(xp, yp,strng(i,:),'FontWeight','bold','FontSize',fnt_headg)
end
f = gcf; %f is the handle of the figure you want to export
print(f,fullfile(pwd,'Exp_3_data_mining'),'-djpeg',['-r',num2str(600)]) %save file


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