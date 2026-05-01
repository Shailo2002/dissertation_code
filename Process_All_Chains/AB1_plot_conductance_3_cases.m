clc; clear all; close all;

z_max = log10(5000);

fnt_headg = 14; fnt_ticks = 12;     fnt_label = 14;    marksize = 8;
                        % Nh, Nw, gap,    marg_h,   marg_w
[ha, pos] = tight_subplot(1, 3, [.06 .07],[.17 .04],[.1 .05]);
close all

figure1 = figure;
set(gcf,'PaperPositionMode','auto')
set(figure1, 'Position', [0 0 1000 300])

%%%%%%%%%
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

pv = pos{1};
subplot = axes('Position',pv,'Parent',figure1,...
    'XMinorTick','on',...
    'YMinorTick','on',...
    'TickDir','in','MinorGridLineStyle','none',...
    'FontSize',fnt_ticks);
%     xlim(subplot,[1 5]);      ylim(subplot,[0.5 4]);
box(subplot,'on');        hold(subplot,'all');

h = histogram(10.^layer_thickness(:,2)./(10.^layer_rho(:,2)),'Normalization','pdf');
% set(h, 'edgecolor', 'none');
plot([250/15 250/15],[0 0.5],'--k', 'linewidth', 1);
xlabel('Conductance (S)','FontWeight','bold','FontSize',fnt_headg)
ylabel('Probability','FontWeight','bold','FontSize',fnt_headg)

%%%%%%%%
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
corr(layer_rho(:,2), layer_thickness(:,2))

pv = pos{2};
subplot = axes('Position',pv,'Parent',figure1,...
    'XMinorTick','on',...
    'YMinorTick','on',...
    'TickDir','in','MinorGridLineStyle','none',...
    'FontSize',fnt_ticks);
%     xlim(subplot,[1 5]);      ylim(subplot,[0.5 4]);
box(subplot,'on');        hold(subplot,'all');

histogram(10.^layer_thickness(:,2)./(10.^layer_rho(:,2)),'Normalization','pdf')
plot([250/15 250/15],[0 0.5],'--k', 'linewidth', 1);
xlabel('Conductance (S)','FontWeight','bold','FontSize',fnt_headg)
ylabel('Probability','FontWeight','bold','FontSize',fnt_headg)

%%%%%%
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
corr(layer_rho(:,2), layer_thickness(:,2))

pv = pos{3};
subplot = axes('Position',pv,'Parent',figure1,...
    'XMinorTick','on',...
    'YMinorTick','on',...
    'TickDir','in','MinorGridLineStyle','none',...
    'FontSize',fnt_ticks);
%     xlim(subplot,[1 5]);      ylim(subplot,[0.5 4]);
box(subplot,'on');        hold(subplot,'all');

histogram(10.^layer_thickness(:,2)./(10.^layer_rho(:,2)),'Normalization','pdf')
plot([250/15 250/15],[0 1],'--k', 'linewidth', 1);
xlabel('Conductance (S)','FontWeight','bold','FontSize',fnt_headg)
ylabel('Probability','FontWeight','bold','FontSize',fnt_headg)

f = gcf; %f is the handle of the figure you want to export
print(f,fullfile(pwd,'exp3_cond_hist'),'-djpeg',['-r',num2str(600)]) %save file
