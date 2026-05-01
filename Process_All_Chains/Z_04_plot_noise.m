clc; clear all;

cfile = 'MT_DC_TD_Chain_Combined.mat';

fprintf('%s','Loading .mat data ...');
load(cfile);
fprintf('%s \n','Done');

figure
subplot(2,1,1)
hold on
for j = 1:CData.nChains_atT1
    if strcmpi(CData.inversion_method, 'MT_DC')
        plot(1:size(sigma,1), sigma(:,j,1), '-k'); % set(gca,'XScale','log'); % axis([0 length(Samples.misfit) 0 10]);
        plot(1:size(sigma,1), sigma(:,j,2), '-r')
    else
        plot(1:size(sigma,1), sigma(:,j), '-')
    end
    % plot(1:size(sigma,1), sigma(:,j), '-'); % set(gca,'XScale','log'); % axis([0 length(Samples.misfit) 0 10]);
end
xlabel('# samples')
ylabel('Lambda')

subplot(2,1,2)
hold on
temp = sigma(:,:,1);
histogram(temp(:),'Normalization','pdf')
temp = sigma(:,:,2);
histogram(temp(:),'Normalization','pdf')
xlabel('Lambda')
ylabel('p(lambda)')
legend ('MT','DCR')
f=gcf;
print(f,'Noise','-djpeg',['-r',num2str(600)]);