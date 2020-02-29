% compare observations toz and zonal wind
clear all
SSWwind = 0;
SPVwind = 40;

tozmonth = 3;
lats = [63,90];
todetrend = 1;
days = [day(datetime('1-Jan-2018'),'dayofyear'),day(datetime('1-Feb-2018'),'dayofyear'),...
    day(datetime('1-Mar-2018'),'dayofyear'),day(datetime('1-Apr-2018'),'dayofyear'),...
    day(datetime('1-May-2018'),'dayofyear'),day(datetime('1-Jun-2018'),'dayofyear'),...
    day(datetime('1-Jul-2018'),'dayofyear'),day(datetime('1-Aug-2018'),'dayofyear'),...
    day(datetime('1-Sep-2018'),'dayofyear'),day(datetime('1-Oct-2018'),'dayofyear'),...
    day(datetime('1-Nov-2018'),'dayofyear'),day(datetime('1-Dec-2018'),'dayofyear'),366];

days_shift = [1,days(8)-days(7)+1,days(9)-days(7)+1,days(10)-days(7)+1,days(11)-days(7)+1,...
    days(12)-days(7)+1,366-days(7)+1,366-days(7)+31+1,366-days(7)+31+28+1,366-days(7)+31+28+31+1,...
    366-days(7)+31+28+31+30+1,366-days(7)+31+28+31+30+31+1,366];

cbrew = cbrewer('qual','Paired',10);
plotozoneextremes = 0;

%% Read in Bodeker toz

tcolat = [63,90];
tozmonth = 3;
toz.timeperiod = [1980,2016];

directory = ['/Volumes/ExternalOne/work/data/BodekerScientific/TCO/daily/6090N/'];
files = dir([directory,'*.nc']);

for i = 1:length(files)
    [~,BSdata(i),~] = Read_in_netcdf([directory,files(i).name]);
    if i == 1
        latindex_bs = find(BSdata(i).latitude >= tcolat(1) & BSdata(i).latitude <= tcolat(2));
%         years = CCMI_years(BSdata(i).date,1);      
%         yearsUnique = unique(years);
    end
    varweighted(i,:) = weightedaverage(squeeze(nanmean(BSdata(i).tco(:,latindex_bs,1:365),1)),BSdata(i).latitude(latindex_bs));    
    %varweighted_shift(i,:) = circshift(varweighted(i,:),365-days(7)-1);
    %varweighted_shift(i,end-365:end) = NaN;
    %take monthly average
%     count = 1;
%     for k = 1:length(yearsUnique)
%         tozdaily_years_shift(i,k,:) = varweighted_shift(i,count:count+364);
%         tozdaily_years(i,k,:) = varweighted(i,count:count+364);
%         count = count+365;
        for j = 1:12
            tozdaily_monthlymean(j,i) = nanmean(varweighted(i,days(j):days(j+1)-1),2);            
        end        
%     end
    
end

tozdaily_monthlymean_detrend = [detrend(tozdaily_monthlymean')+nanmean(tozdaily_monthlymean,2)']';

%%
varweighted_shift = varweighted';
varweighted_shift = varweighted_shift(:);
varweighted_shift = circshift(varweighted_shift,365-days(7)-1);
varweighted_shift(1:365) = NaN;

count= 1;
for k = 1:37
    toz_years_shift(k,:) = varweighted_shift(count:count+364);    
    count = count+365;    
end

%%
upperpct = prctile(squeeze(tozdaily_monthlymean(3,:)),80);
lowerpct = prctile(squeeze(tozdaily_monthlymean(3,:)),20);

upperind = squeeze(tozdaily_monthlymean(3,:)) >= upperpct;
lowerind = squeeze(tozdaily_monthlymean(3,:)) <= lowerpct;

% upperpct = prctile(squeeze(tozdaily_monthlymean_detrend(3,:)),80);
% lowerpct = prctile(squeeze(tozdaily_monthlymean_detrend(3,:)),20);
% 
% upperind = squeeze(tozdaily_monthlymean_detrend(3,:)) >= upperpct;
% lowerind = squeeze(tozdaily_monthlymean_detrend(3,:)) <= lowerpct;


%% Read in wind

timeperiod = 1980:2018;

directory = '/Volumes/ExternalOne/work/data/MERRA/U/daily/';   
obvar = 'U';
oblat = 60;

files = dir([directory,'*.nc']);
for i = 1:length(files)-2
    
    [winddata(i).data] = ncread([directory,files(i).name],obvar);
    [time(i).data] = ncread([directory,files(i).name],'time');
    if i == 1
        [latitude] = ncread([directory,files(i).name],'lat');
        [longitude] = ncread([directory,files(i).name],'lon');
    end
    [~,obslatind] = min(abs(latitude-oblat));
    windzonalmean(i,:) = squeeze(nanmean(winddata(i).data(:,obslatind,:,1:365),1));
    
    
    
end

%% 
windzonalmean_shift = windzonalmean';
windzonalmean_shift = windzonalmean_shift(:);
windzonalmean_shift = circshift(windzonalmean_shift,365-days(7)-1);
windzonalmean_shift(1:365) = NaN;

count = 1;
for k = 1:37
    winddaily_years_shift(k,:) = windzonalmean_shift(count:count+364);    
    count = count+365;    
end

%% plotting
cbrew = cbrewer('qual','Paired',10);
fsize = 18;
monslabel = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

figure;
set(gcf,'position',[100 100 700 800],'color','white');

%toz
subplot(2,1,1)
plot(toz_years_shift(lowerind,:)','color',cbrew(1,:),'LineWidth',2);
hold on
plot(toz_years_shift(upperind,:)','color',cbrew(5,:),'LineWidth',2);
plot(nanmean(toz_years_shift(lowerind,:)),'color',cbrew(2,:),'LineWidth',3);
plot(nanmean(toz_years_shift(upperind,:)),'color',cbrew(6,:),'LineWidth',3);
ylim([260,530]);     
xlim([-5 370]);
plot([0,400],[nanmean(upperpct),nanmean(upperpct)],'--','color','k','LineWidth',2)
plot([0,400],[nanmean(lowerpct),nanmean(lowerpct)],'--','color','k','LineWidth',2)

set(gca,'fontsize',fsize-2)
set(gca,'xtick',days_shift(1:end-1),'xticklabel',monslabel);
ylabel('Dobson units','fontsize',fsize);

title('Observed Arctic average TCO','fontsize',fsize+2)
%wind
subplot(2,1,2)
%plot(winddaily_years_shift','color',[.7 .7 .7],'LineWidth',2);
hold on
plot(winddaily_years_shift(lowerind,:)','color',cbrew(1,:),'LineWidth',2);
plot(winddaily_years_shift(upperind,:)','color',cbrew(5,:),'LineWidth',2);

%plot(nanmean(winddaily_years_shift),'color','k','LineWidth',3);
plot(nanmean(winddaily_years_shift(lowerind,:)),'color',cbrew(2,:),'LineWidth',3);
plot(nanmean(winddaily_years_shift(upperind,:)),'color',cbrew(6,:),'LineWidth',3);
ylim([-30,80]); 
xlim([-5 370]);
plot([0,400],[0 0],'--','color','k','LineWidth',2)
plot([0,400],[40 40],'--','color','k','LineWidth',2)

set(gca,'fontsize',fsize)
set(gca,'xtick',days_shift(1:end-1),'xticklabel',monslabel);    
ylabel('m/s','fontsize',fsize);
xlabel('Month','fontsize',fsize);
title(['Zonal wind at 60',char(176),'N and 10hPa'],'fontsize',fsize+2)

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/Observedtozandwind.pdf'];

export_fig(filename,'-pdf');

%% plotting standard deviations
%% plotting
cbrew = cbrewer('qual','Set1',10);
fsize = 18;
monslabel = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

figure;
set(gcf,'position',[100 100 700 800],'color','white');

%toz
subplot(2,1,1)
plot(nanstd(toz_years_shift(lowerind,:)),'color',cbrew(2,:),'LineWidth',3,'LineStyle','-');
hold on
plot(nanstd(toz_years_shift(upperind,:)),'color',cbrew(1,:),'LineWidth',3,'LineStyle','-');
plot(nanstd(toz_years_shift),'color','k','LineWidth',3,'LineStyle','-.');
ylim([0,40]);     
xlim([-5 370]);

set(gca,'fontsize',fsize-2)
set(gca,'xtick',days_shift(1:end-1),'xticklabel',monslabel);
ylabel('Standard deviation (Dobson units)','fontsize',fsize);

title('Observed Arctic average TCO','fontsize',fsize+2)
%wind
subplot(2,1,2)
box on
hold on
plot(nanstd(winddaily_years_shift(lowerind,:)),'color',cbrew(2,:),'LineWidth',3,'LineStyle','-');
plot(nanstd(winddaily_years_shift(upperind,:)),'color',cbrew(1,:),'LineWidth',3,'LineStyle','-');
plot(nanstd(winddaily_years_shift),'color','k','LineWidth',3,'LineStyle','-.');

ylim([0,30]); 
xlim([-5 370]);

set(gca,'fontsize',fsize)
set(gca,'xtick',days_shift(1:end-1),'xticklabel',monslabel);    
ylabel('standard deviation (m/s)','fontsize',fsize);
xlabel('Month','fontsize',fsize);
title(['Zonal wind at 60',char(176),'N and 10hPa'],'fontsize',fsize+2)

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/STD_Observedtozandwind.pdf'];

export_fig(filename,'-pdf');
