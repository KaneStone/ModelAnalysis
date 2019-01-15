%% plot line plots of predictability run low and high ozone
clear all
close all
ClLevel = 'lowCl';
var = 'TOZ';
directory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
TSdirectory = ['/Volumes/ExternalOne/work/data/predruns/','TS','/',ClLevel,'/'];
files = dir([directory,'*.nc']);
TSfiles = dir([TSdirectory,'*.nc']);
lat = -30;
tozlats = [-90,-75];
countlow = 1;
counthigh = 1;
month = 6;
timeperiod = [1955,1975];
noofyears = timeperiod(2)-timeperiod(1)+1;
count = 1;
for i = 1:length(files)    
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    [~,TSdata(i),~] = Read_in_netcdf([TSdirectory,TSfiles(i).name]);
    [~,latind] = min(abs(TSdata(i).lat-lat));
    years(i).y = CCMI_years(TSdata(i).date);
    dateindfirst = find(years(i).y == timeperiod(1),1);
    dateindlast = find(years(i).y == timeperiod(2),1,'last');
    alldata(:,:,count:count+dateindlast-dateindfirst) = data(i).toz(:,:,dateindfirst:dateindlast);
    count = count+size(data(i).toz(:,:,dateindfirst:dateindlast),3);
    
    if i == 1
        for j = 1:length(tozlats)
            [~,tozlatind(j)] = min(abs(tozlats(j)-data(i).lat));
        end
    end
    
    
    % currently not used
    for j = 1:12        
        TSdata_atlat(i,:,j,:) = squeeze(TSdata(i).TS(:,latind,dateindfirst+j-1:12:dateindlast));
        
        for k = 1:length(TSdata(i).lon)
            TS_detrend(i,k,j,:) = detrend(squeeze(TSdata_atlat(i,k,j,:)),'linear');
        end        
        
        data_zonalmean(i,:,j,:) = squeeze(nanmean(data(i).toz(:,1:16,dateindfirst+j-1:12:dateindlast),1));
        data_allmean(i,j,:) = squeeze(nanmean(data_zonalmean(i,:,j,:),2));           
        lowpercentile(i,j) = prctile(data_allmean(i,j,:),10);
        highpercentile(i,j) = prctile(data_allmean(i,j,:),90);
        lowind(i,j).l = find(data_allmean(i,j,:) <= lowpercentile(i,j));
        highind(i,j).h = find(data_allmean(i,j,:) >= highpercentile(i,j));  
               
        composite_low(j,countlow:countlow+length(lowind(i,j).l)-1) = squeeze(data_allmean(i,j,lowind(i,j).l));
        composite_high(j,counthigh:countlow+length(highind(i,j).h)-1) = squeeze(data_allmean(i,j,highind(i,j).h));                
        
    end
    countlow = countlow+length(lowind(i,j).l);
    counthigh = counthigh+length(highind(i,j).h);                      
end

%weightedaverage
alldata_polarmean = weightedaverage(squeeze(nanmean(alldata(:,tozlatind(1):tozlatind(2),:),1)),data(1).lat(tozlatind(1):tozlatind(2)));

for i = 1:12
    alldata_rearrange(i,:) = alldata_polarmean(:,i:12:end);
    alldata_lowpercentile(i) = prctile(alldata_rearrange(i,:),10);
    alldata_highpercentile(i) = prctile(alldata_rearrange(i,:),90);
    alldata_lowind(i,:) = find(alldata_rearrange(i,:) <= alldata_lowpercentile(i));
    alldata_highind(i,:) = find(alldata_rearrange(i,:) >= alldata_highpercentile(i));  
end

composite_low (composite_low == 0) = NaN;
composite_high (composite_high == 0) = NaN;

countlow = 1;
counthigh = 1;


%% singular month
for i = 1:length(files)
    composite_low_singlemonth(:,countlow:countlow+length(lowind(i,month).l)-1) = squeeze(data_allmean(i,:,lowind(i,month).l));
    composite_high_singlemonth(:,counthigh:countlow+length(highind(i,month).h)-1) = squeeze(data_allmean(i,:,highind(i,month).h));
    TScompositelow_singlemonth(:,countlow:countlow+length(lowind(i,month).l)-1) = squeeze(TSdata_atlat(i,:,12,lowind(i,month+2).l));
    TScompositehigh_singlemonth(:,counthigh:countlow+length(highind(i,month).h)-1) = squeeze(TSdata_atlat(i,:,12,highind(i,month+2).h));
    TS_detrend_compositelow_singlemonth(:,countlow:countlow+length(lowind(i,month).l)-1) = squeeze(TS_detrend(i,:,12,lowind(i,month+2).l));
    TS_detrend_compositehigh_singlemonth(:,counthigh:countlow+length(highind(i,month).h)-1) = squeeze(TS_detrend(i,:,12,highind(i,month+2).h));
    countlow = countlow+length(lowind(i,month).l);
    counthigh = counthigh+length(highind(i,month).h);                      
end

composite_low_singlemonth_mean = nanmean(composite_low_singlemonth,2);
composite_high_singlemonth_mean = nanmean(composite_high_singlemonth,2);
TScompositelow_singlemonth_mean = nanmean(TScompositelow_singlemonth,2);
TScompositehigh_singlemonth_mean = nanmean(TScompositehigh_singlemonth,2);
TS_detrend_compositelow_singlemonth_mean = nanmean(TS_detrend_compositelow_singlemonth,2);
TS_detrend_compositehigh_singlemonth_mean = nanmean(TS_detrend_compositehigh_singlemonth,2);

%% plotting
cbrew = cbrewer('qual','Paired',10);

if strcmp(ClLevel,'highCl');
    Cltitle = 'High Chlorine years';
else strcmp(ClLevel,'highCl');
    Cltitle = 'Low Chlorine years';
end
monthtitle = monthnames(month);
% createfig('medium','on');
% fsize = 20;
% phlow = plot(1:12,composite_low_singlemonth,'color',cbrew(7,:),'LineWidth',2);
% hold on
% phhigh = plot(1:12,composite_high_singlemonth,'color',cbrew(9,:),'LineWidth',2);
% phlowmean = plot(1:12,composite_low_singlemonth_mean,'color',cbrew(8,:),'LineWidth',4);
% phhighmean = plot(1:12,composite_high_singlemonth_mean,'color',cbrew(10,:),'LineWidth',4);
% set(gca,'xtick',1:12,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'fontsize',fsize);
% xlabel('Month','fontsize',fsize+2);
% ylabel('Dobson Units','fontsize',fsize+2);
% title([Cltitle,' 10th and 90th percentiles for ',monthtitle],'fontsize',fsize+4);

%% plotting temperature
monthtitle = monthnames(month+2);
createfig('medium','on');
fsize = 20;
% subplot(2,1,1)
% phlow = plot(TSdata(1).lon,TScompositelow_singlemonth,'color',cbrew(7,:),'LineWidth',2);
% hold on
% phhigh = plot(TSdata(1).lon,TScompositehigh_singlemonth,'color',cbrew(9,:),'LineWidth',2);
% phhigh_mean = plot(TSdata(1).lon,TScompositelow_singlemonth_mean,'color',cbrew(8,:),'LineWidth',4);
% phlow_mean = plot(TSdata(1).lon,TScompositehigh_singlemonth_mean,'color',cbrew(10,:),'LineWidth',4);
% set(gca,'xtick',0:30:360,'xticklabel',0:30:360,'fontsize',fsize);
% xlabel('Longitude','fontsize',fsize+2);
% ylabel('Temperature','fontsize',fsize+2);
% title([Cltitle,' 10th and 90th percentiles for ',monthtitle,' at ',num2str(abs(lat)),'S'],'fontsize',fsize+2);
% subplot(2,1,2)


% %% plotting TOZ
% monthtitle = monthnames(month);
% createfig('medium','on');
% fsize = 20;
% 
% phlow = plot(1:12,composite_low_singlemonth,'color',cbrew(7,:),'LineWidth',2);
% hold on
% phhigh = plot(1:12,composite_high_singlemonth,'color',cbrew(9,:),'LineWidth',2);
% phlow_mean = plot(1:12,composite_low_singlemonth_mean,'color',cbrew(8,:),'LineWidth',4);
% phhigh_mean = plot(1:12,composite_high_singlemonth_mean,'color',cbrew(10,:),'LineWidth',4);
% set(gca,'xtick',1:12,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'fontsize',fsize);
% xlabel('Month','fontsize',fsize+2);
% ylabel('DU','fontsize',fsize+2);
% title([monthtitle,' TOZ 10th and 90th percentiles between 60-90S for ',num2str(timeperiod(1)),'-',num2str(timeperiod(2))],'fontsize',fsize+2);
% lh = legend([phhigh_mean,phlow_mean],'High ozone years','Low ozone years');
% set(lh,'location','SouthWest','box','off')
% 
% filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','compositedifferences/Lineplot_',var,...
%             num2str(timeperiod(1)),'-',num2str(timeperiod(2)),'_','10','and','90','percentile_','6090S.pdf'];
% export_fig(filename,'-pdf');

%%
monthtitle = monthnames(month);
createfig('large','on');
fsize = 20;
count = 1;
subplot(2,1,1)
for i = 1:length(files)
    plot(squeeze(alldata_rearrange(month,count:count+dateindlast/12-1))','color',[.6 .6 .6],...
        'LineWidth',2,'Marker','o','MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor',[.6 .6 .6],'MarkerSize',8);
    hold on
    count = count+dateindlast/12;
    a = find(alldata_highind(month,:) <= noofyears*i & alldata_highind(month,:) >= noofyears*(i-1)+1);
    b = find(alldata_lowind(month,:) <= noofyears*i & alldata_lowind(month,:) >= noofyears*(i-1)+1);
    p1 = plot(alldata_highind(month,a)-(i-1)*noofyears,alldata_rearrange(month,alldata_highind(month,a)),...
        'LineStyle','none','Marker','o','MarkerFaceColor',cbrew(10,:),'MarkerEdgeColor',cbrew(10,:),'MarkerSize',10);
    p2 = plot(alldata_lowind(month,b)-(i-1)*noofyears,alldata_rearrange(month,alldata_lowind(month,b)),...
        'LineStyle','none','Marker','o','MarkerFaceColor',cbrew(8,:),'MarkerEdgeColor',cbrew(8,:),'MarkerSize',10);
    
%     hold on
%     plot(lowind(i,month).l,squeeze(data_allmean(i,month,[lowind(i,month).l])),'LineStyle','none',...
%         'Marker','o','MarkerFaceColor',cbrew(8,:),'MarkerEdgeColor',cbrew(8,:),'MarkerSize',10);
%     plot(highind(i,month).h,squeeze(data_allmean(i,month,[highind(i,month).h])),'LineStyle','none',...
%         'Marker','o','MarkerFaceColor',cbrew(10,:),'MarkerEdgeColor',cbrew(10,:),'MarkerSize',10);
end
set(gca,'xtick',1:2:noofyears,'xticklabel',timeperiod(1):2:timeperiod(2),'fontsize',fsize);
xlabel('Year','fontsize',fsize+2);
xlim([0,noofyears+1])
ylabel('Dobson units','fontsize',fsize+2);


title(['October Average 10th and 90th percentiles between ',num2str(abs(tozlats(1))),'S and ',num2str(abs(tozlats(2))),'S'],'fontsize',fsize+2);

subplot(2,1,2)
plot(alldata_rearrange,'color',[.6 .6 .6],'LineWidth',2);
hold on
for i = 1:size(alldata_highind,2)
    
    ph1 = plot(alldata_rearrange(:,alldata_highind(month,i)),'color',cbrew(10,:),'LineWidth',3);
    ph2 = plot(alldata_rearrange(:,alldata_lowind(month,i)),'color',cbrew(8,:),'LineWidth',3);
    
end
lh = legend([ph1,ph2],'90th percentile','10th percentile');
set(lh,'fontsize',fsize+2,'location','NorthWest','box','off')
set(gca,'xtick',1:1:12,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'fontsize',fsize);
xlabel('Month','fontsize',fsize+2);
xlim([0,13]);
ylabel('Dobson units','fontsize',fsize+2);

title(['October Average 10th and 90th percentiles between ',num2str(abs(tozlats(1))),'S and ',num2str(abs(tozlats(2))),'S'],'fontsize',fsize+2);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','compositedifferences/Lineplot_',var,...
            num2str(timeperiod(1)),'-',num2str(timeperiod(2)),'_','10','and','90',...
            'percentile_',num2str(abs(tozlats(1))),num2str(abs(tozlats(2))),'S.pdf'];
export_fig(filename,'-pdf');

% for i = 1:length(files)
%     
% end
% createfig('medium','on');
% subplot(2,1,1);
% phlow = plot(TSdata(1).lon,TScompositelow_singlemonth_mean-TScompositehigh_singlemonth_mean,'color',[117,112,179]./255,'LineWidth',2);
% set(gca,'xtick',0:30:360,'xticklabel',0:30:360,'fontsize',fsize);
% xlabel('Longitude','fontsize',fsize+2);
% ylabel('Temperature','fontsize',fsize+2);
% subplot(2,1,2);
% phlow = plot(1:12,composite_low_singlemonth_mean-composite_high_singlemonth_mean,'color',[117,112,179]./255,'LineWidth',2);
% set(gca,'xtick',1:12,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'fontsize',fsize);
% xlabel('Month','fontsize',fsize+2);
% ylabel('Dobson Units','fontsize',fsize+2);

% %% plotting detrended temperature
% monthtitle = monthnames(month+2);
% createfig('medium','on');
% fsize = 20;
% phlow = plot(TSdata(1).lon,TS_detrend_compositelow_singlemonth,'color',cbrew(7,:),'LineWidth',2);
% hold on
% phhigh = plot(TSdata(1).lon,TS_detrend_compositehigh_singlemonth,'color',cbrew(9,:),'LineWidth',2);
% phhigh_mean = plot(TSdata(1).lon,TS_detrend_compositelow_singlemonth_mean,'color',cbrew(8,:),'LineWidth',4);
% phlow_mean = plot(TSdata(1).lon,TS_detrend_compositehigh_singlemonth_mean,'color',cbrew(10,:),'LineWidth',4);
% set(gca,'xtick',0:30:360,'xticklabel',0:30:360,'fontsize',fsize);
% xlabel('Longitude','fontsize',fsize+2);
% ylabel('Temperature','fontsize',fsize+2);
% title([Cltitle,' 10th and 90th percentiles for ',monthtitle,' (detrend)'],'fontsize',fsize+2);
% 
% createfig('medium','on');
% subplot(2,1,1);
% phlow = plot(TSdata(1).lon,TS_detrend_compositelow_singlemonth_mean-TS_detrend_compositehigh_singlemonth_mean,'color',[117,112,179]./255,'LineWidth',2);
% set(gca,'xtick',0:30:360,'xticklabel',0:30:360,'fontsize',fsize);
% xlabel('Longitude','fontsize',fsize+2);
% ylabel('Temperature','fontsize',fsize+2);
% subplot(2,1,2);
% phlow = plot(1:12,composite_low_singlemonth_mean-composite_high_singlemonth_mean,'color',[117,112,179]./255,'LineWidth',2);
% set(gca,'xtick',1:12,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'fontsize',fsize);
% xlabel('Month','fontsize',fsize+2);
% ylabel('Dobson Units','fontsize',fsize+2);
