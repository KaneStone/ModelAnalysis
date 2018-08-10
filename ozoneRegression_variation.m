%% plot tropics variations fro Ensemble members and SWOOSH.
clear all
%% import SWOOSH
Stimeperiod = [1995 2016];
%[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201611-latpress-10deg-L31.nc');
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedanomfillo3q_swoosh-v02.6-198401-201712-latpress-2.5deg-L31.nc');
sfields = fieldnames(SWOOSH);
SWOOSH.(sfields{1}) = cat(3,SWOOSH.(sfields{1}),zeros(size(SWOOSH.(sfields{1}),1),size(SWOOSH.(sfields{1}),2),1));
SWOOSH.(sfields{1}) (SWOOSH.(sfields{1}) == 0) = NaN;
%read in SWOOSH regression functions
[SWOOSHregfun] = ozoneRegression_SWOOSHregfun(Stimeperiod);

SWOOSHyears = repmat(1984:2017,[12,1]);
%SWOOSHyears = [SWOOSHyears(:);ones(12,1)*2016];
SWOOSHextract = permute(SWOOSH.(sfields{1})(:,:,SWOOSHyears >= Stimeperiod(1) & SWOOSHyears <= Stimeperiod(2)),[2,1,3]);

%% import highCl
highCllevel = [SWOOSH.level;[.9,.8,.7,.6,.5,.4,.3,.2,.1]'];
highClTimperiod = [1995 2024];
[highClData,WAClat] = ReadInHighClRegression(highCllevel,highClTimperiod,'highCl','O3');

%% WACCM pressure levels

% %% Read in highcl runs
% directory = ['/Volumes/MyBook/work/data/predruns/O3/','highCl','/zonalmean/'];
% files = dir([directory,'*.nc']);
% for i = 1:length(files)
%     [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);    
%     if i == 1
%         highclyears = 1995:2024;
%         highclyears = repmat(highclyears,[12,1]);
%         highclyears = highclyears(:);        
%     end
%     data(i).O3 = data(i).O3(:,:,highclyears >= highClTimperiod(1) & highclyears <= highClTimperiod(2));
%     data(i).PS = data(i).PS(:,highclyears >= highClTimperiod(1) & highclyears <= highClTimperiod(2));
%     
%     pressure(i).p =  permute(repmat(data(i).hyam*100000,[1,size(data(i).PS)]),[2,1,3]) + ...
%         permute(repmat(data(i).hybm,[1,size(data(i).PS)]),[2,1,3]) .* ...
%         double(permute(repmat(data(i).PS,[1,1,length(data(i).lev)]),[1,3,2])); 
%     
%     % integrate to waccm pressure
%     for j = 1:size(data(i).O3,1)
%         [highClData_wacpres(i).O3(:,j,:),~] = intRegPres(squeeze(data(i).O3(j,:,:)),...
%             squeeze(pressure(i).p(j,:,:))./100,data(i).lev);
%     end
%     
% end

%% testing something

plot(squeeze(SWOOSHextract(4,40,:)))

testdata = squeeze(SWOOSHextract(4,40,:));
testdatamean1 = testdata(1:115) - nanmean(testdata(1:115));
testdatamean2 = testdata(116:end) - nanmean(testdata(116:end));
testdatamean3 = [testdatamean1;testdatamean2];
testdatamean4 = testdata - nanmean(testdata);

%%
lats = [-30 30];
%lats = [-80 -50]; %
lev = [25 100];%lev = [1,3]; % 
latind = WAClat >= lats(1) & WAClat <= lats(2);
latindSWOOSH = SWOOSH.lat >= lats(1) & SWOOSH.lat <= lats(2);
%[~,levind] = min(abs(highCllevel - lev));
levind = highCllevel >= lev(1) & highCllevel <= lev(2);
%slevind = highCllevel >= lev(1) & highCllevel <= lev(2);
levind2 = find(levind);
%% 
clearvars highclO3 highclO3_Vertave highclO3vert
for i = 1:9       
    for j = 1:sum(levind)        
        highclO3(j,i,:) = weightedaverage(squeeze(highClData(i).O3(levind2(j),latind,:)),WAClat(latind));    
    end
    highclO3_Vertave(i,:) = squeeze(nanmean(highclO3(:,i,:),1));
    for j = 1:size(highClData(i).O3,1)
        highclO3vert(i,j,:) = weightedaverage(squeeze(highClData(i).O3(j,latind,:)),WAClat(latind));            
    end
%     for j = 1:size(data(i).lev)
%         highclO3vert2(i,j,:) = weightedaverage(squeeze(highClData_wacpres(i).O3(j,latind,:)),WAClat(latind));    
%     end
end

%%

for i = 1:12
    highclO3_anom(:,i,:) = (highclO3_Vertave(:,i:12:end) - nanmean(highclO3_Vertave(:,i:12:end),2))./nanmean(highclO3_Vertave(:,i:12:end),2)*100;
end

highclO3_anom_ts = highclO3_anom(:,:);

highclO3_anom_ts (highclO3_anom_ts >= 30) = NaN;
highclO3_anom_ts (highclO3_anom_ts <= -30) = NaN;
%% SWOOSH
clearvars SWOOSHwa SWOOSHwa_Vertave SWOOSHwa_vert

slatind = SWOOSH.lat >= lats(1) & SWOOSH.lat <= lats(2);

for j = 1:sum(levind)    
    SWOOSHwa(j,:) = weightedaverage(squeeze(SWOOSHextract(levind2(j),slatind,:)),SWOOSH.lat(slatind));
end
SWOOSHwa_Vertave = squeeze(nanmean(SWOOSHwa,1));

for i = 1:size(SWOOSHextract,1)
    SWOOSHwa_vert(i,:) = weightedaverage(squeeze(SWOOSHextract(i,slatind,:)),SWOOSH.lat(slatind));
end

for i = 1:12
    SWOOSH_anom(:,i,:) = (SWOOSHwa_Vertave(:,i:12:end) - nanmean(SWOOSHwa_Vertave(:,i:12:end),2))./nanmean(SWOOSHwa_Vertave(:,i:12:end),2)*100;
end

SWOOSH_anom_ts = SWOOSH_anom(:,:);

%%
plot_vert = 1;
if plot_vert
    
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);
    
    [fig] = createfig('medium','on');
    plot(nanmean(SWOOSHwa_vert,2),log(SWOOSH.level),'-o','LineWidth',3);
    hold on
    highclave = squeeze(nanmean(highclO3vert,1));
    highclave2 = nanmean(highclave,2);

    % highcl_waclev_ave = squeeze(nanmean(highclO3vert2,1));
    % highcl_waclev_ave2 = nanmean(highcl_waclev_ave,2);

    plot(highclave2*1e6,log(highCllevel),'-o','LineWidth',3);
    % plot(highcl_waclev_ave2*1e6,log(data(1).lev),'-o','LineWidth',3);
    set(gca,'ytick',fliplr(logprestick),'yticklabel',fliplr(presticklabel))
    set(gca,'YDir','reverse','fontsize',18)
    ylabel('Pressure (hPa)','fontsize',20);
    xlabel('ppmv','fontsize',20);
    lh = legend('SWOOSH','Ensemble');
    set(lh,'fontsize',20,'box','off');
    title(['Ozone profile averaged between ',num2str(abs(lats(1))),'S and ',num2str(abs(lats(2))),...
    'N, and ',num2str(lev(2)),' and ', num2str(lev(1)) ,' hPa'],'fontsize',22);
    ylim([log(.1),log(300)])
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/Variation/',...
    'Compare_',num2str(abs(lats(1))),'S-',num2str(abs(lats(2))),'N_',num2str(lev(2)),'-',num2str(lev(1)),'_',num2str(Stimeperiod(1)),'-',num2str(Stimeperiod(2))];

    export_fig(filename,'-pdf');
end
%%
cbrew = cbrewer('qual','Set1',10);

fig = figure;
set(fig,'color','white','position',[100 100 1000 500]);
pe = plot(highclO3_anom_ts(:,:)','color',[.7 .7 .7],'LineWidth',2);
hold on

ps = plot(SWOOSH_anom_ts','k','LineWidth',3);

xlim([-11 370]);
%ylim([-11 11]);

set(gca,'fontsize',18,'xtick',1:36:400,'xticklabel',1995:3:2025);
xlabel('Year','fontsize',20);
ylabel('Ozone mixing ratio anomaly (%)','fontsize',20);
title(['Time series averaged between ',num2str(abs(lats(1))),'{\circ}S and ',num2str(abs(lats(2))),...
    '{\circ}N, and ',num2str(lev(2)),' and ', num2str(lev(1)) ,' hPa'],'fontsize',22);

lh = legend([pe(1),ps],'Ensemble members','SWOOSH');
set(lh,'box','off','fontsize',20);

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/Variation/',...
    'TimeSeriesCompare_',num2str(abs(lats(1))),'S-',num2str(abs(lats(2))),'N_',num2str(lev(2)),'-',num2str(lev(1)),'_',num2str(Stimeperiod(1)),'-',num2str(Stimeperiod(2))];

    export_fig(filename,'-pdf');
%%

for  i = 1:9
    bhcl(i,:) = regress(highclO3_anom_ts(i,37:end)',[ones(length(highclO3_anom_ts(i,37:end)),1),[1:length(highclO3_anom_ts(1,37:end))]']);
end
bswoosh = regress(SWOOSH_anom_ts(37:end)',[ones(length(SWOOSH_anom_ts(37:end)),1),[1:length(SWOOSH_anom_ts(37:end))]']);

bhcl_ppd = bhcl(:,2)*120;
bswoosh_ppd = bswoosh(2)*120;

%%

std(SWOOSH_anom_ts(1:end-1))
std(highclO3_anom_ts,0,2)

