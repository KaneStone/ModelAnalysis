

%% import SWOOSH
SStimeperiod = [1984 1999];

[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedeqfillo3q_swoosh-v02.6-198401-201611-latpress-2.5deg-L31.nc');
sfields = fieldnames(SWOOSH);
SWOOSH.(sfields{1}) = cat(3,SWOOSH.(sfields{1}),zeros(size(SWOOSH.(sfields{1}),1),size(SWOOSH.(sfields{1}),2),1));
SWOOSH.(sfields{1}) (SWOOSH.(sfields{1}) == 0) = NaN;
SWOOSHyears = repmat(1984:2015,[12,1]);
SWOOSHyears = [SWOOSHyears(:);ones(12,1)*2016];
SWOOSHextract = permute(SWOOSH.(sfields{1})(:,:,SWOOSHyears >= SStimeperiod(1) & SWOOSHyears <= SStimeperiod(2)),[2,1,3]);

% removing missing lats
SWOOSHextract = SWOOSHextract(:,4:end-3,:);
SWOOSH.lat = SWOOSH.lat(4:end-3);

dateind = SWOOSHyears(1:12:end) >= SStimeperiod(1) & SWOOSHyears(1:12:end)  <= SStimeperiod(2);
dateind_months = SWOOSHyears >= SStimeperiod(1) & SWOOSHyears <= SStimeperiod(2);

%% construct yearly average and extract time period
for i = 1:12
    O3montharrange(:,:,i,:) = SWOOSHextract(:,:,i:12:end);
    O3montharrange_percent(:,:,i,:) = (squeeze(O3montharrange(:,:,i,:)) - nanmean(O3montharrange(:,:,i,:),4))./nanmean(O3montharrange(:,:,i,:),4)*100;
end

O3_montharrange_forreg = permute(O3montharrange_percent(:,:,:),[3,1,2]);
O3yearlyaverage = squeeze(nanmean(O3montharrange,3));
O3yearlyaverage_percent = (O3yearlyaverage - nanmean(O3yearlyaverage,3))./nanmean(O3yearlyaverage,3)*100;
O3yearlyaverage_percent_3y = (O3yearlyaverage - nanmean(O3yearlyaverage(:,:,1:3),3))./nanmean(O3yearlyaverage(:,:,1:3),3)*100;

O3yearlyaverage_foreg = permute(O3yearlyaverage_percent,[3,1,2]);
O3yearlyaverage_foreg_3y = permute(O3yearlyaverage_percent_3y,[3,1,2]);

%% constructing percent anomalies

%% take linear trends
regressfunctions = [ones(sum(dateind),1),[1:sum(dateind)]'];
regressfunctions_months = [ones(sum(dateind_months),1),[1:sum(dateind_months)]'];
for i = 1:size(O3yearlyaverage_foreg,2)
    for j = 1:size(O3yearlyaverage_foreg,3)
        
        b.yearlyaverage(i,j,:) = regress(O3yearlyaverage_foreg(:,i,j),regressfunctions);
        b.monthlytimeseries(i,j,:) = regress(O3_montharrange_forreg(:,i,j),regressfunctions_months);
        b.yearlyaverage_3y(i,j,:) = regress(O3yearlyaverage_foreg_3y(:,i,j),regressfunctions);
        
    end
end

%% plot linear regression
btoplot = permute(cat(3,b.yearlyaverage(:,:,2)*10,b.yearlyaverage_3y(:,:,2)*10,b.monthlytimeseries(:,:,2)*120),[3,2,1]);

prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
    5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
logprestick = log(prestick);

titles = {['Yearly average,',' normalized to entire time series, ', num2str(Stimeperiod(1)),char(8211),num2str(Stimeperiod(2))],...
    ['Yearly average,',' normalized to first 3 years, ', num2str(Stimeperiod(1)),char(8211),num2str(Stimeperiod(2))],...
    ['Monthly time series, normalized to entire time series, ', num2str(Stimeperiod(1)),char(8211),num2str(Stimeperiod(2))]}; 

mtit = {'SWOOSH - combinedanomfillo3q'};

subplotmaps(btoplot,SWOOSH.lat,log(SWOOSH.level),{'div','RdBu'},1,[],16,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
    [-17.5 17.5],22,-90:30:90,-90:30:90,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(1) log(200)],1,'-',0,'');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/WilkaShahSolomonPaper/','combinedanomfillo3q_SWOOSH_trends_',num2str(Stimeperiod(1)),'-',num2str(Stimeperiod(2))];

export_fig(filename,'-pdf');



%% take regression of SD WACCM (chem-only and MAM)
% 
% [bchemonly,predictorschemonly,O3Anomalychemonly] = ...
%     ozoneRegressionTrends(SDWaccmData(2).O3,SDWaccmData(2).U,SDWaccmData(2).solar,SDWaccmData(2).SPE,...
%     SDWaccmData(2).NINO34,SDWaccmData(2).HFS,SDWaccmData(2).HFN,SDWaccmData(2).NO2,SDStimeperiod,1,0,1);%SDWaccmData(2).SPEintpres 

%% 