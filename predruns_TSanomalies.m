% plot average El-nino surface temperature anomalies

clear all

%% Read in TOZ highcl and take percentiles
% years to estimate = 2016:2023

ClLevel = 'highCl';
tozdates = [1995,2024];
directory = ['/Volumes/MyBook/work/data/predruns/','TOZ','/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,'toz',tozdates,[-30,-10],0);

[toz_data2.highcl,toz_years2.highcl,toz_varweighted2.highcl,toz_composite2.highcl,toz_dataMonthArrange2.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,'toz',tozdates,[-90,-75],0);
%% Read in highcl Surface temperature
ClLevel = 'highCl';
timeperiodhigh = [1995,2024];%[1955,1975]

vardirectory = ['/Volumes/MyBook/work/data/predruns/','TS','/',ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl]...
    = predruns_ReadInlayer(vardirectory,varfiles,'TS',timeperiodhigh,[-90,-75],1);


lats = data.highcl.lat;
lons = data.highcl.lon;
%% calculate ENSO
tozmonth = 11;
latlimits = [-5 5];
lonlimits = [190 240];

latindex = lats >= latlimits(1) & lats <= latlimits(2);
lonindex = lons >= lonlimits(1) & lons <= lonlimits(2);


NINO_mn = squeeze(nanmean(dataMonthArrange.highcl(:,tozmonth,:,latindex,lonindex),5));
NINO_mn = squeeze(nanmean(NINO_mn,3));
NINOallmean = nanmean(NINO_mn,2);
NINOstd = std(NINO_mn,1,2);
NINO34 = detrend(((NINO_mn - NINOallmean)./NINOstd)');    

NINOnegative = NINO34 <= -.75;
NINOpositive = NINO34 >= .75;

%% correlate QBO with ENSO

for i = 1:9
    qboesbo_r(i,:) = corr(NINO34(:,i),squeeze(toz_dataMonthArrange.highcl(i,11,:)));
    qboesbo_r2(i,:) = corr(NINO34(:,i),squeeze(toz_dataMonthArrange2.highcl(i,11,:)));
end
figure
%plot(qboesbo_r)
hold on
plot(qboesbo_r2);
%% mean surface temperature



for i = 1:length(varfiles)
    for j = 1:size(data.highcl(i).TS,1)
        for k = 1:size(data.highcl(i).TS,2)
            
            meanTS(j,k,:,i) = detrend(squeeze(nanmean(cat(4,data.highcl(i).TS(j,k,12:12:end-12),data.highcl(i).TS(j,k,13:12:end),data.highcl(i).TS(j,k,14:12:end)),4))) + ...
                nanmean(nanmean(cat(4,data.highcl(i).TS(j,k,12:12:end-12),data.highcl(i).TS(j,k,13:12:end),data.highcl(i).TS(j,k,14:12:end)),4),3);

            %meanTS(j,k,:,i) = detrend(squeeze(data.highcl(i).TS(j,k,12:12:end-12))) + ...
            %    nanmean(data.highcl(i).TS(j,k,12:12:end-12),3);
        
        end
    end
    meanTS_negative(i,:,:) = nanmean(meanTS(:,:,NINOnegative(1:29,i),i),3) - nanmean(meanTS(:,:,:,i),3);
    meanTS_positive(i,:,:) = nanmean(meanTS(:,:,NINOpositive(1:29,i),i),3) - nanmean(meanTS(:,:,:,i),3);       
end

% meanTS_negative = permute(meanTS_negative,[3,1,2]);
% meanTS_positive = permute(meanTS_positive,[3,1,2]);

% %% plotting
% for i = 1:size(meanTS_negative)
% 
%     subplotmaps(meanTS_negative(i,:,:),lons,lats,{'div','RdBu'},1,[],16,{['Negative NINO anomalies - No.',num2str(i)]},'Longitude','latitude','Correlation','on',...
%         [-1 1],18,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],'',1,[0 360],[-90 90],0,'none',1,'Miller Cylindrical');
% 
% 
% %     filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/Ind_','No.',sprintf('%02d',i),monthnames(tozmonth,0,0),'toz','_TS_detrend',...
% %         num2str(thigh(1)),'-',num2str(thigh(2)),'_',num2str(abs(lats2(1))),'-',...
% %         num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,1),'_',ENSOext];
%         
%     %export_fig(filename,'-png')
% end


%%
for i = [1:9]%1:size(meanTS_negative)

    subplotmaps(meanTS_negative(i,:,:),lons,lats,{'div','RdBu'},1,[],16,{['Positive NINO anomalies - No.',num2str(i)]},'Longitude','latitude','Correlation','on',...
        [-1 1],18,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],'',1,[0 360],[-90 90],0,'none',1,'Miller Cylindrical');


%     filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/Ind_','No.',sprintf('%02d',i),monthnames(tozmonth,0,0),'toz','_TS_detrend',...
%         num2str(thigh(1)),'-',num2str(thigh(2)),'_',num2str(abs(lats2(1))),'-',...
%         num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,1),'_',ENSOext];
        
    %export_fig(filename,'-png')
end
