%% predruns investigating sea surface temperatures
clear variables
% read in sea ice
inputs.var = 'ICEFRAC';
inputs.ClLevel = 'highCl';
inputs.timeperiodvar = [1995,2024];
inputs.detrend = 0;
inputs.varmonth = 3;
% Read in surface temperature or other similar variable
vardirectory = ['/Volumes/MyBook/work/data/predruns/',inputs.var,'/',inputs.ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[surfacedata.(inputs.ClLevel).ice.data,surfacedata.(inputs.ClLevel).years,surfacedata.(inputs.ClLevel).ice.composite,...
    surfacedata.(inputs.ClLevel).ice.dataMonthArrange,surfacedata.(inputs.ClLevel).ice.dataMonthArrangeMean]...
    = predruns_ReadInlayer(vardirectory,varfiles,inputs.var,inputs.timeperiodvar,inputs.detrend);

inputs.var = 'TS';

vardirectory = ['/Volumes/MyBook/work/data/predruns/',inputs.var,'/',inputs.ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[surfacedata.(inputs.ClLevel).TS.data,surfacedata.(inputs.ClLevel).years,surfacedata.(inputs.ClLevel).TS.composite,...
    surfacedata.(inputs.ClLevel).TS.dataMonthArrange,surfacedata.(inputs.ClLevel).TS.dataMonthArrangeMean]...
    = predruns_ReadInlayer(vardirectory,varfiles,inputs.var,inputs.timeperiodvar,inputs.detrend);

%%

TSanomaly = surfacedata.highCl.TS.dataMonthArrange - nanmean(surfacedata.highCl.TS.dataMonthArrange,1);
ICEanomaly = surfacedata.highCl.ice.dataMonthArrange - nanmean(surfacedata.highCl.ice.dataMonthArrange,1);

TSanomalyatmonth = squeeze(nanmean(TSanomaly(:,inputs.varmonth,:,:,:),3));
ICEanomalyatmonth = squeeze(nanmean(ICEanomaly(:,inputs.varmonth,:,:,:),3));

%% plotting
longitude = surfacedata.highCl.ice.data(1).lon;
latitude = surfacedata.highCl.ice.data(1).lat;

for i = 1:9
    subplotmaps(permute(TSanomalyatmonth(i,:,:),[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,{['No. ',num2str(i)]},'Longitude','Latitude','','on',...
         [-3,3],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/anomalies/SST/',['No.',sprintf('%02d',i)],monthnames(inputs.varmonth,1,1)];
    export_fig(filename,'-png');
    
    subplotmaps(permute(ICEanomalyatmonth(i,:,:),[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,{['No. ',num2str(i)]},'Longitude','Latitude','','on',...
         [-.1,.1],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/anomalies/ICE/',['No.',sprintf('%02d',i)],monthnames(inputs.varmonth,1,1)];
    export_fig(filename,'-png');
    
end
