% plot NION34

ERA = load('/Volumes/ExternalOne/work/data/predruns/output/NINO34/ERA-Interim_1980_2016');
HighCL = load('/Volumes/ExternalOne/work/data/predruns/output/NINO34/highcl_1995_2014');

%%
comparelines = 0;
if comparelines
    figure
    plot(nanmean(HighCL.NINO34all.highcl([1,2,7],:),1)','LineWidth',2);
    hold on
    plot(ERA.NINO34all(181:end),'k','LineWidth',2);
    figure
    %plot(HighCL.NINO34all.highcl([3,5,9],:)')
    plot(nanmean(HighCL.NINO34all.highcl([3,5,9],:),1)','LineWidth',2);
    hold on
    plot(ERA.NINO34all(181:end),'k','LineWidth',2);
    figure
    plot(nanmean(HighCL.NINO34all.highcl([4,6,8],:),1)','LineWidth',2);
    %plot(HighCL.NINO34all.highcl([4,6,8],:)')
    hold on
    plot(ERA.NINO34all(181:end),'k','LineWidth',2);
end

%% Read in Bodeker
tozmonth = 11;
tcolat = [-90,-75];
timeperiodbs = [1980,2016];
tcolon = 26;
tcolat_QBO = [-30,-10];
[~,BSdata,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/BodekerScientific/TCO/Bodeker_TCO_monavg.nc');

BSyears = 1980:2016;
BSyear_vector = repmat(BSyears,12,1);
BSdata.years = [BSyear_vector(:);ones(length(BSdata.time)-length(BSyear_vector(:)),1)*max(BSyear_vector(:))+1];

BSdateindex(1) = find(BSdata.years == timeperiodbs(1),1,'first');
BSdateindex(2) = find(BSdata.years == timeperiodbs(2),1,'last');

latindex_bsQBO = find(BSdata.lat >= tcolat_QBO(1) & BSdata.lat <= tcolat_QBO(2));
latindex_bs = find(BSdata.lat >= tcolat(1) & BSdata.lat <= tcolat(2));
%toz_zm = squeeze(nanmean(BSdata.tco(:,:,BSdateindex(1)+tozmonth:12:BSdateindex(2)),1));
toz_zm = detrend(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,BSdateindex(1)+tozmonth-1:12:BSdateindex(2)),1)),BSdata.lat(latindex_bs)));
toz_zmQBO = detrend(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bsQBO,BSdateindex(1)+tozmonth-1:12:BSdateindex(2)),1)),BSdata.lat(latindex_bsQBO)));

%% Read in TOZ highcl and take percentiles
tozvar = 'toz';

detrend_ozone = 1;

lats = [-90,-75];

ClLevel = 'highCl';
tozdates = [1995,2016];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);

longitude = toz_data(1).highcl.lon;
latitude = toz_data(1).highcl.lat;

%% equator latitudes
lats = [-30,-10];
[~,~,~,~,eq_toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);

%%
% take the mean of all 12,1,2 elnino
ENS = nanmean(cat(3,squeeze(HighCL.NINOmonth.highcl(:,12,1:end-1)),squeeze(HighCL.NINOmonth.highcl(:,1,2:end)),...
    squeeze(HighCL.NINOmonth.highcl(:,2,2:end))),3);

ENSqbo = squeeze(HighCL.NINOmonth.highcl(:,11,1:end));
    
for i = 1:9
    rQBO(i) = corr(detrend(squeeze(eq_toz_dataMonthArrange.highcl(i,11,1:end))),detrend(ENSqbo(i,:))');
    r(i) = corr(detrend(squeeze(toz_dataMonthArrange.highcl(i,11,1:end-1))),detrend(ENS(i,:))');
    r2(i) = corr(detrend(squeeze(toz_dataMonthArrange.highcl(i,11,1:end))),detrend(ENSqbo(i,:))');
end

ErQBO = corr(ERA.NINO34all(11:12:end),toz_zmQBO');
Er = corr(ERA.NINO34all(13:12:end),toz_zm(1:end-1)');

plot([r,Er])
hold on
plot(r2)
plot(rQBO)

