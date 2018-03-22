%predruns_SAMaccounting

clear all 
tozvar = 'toz';
var = 'PSL';
%area = '6090S';
detrend = 1;
detrend_ozone = 1;
contourplots = 1;
individual_contourplots = 0;
lineplots = 1;
percentile = 50;
lats = [-90,-75];
tozmonth = 9;
indmonth = 11;
varmonth = [12,1,2]; % Can be any number of months in the year (e.g. [12,1,2] for austral summer)
shortnames = 0;
sig = .05;
ENSO = 0;
removeENSO = 0;

%% Read in lowcl variable

ClLevel = 'lowCl';
timeperiodlow = [1955,1976];%[1955,1975]

vardirectory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/'];
varfilespast = dir([vardirectory,'*.nc']);

[data.lowcl,years.lowcl,composite.lowcl,dataMonthArrange.lowcl]...
    = predruns_ReadInlayer(vardirectory,varfilespast,var,timeperiodlow,lats,detrend);

%% Read in highcl variable

ClLevel = 'highCl';
timeperiodhigh = [1995,2016];%[1955,1975]

vardirectory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl]...
    = predruns_ReadInlayer(vardirectory,varfiles,var,timeperiodhigh,lats,detrend);

%% Read in lowGHG variable

ClLevel = 'lowGHG';
timeperiodhigh = [1995,2016];%[1955,1975]

vardirectory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/'];
varfileslowGHG = dir([vardirectory,'*.nc']);

[data.lowGHG,years.lowGHG,composite.lowGHG,dataMonthArrange.lowGHG]...
    = predruns_ReadInlayer(vardirectory,varfileslowGHG,var,timeperiodhigh,lats,detrend);

%% Read in TOZ highcl and take percentiles
ClLevel = 'highCl';
tozdates = [1995,2015];
directory = ['/Volumes/MyBook/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);

[pct_highcl,tozextract.highcl] = predruns_varPercentiles(toz_composite.highcl.montharrange,toz_dataMonthArrange.highcl,...
    tozmonth,percentile,length(tozfiles));

%% Read in TOZ lowcl and take percentiles
ClLevel = 'lowCl';
tozpastdates = [1955,1975];
directory = ['/Volumes/MyBook/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfilespast = dir([directory,'*.nc']);
[toz_data.lowcl,toz_years.lowcl,toz_varweighted.lowcl,toz_composite.lowcl,toz_dataMonthArrange.lowcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfilespast,tozvar,tozpastdates,lats,detrend_ozone);

[pct_lowcl,tozextract.lowcl] = predruns_varPercentiles(toz_composite.lowcl.montharrange,toz_dataMonthArrange.lowcl...
    ,tozmonth,percentile,length(tozfilespast));

%% Read in TOZ lowGHG and take percentiles
ClLevel = 'lowGHG';
tozpastdates = [1995,2015];
directory = ['/Volumes/MyBook/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfileslowGHG = dir([directory,'*.nc']);
[toz_data.lowGHG,toz_years.lowGHG,toz_varweighted.lowGHG,toz_composite.lowGHG,toz_dataMonthArrange.lowGHG] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfileslowGHG,tozvar,tozpastdates,lats,detrend_ozone);

[pct_lowGHG,tozextract.lowGHG] = predruns_varPercentiles(toz_composite.lowGHG.montharrange,toz_dataMonthArrange.lowGHG...
    ,tozmonth,percentile,length(tozfilespast));

%% calculateSAM

fields = fieldnames(dataMonthArrange); 
SAMlats = [-65, -40];

longitude = data(1).lowcl.lon;
latitude = data(1).lowcl.lat;

[~,SAMlat_index(1)] = min(abs(SAMlats(1) - latitude));
[~,SAMlat_index(2)] = min(abs(SAMlats(2) - latitude));

for i = 1:size(fields)
    zonalmean.(fields{i}) = nanmean(dataMonthArrange.(fields{i}),5);
    noyears = size(zonalmean.(fields{i}),3);
    for j = 1:12
        SAM1.(fields{i})(:,j,:) = (zonalmean.(fields{i})(:,j,:,SAMlat_index(1)) - ...
            repmat(nanmean(zonalmean.(fields{i})(:,j,:,SAMlat_index(1)),3),[1,1,noyears]))./...
            repmat(std(zonalmean.(fields{i})(:,j,:,SAMlat_index(1)),1,3),[1,1,noyears]);
        
        SAM2.(fields{i})(:,j,:) = (zonalmean.(fields{i})(:,j,:,SAMlat_index(2)) - ...
            repmat(nanmean(zonalmean.(fields{i})(:,j,:,SAMlat_index(2)),3),[1,1,noyears]))./...
            repmat(std(zonalmean.(fields{i})(:,j,:,SAMlat_index(2)),1,3),[1,1,noyears]);
        
        SAMindex.(fields{i})(:,j,:) = SAM2.(fields{i})(:,j,:) - SAM1.(fields{i})(:,j,:);
        
    end
end


%% SAM toz corr
tozlat = -65;
[~,tozlatindex] = min(abs(tozlat - latitude));
for i = 1:size(SAMindex.highcl,1)
    toz_normanom = (squeeze(toz_dataMonthArrange.highcl(i,10,:)) - nanmean(squeeze(toz_dataMonthArrange.highcl(i,10,:))))./...
        std(squeeze(toz_dataMonthArrange.highcl(i,10,:)),1);
    r(i) = corr(nanmean([squeeze(SAMindex.highcl(i,12,1:end-1)),squeeze(SAMindex.highcl(i,1,2:end)),squeeze(SAMindex.highcl(i,2,2:end))]')',toz_normanom);
    tozsw(i,:,:) = squeeze(toz_data.highcl(i).toz(:,tozlatindex,tozmonth :12:end));
    toz2(i,:,:,:) = squeeze(toz_data.highcl(i).toz(:,1:16,tozmonth :12:end));
end
    

%% read in Bodeker

%% READ in Bodeker 
%tcolat = [75,90];
tcolat = -65;

[~,BSdata,~] = Read_in_netcdf('/Volumes/MyBook/work/data/BodekerScientific/TCO/Bodeker_TCO_monavg.nc');

BSyears = 1980:2016;
BSyear_vector = repmat(BSyears,12,1);
BSdata.years = [BSyear_vector(:);ones(length(BSdata.time)-length(BSyear_vector(:)),1)*max(BSyear_vector(:))+1];

BSdateindex(1) = find(BSdata.years == 1980,1,'first');
BSdateindex(2) = find(BSdata.years == 2016,1,'last');


%latindex_bs = find(BSdata.lat >= tcolat(1) & BSdata.lat <= tcolat(2));
[~,latindex_bs] = min(abs(BSdata.lat - tcolat));
%toz_zm = squeeze(nanmean(BSdata.tco(:,:,BSdateindex(1)+tozmonth:12:BSdateindex(2)),1));
BStoz = squeeze(BSdata.tco(:,latindex_bs,BSdateindex(1)+tozmonth -1:12:BSdateindex(2)));
%toz_am = squeeze(nanmean(toz_zm(latindex_bs,:),1));
%toz_am_QBO = squeeze(nanmean(toz_zm(latindex_bsQBO,:),1));
[BSminlongvalue,BSminloc] = min(BStoz,[],1);
BSminlongvalue (BSminlongvalue < 100) = NaN;
[BSmaxlongvalue,BSmaxloc] = max(BStoz,[],1);
BSmaxlongvalue (BSmaxlongvalue > 400) = NaN;
BSminlong = BSdata.lon(BSminloc);
BSmaxlong = BSdata.lon(BSmaxloc);
%% 

cbrew = cbrewer('qual','Paired',10);
[minlongvalue,minloc] = min(tozsw,[],2);
[maxlongvalue,maxloc] = max(tozsw,[],2);

[minlongvalue2,minloc2] = min(squeeze(min(toz2,[],2)),[],2);
[maxlongvalue2,maxloc2] = max(squeeze(max(toz2,[],2)),[],2);

minlongvalue = squeeze(minlongvalue);
maxlongvalue = squeeze(maxlongvalue);
minloc = squeeze(minloc);
maxloc = squeeze(maxloc);

minlongvalue2 = squeeze(minlongvalue2);
maxlongvalue2 = squeeze(maxlongvalue2);
minloc2 = squeeze(minloc2);
maxloc2 = squeeze(maxloc2);


minlong = longitude(minloc);
maxlong = longitude(maxloc);

for i = 1:size(minlong,1)
    for j = 1:size(minlong,2)
        if minlong(i,j) < 180
            minlong(i,j) = minlong(i,j)+360;
        end
%         if maxlong(i,j) < 180
%             maxlong(i,j) = maxlong(i,j)+360;
%         end
    end
end


difflongvalue = maxlongvalue - minlongvalue;

createfig('medium','on');
plot(repmat(1995:2024,size(minlong,1),1)',minlong','color',cbrew(1,:),'LineWidth',2);
hold on
plot(repmat(1995:2024,size(minlong,1),1)',nanmean(minlong)','color',cbrew(2,:),'LineWidth',4);
plot(repmat(1995:2024,size(minlong,1),1)',maxlong','color',cbrew(5,:),'LineWidth',2);
plot(repmat(1995:2024,size(minlong,1),1)',nanmean(maxlong)','color',cbrew(6,:),'LineWidth',4);

createfig('medium','on');
hold on

pb(1) = plot(1980:2016,BSminlongvalue','color','k','LineWidth',2);
plot(1980:2016,BSmaxlongvalue','color','k','LineWidth',2);
pb(2) = plot(1980:2016,min(squeeze(min(BSdata.tco(:,1:30,10:12:end)))),'--k','LineWidth',2);
plot(1980:2016,max(squeeze(max(BSdata.tco(:,1:30,10:12:end)))),'--k','LineWidth',2);

plot(repmat(1995:2024,size(minlong,1),1)',minlongvalue','color',cbrew(1,:),'LineWidth',2);
ph(1) = plot(1995:2024,nanmean(minlongvalue)','color',cbrew(2,:),'LineWidth',4);

plot(repmat(1995:2024,size(minlong,1),1)',minlongvalue2','color',cbrew(5,:),'LineWidth',2);
ph(2) = plot(1995:2024,nanmean(minlongvalue2)','color',cbrew(6,:),'LineWidth',4);

plot(repmat(1995:2024,size(minlong,1),1)',maxlongvalue','color',cbrew(1,:),'LineStyle','--','LineWidth',2);
ph(3) = plot(1995:2024,nanmean(maxlongvalue)','color',cbrew(2,:),'LineStyle','--','LineWidth',4);

plot(repmat(1995:2024,size(minlong,1),1)',maxlongvalue2','color',cbrew(5,:),'LineStyle','--','LineWidth',2);
ph(4) = plot(1995:2024,nanmean(maxlongvalue2)','color',cbrew(6,:),'LineStyle','--','LineWidth',4);

legend([pb,ph],['Bodeker-NIWA min at ',num2str(tcolat)],'Bodeker-NIWA min over polar cap',['WACCM min at ',num2str(tcolat)],...
    'WACCM min over polar cap',['WACCM max at ',num2str(tcolat)],'WACCM max over polar cap')

xlabel('Year','fontsize',20);
ylabel('DU','fontsize',20);
box on
export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/DifferenceInLonMaxMinTurnaroundTimes',monthnames(tozmonth,1,0),'.pdf'],'-pdf');

createfig('medium','on');
plot(repmat(1995:2024,size(minlong,1),1)',difflongvalue','color',cbrew(1,:),'LineWidth',4);
hold on
plot(repmat(1995:2024,size(minlong,1),1)',nanmean(difflongvalue)','color',cbrew(2,:),'LineWidth',4);
