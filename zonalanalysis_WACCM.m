%% Investigate predruns and WACCM ccmi zonal wave-1 and phase

tozvar = 'toz';
lats = [-90 -60];
detrend_ozone = 0;
sig = .05;
%% Read in predruns TOZ

%% Read in TOZ highcl and take percentiles
ClLevel = 'highCl';
tozdates = [1995,2024];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);


%% Read in TOZ lowcl and take percentiles
ClLevel = 'lowCl';
tozpastdates = [1955,1979];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfilespast = dir([directory,'*.nc']);
[toz_data.lowcl,toz_years.lowcl,toz_varweighted.lowcl,toz_composite.lowcl,toz_dataMonthArrange.lowcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfilespast,tozvar,tozpastdates,lats,detrend_ozone);

latitudes = toz_data.lowcl.lat;
longitudes = toz_data.lowcl.lon;

%% plot ozone at latitude lines

zonallat = -60;
month = 10;
[~,latind] = min(abs(zonallat - latitudes));
fieldnames = fields(toz_data);
colors = {'r','b'};
plotlong = 0;
for j = 1:length(fieldnames)
    for i = 1:length(toz_data.(fieldnames{j}))
        toz_data.(fieldnames{j})(i).toz_atzonallat = squeeze(toz_data.(fieldnames{j})(i).toz(:,latind,month:12:end));
        [toz_data.(fieldnames{j})(i).minvalues, minind] = min(toz_data.(fieldnames{j})(i).toz_atzonallat);
        [toz_data.(fieldnames{j})(i).maxvalues, maxind] = max(toz_data.(fieldnames{j})(i).toz_atzonallat);
        toz_data.(fieldnames{j})(i).minlongitude = longitudes(minind)';
        toz_data.(fieldnames{j})(i).maxlongitude = longitudes(maxind)';
        toz_data.(fieldnames{j})(i).amplitude = toz_data.(fieldnames{j})(i).maxvalues - toz_data.(fieldnames{j})(i).minvalues;
        if plotlong
            plot(toz_data.(fieldnames{j})(i).toz_atzonallat,colors{j});
            hold on
        end
    end
end


%% look at correlations to surface temperature

var = 'T';
level = '400hPa';
varlayer = 0;

%%  Read in lowcl variable

ClLevel = 'lowcl';
timeperiodhigh = [1955,1979];%[1955,1975]

if varlayer
    vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
    varfilespast = dir([vardirectory,'*.nc']);
    [data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl]...
        = predruns_ReadInlayer(vardirectory,varfilespast,var,timeperiodhigh,lats,1);
else
    vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/','500hPa','/'];
    varfilespast = dir([vardirectory,'*.nc']);
    [data.highcl,years.highcl,data_at_level.highcl,composite.highcl,dataMonthArrange.highcl]...
        = predruns_ReadIn4d(vardirectory,varfilespast,var,timeperiodhigh,400,1);
end

%% Read in highcl variable

ClLevel = 'lowCl';
timeperiodlow = [1955,1979];%[1955,1975]

if varlayer
    vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
    varfilespast = dir([vardirectory,'*.nc']);
    [data.lowcl,years.lowcl,composite.lowcl,dataMonthArrange.lowcl]...
        = predruns_ReadInlayer(vardirectory,varfilespast,var,timeperiodlow,lats,1);
else
    vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/',level,'/'];
    varfilespast = dir([vardirectory,'*.nc']);
    [data.lowcl,years.lowcl,data_at_level.lowcl,composite.lowcl,dataMonthArrange.lowcl]...
        = predruns_ReadIn4d(vardirectory,varfilespast,var,timeperiodlow,400,1);
end

%% take correlations 
varmonth = 10;
%abc = squeeze(dataMonthArrange.lowcl(:,12,1,1,:))'

for i = 1:length(latitudes)    
    [rmaxlon.highcl(i,:),pmaxlon.highcl(i,:)] = corr(squeeze(composite.highcl.montharrange(varmonth,:,i,:))',[toz_data.highcl(:).maxlongitude]');
    [rmaxlon.lowcl(i,:),pmaxlon.lowcl(i,:)] = corr(squeeze(composite.lowcl.montharrange(varmonth,:,i,:))',[toz_data.lowcl(:).maxlongitude]');
    [rminlon.highcl(i,:),pminlon.highcl(i,:)] = corr(squeeze(composite.highcl.montharrange(varmonth,:,i,:))',[toz_data.highcl(:).minlongitude]');
    [rminlon.lowcl(i,:),pminlon.lowcl(i,:)] = corr(squeeze(composite.lowcl.montharrange(varmonth,:,i,:))',[toz_data.lowcl(:).minlongitude]');
    [ramp.highcl(i,:),pamp.highcl(i,:)] = corr(squeeze(composite.highcl.montharrange(varmonth,:,i,:))',[toz_data.highcl(:).amplitude]');
    [ramp.lowcl(i,:),pamp.lowcl(i,:)] = corr(squeeze(composite.lowcl.montharrange(varmonth,:,i,:))',[toz_data.lowcl(:).amplitude]');
end

%% take individual correlations
for k = 1:length(fieldnames)
    for j = 1:length(toz_data.(fieldnames{k}))
        for i = 1:length(latitudes)    
            [rmaxlon_ind.(fieldnames{k})(j,:,i),pmaxlon_ind.(fieldnames{k})(j,:,i)] = corr(squeeze(dataMonthArrange.(fieldnames{k})(j,varmonth,:,i,:)),[toz_data.(fieldnames{k})(j).maxlongitude]');            
            [rminlon_ind.(fieldnames{k})(j,:,i),pminlon_ind.(fieldnames{k})(j,:,i)] = corr(squeeze(dataMonthArrange.(fieldnames{k})(j,varmonth,:,i,:)),[toz_data.(fieldnames{k})(j).minlongitude]');            
            [ramp_ind.(fieldnames{k})(j,:,i),pamp_ind.(fieldnames{k})(j,:,i)] = corr(squeeze(dataMonthArrange.(fieldnames{k})(j,varmonth,:,i,:)),[toz_data.(fieldnames{k})(j).amplitude]');
        end
    end
    pmaxlon_ind.(fieldnames{k}) (pmaxlon_ind.(fieldnames{k}) <= sig) = 0;
    pmaxlon_ind.(fieldnames{k}) (pmaxlon_ind.(fieldnames{k}) > sig) = 1;
    pminlon_ind.(fieldnames{k}) (pminlon_ind.(fieldnames{k}) <= sig) = 0;
    pminlon_ind.(fieldnames{k}) (pminlon_ind.(fieldnames{k}) > sig) = 1;
    pamp_ind.(fieldnames{k}) (pamp_ind.(fieldnames{k}) <= sig) = 0;
    pamp_ind.(fieldnames{k}) (pamp_ind.(fieldnames{k}) > sig) = 1;
end

%% plotting

toplot = permute(cat(3,rmaxlon.highcl,rmaxlon.lowcl,rminlon.highcl,rminlon.lowcl,...
    ramp.highcl,ramp.lowcl),[3,2,1]);
toplot_p = permute(cat(3,pmaxlon.highcl,pmaxlon.lowcl,pminlon.highcl,pminlon.lowcl,...
    pamp.highcl,pamp.lowcl),[3,2,1]);
toplot_p (toplot_p <= sig) = 0;
toplot_p (toplot_p > sig) = 1;

cbrew = cbrewer('div','RdBu',16);         

contourtitle = {'Max longitude 1995-2024','Max longitude 1955-1979',...
    'Min longitude 1995-2024','Min longitude 1955-1979',...
    'Amplitude 1995-2024','Amplitude 1955-1979'};       

subplotmaps(toplot,longitudes,latitudes,{'div','RdBu'},1,toplot_p,12,contourtitle,'Longitude','','Correlation','on',...
    [-.5,.5],18,[-90:10:20],[-90:10:20],[longitudes(1:10:end)],[longitudes(1:10:end)],'',1,[0 360],[-90 45],0,'none',1,'Miller Cylindrical');

%% plotting individuals
cbrew = cbrewer('div','RdBu',16);         

%% maxlon

contourtitle = {'1','2','3','4','5','6','7','8','9'};       

subplotmaps(rmaxlon_ind.highcl,longitudes,latitudes,{'div','RdBu'},1,pmaxlon_ind.highcl,12,contourtitle,'Longitude','','Correlation','on',...
    [-.5,.5],18,[-90:10:20],[-90:10:20],[longitudes(1:10:end)],[longitudes(1:10:end)],'',1,[0 360],[-90 45],0,'none',1,'Miller Cylindrical');

contourtitle = {'1','2','3','4','5','6','7','8','9','10'};       

subplotmaps(rmaxlon_ind.lowcl,longitudes,latitudes,{'div','RdBu'},1,pmaxlon_ind.lowcl,12,contourtitle,'Longitude','','Correlation','on',...
    [-.5,.5],18,[-90:10:20],[-90:10:20],[longitudes(1:10:end)],[longitudes(1:10:end)],'',1,[0 360],[-90 45],0,'none',1,'Miller Cylindrical');

%% minlon

contourtitle = {'1','2','3','4','5','6','7','8','9','10'};       

subplotmaps(rminlon_ind.highcl,longitudes,latitudes,{'div','RdBu'},1,pminlon_ind.highcl,12,contourtitle,'Longitude','','Correlation','on',...
    [-.5,.5],18,[-90:10:20],[-90:10:20],[longitudes(1:10:end)],[longitudes(1:10:end)],'',1,[0 360],[-90 45],0,'none',1,'Miller Cylindrical');

contourtitle = {'1','2','3','4','5','6','7','8','9','10'};       

subplotmaps(rminlon_ind.lowcl,longitudes,latitudes,{'div','RdBu'},1,pminlon_ind.lowcl,12,contourtitle,'Longitude','','Correlation','on',...
    [-.5,.5],18,[-90:10:20],[-90:10:20],[longitudes(1:10:end)],[longitudes(1:10:end)],'',1,[0 360],[-90 45],0,'none',1,'Miller Cylindrical');

%% amplitude

contourtitle = {'1','2','3','4','5','6','7','8','9','10'};       

subplotmaps(ramp_ind.highcl,longitudes,latitudes,{'div','RdBu'},1,pamp_ind.highcl,12,contourtitle,'Longitude','','Correlation','on',...
    [-.5,.5],18,[-90:10:20],[-90:10:20],[longitudes(1:10:end)],[longitudes(1:10:end)],'',1,[0 360],[-90 45],0,'none',1,'Miller Cylindrical');

contourtitle = {'1','2','3','4','5','6','7','8','9','10'};       

subplotmaps(ramp_ind.lowcl,longitudes,latitudes,{'div','RdBu'},1,pamp_ind.lowcl,12,contourtitle,'Longitude','','Correlation','on',...
    [-.5,.5],18,[-90:10:20],[-90:10:20],[longitudes(1:10:end)],[longitudes(1:10:end)],'',1,[0 360],[-90 45],0,'none',1,'Miller Cylindrical');

% filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',var,'_TS_',...
%     num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2)),'_',num2str(100-percentile),...
%     'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S.png'];
% 
% export_fig(filename,'-png');
