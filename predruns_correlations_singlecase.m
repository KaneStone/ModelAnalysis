% prediction runs correlations to surface temperature

clear all 
close all
tozvar = 'toz';
var = 'TS';
%area = '6090S';
detrend = 1;
detrend_ozone = 0;
contourplots = 1;
individual_contourplots = 0;
lineplots = 1;
percentile = 20;
lats = [63,90];
%lats = [-90, -63];
tozmonth = 3;
indmonth = 3;
varmonth = [3,4]; % Can be any number of months in the year (e.g. [12,1,2] for austral summer)
shortnames = 0;
sig = .05;
strapme = 0;
ENSO = 0;
removeENSO = 1;
ClLevel = 'highCl';
% read in data 
%directqory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
%files = dir([directory,'*.nc']);

%% Read in highcl variable


timeperiodhigh = [1995,2016];%[1955,1975]
%timeperiodhigh = [1955,1975];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl,dataMonthArrangeMean.highcl]...
    = predruns_ReadInlayer(vardirectory,varfiles,var,timeperiodhigh,lats,detrend);

%% Read in TOZ highcl and take percentiles

tozdates = [1995,2016];
%tozdates = [1955,1974];%[1955,1975]
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);

%%

[pct_highcl,tozextract.highcl,Eachyear] = predruns_varPercentiles(toz_composite.highcl.montharrange,toz_dataMonthArrange.highcl,...
    tozmonth,percentile,length(tozfiles));

longitude = data(1).highcl.lon;
latitude = data(1).highcl.lat;

%%

predruns_diffandcorr_percentiles(dataMonthArrange,dataMonthArrangeMean,toz_dataMonthArrange,varmonth,tozmonth,var,latitude,...
    longitude,timeperiodhigh,lats,removeENSO,Eachyear,ClLevel,1);

%% take individual months
plotdiffmonths = 0;
if plotdiffmonths
    
    varmonthind = [3,4,5,6];
    
    for i = 1:length(varmonthind)
        [~,diff(i,:,:)] = predruns_diffandcorr_percentiles(dataMonthArrange,dataMonthArrangeMean,toz_dataMonthArrange,varmonthind(i),tozmonth,var,latitude,...
            longitude,timeperiodhigh,lats,removeENSO,Eachyear,ClLevel,1);
    end
    
    %%
    
    predruns_plotdiffmonths(diff,longitude,latitude,length(varmonthind));
    
end
    %% plot diffmonths

%% predict regression
TSdataout = load(['/Volumes/ExternalOne/work/data/predruns/output/data/TS_ninoremoved_',num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),num2str(abs(lats(1))),'-',...
num2str(abs(lats(2)))]);
predruns_regmodel(Eachyear,TSdataout,toz_dataMonthArrange,tozmonth,longitude,latitude,1,timeperiodhigh);

%% take correlations of composite and ensemble mean

predruns_enscorr(dataMonthArrange,toz_dataMonthArrange,varmonth,tozmonth,var,latitude,...
    longitude,timeperiodhigh,lats,removeENSO,ClLevel);


