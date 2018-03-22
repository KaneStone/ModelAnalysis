% prediction runs correlations to surface temperature

clear all 
close all
tozvar = 'toz';
var = 'TS';
%area = '6090S';
detrend = 1;
detrend_ozone = 1;
contourplots = 1;
individual_contourplots = 0;
lineplots = 1;
percentile = 50;
lats = [-90,-75];
tozmonth = 11;
indmonth = 11;
varmonth = [12,1,2]; % Can be any number of months in the year (e.g. [12,1,2] for austral summer)
shortnames = 0;
sig = .05;
strapme = 0;
ENSO = 0;
removeENSO = 0;
% read in data 
%directqory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/'];
%files = dir([directory,'*.nc']);
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

%% extracting percentiles

tozextract.lowcl.combine = [toz_composite.lowcl.montharrange(tozmonth,pct_lowcl.highind.a),...
    toz_composite.lowcl.montharrange(tozmonth,pct_lowcl.lowind.a)]';
tozextract.highcl.combine = [toz_composite.highcl.montharrange(tozmonth,pct_highcl.highind.a),...
    toz_composite.highcl.montharrange(tozmonth,pct_highcl.lowind.a)]';
tozextract.lowGHG.combine = [toz_composite.lowGHG.montharrange(tozmonth,pct_lowGHG.highind.a),...
    toz_composite.lowGHG.montharrange(tozmonth,pct_lowGHG.lowind.a)]';

longitude = data(1).lowcl.lon;
latitude = data(1).lowcl.lat;

%% NINO34

[indvarextract,indvarextractmean,indvardifference] = predruns_extractpct(dataMonthArrange,pct_highcl,...
    pct_lowcl,pct_lowGHG,indmonth,length(tozfiles),length(tozfilespast),length(tozfileslowGHG));

%% Calculate NINO

[NINO34,IOD,NH] = predruns_calculateENSO(indvarextract,composite,latitude,longitude,indmonth);

%% take correlations of composite and ensemble mean

predruns_enscorr(dataMonthArrange,toz_dataMonthArrange,varmonth,tozmonth,var,latitude,...
    longitude,timeperiodlow,timeperiodhigh,lats,1);

%% taking correlations

[varextract,varextractmean,vardifference] = predruns_extractpct(dataMonthArrange,pct_highcl,...
    pct_lowcl,pct_lowGHG,varmonth,length(tozfiles),length(tozfilespast),length(tozfileslowGHG));

for i = 1:length(longitude)
    for j = 1:length(latitude)     
        
        if removeENSO
            [b.highcl(i,j,:),bint.highcl(i,j,:,:),r.highcl(i,j,:),rint.highcl(i,j,:,:),stats(i,j).highcl] = ...
                regress([varextract.highcl.highind(:,j,i);varextract.highcl.lowind(:,j,i)],...
                [ones(length(NINO34.highcl),1),NINO34.highcl]);
        
            [b.lowcl(i,j,:),bint.lowcl(i,j,:,:),r.lowcl(i,j,:),rint.lowcl(i,j,:,:),stats(i,j).lowcl] = ...
                regress([varextract.lowcl.highind(:,j,i);varextract.lowcl.lowind(:,j,i)],...
                [ones(length(NINO34.lowcl),1),NINO34.lowcl]);
            
            [b.lowGHG(i,j,:),bint.lowGHG(i,j,:,:),r.lowGHG(i,j,:),rint.lowGHG(i,j,:,:),stats(i,j).lowGHG] = ...
                regress([varextract.lowGHG.highind(:,j,i);varextract.lowGHG.lowind(:,j,i)],...
                [ones(length(NINO34.lowGHG),1),NINO34.lowGHG]);
            
            varforcorr.highcl = r.highcl(i,j,:); 
            varforcorr.lowGHG = r.lowGHG(i,j,:); 
            varforcorr.lowcl = r.lowcl(i,j,:);
        else
            varforcorr.highcl = [varextract.highcl.highind(:,j,i);varextract.highcl.lowind(:,j,i)];
            varforcorr.lowGHG = [varextract.lowGHG.highind(:,j,i);varextract.lowGHG.lowind(:,j,i)];
            varforcorr.lowcl = [varextract.lowcl.highind(:,j,i);varextract.lowcl.lowind(:,j,i)];
        end
        
        [var_toz_correlations.lowcl(i,j),var_toz_correlations.lowcl_p(i,j)] = ...
            corr(squeeze(varforcorr.lowcl),tozextract.lowcl.combine);
        [var_toz_correlations.highcl(i,j),var_toz_correlations.highcl_p(i,j)] = ...
            corr(squeeze(varforcorr.highcl),tozextract.highcl.combine);
        [var_toz_correlations.lowGHG(i,j),var_toz_correlations.lowGHG_p(i,j)] = ...
            corr(squeeze(varforcorr.lowGHG),tozextract.lowGHG.combine);
        
        if ENSO
            [var_NINO34_correlations.lowcl(i,j),var_NINO34_correlations.lowcl_p(i,j)] = ...
                corr([varextract.lowcl.highind(:,j,i);varextract.lowcl.lowind(:,j,i)],NINO34.lowcl);
            [var_NINO34_correlations.highcl(i,j),var_NINO34_correlations.highcl_p(i,j)] = ...
                corr([varextract.highcl.highind(:,j,i);varextract.highcl.lowind(:,j,i)],NINO34.highcl);
            [var_NINO34_correlations.lowGHG(i,j),var_NINO34_correlations.lowGHG_p(i,j)] = ...
                corr([varextract.lowGHG.highind(:,j,i);varextract.lowGHG.lowind(:,j,i)],NINO34.lowGHG);
        end
        
        if strapme
            [bootstat.lowcl(i,j,:)] = bootstrp(500,@corr,...
                [varextract.lowcl.highind(:,j,i);varextract.lowcl.lowind(:,j,i)],...
                tozextract.lowcl.combine);
            [bootstat.highcl(i,j,:)] = bootstrp(500,@corr,...
                [varextract.highcl.highind(:,j,i);varextract.highcl.lowind(:,j,i)],...
                tozextract.highcl.combine);
        end
    end
    
end 

% find location of highest correlation
[maxnum,maxind] = max(var_toz_correlations.lowcl(:));
[maxrow,maxcol] = ind2sub(size(var_toz_correlations.lowcl),maxind);

%% Linear regression prediction

% construct trained regression model
[Mdl,fitInfo] = fitrlinear(tozextract.lowcl.combine(1:end-10),...
    [varextract.lowcl.highind(:,maxcol,maxrow);varextract.lowcl.lowind(1:end-10,maxcol,maxrow)],'Learner','LeastSquares');
YHat = predict(Mdl,tozextract.lowcl.combine(end-10+1:end));
Yact = varextract.lowcl.lowind(end-10+1:end,maxcol,maxrow);
plot(YHat);
figure;
plot(Yact);
figure;
plot(tozextract.lowcl.combine(end-10+1:end));
%%

if contourplots
    addlon = 0;
    if addlon
        lowcltoplot = [var_toz_correlations.lowcl(end,:);var_toz_correlations.lowcl];
        lowcl_h_toplot = [var_toz_correlations.lowcl_p(end,:);var_toz_correlations.lowcl_p];
        highcltoplot = [var_toz_correlations.highcl(end,:);var_toz_correlations.highcl];
        highcl_h_toplot = [var_toz_correlations.highcl_p(end,:);var_toz_correlations.highcl_p];
        lowGHGtoplot = [var_toz_correlations.lowGHG(end,:);var_toz_correlations.lowGHG];
        lowGHG_h_toplot = [var_toz_correlations.lowGHG_p(end,:);var_toz_correlations.lowGHG_p];
        longitudetoplot = [longitude(end);longitude];
    else
        lowcltoplot = var_toz_correlations.lowcl;
        lowcl_h_toplot = var_toz_correlations.lowcl_p;
        highcltoplot = var_toz_correlations.highcl;
        highcl_h_toplot = var_toz_correlations.highcl_p;
        lowGHGtoplot = var_toz_correlations.lowGHG;
        lowGHG_h_toplot = var_toz_correlations.lowGHG_p;
        longitudetoplot = longitude;
    end
    lowcltoplot = reshape(lowcltoplot,[1,size(lowcltoplot)]);
    lowcl_h_toplot (lowcl_h_toplot <= sig) = 0;
    lowcl_h_toplot (lowcl_h_toplot > sig) = 1;
    lowcl_h_toplot = reshape(lowcl_h_toplot,[1,size(lowcl_h_toplot)]);
        
    highcltoplot = reshape(highcltoplot,[1,size(highcltoplot)]);    
    highcl_h_toplot (highcl_h_toplot <= sig) = 0;
    highcl_h_toplot (highcl_h_toplot > sig) = 1;
    highcl_h_toplot = reshape(highcl_h_toplot,[1,size(highcl_h_toplot)]);
    
    lowGHGtoplot = reshape(lowGHGtoplot,[1,size(lowGHGtoplot)]);    
    lowGHG_h_toplot (lowGHG_h_toplot <= sig) = 0;
    lowGHG_h_toplot (lowGHG_h_toplot > sig) = 1;
    lowGHG_h_toplot = reshape(lowGHG_h_toplot,[1,size(lowGHG_h_toplot)]);
    
    if length(varmonth) > 1
        shortnames = 1;            
    end
    
    if removeENSO
        rmENSOappend = 'removeENSO';
    else
        rmENSOappend = [];
    end
        
    
%% plotting lowcl
    cbrew = cbrewer('div','RdBu',16);         
    
    contourtitle = {[monthnames(varmonth,1,shortnames),'{ }',var,'{ }',...
        monthnames(tozmonth,0,0),'{ }',tozvar,' correlations ',...
        num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S, ',num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2))]};       
    
    subplotmaps(lowcltoplot,longitudetoplot,latitude,{'div','RdBu'},1,lowcl_h_toplot,12,contourtitle,'Longitude','','Correlation','on',...
        [-.5,.5],18,[-90:10:20],[-90:10:20],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 0],0,'none',1,'Miller Cylindrical');
    
    if detrend
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',monthnames(tozmonth,0,0),tozvar,'_detrend',var,'_',...
            num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),...
            'S_Tperiod-',monthnames(varmonth,1,shortnames),'_',rmENSOappend,'.png'];
    else
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',var,'_TS_',...
            num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S.png'];
    end
    
    export_fig(filename,'-png');


    %% highcl
    cbrew = cbrewer('div','RdBu',16);    
    
    contourtitle = {[monthnames(varmonth,1,shortnames),'{ }',var,'{ }',...
        monthnames(tozmonth,0,0),'{ }',tozvar,' correlations ',...
        num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S, ',num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2))]};       
    
    subplotmaps(highcltoplot,longitudetoplot,latitude,{'div','RdBu'},1,highcl_h_toplot,12,contourtitle,'Longitude','latitude','Correlation','on',...
        [-.5,.5],22,[-90:10:20],[-90:10:20],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 0],0,'none',1,'Miller Cylindrical');
    
    if detrend
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',monthnames(tozmonth,0,0),tozvar,'_detrend',var,'_',...
            num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),...
            'S_Tperiod-',monthnames(varmonth,1,shortnames),'_',rmENSOappend,'.png'];
    else
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',var,'_TS_',...
            num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),'_',num2str(100-percentile),'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S.png'];
    end
    
    export_fig(filename,'-png');
    
    %% lowGHG
    cbrew = cbrewer('div','RdBu',16);    
    
    contourtitle = {[monthnames(varmonth,1,shortnames),'{ }',var,'{ }',...
        monthnames(tozmonth,0,0),'{ }',tozvar,' correlations ',...
        num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S, ',num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2))]};       
    
    subplotmaps(lowGHGtoplot,longitudetoplot,latitude,{'div','RdBu'},1,lowGHG_h_toplot,12,contourtitle,'Longitude','latitude','Correlation','on',...
        [-.5,.5],22,[-90:10:20],[-90:10:20],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 0],0,'none',1,'Miller Cylindrical');
    
    if detrend
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',monthnames(tozmonth,0,0),tozvar,'_detrend',var,'_',...
            num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),...
            'S_Tperiod-',monthnames(varmonth,1,shortnames),'_',rmENSOappend,'lowGHG.png'];
    else
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',var,'_TS_',...
            num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),'_',num2str(100-percentile),'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S.png'];
    end
    
    export_fig(filename,'-png');
    
end

%% plotting ENSO
if ENSO
    lowcltoplot = [var_NINO34_correlations.lowcl(end,:);var_NINO34_correlations.lowcl];
    lowcltoplot = reshape(lowcltoplot,[1,size(lowcltoplot)]);
    lowcl_h_toplot = [var_NINO34_correlations.lowcl_p(end,:);var_NINO34_correlations.lowcl_p];
    lowcl_h_toplot (lowcl_h_toplot <= sig) = 0;
    lowcl_h_toplot (lowcl_h_toplot > sig) = 1;
    lowcl_h_toplot = reshape(lowcl_h_toplot,[1,size(lowcl_h_toplot)]);
    
    highcltoplot = [var_NINO34_correlations.highcl(end,:);var_NINO34_correlations.highcl];
    highcltoplot = reshape(highcltoplot,[1,size(highcltoplot)]);
    highcl_h_toplot = [var_NINO34_correlations.highcl_p(end,:);var_NINO34_correlations.highcl_p];
    highcl_h_toplot (highcl_h_toplot <= sig) = 0;
    highcl_h_toplot (highcl_h_toplot > sig) = 1;
    highcl_h_toplot = reshape(highcl_h_toplot,[1,size(highcl_h_toplot)]);

    longitudetoplot = [longitude(end);longitude];

    if length(varmonth) > 1
        shortnames = 1;            
    end
    
%% plotting lowcl
    cbrew = cbrewer('div','RdBu',16);         
    
    contourtitle = {[monthnames(varmonth,1,shortnames),'{ }',var,'{ }',...
        monthnames(tozmonth,0,0),'{ }','NINO34',' correlations ',...
        num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S, ',num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2))]};       
    
    subplotmaps(lowcltoplot,longitudetoplot,latitude,{'div','RdBu'},1,lowcl_h_toplot,12,contourtitle,'Longitude','','Correlation','on',...
        [-1,1],18,[-90:10:20],[-90:10:20],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 20],0,'none',1);
    
    if detrend
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',monthnames(tozmonth,0,0),'NINO34','_detrend',var,'_',...
            num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S_Tperiod-',monthnames(varmonth,1,shortnames),'.png'];
    else
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',var,'_TS_',...
            num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S.png'];
    end
    
    export_fig(filename,'-png');


    %% highcl
    cbrew = cbrewer('div','RdBu',16);    
    
    contourtitle = {[monthnames(varmonth,1,shortnames),'{ }',var,'{ }',...
        monthnames(tozmonth,0,0),'{ }','NINO34',' correlations ',...
        num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S, ',num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2))]};       
    
    subplotmaps(highcltoplot,longitudetoplot,latitude,{'div','RdBu'},1,highcl_h_toplot,12,contourtitle,'Longitude','latitude','Correlation','on',...
        [-1,1],18,[-90:10:20],[-90:10:20],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 20],0,'none',1);
    
    if detrend
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',monthnames(tozmonth,0,0),'NINO34','_detrend',var,'_',...
            num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S_Tperiod-',monthnames(varmonth,1,shortnames),'.png'];
    else
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',var,'_TS_',...
            num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),'_',num2str(100-percentile),'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S.png'];
    end
    
    export_fig(filename,'-png');
end
