function [SSWSPVyears] = predruns_SSWSPV(inputs)

% create composite of SSW and SPV events

%% Read in data

directory = '/Volumes/ExternalOne/work/data/predruns/U/highCl/daily/';

files = dir([directory,'*.nc']);

for i = 1:length(files)
    [~,data,~] = Read_in_netcdf([directory,files(i).name]);
    
    %combine into 4d array
    
    dataall(i,:,:,:,:) = data.U;
    if i == 1
        pressure = data.lev;
        latitude = data.lat;
    end
end

[~,latind] = min(abs(inputs.lat - latitude));
[~,presind] = min(abs(inputs.pres - pressure));

preslat_extract = squeeze(nanmean(dataall(:,:,latind,presind,:),2));

%%
nancat = zeros(size(preslat_extract,1),1);
nancat (nancat == 0) = NaN;
count = 1;
for i = 1:30            
    preslat_years(:,i,:) = [preslat_extract(:,count:count+364),nancat];
    count = count+364+1;
end

%% days in years

days = [day(datetime('1-Jan-2018'),'dayofyear'),day(datetime('1-Feb-2018'),'dayofyear'),...
    day(datetime('1-Mar-2018'),'dayofyear'),day(datetime('1-Apr-2018'),'dayofyear'),...
    day(datetime('1-May-2018'),'dayofyear'),day(datetime('1-Jun-2018'),'dayofyear'),...
    day(datetime('1-Jul-2018'),'dayofyear'),day(datetime('1-Aug-2018'),'dayofyear'),...
    day(datetime('1-Sep-2018'),'dayofyear'),day(datetime('1-Oct-2018'),'dayofyear'),...
    day(datetime('1-Nov-2018'),'dayofyear'),day(datetime('1-Dec-2018'),'dayofyear'),365];

daysleap = [day(datetime('1-Jan-2016'),'dayofyear'),day(datetime('1-Feb-2016'),'dayofyear'),...
    day(datetime('1-Mar-2016'),'dayofyear'),day(datetime('1-Apr-2016'),'dayofyear'),...
    day(datetime('1-May-2016'),'dayofyear'),day(datetime('1-Jun-2016'),'dayofyear'),...
    day(datetime('1-Jul-2016'),'dayofyear'),day(datetime('1-Aug-2016'),'dayofyear'),...
    day(datetime('1-Sep-2016'),'dayofyear'),day(datetime('1-Oct-2016'),'dayofyear'),...
    day(datetime('1-Nov-2016'),'dayofyear'),day(datetime('1-Dec-2016'),'dayofyear'),366];

%% Read in ERA data

ERAdirectory = '/Volumes/ExternalOne/work/data/ERA-Interim/U/daily/';
ERAfiles = dir([ERAdirectory,'*.nc']);
for i = 1:length(ERAfiles)
    ERAdata = ncread([ERAdirectory,ERAfiles(i).name],'u');
    ERAlatitude = ncread([ERAdirectory,ERAfiles(i).name],'latitude');
    [~,ERAlatind] = min(abs(ERAlatitude-60));
    ERA_preslat_extract = squeeze(nanmean(ERAdata(:,ERAlatind,:),1));
end

%%
%% accounting for Bodeker time restriction
count = 1;
for i = 1:40
    if i == 2 || i == 6 || i == 10 || i == 14 || i == 18 || i == 22 || i == 26 || i == 30 || i == 34 || i == 38
        ERApreslat_years(i,:) = ERA_preslat_extract(count:count+365);
        count = count+366;
    else
        ERApreslat_years(i,:) = [ERA_preslat_extract(count:count+364);NaN];
        count = count+365;        
    end
end
if inputs.bostimeres
    ERApreslat_years(1,:) = [];
    years = 1980:2016;
else
    years = 1979:2015;
end

% %% plotting averages and compare to ERA
% 
% cbrew = cbrewer('qual','Paired',10);
% 
% createfig('medium','on');
% 
% plot(preslat_years','color',cbrew(1,:),'LineWidth',2);
% hold on
% plot(nanmean(preslat_years),'color',cbrew(2,:),'LineWidth',3);
% 
% plot(ERApreslat_years','color',cbrew(3,:),'LineWidth',2);
% 
% plot(nanmean(ERApreslat_years),'color',cbrew(4,:),'LineWidth',4);

%% 
% find instances that SSWs occur (transition from westerlies to easterlies) 
% and instances that SPVs occur (westerlies greater than 48 ms-1)

for j = 1:size(preslat_years,1)
    for i = 1:size(preslat_years,2)    
        [~,SSWind(j,i).w] = find(squeeze(preslat_years(j,i,1:days(3)-1)) < 0);
        SSWSPVyears.model.SSW(j,i) = ~isempty(SSWind(j,i).w);
        [~,SPVind(j,i).w] = find(squeeze(preslat_years(j,i,1:days(3)-1)) > 42);
        SSWSPVyears.model.SPV(j,i) = ~isempty(SPVind(j,i).w);
    end
end

%%
SPVplotindtemp = SSWSPVyears.model.SPV;
SSWplotindtemp = SSWSPVyears.model.SSW;
SSWSPVyears.model.SSW (SSWSPVyears.model.SSW == SPVplotindtemp) = 0;
SSWSPVyears.model.SPV (SSWSPVyears.model.SPV == SSWplotindtemp) = 0;

SPV_percentage = sum(SSWSPVyears.model.SPV,2)/40;
SSW_percentage = sum(SSWSPVyears.model.SSW,2)/40;

%% same for ERA-interim
if inputs.bostimeres
    endind = 37;
else
    endind = 37;
end
%% January to March -- come back to this
count = 1;
for i = 1:endind
    [~,ERASSWind(i).w] = find(ERApreslat_years(i,1:days(5)-1) < 0);
    % find difference
    diffdates = diff(ERASSWind(i).w);    
    finalwarming = find(diffdates > 10,1,'last');
    
    if ~isempty(diffdates) && isempty(finalwarming) && ERASSWind(i).w(end) < 60
        finalwarming = 1
    end
    
    if ~isempty(ERASSWind(i).w)
        diffdates = [1,diffdates];
    end
    %find central date (first date that winds become easterly);
    
    for j = 1:finalwarming
        if j == 1 || diffdates(j) >= 20 
            Centraldates(count).SSW = datevec(datenum(years(i),1,ERASSWind(i).w(j)));
            count = count+1;
        end
    end
    
    SSWSPVyears.ERA.SSW(i) = ~isempty(ERASSWind(i).w);
    [~,ERASPVind(i).w] = find(ERApreslat_years(i,1:days(4)-1) > 48);
    SSWSPVyears.ERA.SPV(i) = ~isempty(ERASPVind(i).w);
end

%% December
count = 1;
for i = 1:endind
    [~,ERASSWindDec(i).w] = find(ERApreslat_years(i,days(12):end) < 0);    
    % find difference
    diffdatesDec = diff(ERASSWindDec(i).w);    
    finalwarming = find(diffdates > 10,1,'last');
    if ~isempty(ERASSWindDec(i).w)
        diffdatesDec = [1,diffdatesDec];
    end
    %find central date (first date that winds become easterly);
    
    for j = 1:length(diffdatesDec)
        if j == 1 || diffdatesDec(j) >= 20 
            CentraldatesDec(count).SSW = datevec(datenum(years(i),1,ERASSWindDec(i).w(j)));
            count = count+1;
        end
    end
    
    SSWSPVyears.ERA.SSW(i) = ~isempty(ERASSWind(i).w);
    [~,ERASPVind(i).w] = find(ERApreslat_years(i,1:days(4)-1) > 48);
    SSWSPVyears.ERA.SPV(i) = ~isempty(ERASPVind(i).w);
end

%%

ERASPVplotindtemp = SSWSPVyears.ERA.SPV;
ERASSWplotindtemp = SSWSPVyears.ERA.SSW;
SSWSPVyears.ERA.SSW (SSWSPVyears.ERA.SSW == ERASPVplotindtemp) = 0;
SSWSPVyears.ERA.SPV (SSWSPVyears.ERA.SPV == ERASSWplotindtemp) = 0;

ERASPV_percentage = sum(SSWSPVyears.ERA.SPV)/40;
ERASSW_percentage = sum(SSWSPVyears.ERA.SSW)/40;

%%

if inputs.plot
    cbrew = cbrewer('qual','Paired',10);
    createfig('largeportrait','off');
    subplot(2,1,1)
    plot(ERApreslat_years(SSWSPVyears.ERA.SSW,:)','color',cbrew(1,:),'LineWidth',2); 
    hold on; 
    plot(nanmean(ERApreslat_years(SSWSPVyears.ERA.SSW,:))','color',cbrew(2,:),'LineWidth',4);     
    plot(ERApreslat_years(SSWSPVyears.ERA.SPV,:)','color',cbrew(5,:),'LineWidth',2);
    plot(nanmean(ERApreslat_years(SSWSPVyears.ERA.SPV,:))','color',cbrew(6,:),'LineWidth',4);
    ylim([-30,80]);
    subplot(2,1,2)
    for i = 1:size(preslat_years,1)
        plot(squeeze(preslat_years(i,SSWSPVyears.model.SSW(i,:),:))','color',cbrew(1,:),'LineWidth',2); 
        combineSSW(i).w = squeeze(preslat_years(i,SSWSPVyears.model.SSW(i,:),:))';
        hold on; 
    end
    for i = 1:size(preslat_years,1)        
        plot(squeeze(preslat_years(i,SSWSPVyears.model.SPV(i,:),:))','color',cbrew(5,:),'LineWidth',2);        
        combineSPV(i).w = squeeze(preslat_years(i,SSWSPVyears.model.SPV(i,:),:))';
    end
    combineSSWall = [combineSSW(:).w];
    combineSPVall = [combineSPV(:).w];
    plot(nanmean(combineSSWall,2) ,'color',cbrew(2,:),'LineWidth',4);     
    plot(nanmean(combineSPVall,2) ,'color',cbrew(6,:),'LineWidth',4);
    ylim([-30,80]);
    filename = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/CompareSSWsSPVs.pdf';
    export_fig(filename,'-pdf');
    close all
end

end