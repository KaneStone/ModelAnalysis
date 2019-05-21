% replication of Smith et al 2018 NDISC and ERA-Interim plot from 1979-2015
clear all
close all
% Read in NDISC ice data

years = 1979:2015;

%% read in monthly observed sea ice
datatype = 'seaiceindex';
data = predruns_readinNDISC(datatype);
%%
datatype = 'Goddard'; %Goddard %ERA

%% converting to standard grid
if strcmp(datatype,'Goddard')
    if ~exist('/Volumes/ExternalOne/work/data/NSIDC_seaice/output/GoddardMergedStandard_bootstrapped.mat')    
        seaicedata2 = predruns_readinNDISC(datatype);
        seaicedata2_dates = repmat(1979:2017,[12,1]);
        standardlatsbins = 45:1.5:90;
        standardlonsbins = 0:1.5:360;
        standardlats = 45.75:1.5:90;
        standardlons = .75:1.5:360;
        for k = 1:length(seaicedata2)
            testdata = seaicedata2(k).seaiceconcentration(:);
            testlats = seaicedata2(k).latitude(:);            
            testlons = seaicedata2(k).longitude(:);
            testlonsind = (testlons < 0);
            testlons(testlonsind) = testlons(testlonsind) + 360;

            for i = 1:length(standardlatsbins)-1
                for j = 1:length(standardlonsbins)-1
                    Goddardstandardgrid(i,j,k) = nanmean(testdata(testlats >= standardlatsbins(i) & testlats < standardlatsbins(i+1) & testlons >= standardlonsbins(j) & testlons < standardlonsbins(j+1)));
                end
            end
        end
        save('/Volumes/ExternalOne/work/data/NSIDC_seaice/output/GoddardMergedStandard_bootstrapped.mat','Goddardstandardgrid','standardlats','standardlons','seaicedata2_dates');
    else
        load('/Volumes/ExternalOne/work/data/NSIDC_seaice/output/GoddardMergedStandard_bootstrapped.mat');
        seaicedata2.seaiceconcentration = permute(Goddardstandardgrid,[2,1,3]);
        seaicedata2.latitude = standardlats';
        seaicedata2.longitude = standardlons';
    end
else
    seaicedata2 = predruns_readinNDISC(datatype);
    seaicedata2_dates = repmat(1979:2018,[12,1]);
end

seaicedata2_dates = seaicedata2_dates(:);

%%

NSIDC_daily = ReadinNDISCbin;

%%
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

leapyears = [1980,1984,1988,1992,1996,2000,2004,2008,2012,2016];

%% take monthly averages

for i = 1:length(NSIDC_daily.standardgrid)
    for j = 1:12
        if find(NSIDC_daily.years(2) == leapyears)            
            monthave(i).m(:,:,j) = nanmean(NSIDC_daily.standardgrid(i).d(:,:,days(j):days(j+1)-1),3);
        else
            monthave(i).m(:,:,j) = nanmean(NSIDC_daily.standardgrid(i).d(:,:,daysleap(j):days(j+1)-1),3);
        end
    end            
end
monthave_rearrange = cat(4,monthave(:).m);

timeind(1) = find(years(1) == NSIDC_daily.years);
timeind(2) = find(years(end) == NSIDC_daily.years);
monthave_rearrange_extract = monthave_rearrange(:,:,:,timeind(1):timeind(2))./10;

%% detrend data
for i = 1:size(monthave_rearrange_extract,1)
    for j = 1:size(monthave_rearrange_extract,3)
        monthave_re_detrend(i,:,j,:) = detrend(squeeze(monthave_rearrange_extract(i,:,j,:))')+nanmean(squeeze(monthave_rearrange_extract(i,:,j,:))');
    end
end

monthave_re_detrend = permute(monthave_re_detrend,[1,4,3,2]);
monthave_re_detrend = cat(2,monthave_re_detrend,monthave_re_detrend(:,1,:,:));

%% Read in 10 hPa ERA-Interim winds
compareSSW = 0;
if compareSSW
    SSWin.lat = 60;
    SSWin.pres = 10;
    SSWin.plot = 1;
    SSWin.bostimeres = 0;
    [SSWSPVcomposite] = predruns_SSWSPV(SSWin);
else
    SSWSPVcomposite.SSW = logical([0,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,0,1,1,1,1,1,0,0,1,0,0]);
    SSWSPVcomposite.SPV = logical([0,1,1,1,1,0,0,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,0,0,0,0,1,0,1,1,1,1,1,0,0,1,0]);            
    temp1 = SSWSPVcomposite.SSW;
    temp2 = SSWSPVcomposite.SPV;    
    SSWSPVcomposite.SSW (SSWSPVcomposite.SSW == temp2) = 0;
    SSWSPVcomposite.SPV (SSWSPVcomposite.SPV == temp1) = 0;
end

% %% Read in ERA wind
% ERAdirectory = '/Volumes/ExternalOne/work/data/ERA-Interim/U/daily/';
% ERAfiles = dir([ERAdirectory,'*.nc']);
% for i = 1:length(ERAfiles)
%     ERAdata = ncread([ERAdirectory,ERAfiles(i).name],'u');
%     ERAlatitude = ncread([ERAdirectory,ERAfiles(i).name],'latitude');
%     [~,ERAlatind] = min(abs(ERAlatitude-60));
%     ERA_preslat_extract = squeeze(nanmean(ERAdata(:,ERAlatind,:),1));
% end
% 
% count = 1;
% for i = 1:37
%     if i == 2 || i == 6 || i == 10 || i == 14 || i == 18 || i == 22 || i == 26 || i == 30 || i == 34 || i == 38
%         ERApreslat_years(i,:) = ERA_preslat_extract(count:count+365);
%         count = count+366;
%     else
%         ERApreslat_years(i,:) = [ERA_preslat_extract(count:count+364);NaN];
%         count = count+365;        
%     end
% end

%%
% cbrew = cbrewer('qual','Set1',10);
% createfig('medium','on');
% plot(circshift(ERApreslat_years(SSWSPVcomposite.SSW,:),[0,177])','LineWidth',2,'color',cbrew(1,:))
% hold on
% plot(circshift(ERApreslat_years(SSWSPVcomposite.SPV,:),[0,177])','LineWidth',2,'color',cbrew(2,:))

%% extract data associated with SSWs and SPVs

monthave_rearrange_extract_SSW = nanmean(monthave_re_detrend(:,:,:,SSWSPVcomposite.SSW),4);
monthave_rearrange_extract_SPV = nanmean(monthave_re_detrend(:,:,:,SSWSPVcomposite.SPV),4);
monthave_rearrange_extract_clim = nanmean(monthave_re_detrend,4);
SSW_diff = monthave_rearrange_extract_clim - monthave_rearrange_extract_SSW;
SPV_diff = monthave_rearrange_extract_clim - monthave_rearrange_extract_SPV;
SSWSPVdifference = monthave_rearrange_extract_SPV - monthave_rearrange_extract_SSW;
%% plotting
longtouse = [NSIDC_daily.standardlons;360];

monstoave = [1,2,3;4,5,6;7,8,9];
for i = 1:3
    %SSW from climatology
    title = {['Change in ',monthnames(monstoave(i,:),1,0),' sea ice fraction due to winter SSWs',]};        
    SSW_difftouseave = nanmean(permute(SSW_diff(:,:,monstoave(i,:)),[3,2,1]));
    
    [fig,sh] = subplotmaps(SSW_difftouseave,longtouse,NSIDC_daily.standardlats,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',{'Percent'},'on',...
        [-24,24],22,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic',1);
    
    set(get(sh,'Title'),'Position',[.5 1.03 0])

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/recreateSmith',sprintf('%02d',monstoave(i,:)),'_',datatype,'_ChangeSEAICE_SSW','_detrend2'];

    export_fig(filename,'-png');
    
    %SPV from climatology
    title = {['Change in ',monthnames(monstoave(i,:),1,0),' sea ice fraction due to winter SPVs',]};        
    SPV_difftouseave = nanmean(permute(SPV_diff(:,:,monstoave(i,:)),[3,2,1]));    
    [fig,sh] = subplotmaps(SPV_difftouseave,longtouse,NSIDC_daily.standardlats,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',{'Percent'},'on',...
        [-24,24],22,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic',1);
    
    set(get(sh,'Title'),'Position',[.5 1.03 0])

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/recreateSmith',sprintf('%02d',monstoave(i,:)),'_',datatype,'_ChangeSEAICE_SPV','_detrend2'];

    export_fig(filename,'-png');
    
    %SSW difference from SPV
    title = {['Change in ',monthnames(monstoave(i,:),1,0),' sea ice fraction due to wind extremes',]};        
    differencetouse = nanmean(permute(SSWSPVdifference(:,:,monstoave(i,:)),[3,2,1]));    
    [fig,sh] = subplotmaps(differencetouse,longtouse,NSIDC_daily.standardlats,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',{'Percent'},'on',...
        [-24,24],22,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic',1);
    
    set(get(sh,'Title'),'Position',[.5 1.03 0])

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/recreateSmith',sprintf('%02d',monstoave(i,:)),'_',datatype,'_ChangeSEAICE_SSWfromSPV','_detrend2'];

    export_fig(filename,'-png');
end
    
%% Now how does this look due to ozone extremes

%% Read in Bodeker ozone
%% READ in Bodeker 
tcolat = [63,90];
tozmonth = 3;
obstimeperiod = [1980,2016];
[~,BSdata,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/BodekerScientific/TCO/Bodeker_TCO_monavg.nc');

BSyears = 1980:2016;
BSyear_vector = repmat(BSyears,12,1);
BSdata.years = [BSyear_vector(:);ones(length(BSdata.time)-length(BSyear_vector(:)),1)*max(BSyear_vector(:))+1];

ozonedateindex(1) = find(BSdata.years == obstimeperiod(1),1,'first');
ozonedateindex(2) = find(BSdata.years == obstimeperiod(2),1,'last');

latindex_bs = find(BSdata.lat >= tcolat(1) & BSdata.lat <= tcolat(2));

toz_zm = detrend(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+...
    tozmonth-1:12:ozonedateindex(2)),1)),BSdata.lat(latindex_bs))) + ...
    nanmean(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+...
    tozmonth-1:12:ozonedateindex(2)),1)),BSdata.lat(latindex_bs)));
toz_zm_nodetrend = weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+...
    tozmonth-1:12:ozonedateindex(2)),1)),BSdata.lat(latindex_bs));   

toz_zm = toz_zm(1:36); %1980:2015;

lowperc = prctile(toz_zm,20);
highperc = prctile(toz_zm,80);

lowind = find(toz_zm <= lowperc);
highind = find(toz_zm >= highperc);

%%

monthave_re_detrend_matchBS = monthave_re_detrend(:,:,:,2:end);

monthave_rearrange_extract_highozone = nanmean(monthave_re_detrend_matchBS(:,:,:,highind),4);
monthave_rearrange_extract_lowozone = nanmean(monthave_re_detrend_matchBS(:,:,:,lowind),4);
monthave_rearrange_extract_clim2 = nanmean(monthave_re_detrend_matchBS,4);
highozone_diff = monthave_rearrange_extract_clim2 - monthave_rearrange_extract_highozone;
lowozone_diff = monthave_rearrange_extract_clim2 - monthave_rearrange_extract_lowozone;
highlowdifference = monthave_rearrange_extract_lowozone - monthave_rearrange_extract_highozone;

%% plotting
longtouse = [NSIDC_daily.standardlons;360];

monstoave = [1,2,3;4,5,6;7,8,9];
for i = 1:3
    %SSW from climatology
    title = {['Change in ',monthnames(monstoave(i,:),1,0),' sea ice fraction due to high March ozone',]};        
    highozone_difftouseave = nanmean(permute(highozone_diff(:,:,monstoave(i,:)),[3,2,1]));
    
    [fig,sh] = subplotmaps(highozone_difftouseave,longtouse,NSIDC_daily.standardlats,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',{'Percent'},'on',...
        [-24,24],22,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic');
    
    set(get(sh,'Title'),'Position',[.5 1.03 0])

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/recreateSmith',sprintf('%02d',monstoave(i,:)),'_',datatype,'_ChangeSEAICE_highMarchozone','_detrend2'];

    export_fig(filename,'-png');
    
    %SPV from climatology
    title = {['Change in ',monthnames(monstoave(i,:),1,0),' sea ice fraction due to low March ozone',]};        
    lowozone_difftouseave = nanmean(permute(lowozone_diff(:,:,monstoave(i,:)),[3,2,1]));    
    [fig,sh] = subplotmaps(lowozone_difftouseave,longtouse,NSIDC_daily.standardlats,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',{'Percent'},'on',...
        [-24,24],22,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic');
    
    set(get(sh,'Title'),'Position',[.5 1.03 0])

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/recreateSmith',sprintf('%02d',monstoave(i,:)),'_',datatype,'_ChangeSEAICE_lowMarchozone','_detrend2'];

    export_fig(filename,'-png');
    
    %SSW difference from SPV
    title = {['Change in ',monthnames(monstoave(i,:),1,0),' sea ice fraction due to ozone extremes',]};        
    ozonedifferencetouse = nanmean(permute(highlowdifference(:,:,monstoave(i,:)),[3,2,1]));    
    [fig,sh] = subplotmaps(ozonedifferencetouse,longtouse,NSIDC_daily.standardlats,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',{'Percent'},'on',...
        [-24,24],22,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic');
    
    set(get(sh,'Title'),'Position',[.5 1.03 0])

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/recreateSmith',sprintf('%02d',monstoave(i,:)),'_',datatype,'_ChangeSEAICE_Marchozoneextremes','_detrend2'];

    export_fig(filename,'-png');
end
