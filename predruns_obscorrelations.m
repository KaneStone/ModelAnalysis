% using ERA-Interim TS and Bodeker scientific TCO
clear all
tozmonth = 11;
ninomonth = 11;
varmonth = [12,13,14];
varmonth2 = [12,1,2];
timeperiod = [1979,2012];
timeperiodbs = [1980,2012];
shortnames = 0;
removeENSO = 0;
%% Read in ERA-interim surface temperature

ERAdata = ReadinERA('/Volumes/MyBook/work/data/ERA-Interim/TS/TS_ERA-Interim.nc');
ERAyears = 1979:2016;
ERAyear_vector = repmat(ERAyears,12,1);
ERAdata.years = [ERAyear_vector(:);ones(length(ERAdata.time)-length(ERAyear_vector(:)),1)*max(ERAyear_vector(:))+1];
%% READ in Bodeker 
%tcolat = [75,90];
tcolat = -75.5;
tcolon = 26;
tcolat_QBO = [10,30];
[~,BSdata,~] = Read_in_netcdf('/Volumes/MyBook/work/data/BodekerScientific/TCO/Bodeker_TCO_monavg.nc');

BSyears = 1980:2016;
BSyear_vector = repmat(BSyears,12,1);
BSdata.years = [BSyear_vector(:);ones(length(BSdata.time)-length(BSyear_vector(:)),1)*max(BSyear_vector(:))+1];

BSdateindex(1) = find(BSdata.years == timeperiodbs(1),1,'first');
BSdateindex(2) = find(BSdata.years == timeperiodbs(2),1,'last');


%latindex_bs = find(BSdata.lat >= tcolat(1) & BSdata.lat <= tcolat(2));
[~,latindex_bs] = min(abs(BSdata.lat - tcolat));
[~,lonindex_bs] = min(abs(BSdata.lon - tcolon));
latindex_bsQBO = find(BSdata.lat >= tcolat_QBO(1) & BSdata.lat <= tcolat_QBO(2));
%toz_zm = squeeze(nanmean(BSdata.tco(:,:,BSdateindex(1)+tozmonth:12:BSdateindex(2)),1));
toz_zm = squeeze(BSdata.tco(lonindex_bs,latindex_bs,BSdateindex(1)+tozmonth:12:BSdateindex(2)));
%toz_am = squeeze(nanmean(toz_zm(latindex_bs,:),1));
%toz_am_QBO = squeeze(nanmean(toz_zm(latindex_bsQBO,:),1));

%% Read in Halley Station

%% Initialize variables.
filename = '/Volumes/MyBook/work/data/TCOstations/Halley_Total_Ozone_1956-2012.txt';
startRow = 3;

formatSpec = '%4f%7f%7f%7f%7f%7f%7f%7f%7f%7f%f%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);

Halley.Year = dataArray{:, 1};
Halley.Aug = dataArray{:, 2};
Halley.Sep = dataArray{:, 3};
Halley.Oct = dataArray{:, 4};
Halley.Nov = dataArray{:, 5};
Halley.Dec = dataArray{:, 6};
Halley.Jan = dataArray{:, 7};
Halley.Feb = dataArray{:, 8};
Halley.Mar = dataArray{:, 9};
Halley.Apr = dataArray{:, 10};
Halley.Year1 = dataArray{:, 11};
Halley.Nov = detrend(Halley.Nov);
clearvars filename startRow formatSpec fileID dataArray ans;

%% calculate ENSO

ERAdateindex(1) = find(ERAdata.years == timeperiod(1),1,'first');
ERAdateindex(2) = find(ERAdata.years == timeperiod(2),1,'last');

latlimits = [-5 5];
lonlimits = [190 240];

latindex = find(ERAdata.latitude >= latlimits(1) & ERAdata.latitude <= latlimits(2));
lonindex = find(ERAdata.longitude >= lonlimits(1) & ERAdata.longitude <= lonlimits(2));

NINO_mn = squeeze(nanmean(ERAdata.t2m(lonindex,latindex,ERAdateindex(1)+ninomonth:12:ERAdateindex(2)),1));
NINO_mn = squeeze(nanmean(NINO_mn,1));
NINOallmean = nanmean(NINO_mn);
NINOstd = std(NINO_mn,1);
NINO34 = detrend((NINO_mn - NINOallmean)./NINOstd);

%% calculate TS month averages

for i = 1:length(varmonth)
    for j = 1:size(ERAdata.t2m,2)
        for k = 1:size(ERAdata.t2m,1)
            ts(k,j,:,i) = detrend(squeeze(ERAdata.t2m(k,j,ERAdateindex(1)+varmonth(i):12:ERAdateindex(2)+12))) + ...
                repmat(squeeze(nanmean(ERAdata.t2m(k,j,ERAdateindex(1)+varmonth(i):12:ERAdateindex(2)+12),3)),[1,34])';
        end
    end
    ts2(:,:,:,i) = squeeze(ERAdata.t2m(:,:,ERAdateindex(1)+varmonth(i):12:ERAdateindex(2)+12));
end
tsmean = nanmean(ts,4);

%% taking correlations

for i = 1:size(tsmean,1)
    for j = 1:size(tsmean,2)
        
        if removeENSO
            [b(i,j,:),bint(i,j,:,:),r(i,j,:),rint(i,j,:,:),stats(i,j,:)] = regress(squeeze(tsmean(i,j,1:length(Halley.Nov(24:57)))),[ones(length(NINO34),1),NINO34']);
        
            varforcorr = r(i,j,:); 

        else
            varforcorr = tsmean(i,j,1:length(Halley.Nov(24:57)));
        end
        
%         [polar_correlations(i,j),polar_p(i,j)] = ...
%             corr(squeeze(varforcorr),toz_zm);
        
        [polar_correlations(i,j),polar_p(i,j)] = ...
            corr(squeeze(varforcorr),Halley.Nov(24:57));
        
%         [midlat_correlations(i,j),midlat_p(i,j)] = ...
%             corr(squeeze(varforcorr),toz_am_QBO');
%         
%         [NINO_correlations(i,j),NINO_p(i,j)] = ...
%             corr(squeeze(varforcorr),NINO34(1:length(toz_am))');
        
                        
    end
end

%% plotting

lontoplot = [ERAdata.longitude(end);ERAdata.longitude];
polar_correlations = [polar_correlations(end,:);polar_correlations];
polartoplot = reshape(polar_correlations,[1,size(polar_correlations)]);
% midlat_correlations = [midlat_correlations(end,:);midlat_correlations];
% midlattoplot = reshape(midlat_correlations,[1,size(midlat_correlations)]);
% NINO_correlations = [NINO_correlations(end,:);NINO_correlations];
% NINOtoplot = reshape(NINO_correlations,[1,size(NINO_correlations)]);

sig = .05;
polar_p = [polar_p(end,:);polar_p];
polar_p_toplot = reshape(polar_p,[1,size(polar_p)]);
polar_p_toplot (polar_p_toplot <= sig) = 0;
polar_p_toplot (polar_p_toplot > sig) = 1;
% midlat_p = [midlat_p(end,:);midlat_p];
% midlat_p_toplot = reshape(midlat_p,[1,size(midlat_p)]);
% midlat_p_toplot (midlat_p_toplot <= sig) = 0;
% midlat_p_toplot (midlat_p_toplot > sig) = 1;
% NINO_p = [NINO_p(end,:);NINO_p];
% NINO_p_toplot = reshape(NINO_p,[1,size(NINO_p)]);
% NINO_p_toplot (NINO_p_toplot <= sig) = 0;
% NINO_p_toplot (NINO_p_toplot > sig) = 1;

plotting = 1;

if removeENSO
    rmENSOappend = 'removeENSO';
else
    rmENSOappend = [];
end

%% 

if plotting
    %% plotting polar
    cbrew = cbrewer('div','RdBu',16);         

%     contourtitle = {[monthnames(varmonth2,1,shortnames),'{ }','TS','{ }',...
%         monthnames(tozmonth,0,0),'{ }','toz',' correlations ',...
%         num2str(abs(tcolat(1))),'-',num2str(abs(tcolat(2))),'S, ',num2str(timeperiod(1)),'-',num2str(timeperiod(2))]};       

    contourtitle = {'test'};

    subplotmaps(polartoplot,lontoplot,ERAdata.latitude,{'div','RdBu'},1,polar_p_toplot,12,contourtitle,'Longitude','','Correlation','on',...
        [-.5,.5],18,[-90:10:20],[-90:10:20],[lontoplot(2:10:end)],[lontoplot(2:10:end)],'',1,[0 360],[-90 -15],0,'none',1);


%     filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/obs/polar',monthnames(tozmonth,0,0),'toz','_detrendTS_',...
%         num2str(timeperiod(1)),'-',num2str(timeperiod(2)),'_',num2str(abs(tcolat(1))),'-',...
%         num2str(abs(tcolat(2))),'S_Tperiod-',monthnames(varmonth2,1,shortnames),'_',rmENSOappend,'.png'];
% 
%     export_fig(filename,'-png');

    %% plotting midlat
    cbrew = cbrewer('div','RdBu',16);         

    contourtitle = {[monthnames(varmonth2,1,shortnames),'{ }','TS','{ }',...
        monthnames(tozmonth,0,0),'{ }','toz',' correlations ',...
        num2str(abs(tcolat_QBO(1))),'-',num2str(abs(tcolat_QBO(2))),'S, ',num2str(timeperiod(1)),'-',num2str(timeperiod(2))]};       

    subplotmaps(midlattoplot,lontoplot,ERAdata.latitude,{'div','RdBu'},1,midlat_p_toplot,12,contourtitle,'Longitude','','Correlation','on',...
        [-.75,.75],18,[-90:10:20],[-90:10:20],[lontoplot(2:10:end)],[lontoplot(2:10:end)],'',1,[0 360],[90 0],0,'none',1);


    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/obs/midlat',monthnames(tozmonth,0,0),'toz','_detrendTS_',...
        num2str(timeperiod(1)),'-',num2str(timeperiod(2)),'_',num2str(abs(tcolat(1))),'-',...
        num2str(abs(tcolat(2))),'S_Tperiod-',monthnames(varmonth2,1,shortnames),'_',rmENSOappend,'.png'];

    export_fig(filename,'-png');

    %% plotting NINO
    cbrew = cbrewer('div','RdBu',16);         

    contourtitle = {[monthnames(varmonth2,1,shortnames),'{ }','TS','{ }',...
        monthnames(ninomonth,0,0),'{ }','NINO34',' correlations, ',num2str(timeperiod(1)),'-',num2str(timeperiod(2))]};       

    subplotmaps(NINOtoplot,lontoplot,ERAdata.latitude,{'div','RdBu'},1,NINO_p_toplot,12,contourtitle,'Longitude','','Correlation','on',...
        [-.75,.75],18,[-90:10:20],[-90:10:20],[lontoplot(2:10:end)],[lontoplot(2:10:end)],'',1,[0 360],[90 0],0,'none',1);


    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/obs/NINO34',monthnames(tozmonth,0,0),'toz','_detrendTS_',...
        num2str(timeperiod(1)),'-',num2str(timeperiod(2)),'_',num2str(abs(tcolat(1))),'-',...
        num2str(abs(tcolat(2))),'S_Tperiod-',monthnames(varmonth2,1,shortnames),'_',rmENSOappend,'.png'];

    export_fig(filename,'-png');
end