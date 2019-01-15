% recreate Justin's figure

% using ERA-Interim TS and Bodeker scientific TCO
clear all
tozmonth = 11;
ninomonth = 11;
varmonth = [12,13,14];
varmonth2 = [12,1,2];
timeperiod = [1979,2016];
shortnames = 0;
removeENSO = 0;
%% Read in ERA-interim surface temperature

ERAdata = ReadinERA('/Volumes/ExternalOne/work/data/ERA-Interim/TS/TS_ERA-Interim.nc');
ERAyears = 1979:2016;
ERAyear_vector = repmat(ERAyears,12,1);
ERAdata.years = [ERAyear_vector(:);ones(length(ERAdata.time)-length(ERAyear_vector(:)),1)*max(ERAyear_vector(:))+1];

%% Read in Halley Station

%% Initialize variables.
filename = '/Volumes/ExternalOne/work/data/TCOstations/Halley_Total_Ozone_1956-2012.txt';
startRow = 3;

formatSpec = '%4f%7f%7f%7f%7f%7f%7f%7f%7f%7f%f%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fclose(fileID);

Halley.Year = dataArray{:, 1};
Halley.Nov = dataArray{:, 5};
Halley.Nov = detrend(Halley.Nov);
clearvars filename startRow formatSpec fileID dataArray ans;

%% calculate ENSO

ERAdateindex(1) = find(ERAdata.years == timeperiod(1),1,'first');
ERAdateindex(2) = find(ERAdata.years == timeperiod(2),1,'last');

latlimits = [-5 5];
lonlimits = [190 240];

latindex = find(ERAdata.latitude >= latlimits(1) & ERAdata.latitude <= latlimits(2));
lonindex = find(ERAdata.longitude >= lonlimits(1) & ERAdata.longitude <= lonlimits(2));

NINO_mn = squeeze(nanmean(ERAdata.t2m(lonindex,latindex,ERAdateindex(1)+ninomonth-1:12:ERAdateindex(2)),1));
NINO_mn = squeeze(nanmean(NINO_mn,1));
NINOallmean = nanmean(NINO_mn);
NINOstd = std(NINO_mn,1);
NINO34 = detrend((NINO_mn - NINOallmean)./NINOstd);

%% calculate TS month averages

for i = 1:length(varmonth)
    for j = 1:size(ERAdata.t2m,2)
        for k = 1:size(ERAdata.t2m,1)
            ts(k,j,:,i) = detrend(squeeze(ERAdata.t2m(k,j,ERAdateindex(1)+varmonth(i)-1:12:ERAdateindex(2)+2))) + ...
                repmat(squeeze(nanmean(ERAdata.t2m(k,j,ERAdateindex(1)+varmonth(i)-1:12:ERAdateindex(2)+2),3)),[1,length(timeperiod(1):timeperiod(2))])';
        end
    end
    ts2(:,:,:,i) = squeeze(ERAdata.t2m(:,:,ERAdateindex(1)+varmonth(i)-1:12:ERAdateindex(2)+2));
end
tsmean = nanmean(ts,4);

%% taking correlations
Halleydateindex(1) = find(Halley.Year == timeperiod(1));
Halleydateindex(2) = find(Halley.Year == timeperiod(2));
for i = 1:size(tsmean,1)
    for j = 1:size(tsmean,2)
        
        if removeENSO
            [b(i,j,:),bint(i,j,:,:),r(i,j,:),rint(i,j,:,:),stats(i,j,:)] = regress(squeeze(tsmean(i,j,1:length(Halley.Nov(Halleydateindex(1):Halleydateindex(2))))),[ones(length(NINO34),1),NINO34']);        

            varforcorr = [squeeze(tsmean(i,j,1:length(Halley.Nov(Halleydateindex(1):Halleydateindex(2)))))' - squeeze(b(i,j,2))*NINO34]';
            
        else
            varforcorr = tsmean(i,j,1:length(Halley.Nov(Halleydateindex(1):Halleydateindex(2))));
        end
        
        [polar_correlations(i,j),polar_p(i,j)] = ...
               corr(squeeze(varforcorr),Halley.Nov(Halleydateindex(1):Halleydateindex(2)));
    end
end

%% taking rolling correlations

[rollcorr,rolltrend] = rollingCorrelations(squeeze(tsmean),Halley.Nov(Halleydateindex(1):Halleydateindex(2)),15);


%% plotting

% lontoplot = [ERAdata.longitude(end);ERAdata.longitude];
lontoplot = [ERAdata.longitude];
% polar_correlations = [polar_correlations(end,:);polar_correlations];
 polartoplot = reshape(polar_correlations,[1,size(polar_correlations)]);

sig = .05;
% polar_p = [polar_p(end,:);polar_p];
polar_p_toplot = reshape(polar_p,[1,size(polar_p)]);
polar_p_toplot (polar_p_toplot <= sig) = 0;
polar_p_toplot (polar_p_toplot > sig) = 1;

plotting = 1;

if removeENSO
    rmENSOappend = 'removeENSO';
else
    rmENSOappend = [];
end

lat = ERAdata.latitude;
lon = ERAdata.longitude;

save(['/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefs_Justins_',num2str(timeperiod(1)),'-',num2str(timeperiod(2)),'_',rmENSOappend],...
    'polar_correlations','polar_p','rollcorr','rolltrend','lat','lon');

%% 

if plotting
    %% plotting polar
    cbrew = cbrewer('div','RdBu',16);                

    contourtitle = {['Halley Station ',num2str(timeperiod(1)),'-',num2str(timeperiod(2)),' November TCO - ERA-interim DJF surface temperature correlations']};

    subplotmaps(polartoplot,lontoplot,ERAdata.latitude,{'div','RdBu'},1,polar_p_toplot,16,contourtitle,'Longitude','','Correlation','on',...
        [-.5,.5],18,[-90:10:20],[-90:10:20],[lontoplot(2:10:end)],[lontoplot(2:10:end)],'',1,[0 360],[-90 0],0,'none',1,'Miller Cylindrical');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/','correlations/maps/','Justins_figure',num2str(timeperiod(1)),'-',num2str(timeperiod(2))];
    
%     export_fig(filename,'-png','-r300');
    export_fig(filename,'-pdf');

    %% plotting rolling correlation trends
    
    contourtitle = {'Halley Station November TCO - ERA-interim DJF rolling correlation trends (15 year window) 1979-2013'};

    subplotmaps(reshape(rolltrend.r(:,:,2),[1,size(rolltrend.r(:,:,2))])*10,lontoplot,ERAdata.latitude,{'div','RdBu'},1,[],12,contourtitle,'Longitude','','r/decade','on',...
        [-.75,.75],18,[-90:10:20],[-90:10:20],[lontoplot(2:10:end)],[lontoplot(2:10:end)],'',1,[0 360],[-90 0],0,'none',1,'Miller Cylindrical');

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/','Justins_figure19792013_rolling.png'];
    
    export_fig(filename,'-png','-r300');
    
    
end
