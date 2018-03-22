%  ERA Interim differences in low and high percentiles in ozone

% using ERA-Interim TS and Bodeker scientific TCO
clear all
tozmonth = 3; %3
ninomonth = 3; %3
varmonth = [3,4]; %3,4
timeperiod = [1980,2016];
timeperiodbs = [1980,2016];
shortnames = 0;
removeENSO = 1;
percentiles = 0;

useHalley = 0;
useBodeker = 1;
%% Read in ERA-interim surface temperature

ERAdata = ReadinERA('/Volumes/MyBook/work/data/ERA-Interim/T/T_ERA-Interim.nc');
ERAyears = 1979:2016;
ERAyear_vector = repmat(ERAyears,12,1);
ERAdata.years = [ERAyear_vector(:);ones(length(ERAdata.time)-length(ERAyear_vector(:)),1)*max(ERAyear_vector(:))+1];

%% extracting temperatures
ERAlats = [63,90];
ERAlats2 = [30,60];
ERAlatindex = ERAdata.latitude >= ERAlats(1) & ERAdata.latitude <= ERAlats(2);
ERAlatindex2 = ERAdata.latitude >= ERAlats2(1) & ERAdata.latitude <= ERAlats2(2);
if ERAlats(1) < 0
    ERAareamean1 = squeeze(nanmean(ERAdata.t(:,ERAlatindex,:,:),1));
    ERAareamean2 = squeeze(nanmean(ERAdata.t(:,ERAlatindex2,:,:),1));
else
    ERAareamean1 = squeeze(nanmean(ERAdata.t(:,ERAlatindex,:,13:end),1));
    ERAareamean2 = squeeze(nanmean(ERAdata.t(:,ERAlatindex2,:,13:end),1));
end
%%
for i = 1:size(ERAareamean1,2)
    [ERAwa(i,:)] = weightedaverage(squeeze(ERAareamean1(:,i,:)),ERAdata.latitude(ERAlatindex));        
    [ERAwa2(i,:)] = weightedaverage(squeeze(ERAareamean2(:,i,:)),ERAdata.latitude(ERAlatindex2));        
end

for i = 1:12
    ERAmonths(:,:,i) = ERAwa(:,i:12:end-2);
    ERAmonths2(:,:,i) = ERAwa2(:,i:12:end-2);
end

if ERAlats(1) < 0
    ERAmonths = circshift(ERAmonths,[0,0,7]);
    ERAmonths2 = circshift(ERAmonths2,[0,0,7]);
else    
end
%%
if useBodeker
    %% READ in Bodeker 
    tcolat = [63,90];
  
    tcolon = 26;
    tcolat_QBO = [10,30];
    [~,BSdata,~] = Read_in_netcdf('/Volumes/MyBook/work/data/BodekerScientific/TCO/Bodeker_TCO_monavg.nc');

    BSyears = 1980:2016;
    BSyear_vector = repmat(BSyears,12,1);
    BSdata.years = [BSyear_vector(:);ones(length(BSdata.time)-length(BSyear_vector(:)),1)*max(BSyear_vector(:))+1];

    ozonedateindex(1) = find(BSdata.years == timeperiodbs(1),1,'first');
    ozonedateindex(2) = find(BSdata.years == timeperiodbs(2),1,'last');


    %latindex_bs = find(BSdata.lat >= tcolat(1) & BSdata.lat <= tcolat(2));
    latindex_bs = find(BSdata.lat >= tcolat(1) & BSdata.lat <= tcolat(2));
    %toz_zm = squeeze(nanmean(BSdata.tco(:,:,BSdateindex(1)+tozmonth:12:BSdateindex(2)),1));
    toz_zm = detrend(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+tozmonth-1:12:ozonedateindex(2)),1)),BSdata.lat(latindex_bs))) + ...
        nanmean(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+tozmonth-1:12:ozonedateindex(2)),1)),BSdata.lat(latindex_bs)));
    toz_zm_nodetrend = weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+tozmonth-1:12:ozonedateindex(2)),1)),BSdata.lat(latindex_bs));        
    %toz_am = squeeze(nanmean(toz_zm(latindex_bs,:),1));
    %toz_am_QBO = squeeze(nanmean(toz_zm(latindex_bsQBO,:),1));

end

%%
if useHalley
    %% Read in Halley Station

    %% Initialize variables.
      tcolat = [-90,-63];
    filename = '/Volumes/MyBook/work/data/TCOstations/Halley_Total_Ozone_1956-2012.txt';
    startRow = 3;

    formatSpec = '%4f%7f%7f%7f%7f%7f%7f%7f%7f%7f%f%[^\n\r]';

    fileID = fopen(filename,'r');

    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

    fclose(fileID);

    Halley.Year = dataArray{:, 1};
    Halley.Nov = dataArray{:, 5};
    Halley.Nov = detrend(Halley.Nov);
        
    ozonedateindex(1) = find(Halley.Year == timeperiod(1));
    ozonedateindex(2) = find(Halley.Year == timeperiod(2));
    
    ozonetimeextract = Halley.Nov(ozonedateindex(1):ozonedateindex(2)-1)';
    toz_zm = detrend(Halley.Nov(ozonedateindex(1):ozonedateindex(2)-1)');
    toz_zm_nodetrend = Halley.Nov(ozonedateindex(1):ozonedateindex(2)-1)';
            
    clearvars filename startRow formatSpec fileID dataArray ans;
end
percentile = 20;
pct.lowerpercentile = prctile(toz_zm,percentile);
pct.upperpercentile = prctile(toz_zm,100-percentile);
pct.lowerind = find(toz_zm <= pct.lowerpercentile);
pct.upperind = find(toz_zm >= pct.upperpercentile);

%%

tempdiff(1,:,:) = squeeze(nanmean(ERAmonths(:,pct.lowerind,:),2)) - squeeze(nanmean(ERAmonths(:,pct.upperind,:),2));
%tempdiff(2,:,:) = squeeze(nanmean(ERAmonths2(:,pct.lowerind,:),2)) - squeeze(nanmean(ERAmonths2(:,pct.upperind,:),2));
tempdiff = permute(tempdiff,[1,3,2]);
%% plot
if tcolat(1) < 0
    xticklab = {'Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May'};
else
    tempdiff = circshift(tempdiff,3,2);
    xticklab = {'Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'};    
end
subplotmaps(tempdiff,1:12,log(double(ERAdata.level)),{'div','RdBu'},1,[],20,{'ERA-Interim 63-90N temperature difference over 1979-2016'},'Month','Pressure (hPa)','Temperature (K)','on',...
    [-15 15],22,[1:12],xticklab,[log([1,10,100,1000])],[1,10,100,1000],'',1,[1,12],[log(1) log(1000)],1,'-',0,'none');

%rtoplot = circshift((abs(squeeze(r.ind.highcl(i,:,:)))'),size((abs(squeeze(r.ind.highcl(i,:,:)))'),2)/2,2);
%m_contour(x,lats,rtoplot,[.6 .6],'color',[77,175,74]/255,'Linewidth',3,'LineStyle','-');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/T/','North'];

export_fig(filename,'-png');

