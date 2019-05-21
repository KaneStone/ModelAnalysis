%% create time series of integrated sea ice
clear all

%%inputs
mon = 9;

%% read in observed sea ice
datatype = 'seaiceindex';
data = predruns_readinNDISC(datatype);
%%
datatype = 'ERA'; %Goddard %ERA

%% testing lat long standard grid
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
compareSSW = 1;
if compareSSW
    SSWin.lat = 60;
    SSWin.pres = 10;
    SSWin.plot = 0;
    inputs.bostimseries
    [SSWSPVcomposite] = predruns_SSWSPV(SSWin);
end

%%
% if strcmp(datatype,'source')
%     for i = 1:length(seaicedata2)
%         seaicedataweighted(i,:,:) = seaicedata2(i).seaiceconcentration.*cosd(seaicedata2(i).latitude);    
%         seaicedatanonweighted(i,:,:) = seaicedata2(i).seaiceconcentration;
%     end
%     seaicedataintegration = nansum(seaicedataweighted(:,:),2);
%     seaicedetrend = detrend(seaicedataintegration);
%     seaicedetrendextract = seaicedetrend(seaicedata2_dates >= 1980 & seaicedata2_dates <= 2016);
% elseif strcmp(datatype,'ERA')

    seaicedataweighted = permute(seaicedata2.seaiceconcentration,[2,1,3]).*cosd(seaicedata2.latitude);    
	seaicedatanonweighted = seaicedata2.seaiceconcentration;
    lats = [50,90];
    latsindex = find(seaicedata2.latitude > lats(1) & seaicedata2.latitude < lats(2));
    seaicedataintegration = permute(seaicedataweighted(latsindex,:,:),[3,1,2]);
    seaicedataintegration = nansum(seaicedataintegration(:,:),2);
    seaicedetrend = detrend(seaicedataintegration);
    seaicedetrendextract = seaicedetrend(seaicedata2_dates >= 1980 & seaicedata2_dates <= 2016);
    
    
  
    seaicedataweighted = permute(seaicedataweighted,[3,2,1]);
    seaicedatanonweighted = permute(seaicedatanonweighted,[3,1,2]);
    
% end

%%
seaicedatanonweighted_extract = seaicedatanonweighted(seaicedata2_dates >= 1980 & seaicedata2_dates <= 2016,:,:);
seaicedataweighted_extract = seaicedataweighted(seaicedata2_dates >= 1980 & seaicedata2_dates <= 2016,:,:);
for i = 1:12
    seaicedatanonweighted_reshape(i,:,:,:) = seaicedatanonweighted_extract(i:12:end,:,:);
    seaicedataweighted_reshape(i,:,:,:) = seaicedataweighted_extract(i:12:end,:,:);
    for j = 1:size(seaicedatanonweighted_reshape,4)
        %detrend
        seaicedatanonweighted_reshape_detrend(i,:,:,j) = detrend(squeeze(seaicedatanonweighted_reshape(i,:,:,j)));
        seaicedataweighted_reshape_detrend(i,:,:,j) = detrend(squeeze(seaicedataweighted_reshape(i,:,:,j)))+nanmean(squeeze(seaicedataweighted_reshape(i,:,:,j)));
    end
end

%% read in TCO

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

lowperc = prctile(toz_zm,20);
highperc = prctile(toz_zm,80);

lowind = find(toz_zm <= lowperc);
highind = find(toz_zm >= highperc);

%% calculate difference 
%if strcmp(datatype,'ERA')

seasons = 1;
wind = 0;

if wind
    windext = 'SSWSPV';
else
    windext = 'MarchOzone';
end

if seasons 
    mons = [4,5,6;7,8,9;10,11,12];
else
    mons = [1:12]';
end
for i = 1:length(mons) %mon = 1:12
    longtouse = [seaicedata2.longitude;seaicedata2.longitude(end)+(seaicedata2.longitude(2)-seaicedata2.longitude(1))];
    if wind
        title = {['Change in ',monthnames(mons(i,:),1,0),' sea ice fraction due to winter SSWs and SPVs',]};        
    else
        title = {['Change in ',monthnames(mons(i,:),1,0),' sea ice fraction due to ',monthnames(tozmonth,0,0), ' ozone extremes']};        
    end
    ctitle = {'Ice fraction'};
    clims = [-.3 .3];
    titext = 'ICEFRAC';
    toplotp = [];

    %difference (isnan(difference)) = -10;
    if seasons
        if wind
            difference = squeeze(nanmean(nanmean(seaicedatanonweighted_reshape_detrend(mons(i,:),SSWSPVcomposite.ERA.SSW,:,:),1),2) - nanmean(nanmean(seaicedatanonweighted_reshape_detrend(mons(i,:),SSWSPVcomposite.ERA.SPV,:,:),1),2));
        else
            difference = squeeze(nanmean(nanmean(seaicedatanonweighted_reshape_detrend(mons(i,:),highind,:,:),1),2) - nanmean(nanmean(seaicedatanonweighted_reshape_detrend(mons(i,:),lowind,:,:),1),2));
        end
        
    else
        if wind
            difference = squeeze(nanmean(seaicedatanonweighted_reshape_detrend(i,SSWSPVcomposite.ERA.SSW,:,:),2) - nanmean(seaicedatanonweighted_reshape_detrend(i,SSWSPVcomposite.ERA.SPV,:,:),2));
        else            
            difference = squeeze(nanmean(seaicedatanonweighted_reshape_detrend(i,highind,:,:),2) - nanmean(seaicedatanonweighted_reshape_detrend(i,lowind,:,:),2));
        end        
    end
    difference = permute(repmat(difference,[1,1,1]),[3,1,2]);
    difference (isnan(difference)) = 0;
    difference = cat(2,difference,difference(:,1,:));        
    
    [fig,sh] = subplotmaps(difference,longtouse,seaicedata2.latitude,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',ctitle,'on',...
        clims,22,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic');
    
    set(get(sh,'Title'),'Position',[.5 1.03 0])

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/',sprintf('%02d',mons(i,:)),'_',datatype,'_ChangeSEAICE_',windext,'_detrend'];

    export_fig(filename,'-png');
    
end

%% plot sea ice extent
seasons = 1;
if seasons 
    mons = [4,5,6;7,8,9;10,11,12];
else
    mons = [1:12]';
end

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewer('qual','Set2',8);
createfig('medium','on')

clearvars title

for i = 1:length(mons) %mon = 1:12
    longtouse = [seaicedata2.longitude;seaicedata2.longitude(end)+(seaicedata2.longitude(2)-seaicedata2.longitude(1))];
    %ctitle = {'Ice fraction'};
    %clims = [-.3 .3];    
    toplotp = [];

    %difference (isnan(difference)) = -10;
    if seasons
        lowmean = squeeze(nanmean(nanmean(seaicedatanonweighted_reshape(mons(i,:),lowind,:,:),1),2));
        highmean = squeeze(nanmean(nanmean(seaicedatanonweighted_reshape(mons(i,:),highind,:,:),1),2));
    else
        lowmean = squeeze(nanmean(nanmean(seaicedatanonweighted_reshape(i,lowind,:,:),1),2));
        highmean = squeeze(nanmean(nanmean(seaicedatanonweighted_reshape(i,highind,:,:),1),2));
    end
    
    lowmean = cat(1,lowmean,lowmean(1,:));
    highmean = cat(1,highmean,highmean(1,:));
    createfig('medium','on');
    m_proj('Stereographic','lon',-180,'lat',90,'rad',abs(90-45))
    
    %m_proj('miller','lon',-180,'lat',ylimit(1));
    [~,h1] = m_contour(longtouse,seaicedata2.latitude,lowmean',[0, .15]);
    hold on
    [h11] = plot(1,1,'color',cbrewqual(1,:),'LineWidth',3);
    set(h1,'LineWidth',3,'color',cbrewqual(1,:));
    
    [~,h2] = m_contour(longtouse,seaicedata2.latitude,highmean',[0, .15]);
    [h12] = plot(1,1,'color',cbrewqual(2,:),'LineWidth',3);
    set(h2,'LineWidth',3,'color',cbrewqual(2,:));

    %m_coast('color','k','LineWidth',1);
    m_coast('patch',[.5 .5 .5],'edgecolor','none','LineWidth',1);    
    m_grid('ytick',0:15:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',20,'LineWidth',2);    
    %m_grid('ytick',[longtouse(1:24:end-1)] ,'xtick',6,'XaxisLocation','bottom','fontsize',20);            
    
%     [fig,sh] = subplotmaps(difference,longtouse,seaicedata2.latitude,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',ctitle,'on',...
%         clims,22,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic');
    
    tit = {['Change in ',monthnames(mons(i,:),1,0),' sea ice extent due to ',monthnames(tozmonth,0,0), ' ozone extremes']};     
    th = title(tit);
    set(th,'units','Normalize');
    set(th,'Position',[.5 1.05 0],'fontsize',20);
    lh = legend([h11,h12],'Low ozone','High ozone');
    set(lh,'box','off','fontsize',20,'orientation','horizontal','location','south');
    lhpos = get(lh,'position');
    set(lh,'position',[lhpos(1) lhpos(2)-.1 lhpos(3) lhpos(4)]);
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/',sprintf('%02d',mons(i,:)),'_',datatype,'_ChangeSEAICE_EXTENT_MarchOzoneExtremes_detrend'];
    export_fig(filename,'-png');
end

%% seas for accumulation metric

% Barents Kara sea; Laptet, East Siberian, and Chukchi seas, Greenland
% Seas, Sea of Okhotsk and Berling Sea

Seas.lon = [20,80;90,200;330,30;140,200];
Seas.lat = [65,80;65,80;65,85;45,65];

for i = 1:size(Seas.lon,1)
    seaslatind = seaicedata2.latitude >= Seas.lat(i,1) & seaicedata2.latitude < Seas.lat(i,2);
    if i == 3
        seaslonind = seaicedata2.longitude >= Seas.lon(i,1) | seaicedata2.longitude < Seas.lon(i,2);
    else
        seaslonind = seaicedata2.longitude >= Seas.lon(i,1) & seaicedata2.longitude < Seas.lon(i,2);
    end
    %calculate total
    fractionextract = seaicedataweighted_reshape_detrend(:,:,seaslonind,seaslatind);
    fractionInt(i,:,:) = nansum(fractionextract(:,:,:),3);
    
end
    %take difference
fractionint_difference = nanmean(fractionInt(:,:,highind),3) - nanmean(fractionInt(:,:,lowind),3);
createfig('medium','on');
ph = plot(fractionint_difference','Linewidth',2);
lh = legend('Barents-Kara','Laptet, East Siberian, Chukchi','Greenland','Okhotsk and Berling');
set(lh,'fontsize',20,'box','off');
%% extract time period 1980 - 2016 from sea ice data

seaicetimeextract = seaicedetrend(12+mon:12:end-12)';
seaicetimeextract_nodetrend = seaicedataintegration(12+mon:12:end-12)';
%% plotting
createfig('medium','on')
plot(1:length([lowind,highind]),seaicetimeextract([lowind,highind]),'-o','Linewidth',3)
ylabel('Detrended sea ice integration','fontsize',24);
yyaxis right
plot(1:length([lowind,highind]),toz_zm([lowind,highind]),'-o','Linewidth',3);
set(gca,'ydir','reverse','fontsize',20)
ylabel('Total column ozone upper and lower 20th percentiles (DU)','fontsize',20);
xlabel('Composite year index','fontsize',20);

createfig('medium','on')
plot(1:length([lowind,highind]),seaicetimeextract_nodetrend([lowind,highind]),'-o','Linewidth',3)
ylabel('Sea ice integration','fontsize',24);
yyaxis right
plot(1:length([lowind,highind]),toz_zm([lowind,highind]),'-o','Linewidth',3);
set(gca,'ydir','reverse','fontsize',20)
ylabel('Total column ozone upper and lower 20th percentiles (DU)','fontsize',20);
xlabel('Composite year index','fontsize',20);