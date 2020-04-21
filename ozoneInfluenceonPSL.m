% Look at PSL and 500hPa temperature change due to ozone extremes.

%% Read in PSL
directory = ['/Volumes/ExternalOne/work/data/predruns/PSL/highCl/'];
files = dir([directory,'*.nc']);

[PSL.data,PSL.years,PSL.composite,...
            PSL.dataMonthArrange,PSL.dataMonthArrangeMean]...
            = predruns_ReadInlayer(directory,files,'PSL',[1995,2024],1);
        
        
%% Read in 500 hPa temperature
directory = ['/Volumes/ExternalOne/work/data/predruns/T/highCl/500hPa/'];
files = dir([directory,'*.nc']);

[T500.data,T500.years,T500.composite,...
            T500.dataMonthArrange,T500.dataMonthArrangeMean]...
            = predruns_ReadInlayer(directory,files,'T',[1995,2024],1);
        
%% Read in ozone
directory = ['/Volumes/ExternalOne/work/data/predruns/toz/highCl/'];
files = dir([directory,'*.nc']);

[toz.data,toz.years,toz.composite,...
            toz.dataMonthArrange,toz.dataMonthArrangeMean]...
            = predruns_ReadInlayer(directory,files,'toz',[1995,2024],1);
        
%% calculate toz percentiles

for i = 1:9
    tozlonmean(i,:,:,:) = squeeze(nanmean(toz.dataMonthArrange(i+1,3,:,:,:),5));
    tozweighted(i,:) = weightedaverage(squeeze(tozlonmean(i,:,82:96))',toz.data(1).lat(82:96));
    upper(i) = prctile(tozweighted(i,:),80);
    lower(i) = prctile(tozweighted(i,:),80);
    upperpct(i,:) = find(tozweighted(i,:) >= upper(i));
    lowerpct(i,:) = find(tozweighted(i,:) <= lower(i));
end

%% extract PSL
mons = [7,8,9];
clearvars PSLupper PSLlower

for i = 1:9
    PSLupper(i,:,:,:) = squeeze(nanmean(PSL.dataMonthArrange(i,mons,upperpct(i,:),:,:),2));
    PSLlower(i,:,:,:) = squeeze(nanmean(PSL.dataMonthArrange(i,mons,lowerpct(i,:),:,:),2));
end

PSLupper = permute(PSLupper,[3,4,1,2]);
PSLlower = permute(PSLlower,[3,4,1,2]);

PSLuppermean = nanmean(PSLupper(:,:,:),3);
PSLlowermean = nanmean(PSLlower(:,:,:),3);

PSLdifference = (PSLuppermean - PSLlowermean)./100;

%% extract T500
mons = [7,8,9];
clearvars T500upper T500lower

for i = 1:9
    T500upper(i,:,:,:) = squeeze(nanmean(T500.dataMonthArrange(i,mons,upperpct(i,:),:,:),2));
    T500lower(i,:,:,:) = squeeze(nanmean(T500.dataMonthArrange(i,mons,lowerpct(i,:),:,:),2));
end

T500upper = permute(T500upper,[3,4,1,2]);
T500lower = permute(T500lower,[3,4,1,2]);

T500uppermean = nanmean(T500upper(:,:,:),3);
T500lowermean = nanmean(T500lower(:,:,:),3);

T500difference = (T500uppermean - T500lowermean);

%% plotting
toplot = cat(2,PSLdifference,PSLdifference(:,2));
lontoplot = [PSL.data(1).lon;360];
createfig('largesquare','on')
m_proj('Stereographic','lon',0,'lat',90,'rad',60);    
axm = axesm('MapProjection','stereo','MapLatLimit',[30 90]);%,'Grid','On','MLineLocation',60,...
contourfm(PSL.data(1).lat,lontoplot,toplot,-2:.2:2);
%caxis([-10,10])
framem; axis off; tightmap; 
load coastlines
geoshow(coastlat, coastlon, 'Color', 'black','LineWidth',3)
%patchm(coastlat,coastlon,[.5 .5 .5])  
colorbar
%% plotting
toplot = cat(2,T500difference,T500difference(:,2));
lontoplot = [PSL.data(1).lon;360];
createfig('largesquare','on')
m_proj('Stereographic','lon',0,'lat',90,'rad',60);    
axm = axesm('MapProjection','stereo','MapLatLimit',[30 90]);%,'Grid','On','MLineLocation',60,...
contourfm(PSL.data(1).lat,lontoplot,toplot);
%caxis([-10,10])
framem; axis off; tightmap; 
load coastlines
geoshow(coastlat, coastlon, 'Color', 'black','LineWidth',3)
%patchm(coastlat,coastlon,[.5 .5 .5])  
colorbar
