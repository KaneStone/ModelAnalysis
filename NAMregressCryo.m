% Read in and calculate NAM and regress onto sea ice and snow depth

%% Read in PSL
directory = ['/Volumes/ExternalOne/work/data/predruns/PSL/highCl/'];
files = dir([directory,'*.nc']);

[PSL.data,PSL.years,PSL.composite,...
            PSL.dataMonthArrange,PSL.dataMonthArrangeMean]...
            = predruns_ReadInlayer(directory,files,'PSL',[1995,2024],1);
        
        
%% Read in sea ice
directory = ['/Volumes/ExternalOne/work/data/predruns/ICEFRAC/highCl/'];
files = dir([directory,'*.nc']);

[ICEFRAC.data,ICEFRAC.years,ICEFRAC.composite,...
            ICEFRAC.dataMonthArrange,ICEFRAC.dataMonthArrangeMean]...
            = predruns_ReadInlayer(directory,files,'ICEFRAC',[1995,2024],1);
        
%% Read in snow

directory = ['/Volumes/ExternalOne/work/data/predruns/SNOWHICE/highCl/'];
files = dir([directory,'*.nc']);

[SNOWHICE.data,SNOWHICE.years,SNOWHICE.composite,...
            SNOWHICE.dataMonthArrange,SNOWHICE.dataMonthArrangeMean]...
            = predruns_ReadInlayer(directory,files,'SNOWHICE',[1995,2024],1);
        
directory = ['/Volumes/ExternalOne/work/data/predruns/SNOWHLND/highCl/'];
files = dir([directory,'*.nc']);

[SNOWHLND.data,SNOWHLND.years,SNOWHLND.composite,...
            SNOWHLND.dataMonthArrange,SNOWHLND.dataMonthArrangeMean]...
            = predruns_ReadInlayer(directory,files,'SNOWHLND',[1995,2024],1);
for i = 1:9
    SNOWDEPTH(i).data = SNOWHLND.data(i+1).SNOWHLND;% + SNOWHICE.data(i+1).SNOWHICE;                
end

%%

[~,lf,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/predruns/landfrac/LANDFRAC_b.e11.BWTREFC2.f19_g16.ccmi32.001.cam.h0.2097-08.nc');
lfrep = repmat(lf.LANDFRAC,[1,1,360]);     

%% calculate NAM through EOF

latitude = PSL.data(1).lat;
longitude = PSL.data(1).lon;

lats = [30,90];
latind = PSL.data(1).lat >= lats(1) & PSL.data(1).lat <= lats(2);

eofmonth = [1];

for i = 1:9
    
% month arrange
for j = 1:12
    PSL.montharrange(i,:,:,j,:) = PSL.data(i).PSL(:,:,j:12:end)-nanmean(PSL.data(i).PSL(:,:,j:12:end),3);
    ICEFRAC.montharrange(i,:,:,j,:) = ICEFRAC.data(i).ICEFRAC(:,:,j:12:end)-nanmean(ICEFRAC.data(i).ICEFRAC(:,:,j:12:end),3);
    SNOW.montharrange(i,:,:,j,:) = SNOWDEPTH(i).data(:,:,j:12:end) - nanmean(SNOWDEPTH(i).data(:,:,j:12:end),3);
end
PSL.latextract(i,:,:,:,:) = PSL.montharrange(i,:,latind,:,:);
ICEFRAC.latextract(i,:,:,:,:) = ICEFRAC.montharrange(i,:,latind,:,:);
SNOW.latextract(i,:,:,:,:) = SNOW.montharrange(i,:,latind,:,:);
    
[eof_maps(i,:,:,:),pc(i,:,:),expvar(i,:)]=eof(squeeze(nanmean(PSL.latextract(i,:,:,eofmonth,:),4)),4); %maybe need to weight here?

% normalize 


%pcout(i,:) = squeeze(pc(i,1,:))'./std(squeeze(pc(i,1,:))',0,2);  
pcout(i,:) = (squeeze(pc(i,1,:))'-nanmean(squeeze(pc(i,1,:))'))./std(squeeze(pc(i,1,:))',0,2);  
%pcout(i,:) = squeeze(pc(i,1,:))'./max(squeeze(pc(i,1,:))',[],2);


if eof_maps(i,10,end,1) >= 0
    eof_maps(i,:,:,1) = eof_maps(i,:,:,1).*-1;
    pcout = pcout .* -1;
end
end

% eof_maps([2,3,9],:,:,:) = eof_maps([2,3,9],:,:,:)*-1;
% pcout([2,3,9],:) = pcout([2,3,9],:)*-1;
% pcout = pcout .* -1;
% eof_maps = eof_maps.*-1;

%%
createfig('largesquare','on')
for i = 1:9
    subplot(3,3,i)
    m_proj('Stereographic','lon',0,'lat',90,'rad',60);    
    axm = axesm('MapProjection','stereo','MapLatLimit',[30 90]);%,'Grid','On','MLineLocation',60,...
        %'PlineLocation',15,'MLabelLocation',60,'PLabelLocation',[60,75],'MeridianLabel',...
        %'on','ParallelLabel','on','MLabelParallel',45,'LabelRotation','on','fontsize',18); 

    contourfm(latitude(latind),longitude,squeeze(eof_maps(i,:,:,1))');

    framem; axis off; tightmap; 
    load coastlines
    geoshow(coastlat, coastlon, 'Color', 'black')
    patchm(coastlat,coastlon,[.5 .5 .5])  
end

%%
toplot = squeeze(nanmean(eof_maps(:,:,:,1),1))';
toplot = cat(2,toplot,toplot(:,1));
lontoplot = [longitude;360];
createfig('largesquare','on')
m_proj('Stereographic','lon',0,'lat',90,'rad',60);    
axm = axesm('MapProjection','stereo','MapLatLimit',[30 90]);%,'Grid','On','MLineLocation',60,...
contourfm(latitude(latind),lontoplot,toplot);
framem; axis off; tightmap; 
load coastlines
geoshow(coastlat, coastlon, 'Color', 'black')
%patchm(coastlat,coastlon,[.5 .5 .5])  
%% now regress onto ICEFRAC and SNOWDEPTH
cryomonth = [4,5,6];
monthextract = squeeze(nanmean(ICEFRAC.latextract(:,:,:,cryomonth,:),4));
monthextractsnow = squeeze(nanmean(SNOW.latextract(:,:,:,cryomonth,:),4));
lfrep2 = lfrep(:,96-31:96,1);
for i = 1:9
    predictors = [ones(30,1),pcout(i,:)'];
    for j = 1:size(ICEFRAC.latextract,2)
        for k = 1:size(ICEFRAC.latextract,3)
            b.ICEFRAC(i,j,k,:) = regress(squeeze(monthextract(i,j,k,:)),predictors);
            model.ICEFRAC(i,j,k,:) = b.ICEFRAC(i,j,k,2)*predictors(:,2)+b.ICEFRAC(i,j,k,1);
            b.SNOW(i,j,k,:) = regress(squeeze(monthextractsnow(i,j,k,:)),predictors);
            model.SNOW(i,j,k,:) = b.SNOW(i,j,k,2)*predictors(:,2)+b.SNOW(i,j,k,1);
        end
    end
end
    

b.ICEFRACmean = squeeze(nanmean(b.ICEFRAC(:,:,:,2),1));
b.SNOWmean = squeeze(nanmean(b.SNOW(:,:,:,2),1));
b.SNOWmean (lfrep2 < .2) = NaN; 
cbrew = cbrewer('div','RdBu',20);
toplot = cat(1,b.ICEFRACmean,b.ICEFRACmean(1,:));
lontoplot = [longitude;360];
createfig('largesquare','on')
m_proj('Stereographic','lon',0,'lat',90,'rad',45);    
axm = axesm('MapProjection','stereo','MapLatLimit',[45 90]);%,'Grid','On','MLineLocation',60,...
contourfm(latitude(latind),lontoplot,toplot',-.03:.003:.03,'LineStyle','none');
caxis([-.03,.03])
colormap(cbrew);
ch = colorbar;
set(get(ch,'ylabel'),'string','Sea ice fraction (s.d.)^-^{1}','fontsize',20+2)
framem; axis off; tightmap; 
load coastlines
geoshow(coastlat, coastlon, 'Color', 'k')
patchm(coastlat,coastlon,[.5 .5 .5])  

set(gca,'fontsize',20)
title('April-May-June sea ice fraction response to January NAM');

print(gcf,'/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePredCryosphere/Draft/Revision/figuresForResponse/JanNAMicefrac','-depsc','-painters');

%%
cbrew = cbrewer('div','RdBu',20);
toplot = cat(1,b.SNOWmean,b.SNOWmean(1,:));
toplot (toplot < -.005) = NaN;
createfig('largesquare','on')
m_proj('Stereographic','lon',0,'lat',90,'rad',45);    
axm = axesm('MapProjection','stereo','MapLatLimit',[45 90]);%,'Grid','On','MLineLocation',60,...
contourfm(latitude(latind),lontoplot,toplot',-.01:.001:.01,'LineStyle','none');
caxis([-.01,.01])
colormap(cbrew);
ch = colorbar;
set(get(ch,'ylabel'),'string','Snow Depth (m) (s.d.)^-^{1}','fontsize',20+2)
framem; axis off; tightmap; 
load coastlines
geoshow(coastlat, coastlon, 'Color', 'k','LineWidth',2)
%patchm(coastlat,coastlon,[.5 .5 .5])  

set(gca,'fontsize',20)
title('April-May-June snow depth response to January NAM');

print(gcf,'/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePredCryosphere/Draft/Revision/figuresForResponse/JanNAMsnowdepth','-depsc','-painters');
