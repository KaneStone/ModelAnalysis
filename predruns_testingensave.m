% testing ensave

clear all

[~,ENSAVE_TS,~] = Read_in_netcdf('/Volumes/MyBook/work/data/predruns/TS/highcl/ensave/TS_b.e11.BWTREFC2.f19_g16.ccmi34.HighCl.ENSAVE.cam.h0.nc');

[~,ENSAVE_TOZ,~] = Read_in_netcdf('/Volumes/MyBook/work/data/predruns/TOZ/highcl/ensave/TOZ_b.e11.BWTREFC2.f19_g16.ccmi34.HighCl.ENSAVE.nc');

%%
lats = [75,90];
TSmean = nanmean(cat(4,ENSAVE_TS.TS(:,:,3:12:end-120),ENSAVE_TS.TS(:,:,4:12:end-120)),4);
Tozextract = ENSAVE_TOZ.toz(:,ENSAVE_TOZ.lat >= lats(1) & ENSAVE_TOZ.lat <= lats(2),3:12:end-120);

for i = 1:size(TSmean,1)
    Tozmean = weightedaverage(squeeze(Tozextract(i,:,:)),ENSAVE_TOZ.lat(ENSAVE_TOZ.lat >= lats(1) & ENSAVE_TOZ.lat <= lats(2)));
end

for i = 1:size(TSmean,1)
    for j = 1:size(TSmean,2)
        r(1,i,j) = corr(detrend(squeeze(TSmean(i,j,:))),detrend(Tozmean'));
    end
end


%%
%contourf(ENSAVE_TOZ.lon,ENSAVE_TOZ.lat(49:end),r(:,49:end)',18);

cbrew = cbrewer('div','RdBu',16);                

contourtitle = {'testing'};

subplotmaps(r,ENSAVE_TOZ.lon,ENSAVE_TOZ.lat,{'div','RdBu'},1,[],12,contourtitle,'Longitude','','Correlation','on',...
    [-.75,.75],18,[-90:10:20],[-90:10:20],[ENSAVE_TOZ.lon(2:10:end)],[ENSAVE_TOZ.lon(2:10:end)],'',1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

%filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/','Justins_figure19792013.png'];

%export_fig(filename,'-png','-r300');