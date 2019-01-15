% compare difference between ENSO and no ENSO signals in correlations

clear all
SHreg = load('/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefsNovembertoz_TS_detrend1995-2024_90-75S_Tperiod-DecJanFeb_.mat');
SHregshort = load('/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefsNovembertoz_TS_detrend1995-2016_90-75S_Tperiod-DecJanFeb_.mat');

obsreg = load('/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefs_Justins_1979-2013_.mat'); 
obsshort = load('/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefs_Justins_1979-2000_.mat');

tozmonth = 11;
varmonth = [12,1,2];

rolltrends = cat(4,SHreg.rolltrend.ind.highcl(:).r);
rolltrends = permute(squeeze(rolltrends(:,:,2,:)),[3,1,2]);
%diff_highcl = abs(SHreg.r.ind.highcl) - abs(SHregshort.r.ind.highcl);
diff =  abs(rolltrends*15) - abs(SHregshort.r.ind.highcl);
diffobs = abs(obsreg.rolltrend.r(:,:,2)*15) - abs(obsshort.polar_correlations);  
diffobs = reshape(diffobs,[1,size(diffobs)]);
%%
model = 1;
if model

    for i = 1:size(SHreg.r.ind.highcl,1)

        subplotmaps(diff(i,:,:),SHreg.lons,SHreg.lats,{'div','RdBu'},1,[],16,{['With ENSO - without ENSO ', 'No.',num2str(i)]},'Longitude','latitude','r/year','on',...
            [-1 1],18,[-90:10:90],[-90:10:90],[SHreg.lons(1:10:end)],[SHreg.lons(1:10:end)],'',1,[0 360],[-90 0],0,'none',1,'Miller Cylindrical');

        %rtoplot = circshift((abs(squeeze(r.ind.highcl(i,:,:)))'),size((abs(squeeze(r.ind.highcl(i,:,:)))'),2)/2,2);
        %m_contour(x,lats,rtoplot,[.6 .6],'color',[77,175,74]/255,'Linewidth',3,'LineStyle','-');

        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/consistentcorr/cc_','No.',sprintf('%02d',i),monthnames(tozmonth,0,0),'toz','_TS_detrend',...
            'South_1995-2024_','Tperiod-',monthnames(varmonth,1,1),'_','15year_rolltrend_withENSO_diff'];

        export_fig(filename,'-png');
    end
end
%% 
obs = 1;
if obs
    subplotmaps(diffobs,obsreg.lon,obsreg.lat,{'div','RdBu'},1,[],16,{'Correlation trend differences'},'Longitude','latitude','r/year','on',...
            [-1 1],18,[-90:10:90],[-90:10:90],[obsreg.lon(1:10:end)],[obsreg.lon(1:10:end)],'',1,[0 360],[-90 0],0,'none',1,'Miller Cylindrical');

        %rtoplot = circshift((abs(squeeze(r.ind.highcl(i,:,:)))'),size((abs(squeeze(r.ind.highcl(i,:,:)))'),2)/2,2);
        %m_contour(x,lats,rtoplot,[.6 .6],'color',[77,175,74]/255,'Linewidth',3,'LineStyle','-');

        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/consistentcorr/Obs_','No.',sprintf('%02d',i),monthnames(tozmonth,0,0),'toz','_TS_detrend',...
            'South_1979-2013_','Tperiod-',monthnames(varmonth,1,1),'_','15year_rolltrend_withENSO_diff'];

        export_fig(filename,'-png');
end
