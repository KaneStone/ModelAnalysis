% compare difference between ENSO and no ENSO signals in correlations

clear all
% SHreg = load('/Volumes/MyBook/work/data/predruns/output/regression/regcoefsNovembertoz_TS_detrend1995-2016_90-75S_Tperiod-DecJanFeb_.mat');
% SHreg_rmENSO = load('/Volumes/MyBook/work/data/predruns/output/regression/regcoefsNovembertoz_TS_detrend1995-2016_90-75S_Tperiod-DecJanFeb_rmENSO_lagged.mat');

SHreg = load('/Volumes/MyBook/work/data/predruns/output/regression/regcoefsMarchtoz_TS_detrend1995-2016_63-90S_Tperiod-MarApr_.mat');
SHreg_rmENSO = load('/Volumes/MyBook/work/data/predruns/output/regression/regcoefsMarchtoz_TS_detrend1995-2016_63-90S_Tperiod-MarApr_rmENSO_lagged.mat');

tozmonth = 3;
varmonth = [3,4];

%ENSOdiff_highcl = SHreg.r.ind.highcl - SHreg_rmENSO.r.ind.highcl;
ENSOdiff_highcl = abs(SHreg.r.ind.highcl) - abs(SHreg_rmENSO.r.ind.highcl);

for i = 1:size(SHreg.r.ind.highcl,1)

    subplotmaps(ENSOdiff_highcl(i,:,:),SHreg.lons,SHreg.lats,{'div','RdBu'},1,[],16,{['With ENSO - without ENSO ', 'No.',num2str(i)]},'Longitude','latitude','Correlation','on',...
        [-.5 .5],18,[-90:10:90],[-90:10:90],[SHreg.lons(1:10:end)],[SHreg.lons(1:10:end)],'',1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

    %rtoplot = circshift((abs(squeeze(r.ind.highcl(i,:,:)))'),size((abs(squeeze(r.ind.highcl(i,:,:)))'),2)/2,2);
    %m_contour(x,lats,rtoplot,[.6 .6],'color',[77,175,74]/255,'Linewidth',3,'LineStyle','-');

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/ENSOsignal/ENSO_Ind_','No.',sprintf('%02d',i),monthnames(tozmonth,0,0),'toz','_TS_detrend',...
        'South_1995-2016_abs','Tperiod-',monthnames(varmonth,1,1)];

    export_fig(filename,'-png');
end

%% Observations

obsreg = load('/Volumes/MyBook/work/data/predruns/output/regression/regcoefs_Justins_.mat');
obsreg_rmENSO = load('/Volumes/MyBook/work/data/predruns/output/regression/regcoefs_Justins_rmENSO.mat');

obsENSOdiff = reshape(abs(obsreg.polar_correlations) - abs(obsreg_rmENSO.polar_correlations),[1,size(obsreg.polar_correlations)]);

subplotmaps(obsENSOdiff,obsreg.lons,obsreg.lats,{'div','RdBu'},1,[],16,{'With ENSO - without ENSO'},'Longitude','latitude','Correlation','on',...
    [-.5 .5],18,[-90:10:90],[-90:10:90],[obsreg.lons(1:10:end)],[obsreg.lons(1:10:end)],'',1,[0 360],[-90 0],0,'none',1,'Miller Cylindrical');

    %rtoplot = circshift((abs(squeeze(r.ind.highcl(i,:,:)))'),size((abs(squeeze(r.ind.highcl(i,:,:)))'),2)/2,2);
    %m_contour(x,lats,rtoplot,[.6 .6],'color',[77,175,74]/255,'Linewidth',3,'LineStyle','-');

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/ENSOsignal/OBS_ENSO_Ind_',monthnames(tozmonth,0,0),'toz','_TS_detrend',...
        'South_1995-2016_abs','Tperiod-',monthnames(varmonth,1,1),'lagged'];

    export_fig(filename,'-png');