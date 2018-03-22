% predruns_plotrollingcorr

clear all
%SHreg = load('/Volumes/MyBook/work/data/predruns/output/regression/regcoefsNovembertoz_TS_detrend1995-2024_90-63S_Tperiod-DecJanFeb_rmENSO_lagged_roll20.mat');
SHreg = load('/Volumes/MyBook/work/data/predruns/output/regression/regcoefsMarchtoz_TS_detrend1995-2024_63-90S_Tperiod-MarApr_rmENSO_lagged_roll20.mat');

%% take rollcorr ensemble mean

composite = nanmean(cat(4,SHreg.rollcorr.ind.highcl(:).r),4);
compositep = nanmean(cat(4,SHreg.rollcorr.ind.highcl(:).p),4);

%% Read in ozone
lats = [63,90];
tozmonth = 3;

%% Read in TOZ highcl and take percentiles
tozvar = 'toz';
detrend_ozone = 0;
ClLevel = 'highCl';
tozdates = [1995,2023];
directory = ['/Volumes/MyBook/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);




%%
rollwindow = 20;
ozonetoplot = squeeze(nanmean(toz_dataMonthArrange.highcl(:,11,:),1));

%% plot rolling correlations

    
if lats(1) < 0
    ylims = [-90 0];
    hemext = 'south';
    if rollwindow == 10
        rollyears = 2000:2018;
    elseif rollwindow == 20
        rollyears = 2006:2014;
    end
else
    ylims = [0 90];
    hemext = 'north';
    if rollwindow == 10
        rollyears = 2000:2019;
    elseif rollwindow == 20
        rollyears = 2006:2015;
    end
end

%% First plot maps of 10 year correlations
sig = .05;
toplotrollcorr = permute(composite,[3,1,2]);
toplotp = permute(compositep,[3,1,2]);    
toplotp (toplotp <= sig) = 0;
toplotp (toplotp > sig) = 1;

for i = 1:size(toplotrollcorr,1)

    subplotmaps(toplotrollcorr(i,:,:),SHreg.lons,SHreg.lats,{'div','RdBu'},1,[],16,{num2str(rollyears(i))},'Longitude','','Correlation','on',...
        [-.5 .5],22,[-90:10:20],[-90:10:20],[SHreg.lons(2:10:end)],[SHreg.lons(2:10:end)],'',1,[0 360],ylims,0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/Rolling/maps/','ModelComposite_figure_rolling',num2str(rollwindow),'_',sprintf('%02d',i),'_',hemext,'_',num2str(tozdates(1)),'-',num2str(tozdates(2))];    
    
    export_fig(filename,'-png');
end

%% plotting individual members

for j = 1:length(SHreg.rollcorr.ind.highcl)
    toplotrollcorr = permute(SHreg.rollcorr.ind.highcl(j).r,[3,1,2]);
    for i = 1:size(toplotrollcorr,1)
        subplotmaps(toplotrollcorr(i,:,:),SHreg.lons,SHreg.lats,{'div','RdBu'},1,[],16,{num2str(rollyears(i))},'Longitude','','Correlation','off',...
            [-1 1],22,[-90:10:20],[-90:10:20],[SHreg.lons(2:10:end)],[SHreg.lons(2:10:end)],'',1,[0 360],ylims,0,'none',1,'Miller Cylindrical');
        if ~exist(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/Rolling/maps/Model',sprintf('%02d',j),'/'],'dir')
            mkdir(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/Rolling/maps/Model',sprintf('%02d',j),'/']);
        end            
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/Rolling/maps/Model',sprintf('%02d',j),'/','ModelComposite_figure_rolling',num2str(rollwindow),'_',sprintf('%02d',i),'_',hemext,'_',num2str(tozdates(1)),'-',num2str(tozdates(2))];    

        export_fig(filename,'-png');
    end
end

%% 
plotroll = 1;
if plotroll
    if lats(1) > 0
        loc(1).lats = [28,16,15]; loc(2).lats = [30,40]; loc(3).lats = [47,56,71]; loc(4).lats = [64,53];
        loc(1).lons = [40,92,125]; loc(2).lons = [228,352]; loc(3).lons = [262,295,309]; loc(4).lons = [100,167];   
    else
        loc(1).lats = [-35,-30,-22]; loc(2).lats = [-29,-21]; loc(3).lats = [-61,-55,-60]; loc(4).lats = [-50,-48];
        loc(1).lons = [130,295,27]; loc(2).lons = [210,165]; loc(3).lons = [62,140,242]; loc(4).lons = [327,57];   
    end
    titles = {'Positive correlations - tropical and mid-latitudes','Negative correlations - tropical and mid-latitudes',...
        'Positive correlations - polar-latitudes','Negative correlations - polar-latitudes'};
    
    plotRollingCorrelations(composite,ozonetoplot,loc,length(loc),[2,2],...
    titles,18,rollwindow,rollyears,SHreg.lats,SHreg.lons,'highcl')    
end


%%
% %%
% rollyears = 2003:2016;
% composite = nanmean(cat(4,SHreg.rollcorr.ind.highcl(:).r),4);   
% createfig('largeportrait','on')
% cbrewqual = cbrewer('qual','Set1',10);
% %-----------midlat-positive----------------------
% 
% fsize = 14;
% 
% latstoplot = [-35,-35,-30];
% lonstoplot = [135,20,265];
% 
% subplot(2,2,1);
% for i = 1:length(latstoplot)
%     [~,latindex] = min(abs(SHreg.lats - latstoplot(i))); 
%     [~,lonindex] = min(abs(SHreg.lons - lonstoplot(i)));
%             
%     ph(i) = plot(rollyears,squeeze(composite(lonindex,latindex,:)),'LineWidth',2,'color',cbrewqual(i,:));
%     hold on
%     forleg{i} = [num2str(latstoplot(i)),'N, ',num2str(lonstoplot(i)),'E'];
% end
% ylabel('Rolling correlation','fontsize',fsize+2);
% yyaxis right
% set(gca,'YDir','reverse')
% plot(rollyears,ozonesmooth(ceil(rollwindow/2):end-ceil(rollwindow/2)),'LineWidth',3,'color','k');
% forleg{i+1} = 'Normalized ozone';
% lh = legend(forleg);    
% set(lh,'box','off','location','SouthEast','fontsize',fsize)
% 
% set(gca,'fontsize',fsize);
% 
% xlabel('Year','fontsize',fsize+2);
% title('Positive correlations - mid-latitudes','fontsize',fsize+2);
% 
% %------------midlat-negative---------------------
% latstoplotneg = [-30,-22];
% lonstoplotneg = [210,165];
% subplot(2,2,2);
% for i = 1:length(latstoplotneg)
%     [~,latindex] = min(abs(SHreg.lats - latstoplotneg(i))); 
%     [~,lonindex] = min(abs(SHreg.lons - lonstoplotneg(i)));
%             
%     ph(i) = plot(rollyears,squeeze(composite(lonindex,latindex,:)),'LineWidth',2,'color',cbrewqual(i,:));
%     hold on
%     forleg{i} = [num2str(latstoplotneg(i)),'N, ',num2str(lonstoplotneg(i)),'E'];
% end
% yyaxis right
% %set(gca,'YDir','reverse')
% plot(rollyears,ozonesmooth(ceil(rollwindow/2):end-ceil(rollwindow/2)),'LineWidth',3,'color','k');
% forleg{i+1} = 'Normalized ozone';
% lh = legend(forleg);    
% set(lh,'box','off','location','NorthEast','fontsize',fsize)
% 
% ylabel('Normalized TOZ','fontsize',fsize+2);
% title('Negative correlations - mid-latitudes','fontsize',fsize+2);
% %------------polar-positive---------------------
% latstoplot = [-63,-55,-61];
% lonstoplot = [25,140,240];
% 
% subplot(2,2,3);
% for i = 1:length(latstoplot)
%     [~,latindex] = min(abs(SHreg.lats - latstoplot(i))); 
%     [~,lonindex] = min(abs(SHreg.lons - lonstoplot(i)));
%         
%     yyaxis left
%     ph(i) = plot(rollyears,squeeze(composite(lonindex,latindex,:)),'LineWidth',2,'color',cbrewqual(i,:));
%     hold on
%     forleg{i} = [num2str(latstoplot(i)),'N, ',num2str(lonstoplot(i)),'E'];
% end
% ylabel('Rolling correlation','fontsize',fsize+2);
% yyaxis right
% set(gca,'YDir','reverse')
% plot(rollyears,ozonesmooth(ceil(rollwindow/2):end-ceil(rollwindow/2)),'LineWidth',3,'color','k');
% forleg{i+1} = 'Normalized ozone';
% lh = legend(forleg);    
% set(lh,'box','off','location','SouthEast','fontsize',fsize)
% xlabel('Year','fontsize',fsize+2);
% title('Positive correlations - polar-latitudes','fontsize',fsize+2);
% 
% %------------polar-negative---------------------
% latstoplotneg = [-55,-43];
% lonstoplotneg = [310,105];
% subplot(2,2,4);
% for i = 1:length(latstoplotneg)
%     [~,latindex] = min(abs(SHreg.lats - latstoplotneg(i))); 
%     [~,lonindex] = min(abs(SHreg.lons - lonstoplotneg(i)));
%         
%     yyaxis left
%     ph(i) = plot(rollyears,squeeze(composite(lonindex,latindex,:)),'LineWidth',2,'color',cbrewqual(i,:));
%     hold on
%     forleg{i} = [num2str(latstoplotneg(i)),'N, ',num2str(lonstoplotneg(i)),'E'];
% end
% yyaxis right
% plot(rollyears,ozonesmooth(ceil(rollwindow/2):end-ceil(rollwindow/2)),'LineWidth',3,'color','k');
% forleg{i+1} = 'Normalized ozone';
% lh = legend(forleg);    
% set(lh,'box','off','location','NorthEast','fontsize',fsize)
% 
% ylabel('Normalized TOZ','fontsize',fsize+2);
% xlabel('Year','fontsize',fsize+2);
% title('Negative correlations - polar-latitudes','fontsize',fsize+2);
% 
% annotation('textbox',[0 .905 1 .1],'String','Highcl Southern Hemisphere - selected rolling correlations 15 year window, 1995-2025',...
%     'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,'EdgeColor','none','fontweight','bold')
% 
% filename = '/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/correlations/Rolling/Highncl_SouthRolling_15year_1979-2016';
% export_fig(filename,'-pdf');
