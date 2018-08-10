%Observations and model Figure 1

%% Model --------------------------------------------------------------------------------

%% prediction runs correlations to surface temperature

clear all 
close all
tozvar = 'toz';
var = 'TS';
%area = '6090S';
detrend = 1;
detrend_ozone = 0;
contourplots = 1;
individual_contourplots = 0;
lineplots = 1;
percentile = 20;
lats = [63,90];
%lats = [-90, -63];
tozmonth = 3;
indmonth = 3;
varmonth = [3,4]; % Can be any number of months in the year (e.g. [12,1,2] for austral summer)
shortnames = 0;
sig = .05;
strapme = 0;
ENSO = 0;
removeENSO = 1;
ClLevel = 'highCl';
% read in data 
%directqory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/'];
%files = dir([directory,'*.nc']);

%% Read in highcl variable


timeperiodhigh = [1995,2016];%[1955,1975]
%timeperiodhigh = [1955,1975];%[1955,1975]

vardirectory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl,dataMonthArrangeMean.highcl]...
    = predruns_ReadInlayer(vardirectory,varfiles,var,timeperiodhigh,lats,detrend);

%% Read in TOZ highcl and take percentiles

tozdates = [1995,2016];
%tozdates = [1955,1974];%[1955,1975]
directory = ['/Volumes/MyBook/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);

%%

[pct_highcl,tozextract.highcl,Eachyear] = predruns_varPercentiles(toz_composite.highcl.montharrange,toz_dataMonthArrange.highcl,...
    tozmonth,percentile,length(tozfiles));

longitude = data(1).highcl.lon;
latitude = data(1).highcl.lat;

%%

[~,diffcomp,diffp,compcorr,corrp] = predruns_diffandcorr_percentiles(dataMonthArrange,...
    dataMonthArrangeMean,toz_dataMonthArrange,varmonth,tozmonth,var,latitude,longitude,...
    timeperiodhigh,lats,removeENSO,Eachyear,ClLevel,0);

%% Observations -----------------------------------------------------------------------------

% prediction runs correlations to surface temperature

% recreate Justin's figure

% using ERA-Interim TS and Bodeker scientific TCO
clearvars -except diffcomp diffp compcorr corrp latitude longitude
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

ERAdata = ReadinERA('/Volumes/MyBook/work/data/ERA-Interim/TS/TS_ERA-Interim.nc');
ERAyears = 1979:2016;
ERAyear_vector = repmat(ERAyears,12,1);
ERAdata.years = [ERAyear_vector(:);ones(length(ERAdata.time)-length(ERAyear_vector(:)),1)*max(ERAyear_vector(:))+1];

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
takediff = 1;
if takediff

    %% calculate observed differences
    if tcolat(1) < 0
        areaslon = [140,160;15,30;280,305];
        areaslat = [-38 -25;-35,-15;-40,-22];
        hem = 0;
    else
        areaslon = [70,90;300,320;245,270;70,90];
        areaslat = [65,75;65,78;43,65;25,35];
        hem = 1;
    end    
    monin = [3,4];
    
    [differences,differencesp] = obsPercentileDifferences(ERAdata.t2m,toz_zm_nodetrend,20,monin,tozmonth,ERAdata.longitude,ERAdata.latitude,hem,[]); % 1 = north totaldiff(i)            
    
%     for i = 1:length(monin)
%         [differences2(i).m] = obsPercentileDifferences(ERAdata.t2m,toz_zm_nodetrend,20,monin(i),tozmonth,ERAdata.longitude,ERAdata.latitude,hem,[]); % 1 = north totaldiff(i)            
%     end
    
    %%
%     for i = 1:length(monin)
%         for j = 1:size(areaslon,1)
%             diffarea = squeeze(differences2(i).m(1,ERAdata.longitude >= areaslon(j,1) & ERAdata.longitude <= areaslon(j,2),...
%                 ERAdata.latitude >= areaslat(j,1) & ERAdata.latitude <= areaslat(j,2)));
%             meandiff(i,j) = nanmean(diffarea(:));
%         end
%     end

    %%    
    
%     lnames = {'Russia','Greenland','North America','Himalayas'};
%     filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/compositedifferences/observations/North_',num2str(20),'th_percentile_lines'];
%     xticklab = {'Mar','Apr','May','Jun'};    
%     
%     createfig('medium','on')
%     plot(meandiff,'linewidth',3);
%     hold on
%     plot([0 5],[0,0],'--k','lineWidth',2);
%     xlim([.9,4.1]);
%     set(gca,'fontsize',20,'xtick',1:5,'xticklabel',xticklab)
%     xlabel('Month','fontsize',22);
%     ylabel('Temperature difference (K)','fontsize',22);
%     title('Surface temperature difference between high and low ozone years','fontsize',24)
%     lh = legend(lnames);
%     set(lh,'box','off','fontsize',24)
%     
%     export_fig(filename,'-png');
end
%% calculate ENSO
if ninomonth <= 2
    yearend = 0;
else
    yearend = 0;
end
ERAdateindex(1) = find(ERAdata.years == timeperiod(1),1,'first');
ERAdateindex(2) = find(ERAdata.years == timeperiod(2),1,'last');

latlimits = [-5 5];
lonlimits = [190 240];

latindex = find(ERAdata.latitude >= latlimits(1) & ERAdata.latitude <= latlimits(2));
lonindex = find(ERAdata.longitude >= lonlimits(1) & ERAdata.longitude <= lonlimits(2));

for j = 1:12
    NINO_mn(:,j,:) = squeeze(nanmean(ERAdata.t2m(lonindex,latindex,ERAdateindex(1)+j-1:12:ERAdateindex(2)),1));
    NINO_mn2(j,:) = squeeze(nanmean(NINO_mn(:,j,:),1));
    NINOmonth(j,:) = (NINO_mn2(j,:) - nanmean(NINO_mn2(j,:),2))./std(NINO_mn2(j,:),1,2);
end   

NINO34all = NINOmonth(:);
NINO34 = NINOmonth(ninomonth,:);

%save('/Volumes/MyBook/work/data/predruns/output/NINO34/ERA-Interim_1980_2016','NINO34all','NINO34');
    
for i = 1:12    
    for j = 1:size(ERAdata.t2m,2)
        for k = 1:size(ERAdata.t2m,1)
            ts(k,j,:,i) = detrend(squeeze(ERAdata.t2m(k,j,ERAdateindex(1)+i-1:12:ERAdateindex(2)))) + ...
                repmat(squeeze(nanmean(ERAdata.t2m(k,j,ERAdateindex(1)+i-1:12:ERAdateindex(2)),3)),[1,timeperiod(2)-timeperiod(1)+1])';
        end
    end
    %ts2(:,:,:,i) = squeeze(ERAdata.t2m(:,:,ERAdateindex(1)+varmonth(i):12:ERAdateindex(2)+12));
end
ts = permute(ts,[1,2,4,3]);
ts = ts(:,:,:);
%tsmean = nanmean(ts,4);

%% taking correlations
if tcolat(1) > 0
    endregress = 0;   
else
    endregress = 24-varmonth(end);
end
laglength = 3;
count = 1;
for i = 1:size(ts,1) % lon
    for j = 1:size(ts,2) % lat
        for m = 1:length(varmonth)
            if removeENSO            
                for k = 1:laglength
                    regressors = [ones(length(NINO34all(varmonth(m)-k+1:12:end-k)),1),NINO34all(varmonth(m)-k+1:12:end-k)];
                    [b(i,j,m,k,:),bint(i,j,m,k,:,:),r(i,j,m,k,:),rint(i,j,m,k,:,:),stats(i,j,m,k,:)] = ...
                        regress(squeeze(ts(i,j,varmonth(m):12:end-endregress)),regressors);
                    %[b(i,j,k,:),bint(i,j,:,:),r(i,j,:),rint(i,j,:,:),stats(i,j,:)] = regress(squeeze(tsmean(i,j,1:length(toz_zm))),[ones(length(NINO34),1),NINO34']);
                end
            
        
                %varforcorr = r(i,j,:); 
                [~,llc(i,j,m)] = max(abs(squeeze(b(i,j,m,:,2))));
                blag(i,j,m,:) = b(i,j,m,llc(i,j,m),:);
                varforcorr(i,j,m,:) = squeeze(ts(i,j,varmonth(m):12:end-laglength))' - squeeze(blag(i,j,m,2))*squeeze(NINO34all(varmonth(m)-llc(i,j,m)+1:12:end-laglength))';                        
            
            else
                varforcorr(i,j,m,:) = squeeze(ts(i,j,varmonth(m):12:end-endregress));
            end
        end
        
        if percentiles
            if count == 1
                num = ceil(length(toz_zm_nodetrend)/100*20);
                [toz_sort,toz_sortind] = sort(toz_zm_nodetrend);
                 toz_zmforcorr = toz_zm([toz_sortind(1:num),toz_sortind(end-num+1:end)]);
                count = -1;
            end
            varForCorrMonthMean(i,j,:) = squeeze(nanmean(varforcorr(i,j,:,[toz_sortind(1:num),toz_sortind(end-num+1:end)]),3));
            
        else    
            toz_zmforcorr = toz_zm;
            varForCorrMonthMean(i,j,:) = squeeze(nanmean(varforcorr(i,j,:,:),3));
        end
        [polar_correlations(i,j),polar_p(i,j)] = ...
           corr(squeeze(varForCorrMonthMean(i,j,:)),toz_zmforcorr(1:size(varForCorrMonthMean(i,j,:),3))');        
    end
end

%% interpolate observations onto WACCM grid

for i = 1:length(ERAdata.latitude)
    dataloninterp(:,i) = interp1(ERAdata.longitude,polar_correlations(:,i),longitude);
    dataploninterp(:,i) = interp1(ERAdata.longitude,polar_p(:,i),longitude);
    diffloninterp(:,i) = interp1(ERAdata.longitude,squeeze(differences(1,:,i)),longitude);
    diffploninterp(:,i) = interp1(ERAdata.longitude,squeeze(differencesp(1,:,i)),longitude);
end
%%
for i = 1:length(longitude)
    polar_correlationsinterp(1,i,:) = interp1(ERAdata.latitude,dataloninterp(i,:),latitude);
    polar_correlationsinterp_p(1,i,:) = interp1(ERAdata.latitude,dataploninterp(i,:),latitude);
    differences_interp(1,i,:) = interp1(ERAdata.latitude,diffloninterp(i,:),latitude);
    differences_interp_p(1,i,:) = interp1(ERAdata.latitude,diffploninterp(i,:),latitude);
end

differences_interp_p (differences_interp_p ~= 0) = -1;

%% plotting
polar_correlationsinterp_p (polar_correlationsinterp_p <=.05) = 0;
polar_correlationsinterp_p (polar_correlationsinterp_p >.05) = -1;
%lontoplot = ERAdata.longitude;
%polar_correlations = [polar_correlations(end,:);polar_correlations];
polartoplot = cat(1,polar_correlationsinterp,compcorr);%,differences_interp,diffcomp);
polar_p_toplot = cat(1,polar_correlationsinterp_p,corrp);%,differences_interp_p,diffp);

difftoplot = cat(1,differences_interp,diffcomp);
diff_p_toplot = cat(1,differences_interp_p,diffp);

% sig = .05;
% polar_p_toplot = reshape(polar_p,[1,size(polar_p)]);
% polar_p_toplot (polar_p_toplot <= sig) = 0;
% polar_p_toplot (polar_p_toplot > sig) = -1;


%% plotting polar    
plotting = 1;
fsize = 16;
contourtitle = {'Observed correlations','Ensemble composite correlations'};
contourtitle2 = {'Observed temperature difference','Ensemble composite temperature difference'};
corrtitles = {'Correlation','Temperature'};
intervals = [-1 1;-1 1; -5 5;-5 5];

subplotmaps(polartoplot,longitude,latitude,{'div','RdBu'},1,polar_p_toplot,16,contourtitle,'Longitude','Latitude','Correlation','on',...
        [-1 1],22,[-90:10:90],[-90:10:90],[longitude(2:10:end)],[longitude(2:10:end)],'',1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
    
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/Correlation_comparison'];
 
print(filename,'-depsc');

%export_fig(filename,'-pdf');
%%    
[fig,~] = subplotmaps(difftoplot,longitude,latitude,{'div','RdBu'},1,diff_p_toplot,16,contourtitle2,'Longitude','Latitude','Temperature','on',...
    [-5 5],22,[-90:10:90],[-90:10:90],[longitude(2:10:end)],[longitude(2:10:end)],'',1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

child = get(fig,'children');

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([3,4,10,8],:);
cbrewqual3 = cbrewqual2./2;

axes1 = child(4);
axes2 = child(2);

axes(axes1);

obslats = [70,67,50,30];
obslons = [310-360,80,260-360,80];

Mks = ['p','s','o','v'];

for i = 1:length(obslats)
    m_plot(obslons(i),obslats(i),'LineStyle','none','Marker',Mks(i),'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
end

complats = [66,70,52,35];
complons = [310-360,60,255-360,60];

axes(axes2);

for i = 1:length(obslats)
    m_plot(complons(i),complats(i),'LineStyle','none','Marker',Mks(i),'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
end

filename2 = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/Diff_comparison'];

print(filename2,'-depsc');

%%

%export_fig(filename2,'-pdf')

% 
% if plotting
%     fig = figure;
%     set(fig,'position',[100 100 1000 700],'color','white');       
%     
%     
%     
%     for i = 1:size(polartoplot,1)./2
%         
%         maxcontourinterval = intervals(i,2);
%         mincontourinterval = intervals(i,1);
%     
%         contourstep = (maxcontourinterval - mincontourinterval) / (22-2);
%         
%         contourinterval = mincontourinterval-contourstep:contourstep:maxcontourinterval+contourstep;
% 
%         lcontour = length(contourinterval)-2;
%         
%         colourmap = flipud(cbrewer('div','RdBu',lcontour));
%     colourmap(ceil(lcontour/2),:) = [];
%         
%         sp(i) = subplot(4,1,i);
%         sp_pos(i,:) = get(sp(i),'position');
%         
% %         if i == 1 || i == 3
% %             set(sp(i),'position',[sp_pos(i,1)-.06 sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)]);
% %         else
% %             set(sp(i),'position',[sp_pos(i,1)-.02 sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)]);
% %         end
%                 
%         m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
%         if i == 1
%             x(size(x,1)/2+1:size(x,1)) = x(size(x,1)/2+1:size(x,1)) - 360;
%             x = [x(size(x,1)/2+1:size(x,1));x(1:size(x,1)/2)];
%             data = circshift(polartoplot,size(polartoplot,2)/2,2);
%             pdata = circshift(polar_p_toplot,size(polar_p_toplot,2)/2,2);
%         end
%            
%         [~,h] = m_contourf(x,y,squeeze(data(i,:,:))',contourinterval,'LineStyle','-','LineColor',[.6 .6 .6]);
%         m_coast('color','k','LineWidth',2);
%         m_grid('ytick',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90],'XaxisLocation','bottom','fontsize',fsize);        
%         
%         colormap(colourmap);  
%         
%         caxis([intervals(i,1) intervals(i,2)]);
%         
%         title(contourtitle{i},'fontsize',fsize+2);
%         
%         yl = ylabel('Latitude','fontsize',fsize);
%         ylpos = get(yl,'position');
%         set(yl,'position',[ylpos(1)-.4,ylpos(2),ylpos(3)])
%         if i == 4
%             xlabel('Longitude','fontsize',fsize);
%         end
%         
%         sp_pos(i,:) = get(sp(i),'position');
%         if i == 2
%             ch = colorbar;
%             cbaxloc = get(ch,'Position');
%             
%             
%             set(ch,'Position',[cbaxloc(1)+cbaxloc(1)/13.5, cbaxloc(2),...
%                 cbaxloc(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(2,2)-.005],...
%                 'Box','on','YAxisLocation','right','fontsize',fsize);        
%             
%             set(ch,'YTick',[intervals(i,1):contourstep*2:intervals(i,2)])
%             
%         elseif i == 4
%             ch = colorbar;
%             cbaxloc = get(ch,'Position');
%             
%         
%             set(ch,'Position',[cbaxloc(1)+cbaxloc(1)/13.5, cbaxloc(2),...
%                 cbaxloc(3), sp_pos(3,2)+sp_pos(3,4)-sp_pos(4,2)-.005],...
%                 'Box','on','YAxisLocation','right','fontsize',fsize);  
%             
%             set(ch,'YTick',[intervals(i,1):contourstep*2:intervals(i,2)])%,'YtickLabel',{'',minmax(1):contourstep*2:minmax(2),''});
%             
%         end
%         
%         hold on
%         lathatch = repmat(latitude,[1,size(pdata(i,:,:),2)])';
%         lonhatch = repmat(x,[1,size(pdata(i,:,:),3)]);
%         for k = size(polar_p_toplot,3)                
%             latvec = lathatch(squeeze(pdata(i,:,:)) == 0);
%             lonvec = lonhatch(squeeze(pdata(i,:,:)) == 0);                
% 
%             %latvec = lathatch(squeeze(pdata(1,:,k)) == 0);
%             %lonvec = lonhatch(squeeze(pdata) == 0);                
%         end            
% 
%          m_line(lonvec(2:4:end),latvec(2:4:end),'LineStyle','none','marker','.','color','k')
%         
%     end
%     
% %     subplotmaps(polartoplot,longitude,latitude,{'div','RdBu'},1,polar_p_toplot,16,contourtitle,'Longitude','','Correlation','on',...
% %         clims,22,[-90:10:20],[-90:10:20],[lontoplot(2:10:end)],[lontoplot(2:10:end)],'',1,[0 360],ylims,0,'none',1,'Miller Cylindrical');
%    
% % filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/','correlations/maps/','Dianes_figure_',percext,'_',hemext,'_',locext,'_',rmENSOappend,...
% %     num2str(timeperiod(1)),'-',num2str(timeperiod(2))];    
% 
% % export_fig(filename,'-pdf');
% % export_fig(filename,'-png');       

    
%end