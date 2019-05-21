function [corr2corr,corr2diff,ERAdata,toz_zm_nodetrend] = predruns_compareERA(obstimeperiod,modelcorrelations,modeldifferences,inputs,waclon,waclat)


if inputs.removeENSO
    ensoext = 'removeENSO';
else
    ensoext = 'noremoveENSO';
end

ifmore = inputs.varmonth > 12;

if sum(ifmore) > 0
    ext = 12;
    extts = 4;
else
    ext = 0;
    extts = 4;
end

%% Observations -----------------------------------------------------------------------------
% using ERA-Interim TS and Bodeker scientific TCO
if inputs.timeperiodvar(2) < 2016
    timeperiod = [1995,2015];
    timeperiodbs = [1995,2015];
elseif inputs.timeperiodvar(1) > 1993
    timeperiod = [1990,2016];
    timeperiodbs = [1990,2016];
else
    timeperiod = [1980,2016];
    timeperiodbs = [1980,2016];
end

useBodeker = 1;
%% Read in ERA-interim surface temperature

ERAdata = ReadinERA('/Volumes/ExternalOne/work/data/ERA-Interim/TS/TS_ERA-Interim.nc');
ERAyears = 1979:2016;
ERAyear_vector = repmat(ERAyears,12,1);
ERAdata.years = [ERAyear_vector(:);ones(length(ERAdata.time)-length(ERAyear_vector(:)),1)*max(ERAyear_vector(:))+1];


if useBodeker
    %% READ in Bodeker 
    tcolat = inputs.lats;
  
    tcolon = 26;
    tcolat_QBO = [10,30];
    [~,BSdata,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/BodekerScientific/TCO/Bodeker_TCO_monavg.nc');

    BSyears = 1980:2016;
    BSyear_vector = repmat(BSyears,12,1);
    BSdata.years = [BSyear_vector(:);ones(length(BSdata.time)-length(BSyear_vector(:)),1)*max(BSyear_vector(:))+1];

    ozonedateindex(1) = find(BSdata.years == obstimeperiod(1),1,'first');
    ozonedateindex(2) = find(BSdata.years == obstimeperiod(2),1,'last');
    
    latindex_bs = find(BSdata.lat >= tcolat(1) & BSdata.lat <= tcolat(2));
    
    toz_zm = detrend(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+...
        inputs.tozmonth-1:12:ozonedateindex(2))-ext,1)),BSdata.lat(latindex_bs))) + ...
        nanmean(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+...
        inputs.tozmonth-1:12:ozonedateindex(2)-ext),1)),BSdata.lat(latindex_bs)));
    toz_zm_nodetrend = weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+...
        inputs.tozmonth-1:12:ozonedateindex(2)-ext),1)),BSdata.lat(latindex_bs));        

    
end
ERAdateindex(1) = find(ERAdata.years == obstimeperiod(1),1,'first');
ERAdateindex(2) = find(ERAdata.years == obstimeperiod(2),1,'last');
ERAextract = ERAdata.t2m(:,:,ERAdateindex(1):ERAdateindex(2));

if inputs.justextractobs
    corr2corr = [];
    corr2diff = [];
    return
end


if inputs.takediff

    %% calculate ENSO

    %ERAdateindex(1) = find(ERAdata.years == obstimeperiod(1),1,'first');
    %ERAdateindex(2) = find(ERAdata.years == obstimeperiod(2),1,'last');

    latlimits = [-5 5];
    lonlimits = [190 240];

    latindex = find(ERAdata.latitude >= latlimits(1) & ERAdata.latitude <= latlimits(2));
    lonindex = find(ERAdata.longitude >= lonlimits(1) & ERAdata.longitude <= lonlimits(2));

    for j = 1:12
        NINO_mn(:,j,:) = squeeze(nanmean(ERAextract(lonindex,latindex,j:12:end),1));
        NINO_mn2(j,:) = squeeze(nanmean(NINO_mn(:,j,:),1));
        NINOmonth(j,:) = (NINO_mn2(j,:) - nanmean(NINO_mn2(j,:),2))./std(NINO_mn2(j,:),1,2);
    end   

    NINO34all = NINOmonth(:);

    for i = 1:12    
        for j = 1:size(ERAdata.t2m,2)
            for k = 1:size(ERAdata.t2m,1)
                ts(k,j,:,i) = detrend(squeeze(ERAextract(k,j,i:12:end))) + ...
                    repmat(squeeze(nanmean(ERAextract(k,j,i:12:end),3)),...
                    [1,obstimeperiod(2)-obstimeperiod(1)+1])';
            end
        end    
    end
    ts = permute(ts,[1,2,4,3]);
    ts = ts(:,:,:);
    
    
    %% calculate observed differences  
    
    ifmore = inputs.varmonth > 12;

    if sum(ifmore) > 0
        ext = 6;
    else
        ext = 0;
    end
    
    monin = inputs.varmonthtomean;
    if ~exist(['/Volumes/ExternalOne/work/data/predruns/output/data/obs/','obs_perc_diff',monthnames(monin,1,1),'_and_',monthnames(inputs.varmonth,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',ensoext,'.mat'])
        [differences,differencesp,obspct,obstemp] = obsPercentileDifferences(ts,toz_zm_nodetrend,20,...
            monin,inputs.tozmonth,ERAdata.longitude,ERAdata.latitude,[],ext); 
        for i = 1:length(inputs.varmonth)
            [differences_ind(i,:,:),differencesp_ind(i,:,:),~,~] = obsPercentileDifferences(ts,toz_zm_nodetrend,20,...
                inputs.varmonth(i),inputs.tozmonth,ERAdata.longitude,ERAdata.latitude,[],ext); 
        end
        save(['/Volumes/ExternalOne/work/data/predruns/output/data/obs/','obs_perc_diff',monthnames(monin,1,1),'_and_',...
            monthnames(inputs.varmonth,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),...
            '_',ensoext,'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))],...
            'differences','differencesp','differences_ind','differencesp_ind','obspct','obstemp','toz_zm_nodetrend');
    else
        load(['/Volumes/ExternalOne/work/data/predruns/output/data/obs/','obs_perc_diff',...
            monthnames(monin,1,1),'_and_',monthnames(inputs.varmonth,1,1),'_',...
            num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',ensoext,'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'.mat']);
    end



%% taking correlations
if tcolat(1) > 0
    endregress = 0;   
else
    endregress = 24-inputs.varmonthtomean(end);
end
laglength = 3;
count = 1;
for i = 1:size(ts,1) % lon
    for j = 1:size(ts,2) % lat
        for m = 1:length(inputs.varmonthtomean)
            if inputs.removeENSO            
                for k = 1:laglength
                    regressors = [ones(length(NINO34all(inputs.varmonthtomean(m)-k+1:12:end-k)),1),...
                        NINO34all(inputs.varmonthtomean(m)-k+1:12:end-k)];
                    [b(i,j,m,k,:),bint(i,j,m,k,:,:),r(i,j,m,k,:),rint(i,j,m,k,:,:),stats(i,j,m,k,:)] = ...
                        regress(squeeze(ts(i,j,inputs.varmonthtomean(m):12:end-endregress)),regressors);                   
                end
                                    
                [~,llc(i,j,m)] = max(abs(squeeze(b(i,j,m,:,2))));
                blag(i,j,m,:) = b(i,j,m,llc(i,j,m),:);
                varforcorr(i,j,m,:) = squeeze(ts(i,j,inputs.varmonthtomean(m):12:end-laglength))' - ...
                    squeeze(blag(i,j,m,2))*squeeze(NINO34all(inputs.varmonthtomean(m)-llc(i,j,m)+1:12:end-laglength))';                        
            
            else
                varforcorr(i,j,m,:) = squeeze(ts(i,j,inputs.varmonthtomean(m):12:end-endregress));
            end
        end
        
        if inputs.obsperc
            if count == 1
                num = ceil(length(toz_zm_nodetrend)/100*inputs.percentile);
                [toz_sort,toz_sortind] = sort(toz_zm_nodetrend);
                 toz_zmforcorr = toz_zm([toz_sortind(1:num),toz_sortind(end-num+1:end)]);
                count = -1;
            end
            varForCorrMonthMean(i,j,:) = squeeze(nanmean(varforcorr(i,j,:,[toz_sortind(1:num),...
                toz_sortind(end-num+1:end)]),3));
            
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
    dataloninterp(:,i) = interp1(ERAdata.longitude,polar_correlations(:,i),waclon);
    dataploninterp(:,i) = interp1(ERAdata.longitude,polar_p(:,i),waclon);
    diffloninterp(:,i) = interp1(ERAdata.longitude,squeeze(differences(1,:,i)),waclon);
    diffploninterp(:,i) = interp1(ERAdata.longitude,squeeze(differencesp(1,:,i)),waclon);
end
%%
for i = 1:length(waclon)
    polar_correlationsinterp(1,i,:) = interp1(ERAdata.latitude,dataloninterp(i,:),waclat);
    polar_correlationsinterp_p(1,i,:) = interp1(ERAdata.latitude,dataploninterp(i,:),waclat);
    differences_interp(1,i,:) = interp1(ERAdata.latitude,diffloninterp(i,:),waclat);
    differences_interp_p(1,i,:) = interp1(ERAdata.latitude,diffploninterp(i,:),waclat);
end

%differences_interp_p (differences_interp_p ~= 0) = -;

%% taking 2d correlations
corr2corr = corr2(squeeze(polar_correlationsinterp),squeeze(modelcorrelations.corrmean));
corr2diff = corr2(squeeze(differences_interp),squeeze(modeldifferences.ensmean));

%% plotting

polartoplot = cat(1,polar_correlationsinterp,modelcorrelations.corrmean);
polar_p_toplot = cat(1,polar_correlationsinterp_p,modelcorrelations.sig.composite_all_pct);
polar_p_toplot (polar_p_toplot <=.05) = 0;
polar_p_toplot (polar_p_toplot >.05) = -1;

modeldifferences.ttest.ensmean (modeldifferences.ttest.ensmean == 0) = -1;
modeldifferences.ttest.ensmean (modeldifferences.ttest.ensmean == 1) = 0;

%difftoplot = cat(1,differences_interp,nanmean(modeldifferences.individual,1));
%difftoplot = cat(1,differences_interp,nanmean(modeldifferences.indmonths.individual(:,:,:,2),1));
difftoplot = cat(1,differences_interp,modeldifferences.composite);
%difftoplot = cat(1,differences_interp,modeldifferences.ensmean);
diff_p_toplot = cat(1,differences_interp_p,modeldifferences.ttest.ensmean);

%% plotting polar    

if inputs.lats(1) < 0
    ylims = [-90 0];
    diffclims = [-3 3];
else
    ylims = [0 90];
    diffclims = [-5 5];
end

%%
plotting = 1;
fsize = 16;
contourtitle = {'Observed correlations','Ensemble composite correlations'};
contourtitle2 = {'Observed temperature difference during ozone extremes','Ensemble composite temperature difference during ozone extremes'};
corrtitles = {'Correlation','Temperature'};
intervals = [-1 1;-1 1; -5 5;-5 5];

[fig,fh] = subplotmaps(polartoplot,waclon,waclat,{'div','RdBu'},1,polar_p_toplot,16,contourtitle,'Longitude','Latitude','Correlation','on',...
    [-.75 .75],22,[waclon(1:24:end)]-180,[waclon(1:24:end)]-180,[ylims(1):15:ylims(2)],[ylims(1):15:ylims(2)],'',1,[0 360],ylims,0,'none',1,'Miller Cylindrical');
    
child = get(fig,'children');
axes1 = child(4);
axes2 = child(2);       

axes(axes1);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','a','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',18,... % 
    'EdgeColor','none','fontweight','bold');    

axes(axes2);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','b','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',18,... % 
    'EdgeColor','none','fontweight','bold');    


filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/Correlation_comparison_',inputs.ClLevel{1},...
    '_',monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'_',ensoext];
 
print(filename,'-depsc');

%export_fig(filename,'-pdf');

%% calculate areas of largest difference
%ERAinterim
if inputs.lats(1) < 0
    areas_lons = [115,155;10,40;285,320;280,300]; %lons (East Russia, West Russia,America,Asia)
    areas_lats = [-40,-15;-35,-15;-40,-15;-55,-40]; %lats
else
    areas_lons = [85,180;30,85;240,290;30,120]; %lons (East Russia, West Russia,America,Asia)
    areas_lats = [55,80;55,80;40,65;20,45]; %lats
end


for i = 1:size(areas_lons,1)
    lats = waclat > areas_lats(i,1) & waclat < areas_lats(i,2);
    lons = waclon > areas_lons(i,1) & waclon < areas_lons(i,2);
    latextract = waclat(lats);
    lonextract = waclon(lons);
    [latmesh,lonmesh] = meshgrid(latextract,lonextract);

    latmesh = latmesh';
    lonmesh = lonmesh';
    for k = 1:size(difftoplot,1)           
        
        diff = permute(difftoplot,[1,3,2]);                   
        mult = squeeze(diff(k,lats,lons));          
        if i == 1 && inputs.lats(1) < 0
            [maxval(i,k),maxind] = max(mult(:));        
        else
            [maxval(i,k),maxind] = max(abs(mult(:)));        
        end
        lattoplot(k,i) = latmesh(maxind);
        lontoplot(k,i) = lonmesh(maxind);
    end
end

for i = 1:numel(areas_lons)
    if areas_lons(i) > 180
        areas_lons(i) = areas_lons(i) - 360;
    end
end
lontoplot2 = lontoplot;
for i = 1:numel(lontoplot)
    if lontoplot(i) > 180
        lontoplot(i) = lontoplot(i) - 360;
    end
end
% WACCM

%%    

cbrewdiff = cbrewer('div','RdBu',21);
cbrewdiff = flipud(cbrewdiff);
cbrewdiff = [cbrewdiff(1:10,:);cbrewdiff(12:end,:)];

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([4,1,10,8],:);
cbrewqual3 = cbrewqual2./2;

lstyle = {':','-','-.','--','-'};    
Mks = {'d','s','o','v'};

[fig,~] = subplotmaps(difftoplot,waclon,waclat,{'div','RdBu'},1,diff_p_toplot,16,contourtitle2,'Longitude','Latitude','Temperature difference (K)','on',...
    diffclims,22,[waclon(1:24:end)]-180,[waclon(1:24:end)]-180,[ylims(1):15:ylims(2)],[ylims(1):15:ylims(2)],'',1,[0 360],ylims,0,'none',1,'Miller Cylindrical');


child = get(fig,'children');
axes1 = child(4);
axes2 = child(2);

axes(axes1);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','a','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',18,... % 
    'EdgeColor','none','fontweight','bold');    

axes(axes2);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','b','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',18,... % 
    'EdgeColor','none','fontweight','bold');  


if inputs.includemarkers    

    axes(axes1);

    obslats = [70,67,57,30];%Greenland,Russia,USA,India
    obslons = [310-360,70,266-360,78];    

    for i = 1:length(obslats)
        %m_plot(obslons(i),obslats(i),'LineStyle','none','Marker',Mks(i),'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
        m_plot(lontoplot(1,i),lattoplot(1,i),'LineStyle','none','Marker',Mks{i},'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
        m_plot([areas_lons(i,1),areas_lons(i,2),areas_lons(i,2),areas_lons(i,1),areas_lons(i,1)],[areas_lats(i,1),areas_lats(i,1),areas_lats(i,2),areas_lats(i,2),areas_lats(i,1)],...
            'LineStyle','-','LineWidth',4,'color',cbrewqual2(i,:),'LineStyle',lstyle{i});
    end

    complats = [67,66,52,36];
    complons = [310-360,75,255-360,75];

    axes(axes2);

    for i = 1:length(obslats)
        %m_plot(complons(i),complats(i),'LineStyle','none','Marker',Mks(i),'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
        m_plot(lontoplot(2,i),lattoplot(2,i),'LineStyle','none','Marker',Mks{i},'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2);
        m_plot([areas_lons(i,1),areas_lons(i,2),areas_lons(i,2),areas_lons(i,1),areas_lons(i,1)],[areas_lats(i,1),areas_lats(i,1),areas_lats(i,2),areas_lats(i,2),areas_lats(i,1)],...
            'LineStyle','-','LineWidth',4,'color',cbrewqual2(i,:),'LineStyle',lstyle{i});
    end        
    
end

filename2 = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/Diff_comparison_',...
    inputs.ClLevel{1},'_',monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))...
    ,num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),ensoext,'composite'];

print(filename2,'-depsc');

%% plot line plots

fsize = 18;
obslats = [70,67,57,30];%Greenland,Russia,USA,India
obslons = [310,70,266,78];
complats = [68,66,52,36];
complons = [310,75,255,75];
for i = 1:length(obslats)
    [~,latind(i)] = min(abs(lattoplot(1,i)-ERAdata.latitude));
    [~,lonind(i)] = min(abs(lontoplot2(1,i)-ERAdata.longitude));
end
for i = 1:length(complats)
    [~,modlatind(i)] = min(abs(lattoplot(2,i)-waclat));
    [~,modlonind(i)] = min(abs(lontoplot2(2,i)-waclon));
end

if inputs.lats(1) < 0
    ylimax = [-2.5 2.5];
    xticklab = {'November','December','January','February','March'};
    lnames = {['Australia, (',num2str(round(abs(lontoplot(1,1)))),'{\circ}','W, ',num2str(round(abs(lattoplot(1,1)))),'{\circ}','N)'],...
        ['Southern Africa, (',num2str(round(abs(lontoplot(1,2)))),'{\circ}','W, ',num2str(round(abs(lattoplot(1,2)))),'{\circ}','N)'],...
        ['South America, (',num2str(round(abs(lontoplot(1,3)))),'{\circ}','W, ',num2str(round(abs(lattoplot(1,3)))),'{\circ}','N)'],...
        ['SA peninsula , (',num2str(round(abs(lontoplot(1,4)))),'{\circ}','W, ',num2str(round(abs(lattoplot(1,4)))),'{\circ}','N)']};
    lnamesmod = {['Australia, (',num2str(round(abs(lontoplot(2,1)))),'{\circ}','W, ',num2str(round(abs(lattoplot(2,1)))),'{\circ}','N)'],...
        ['Southern Africa, (',num2str(round(abs(lontoplot(2,2)))),'{\circ}','W, ',num2str(round(abs(lattoplot(2,2)))),'{\circ}','N)'],...
        ['South America, (',num2str(round(abs(lontoplot(2,3)))),'{\circ}','W, ',num2str(round(abs(lattoplot(2,3)))),'{\circ}','N)'],...
        ['SA peninsula, (',num2str(round(abs(lontoplot(2,4)))),'{\circ}','W, ',num2str(round(abs(lattoplot(2,4)))),'{\circ}','N)']};
else
    ylimax = [-7.5 7.5];
    xticklab = {'March','April','May','June','July'};
    lnames = {['eastern Russia, (',num2str(round(abs(lontoplot(1,1)))),'{\circ}','W, ',num2str(round(abs(lattoplot(1,1)))),'{\circ}','N)'],...
        ['western Russia, (',num2str(round(abs(lontoplot(1,2)))),'{\circ}','W, ',num2str(round(abs(lattoplot(1,2)))),'{\circ}','N)'],...
        ['North America, (',num2str(round(abs(lontoplot(1,3)))),'{\circ}','W, ',num2str(round(abs(lattoplot(1,3)))),'{\circ}','N)'],...
        ['southern Asia, (',num2str(round(abs(lontoplot(1,4)))),'{\circ}','W, ',num2str(round(abs(lattoplot(1,4)))),'{\circ}','N)']};
    lnamesmod = {['eastern Russia, (',num2str(round(abs(lontoplot(2,1)))),'{\circ}','W, ',num2str(round(abs(lattoplot(2,1)))),'{\circ}','N)'],...
        ['western Russia, (',num2str(round(abs(lontoplot(2,2)))),'{\circ}','W, ',num2str(round(abs(lattoplot(2,2)))),'{\circ}','N)'],...
        ['North America, (',num2str(round(abs(lontoplot(2,3)))),'{\circ}','W, ',num2str(round(abs(lattoplot(2,3)))),'{\circ}','N)'],...
        ['southern Asia, (',num2str(round(abs(lontoplot(2,4)))),'{\circ}','W, ',num2str(round(abs(lattoplot(2,4)))),'{\circ}','N)']};
end



% cbrewqual = cbrewer('qual','Set1',10);
% cbrewqual2 = cbrewqual([3,4,10,8],:);

fig = figure;
set(fig,'position',[100 100 1200 450],'color','white');
lstyle = {':','-','-.','--'};        
sp(1) = subplot(1,2,1);
sp_pos(1,:) = get(sp(1),'position');
for i = 1:length(obslats)    
    pho(i) = plot(differences_ind(:,lonind(i),latind(i)),Mks{i},'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:),...
        'MarkerSize',10,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:));
    hold on
end
plot([0 5],[0,0],'--k','lineWidth',2);
xlim([.5,5.5]);
ylim(ylimax);
set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
xlabel('Month','fontsize',fsize+2);
ylabel('Temperature difference (K)','fontsize',fsize+2);
title('Observations','fontsize',fsize+4)
lh = legend(pho,lnames);
set(lh,'box','off','fontsize',fsize-4,'location','NorthEast')

sp(2) = subplot(1,2,2);
totimes = [-3,-1,1,3];
for i = 1:length(complats)      
    %plot(squeeze(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:))','linewidth',1,'LineStyle',lstyle{i},'color',cbrewqual2(i,:)./1.5);
    
    hold on
%     errorbar([1:5]+(totimes(i)./20),squeeze(nanmean(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:),1)),squeeze(std(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:),0,1)),...
%         'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
%     ph(i) = plot([1:5]+(totimes(i)./20),squeeze(nanmean(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:),1)),Mks{i},...
%         'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:),'MarkerSize',10,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:));
    
    errorbar([1:5]+(totimes(i)./20),squeeze(modeldifferences.indcomposite(:,modlonind(i),modlatind(i),:)),squeeze(std(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:),0,1)),... 
        'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
    ph(i) = plot([1:5]+(totimes(i)./20),squeeze(modeldifferences.indcomposite(:,modlonind(i),modlatind(i),:)),Mks{i},...
        'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:),'MarkerSize',10,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:));
    
end
box on
sp_pos(2,:) = get(sp(2),'position');
set(sp(2),'position',[sp_pos(2,1)-.05,sp_pos(2,2),sp_pos(2,3),sp_pos(1,4)]);
set(sp(1),'position',[sp_pos(1,1),sp_pos(1,2),sp_pos(1,3),sp_pos(1,4)]);

sp_pos(1,:) = get(sp(1),'position');
sp_pos(2,:) = get(sp(2),'position');

annotation('textbox',[sp_pos(1,1),sp_pos(1,2),sp_pos(1,3:4)],'String','c','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',18,... % 
    'EdgeColor','none','fontweight','bold');  

annotation('textbox',[sp_pos(2,1),sp_pos(2,2),sp_pos(2,3:4)],'String','d','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',18,... % 
    'EdgeColor','none','fontweight','bold');  

%plot(meandiff,'linewidth',3);
hold on
plot([0 5],[0,0],'--k','lineWidth',2);
xlim([.5,5.5]);
ylim([ylimax]);
set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
xlabel('Month','fontsize',fsize+2);

title('Ensemble composite','fontsize',fsize+4)
lh = legend(ph,lnamesmod);
set(lh,'box','off','fontsize',fsize-4,'location','NorthEast')
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/20th_percentile_lines_',...
    inputs.ClLevel{1},'_',monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'_',ensoext,'composite'];
export_fig(filename,'-pdf');

end
