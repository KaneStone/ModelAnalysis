function [] = Seaice_SnowAverage(differences,Seasplot,...
    modlat,modlon,obslat,obslon,inputs,mons,surfacedata,observations,pct)

% plots sea ice extent for ensemble and observations
colors = cbrewer('qual','Set1',10);
colors2 = cbrewer('qual','Set2',7);

%% take weighted average of different regions
mon = 4;
% for i = 1:size(Seasplot.obs.lat,1)
%     obslatsind = obslat > Seasplot.lat(i,1) & obslat < Seasplot.lat(i,2);
%     obslonsind = obslon > Seasplot.lon(i,1) & obslon < Seasplot.lon(i,2);
%     modlatsind = modlat > Seasplot.lat(i,1) & modlat < Seasplot.lat(i,2);
%     modlonsind = modlon > Seasplot.lon(i,1) & modlon < Seasplot.lon(i,2);
%         
%     obslatextract = obslat(obslatsind);
%     obslonextract = obslon(obslonsind);
%     modlatextract = modlat(modlatsind);
%     modlonextract = modlon(modlonsind);
%     
%     [obslatmesh,obslonmesh] = meshgrid(obslatextract,obslonextract);
%     [modlatmesh,modlonmesh] = meshgrid(modlatextract,modlonextract);
% 
%     obslatmesh = obslatmesh';
%     obslonmesh = obslonmesh';
%     modlatmesh = modlatmesh';
%     modlonmesh = modlonmesh';
%     
%     for k = 1:size(differences.observations.difference,1)           
%         
%         obsdiff = permute(differences.observations.difference,[1,3,2]);                   
%         obsmult = squeeze(obsdiff(k,obslatsind,obslonsind));          
%         if i == 1 || i == 3
%             [obsmaxval(i,k),obsmaxind(i,k)] = min(obsmult(:));        
%         else
%             [obsmaxval(i,k),obsmaxind(i,k)] = max(obsmult(:));        
%         end
%         obslattoplot(k,i) = obslatmesh(obsmaxind(i,k));
%         obslontoplot(k,i) = obslonmesh(obsmaxind(i,k));
%     end
%     obstemp = squeeze(obsdiff(:,obslatsind,obslonsind));        
%     obsmaxfm(:,i) = obstemp(:,obsmaxind(i,mon));  
%     
%     for k = 1:size(differences.difference.ens,3)           
%         
%         moddiff = permute(differences.difference.ens,[3,1,2]);                   
%         modmult = squeeze(moddiff(k,modlatsind,modlonsind));          
%         if i == 1 || i == 3
%             [modmaxval(i,k),modmaxind(i,k)] = min(modmult(:));        
%         else
%             [modmaxval(i,k),modmaxind(i,k)] = max(modmult(:));        
%         end
%         modlattoplot(k,i) = modlatmesh(modmaxind(i,k));
%         modlontoplot(k,i) = modlonmesh(modmaxind(i,k));
%         for j = 1:size(differences.difference.ind,1)
%             moddiffind = permute(differences.difference.ind,[4,3,1,2]);                   
%             modmultind = squeeze(moddiffind(j,k,modlatsind,modlonsind));          
%             if i == 1 || i == 3
%                 [modmaxval(i,k),modmaxind(i,k)] = min(modmult(:));        
%             else
%                 [modmaxval(i,k),modmaxind(i,k)] = max(modmult(:));        
%             end
%             modlattoplot(k,i) = modlatmesh(modmaxind(i,k));
%             modlontoplot(k,i) = modlonmesh(modmaxind(i,k));
%         end
%     end
%     temp = squeeze(moddiff(:,modlatsind,modlonsind));
%     modmaxfm(:,i) = temp(:,modmaxind(i,mon));  
% end

%% remove ENSO

%% load in detrended and ENSO removed differences
load('/Volumes/ExternalOne/work/data/predruns/output/data/detrendedENSOremoved_differences.mat');
contourlines = 0;
modeldata = surfacedata.highCl.dataMonthArrange;
obsdata = observations.montharrange;

for i = 1:size(Seasplot.lon,1)
    latind = modlat >= Seasplot.lat(i,1) & modlat <= Seasplot.lat(i,2);
    latindobs = obslat >= Seasplot.lat(i,1) & obslat <= Seasplot.lat(i,2);
    
    lonind = modlon >= Seasplot.lon(i,1) & modlon < Seasplot .lon(i,2);
    lonindobs = obslon >= Seasplot.lon(i,1) & obslon < Seasplot.lon(i,2);
    
    alldataextent(i).sea = modeldata(:,:,:,latind,lonind);
    modelens = repmat(permute(modelensdifferences(latind,lonind,1),[3,4,5,1,2]),[9,12,30,1,1]);   
    
    
    obsdataextent(i).sea = obsdata(:,:,lonindobs,latindobs);
    obs = repmat(permute(observeddifferences(1,lonindobs,latindobs),[4,1,2,3]),[12,37,1,1,]);
        
    if contourlines
        if i < 3
            alldataextent(i).sea (modelens > 0) = NaN;
            obsdataextent(i).sea (obs > 0) = NaN;
        else
            alldataextent(i).sea (modelens < 0) = NaN;
            obsdataextent(i).sea (obs < 0) = NaN;
        end
    end
    
    snowextenttemp(i).sea = nanmean(alldataextent(i).sea(:,:,:,:),4);
    obssnowextenttemp(i).sea = nanmean(obsdataextent(i).sea(:,:,:),3);
    
    laglength = 3;
    if inputs.removeENSO
        for m = 1:12
            for lag = 1:laglength                                                            
                ENSOyearind = observations.ENSO;             
                if m == 1 || m == 2
                    m2 = 3;
                else
                    m2 = m;
                end
                ensopredictors = [squeeze(ENSOyearind(m2-lag+1:12:end,:)),ones(length(squeeze(ENSOyearind(m2-lag+1:12:end,:))),1)];
                benso(lag,:,:) = ensopredictors\squeeze(obssnowextenttemp(i).sea(m,:))';

                %                     ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonth(m)-lag+1:12:end,:))',ones(size(alldata,3)-ext,1)];
                %                     benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedataind(k,:,m,j,:));                        
                %find max lag
            end
            [~,bensomax_ind] = max(abs(benso),[],1);                
            %[~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   

            bensomax = benso(bensomax_ind(1),1);
            %bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
            obssnowextent(i).sea(m,:) = squeeze(obssnowextenttemp(i).sea(m,:))' - bensomax.*squeeze(ENSOyearind(m2-bensomax_ind(1)+1:12:end,:));   
            
        end
    else
        obssnowextent(i).sea =  obssnowextenttemp(i).sea;
    end
    
    % remove model
    if inputs.removeENSO
        for k = 1:size(surfacedata.highCl.ENSO,1)
            for m = 1:12
                for lag = 1:laglength                                                            
                    ENSOyearind = surfacedata.highCl.ENSO(k,:);             
                    if m == 1 || m == 2
                        m2 = 3;
                    else
                        m2 = m;
                    end
                    ensopredictors = [squeeze(ENSOyearind(m2-lag+1:12:end))',ones(length(squeeze(ENSOyearind(m2-lag+1:12:end))),1)];
                    benso(lag,:,:) = ensopredictors\squeeze(snowextenttemp(i).sea(k,m,:));

                    %                     ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonth(m)-lag+1:12:end,:))',ones(size(alldata,3)-ext,1)];
                    %                     benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedataind(k,:,m,j,:));                        
                    %find max lag
                end
                [~,bensomax_ind] = max(abs(benso),[],1);                
                %[~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   

                bensomax = benso(bensomax_ind(1),1);
                %bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
                snowextent(i).sea(k,m,:) = squeeze(snowextenttemp(i).sea(k,m,:))' - bensomax.*squeeze(ENSOyearind(m2-bensomax_ind(1)+1:12:end));   

            end
        end
    else
        snowextent(i).sea =  snowextenttemp(i).sea;
    end
    
    for j = 1:12
        snowextent(i).detrend(:,j,:) = (detrend(squeeze(snowextent(i).sea(:,j,:))')+nanmean(squeeze(snowextent(i).sea(:,j,:))',1))';
        obssnowextent(i).detrend(j,:) = (detrend(squeeze(obssnowextent(i).sea(j,:))')+nanmean(squeeze(obssnowextent(i).sea(j,:))',1))';
    end
end

%% find differences
for i = 1:length(snowextent)
    for j = 1:size(snowextent(i).detrend,1)    
        snowextent(i).low(j,:,:) = snowextent(i).detrend(j,:,pct.highCl.ind.lowind(j,:));
        snowextent(i).high(j,:,:) = snowextent(i).detrend(j,:,pct.highCl.ind.highind(j,:));        
        snowextent(i).difference(j,:) = nanmean(snowextent(i).high(j,:,:),3) - nanmean(snowextent(i).low(j,:,:),3);
    end
    obssnowextent(i).low = obssnowextent(i).detrend(:,observations.toz.lowind);
    obssnowextent(i).high = obssnowextent(i).detrend(:,observations.toz.highind);
    obssnowextent(i).difference = nanmean(obssnowextent(i).high,2) - nanmean(obssnowextent(i).low,2);
    
    snowextent(i).ensdifference = squeeze(nanmean(snowextent(i).difference,1));    
    snowextent(i).std = squeeze(std(snowextent(i).difference,0,1));
end

%% testing
mont = 5;
loc = 3;
figure; plot([obssnowextent(loc).low(mont,:),obssnowextent(loc).high(mont,:)])
hold on
yyaxis right
plot([observations.toz.zm(observations.toz.lowind),observations.toz.zm(observations.toz.highind)]);

figure; plot(obssnowextent(loc).detrend(mont,1:31));
hold on
yyaxis right
plot(observations.toz.zm(1:31))

corr([obssnowextent(loc).low(mont,:),obssnowextent(loc).high(mont,:)]',[observations.toz.zm(observations.toz.lowind),observations.toz.zm(observations.toz.highind)]')
corr(obssnowextent(loc).detrend(mont,1:31)',observations.toz.zm(1:31)')
%% find specific locations % this isn't used
if size(Seasplot.obs.lat,2) == 1
    for i = 1:length(Seasplot.obs.lat)
        [~,obslatind(i)] = min(abs(obslat - Seasplot.obs.lat(i)));
        [~,modlatind(i)] = min(abs(modlat - Seasplot.mod.lat(i)));
        [~,obslonind(i)] = min(abs(obslon - Seasplot.obs.lon(i)));
        [~,modlonind(i)] = min(abs(modlon - Seasplot.mod.lon(i)));
        obsplotextract(:,i) = differences.observations.difference(:,obslonind(i),obslatind(i));
        modplotextract(:,i) = squeeze(differences.difference.ens(modlatind(i),modlonind(i),:));
        modindplotextract(:,:,i) = squeeze(differences.difference.ind(modlatind(i),modlonind(i),:,:));
    end
    modindstd = squeeze(std(modindplotextract,0,2));
else
    for i = 1:size(Seasplot.obs.lat,1)
        
        obsdiff = differences.observations.difference;
        obsdiffseason = repmat(nanmean(differences.observations.difference([2,3,4],:,:),1),[length(inputs.varmonth),1,1]);
        moddiff = differences.difference.ens;
        moddiffseason = repmat(nanmean(differences.difference.ens(:,:,[2,3,4]),3),[1,1,length(inputs.varmonth)]);
        moddiffind = differences.difference.ind;
        moddiffindseason = repmat(moddiffseason,[1,1,1,9]);
        %moddiffindseason = repmat(nanmean(differences.difference.ind(:,:,[4,5,6],:),3),[1,1,12],1);
        if i < 3
%             obsdiff (obsdiffseason > -.005) = NaN;
%             moddiff (moddiffseason > -.005) = NaN;            
%             moddiffind (moddiffseason > -.005) = NaN;
            obsdiff (obsdiffseason > 0) = NaN;
            moddiff (moddiffseason > 0) = NaN;            
            moddiffind (moddiffseason > 0) = NaN;
        else
             obsdiff (obsdiffseason < 0) = NaN;
             moddiff (moddiffseason < 0) = NaN;
             moddiffind (moddiffseason < 0) = NaN;
        end
        
        obslatsind = obslat > Seasplot.obs.lat(i,1) & obslat < Seasplot.obs.lat(i,2);
        obslonsind = obslon > Seasplot.obs.lon(i,1) & obslon < Seasplot.obs.lon(i,2);
        modlatsind = modlat > Seasplot.mod.lat(i,1) & modlat < Seasplot.mod.lat(i,2);
        modlonsind = modlon > Seasplot.mod.lon(i,1) & modlon < Seasplot.mod.lon(i,2);
        
        obsdiffextract(i).a = obsdiff(:,obslonsind,obslatsind);
        obsplotextract(:,i) = nanmean(squeeze(obsdiffextract(i).a(:,:)),2);
        
        moddiffextract(i).a = permute(moddiff(modlatsind,modlonsind,:),[3,2,1]);
        modplotextract(:,i) = nanmean(squeeze(moddiffextract(i).a(:,:)),2);
        
        moddiffindextract(i).a = permute(moddiffind(modlatsind,modlonsind,:,:),[4,3,2,1]);
        moddiffindmean(:,i,:) = nanmean(squeeze(moddiffindextract(i).a(:,:,:)),3);               
        
    end
    modindstd = squeeze(std(moddiffindmean,0,1))';
end


%%
obsplotextract = [obssnowextent(:).difference];
obsplotextract = obsplotextract(3:end,:);
modplotextract = cat(1,snowextent(:).ensdifference)';
modplotextract = modplotextract(3:end,:);
modindstd = cat(1,snowextent(:).std)';
modindstd = modindstd(3:end,:);
%% plot
filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
createfig('largeportrait','on')
lwidth = 3;
fsize = 18;
colors = cbrewer('qual','Set1',10);
Seanames = {'Northern Russia','Eastern Russia','South East Alaska'};
lims = [[-.075 .05];[-.075 .05];[-.05 .1]];
yticks = [.025,.025,.025];
% sea ice extent (sum up grid point area for all grid points above inputs.fraction)
if strcmp(inputs.var,'ICEFRAC')
    stit  = {'c','d','a','b'};
else
    stit  = {'c','b','a'};
end

for i = 1:length(Seasplot.obs.lon)
    sp(i) = subplot(3,1,i);
    sppos(i,:) = get(sp(i),'position');
    
    %MOVE PLOTS TO MAKE LOOK NICE
    plot([-1 13],[0 0],'--k');
    hold on
    %eh(i) = errorbar([1:size(mons,1)]',seaiceextent.ensemble.difference(:,i),seaiceextent.ind.std(:,i),'color',colors(1,:),'LineWidth',lwidth-1);    
%     ph(i) = plot([1:size(mons,1)]',modmaxfm(:,i),'LineWidth',lwidth,'Color',...
%         colors(1,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
%     pho(i) = plot([1:size(mons,1)]',obsmaxfm(:,i),'LineWidth',lwidth,'Color',...
%         colors(2,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
    
    eh(i) = errorbar(mons,modplotextract(:,i),modindstd(:,i),'color',colors(1,:),'LineWidth',lwidth-1);    
    ph(i) = plot(mons,modplotextract(:,i),'LineWidth',lwidth,'Color',...
        colors(1,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
    pho(i) = plot(mons,obsplotextract(:,i),'LineWidth',lwidth,'Color',...
        colors(2,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
    
    
    set(gca,'xtick',mons,'xticklabel',monthnames(mons,0,'single'),...
        'ytick',-.5:yticks(i):.5,'yticklabel',-.5:yticks(i):.5,'fontsize',fsize);
    
    xlim([2.5 12.5])
    ylim(lims(i,:));
    
    ylabel(['Snow depth difference (m)'],'fontsize',fsize)
    
    if i == 3
        xlabel('Month','fontsize',fsize);        
    end
    title(Seanames{i},'fontsize',fsize);
    if i == 1
        lh = legend([ph(i),pho(i)],'Model Ensemble','Observations');
        set(lh,'box','off','location','NorthEast','fontsize',fsize);
    end
    
    
    
end

axpos = tightplots(gcf,10,5,.05,[]);
for i = 1:size(axpos,1)
    annotation('textbox',[axpos(i,1),axpos(i,2),axpos(i,3:4)],'String',stit{i},...
        'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',18,... % 
    'EdgeColor','none','fontweight','bold');    
end

if inputs.removeENSO
    ensoext = 're';
else
    ensoext = [];
end

filename = [filedir,inputs.var,'_LINEPLOT_from_',...
    monthnames(inputs.tozmonth,1,'long'),'_Arcticozoneextremes_over_',...
    num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',ensoext];
export_fig(filename,'-pdf');

% % sea ice area (sum up fraction of sea ice for all grid point above inputs.fraction)
% createfig('medium','on')
% lims = [[-1.25e5 1.75e5];[-.5e5 4e5];[-1.5e5 1e5];[-1e5 5e5]];
% yticks = [.5e5,.5e5,.5e5,1e5];
% for i = 1:size(Seasplot.lon,1)
%     subplot(2,2,i)
%     plot([-1 13],[0 0],'--k');
%     hold on
%     eh(i) = errorbar(mons,seaicearea.ensemble.difference(:,i),seaicearea.ind.std(:,i),'color',colors(1,:),'LineWidth',lwidth-1);    
%     ph(i) = plot(mons,seaicearea.ensemble.difference(:,i),'LineWidth',lwidth,'Color',...
%         colors(1,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
%     set(gca,'xtick',mons,'xticklabel',mons,'ytick',-8e5:yticks(i):8e5,'yticklabel',-8:yticks(i)/1e5:8,'fontsize',fsize);
%     xlim([.5 12.5])
%     ylim(lims(i,:));
%     if i == 1 || i == 3
%         ylabel(['Sea ice area( ',char(215),'10^5 km^2)'],'fontsize',fsize+2)
%     end
%     if i == 3 || i == 4
%         xlabel('Month','fontsize',fsize+2);        
%     end
%     title(Seanames{i},'fontsize',fsize);
%     
% end

% % both on same graph
% createfig('medium','on')
% lims = [[-2e5 2e5];[-2e5 4e5];[-2.5e5 1e5];[-2e5 8.5e5]];
% yticks = [.5e5,1e5,.5e5,1e5];
% for i = 1:size(Seasplot.lon,1)
%     subplot(2,2,i)
%     plot([-1 13],[0 0],'--k');
%     hold on
%     eh(i) = errorbar(mons,seaiceextent.ensemble.difference(:,i),seaiceextent.ind.std(:,i),'color',colors(1,:),'LineWidth',lwidth-1);    
%     ph(i) = plot(mons,seaiceextent.ensemble.difference(:,i),'LineWidth',lwidth,'Color',...
%         colors(1,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
%     eh(i) = errorbar(mons,seaicearea.ensemble.difference(:,i),seaicearea.ind.std(:,i),'color',colors(2,:),'LineWidth',lwidth-1);    
%     ph(i) = plot(mons,seaicearea.ensemble.difference(:,i),'LineWidth',lwidth,'Color',...
%         colors(2,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
%     
%     
%     set(gca,'xtick',mons,'xticklabel',mons,'ytick',-8e5:yticks(i):8e5,'yticklabel',-8:yticks(i)/1e5:8,'fontsize',fsize);
%     xlim([.5 12.5])
%     ylim(lims(i,:));
%     if i == 1 || i == 3
%         ylabel(['Sea ice area( ',char(215),'10^5 km^2)'],'fontsize',fsize+2)
%     end
%     if i == 3 || i == 4
%         xlabel('Month','fontsize',fsize+2);        
%     end
%     title(Seanames{i},'fontsize',fsize);
    
% end


end