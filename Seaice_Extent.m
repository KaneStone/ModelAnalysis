function [] = Seaice_Extent(differences,Seasplot,...
    modlat,modlon,obslat,obslon,inputs,mons,surfacedata,observations,pct)

% plots sea ice extent for ensemble and observations
colors = cbrewer('qual','Set1',10);
colors2 = cbrewer('qual','Set2',7);
%% Read in area grid cells
areacell = gridcellarea(modlat,modlon);
areacellobs = gridcellarea(obslat,obslon);

% I need to calculate the sea ice extent for all time periods first before
% detrending.

%% calculate sea ice extent for all years before detrending.

modeldata = surfacedata.highCl.dataMonthArrange;
modeldata (modeldata < inputs.fraction) = 0;
modeldata (modeldata >= inputs.fraction) = 1;

obsdata = observations.montharrange;
obsdata (obsdata < inputs.fraction) = 0;
obsdata (obsdata >= inputs.fraction) = 1;

for i = 1:size(Seasplot.lon,1)
    latind = modlat >= Seasplot.lat(i,1) & modlat <= Seasplot.lat(i,2);
    latindobs = obslat >= Seasplot.lat(i,1) & obslat <= Seasplot.lat(i,2);
    if i == 3
        lonind = modlon >= Seasplot.lon(i,1) | modlon < Seasplot.lon(i,2)-360;
        lonindobs = obslon >= Seasplot.lon(i,1) | obslon < Seasplot.lon(i,2)-360;
    else
        lonind = modlon >= Seasplot.lon(i,1) & modlon < Seasplot .lon(i,2);
        lonindobs = obslon >= Seasplot.lon(i,1) & obslon < Seasplot.lon(i,2);
    end
    alldataextent(i).sea = permute(areacell(latind,lonind).*permute(modeldata(:,:,:,latind,lonind),[4,5,1,2,3]),[3,4,5,1,2]);
    seaiceextenttemp(i).sea = nansum(alldataextent(i).sea(:,:,:,:),4);
    
    obsdataextent(i).sea = permute(areacellobs(latindobs,lonindobs).*permute(obsdata(:,:,lonindobs,latindobs),[4,3,1,2]),[3,4,1,2]);
    obsseaiceextenttemp(i).sea = nansum(obsdataextent(i).sea(:,:,:),3);
    
    %% detrend first
    
    for j = 1:12        
        seaiceextent(i).detrend(:,j,:) = (detrend(squeeze(seaiceextenttemp(i).sea(:,j,:))')+nanmean(squeeze(seaiceextenttemp(i).sea(:,j,:))',1))';
        obsseaiceextent(i).detrend(j,:) = (detrend(squeeze(obsseaiceextenttemp(i).sea(j,:))')+nanmean(squeeze(obsseaiceextenttemp(i).sea(j,:))',1))';
    end
    %% remove enso now observations
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
                benso(lag,:,:) = ensopredictors\squeeze(obsseaiceextent(i).detrend(m,:))';

                %                     ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonth(m)-lag+1:12:end,:))',ones(size(alldata,3)-ext,1)];
                %                     benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedataind(k,:,m,j,:));                        
                %find max lag
            end
            [~,bensomax_ind] = max(abs(benso),[],1);                
            %[~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   

            bensomax = benso(bensomax_ind(1),1);
            %bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
            obsseaiceextent(i).sea(m,:) = squeeze(obsseaiceextent(i).detrend(m,:))' - bensomax.*squeeze(ENSOyearind(m2-bensomax_ind(1)+1:12:end,:));   
            
        end
    else
        obsseaiceextent(i).sea =  obsseaiceextent(i).detrend;
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
                    benso(lag,:,:) = ensopredictors\squeeze(seaiceextent(i).detrend(k,m,:));

                    %                     ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonth(m)-lag+1:12:end,:))',ones(size(alldata,3)-ext,1)];
                    %                     benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedataind(k,:,m,j,:));                        
                    %find max lag
                end
                [~,bensomax_ind] = max(abs(benso),[],1);                
                %[~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   

                bensomax = benso(bensomax_ind(1),1);
                %bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
                seaiceextent(i).sea(k,m,:) = squeeze(seaiceextent(i).detrend(k,m,:))' - bensomax.*squeeze(ENSOyearind(m2-bensomax_ind(1)+1:12:end));   

            end
        end
    else
        seaiceextent(i).sea =  seaiceextent(i).detrend;
    end
    
%     for j = 1:12
%         seaiceextent(i).sea(:,j,:) = (detrend(squeeze(seaiceextent(i).sea(:,j,:))')+nanmean(squeeze(seaiceextent(i).sea(:,j,:))',1))';
%         obsseaiceextent(i).sea(j,:) = (detrend(squeeze(obsseaiceextent(i).sea(j,:))')+nanmean(squeeze(obsseaiceextent(i).sea(j,:))',1))';
%     end
end

%% extract ozone percentiles
for i = 1:length(seaiceextent)
    for j = 1:size(seaiceextent(i).detrend,1)    
        seaiceextent(i).low(j,:,:) = seaiceextent(i).sea(j,:,pct.highCl.ind.lowind(j,:));
        seaiceextent(i).high(j,:,:) = seaiceextent(i).sea(j,:,pct.highCl.ind.highind(j,:));        
        seaiceextent(i).difference(j,:) = nanmean(seaiceextent(i).high(j,:,:),3) - nanmean(seaiceextent(i).low(j,:,:),3);
    end
    obsseaiceextent(i).low = obsseaiceextent(i).sea(:,observations.toz.lowind);
    obsseaiceextent(i).high = obsseaiceextent(i).sea(:,observations.toz.highind);
    obsseaiceextent(i).difference = nanmean(obsseaiceextent(i).high,2) - nanmean(obsseaiceextent(i).low,2);
    
    seaiceextent(i).ensdifference = squeeze(nanmean(seaiceextent(i).difference,1));    
    seaiceextent(i).std = squeeze(std(seaiceextent(i).difference,0,1));
end



% %% flag ensemble and individuals
% % ensemble
% 
% below.low = differences.extremesLow.ens;
% below.high = differences.extremesHigh.ens;
% below.low (differences.extremesLow.ens < inputs.fraction) = 0;
% below.low (differences.extremesLow.ens >= inputs.fraction) = 1;
% below.high (differences.extremesHigh.ens < inputs.fraction) = 0;
% below.high (differences.extremesHigh.ens >= inputs.fraction) = 1;
% 
% % individual
% lowmeanind = differences.extremesLow.ind; %[months,lat,lon,fileno.]
% highmeanind = differences.extremesHigh.ind;
% indbelow.low = differences.extremesLow.ind;
% indbelow.high = differences.extremesHigh.ind;
% indbelow.low (differences.extremesLow.ind < inputs.fraction) = 0;
% indbelow.low (differences.extremesLow.ind >= inputs.fraction) = 1;
% indbelow.high (differences.extremesHigh.ind < inputs.fraction) = 0;
% indbelow.high (differences.extremesHigh.ind >= inputs.fraction) = 1;
% 
% %% flag observations
% %% flag ensemble and individuals
% % ensemble
% below.obs.low = differences.observations.lowdata;
% below.obs.high = differences.observations.highdata;
% below.obs.low (differences.observations.lowdata < inputs.fraction) = 0;
% below.obs.low (differences.observations.lowdata >= inputs.fraction) = 1;
% below.obs.high (differences.observations.highdata < inputs.fraction) = 0;
% below.obs.high (differences.observations.highdata >= inputs.fraction) = 1;
% below.obs.low = permute(below.obs.low,[3,2,1]);
% below.obs.high = permute(below.obs.high,[3,2,1]);
% %% difference in sea ice area
% 
% below.area.low = differences.extremesLow.ens;
% below.area.high = differences.extremesHigh.ens;
% below.area.low (differences.extremesLow.ens < inputs.fraction) = 0;
% below.area.high (differences.extremesHigh.ens < inputs.fraction) = 0;
% 
% indbelow.area.low = lowmeanind;
% indbelow.area.high = highmeanind;
% indbelow.area.low (lowmeanind < inputs.fraction) = 0;
% indbelow.area.high (highmeanind < inputs.fraction) = 0;

%% calculating sea ice extent for different regions

% for i = 1:size(Seasplot.lon,1)
%     latind = modlat >= Seasplot.lat(i,1) & modlat <= Seasplot.lat(i,2);
%     latindobs = obslat >= Seasplot.lat(i,1) & obslat <= Seasplot.lat(i,2);
%     if i == 3
%         lonind = modlon >= Seasplot.lon(i,1) | modlon < Seasplot.lon(i,2)-360;
%         lonindobs = obslon >= Seasplot.lon(i,1) | obslon < Seasplot.lon(i,2)-360;
%     else
%         lonind = modlon >= Seasplot.lon(i,1) & modlon < Seasplot.lon(i,2);
%         lonindobs = obslon >= Seasplot.lon(i,1) & obslon < Seasplot.lon(i,2);
%     end
%     % model
%     seaiceextent.ensemble.low(i).a = permute(areacell(latind,lonind).*below.low(latind,lonind,:),[3,1,2]);
%     seaiceextent.ensemble.low(i).a = sum(seaiceextent.ensemble.low(i).a(:,:),2);
%     seaiceextent.ensemble.high(i).a = permute(areacell(latind,lonind).*below.high(latind,lonind,:),[3,1,2]);
%     seaiceextent.ensemble.high(i).a = sum(seaiceextent.ensemble.high(i).a(:,:),2);
%     seaiceextent.ensemble.difference(:,i) = seaiceextent.ensemble.high(i).a - seaiceextent.ensemble.low(i).a;
%     
%     %observations
%     seaiceextent.observations.low(i).a = permute(areacellobs(latindobs,lonindobs).*below.obs.low(latindobs,lonindobs,:),[3,1,2]);
%     seaiceextent.observations.low(i).a = nansum(seaiceextent.observations.low(i).a(:,:),2);
%     seaiceextent.observations.high(i).a = permute(areacellobs(latindobs,lonindobs).*below.obs.high(latindobs,lonindobs,:),[3,1,2]);
%     seaiceextent.observations.high(i).a = nansum(seaiceextent.observations.high(i).a(:,:),2);
%     seaiceextent.observations.difference(:,i) = seaiceextent.observations.high(i).a - seaiceextent.observations.low(i).a;    
%     
%     %individuals
%     seaiceextent.ind.low(i).a = permute(areacell(latind,lonind).*indbelow.low(latind,lonind,:,:),[4,3,1,2]);
%     seaiceextent.ind.low(i).a = sum(seaiceextent.ind.low(i).a(:,:,:),3);
%     seaiceextent.ind.high(i).a = permute(areacell(latind,lonind).*indbelow.high(latind,lonind,:,:),[4,3,1,2]);
%     seaiceextent.ind.high(i).a = sum(seaiceextent.ind.high(i).a(:,:,:),3);
%     seaiceextent.ind.difference(:,:,i) = seaiceextent.ind.high(i).a - seaiceextent.ind.low(i).a;
%     seaiceextent.ensemble2.difference(:,i) = nanmean(seaiceextent.ind.difference(:,:,i));
%     
%     seaicearea.ensemble.low(i).a = permute(areacell(latind,lonind).*below.area.low(latind,lonind,:),[3,1,2]);
%     seaicearea.ensemble.low(i).a = sum(seaicearea.ensemble.low(i).a(:,:),2);
%     seaicearea.ensemble.high(i).a = permute(areacell(latind,lonind).*below.area.high(latind,lonind),[3,1,2]);
%     seaicearea.ensemble.high(i).a = sum(seaicearea.ensemble.high(i).a(:,:),2);
%     seaicearea.ensemble.difference(:,i) = seaicearea.ensemble.high(i).a - seaicearea.ensemble.low(i).a;
%     
%     seaicearea.ind.low(i).a = permute(areacell(latind,lonind).*indbelow.area.low(latind,lonind,:,:),[4,3,1,2]);
%     seaicearea.ind.low(i).a = sum(seaicearea.ind.low(i).a(:,:,:),3);
%     seaicearea.ind.high(i).a = permute(areacell(latind,lonind).*indbelow.area.high(latind,lonind,:,:),[4,3,1,2]);
%     seaicearea.ind.high(i).a = sum(seaicearea.ind.high(i).a(:,:,:),3);
%     seaicearea.ind.difference(:,:,i) = seaicearea.ind.high(i).a - seaicearea.ind.low(i).a;
%     
% end
% 
% seaicearea.ind.std = squeeze(std(seaicearea.ind.difference,0,1));
% seaiceextent.ind.std = squeeze(std(seaiceextent.ind.difference,0,1));


%% plot contour lines
%% Zoom in on different seas
% Barents Kara sea; Laptet, East Siberian, and Chukchi seas, Greenland
% Seas; Sea of Okhotsk and Berling Sea
seanames = {'BKS','LEC','GS','OB'};
obslow = permute(differences.observations.lowdata,[3,2,1]);
obshigh = permute(differences.observations.highdata,[3,2,1]);
 filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
if inputs.clines
    
    %x(size(x,1)/2+1:size(x,1)) = x(size(x,1)/2+1:size(x,1)) - 360;
    %x = [x(size(x,1)/2+1:size(x,1));x(1:size(x,1)/2)];
    for i = 1:size(Seasplot.lon,1)
        for j = 1:size(differences.extremesLow.ens,3)
            createfig('medium','on') 
            if i == 3
                datalow = cat(2,differences.extremesLow.ens(:,:,j),differences.extremesLow.ens(:,1:11,j));%circshift(differences.extremesLow.ens(j,:,:),size(differences.extremesLow.ens(j,:,:),3)/2,2);
                datahigh = cat(2,differences.extremesHigh.ens(:,:,j),differences.extremesHigh.ens(:,1:11,j));%circshift(differences.extremesHigh.ens(j,:,:),size(differences.extremesHigh.ens(j,:,:),3)/2,2);                
                x = [modlon;modlon(1:11)+360];
            else
                datalow = differences.extremesLow.ens(:,:,j);%circshift(differences.extremesLow.ens(j,:,:),size(differences.extremesLow.ens(j,:,:),3)/2,2);
                datahigh = differences.extremesHigh.ens(:,:,j);%circshift(highmeanfinal(j,:,:),size(highmeanfinal(j,:,:),3)/2,2);
                x = modlon;
            end
            m_proj('lambert','lon',[Seasplot.lon(i,1) Seasplot.lon(i,2)],'lat',[Seasplot.lat(i,1) Seasplot.lat(i,2)]);
            [~,h] = m_contour(x,modlat,squeeze(datalow(:,:)),[inputs.fraction,inputs.fraction]);
            hold on
            set(h,'LineWidth',3,'color',colors(1,:));
            hl = plot(1,1,'color',colors(1,:),'LineWidth',3);

            [~,h2] = m_contour(x,modlat,squeeze(datahigh(:,:)),[inputs.fraction,inputs.fraction]);
            set(h2,'LineWidth',3,'color',colors(2,:));
            hl2 = plot(1,1,'color',colors(2,:),'LineWidth',3);
            m_coast('patch',colors2(7,:));
            m_grid('box','fancy','tickdir','in');
            th = title(['Ensemble average ',monthnames(mons(j,:),1,'single'),' sea ice extent relative to ',...
                monthnames(inputs.tozmonth,1,'long'),' Arctic ozone extemes'],'fontsize',20);
            set(th,'units','normalized')
            set(th,'position',[.5 1.045 0])
            filename = [filedir,sprintf('%02d',j),'_iceextenddiff_from_',...
                monthnames(inputs.tozmonth,1,'long'),'_Arcticozoneextremes_over_',seanames{i},...
                num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
            lh = legend([hl,hl2],'relative to TCO lower 20th percentile','relative to TCO upper 20th percentile');
            set(lh,'box','off','fontsize',20,'location','southoutside')
            export_fig(filename,'-png');
            close 1
        end
    end
    
    % observations
    for i = 1:size(Seasplot.lon,1)
        for j = 1:size(obslow,3)
            createfig('medium','on') 
            if i == 3
                datalow = cat(2,obslow(:,:,j),obslow(:,1:21,j));
                datahigh = cat(2,obshigh(:,:,j),obshigh(:,1:21,j));
                x = [obslon;obslon(1:21)+360];
            else
                datalow = obslow(:,:,j);
                datahigh = obshigh(:,:,j);
                x = obslon;
            end
            if i == 4
                a = 1;
            end
            m_proj('lambert','lon',[Seasplot.lon(i,1) Seasplot.lon(i,2)],'lat',[Seasplot.lat(i,1) Seasplot.lat(i,2)]);
            [~,h] = m_contour(x,obslat,datalow,[inputs.fraction,inputs.fraction]);
            hold on
            set(h,'LineWidth',3,'color',colors(1,:));
            hl = plot(1,1,'color',colors(1,:),'LineWidth',3);

            [~,h2] = m_contour(x,obslat,datahigh,[inputs.fraction,inputs.fraction]);
            set(h2,'LineWidth',3,'color',colors(2,:));
            hl2 = plot(1,1,'color',colors(2,:),'LineWidth',3);
            m_coast('patch',colors2(7,:));
            m_grid('box','fancy','tickdir','in');
            th = title(['Ensemble average ',monthnames(mons(j,:),1,'single'),' sea ice extent relative to ',...
                monthnames(inputs.tozmonth,1,'long'),' Arctic ozone extemes'],'fontsize',20);
            set(th,'units','normalized')
            set(th,'position',[.5 1.045 0])
            filename = [filedir,sprintf('%02d',j),'OBS_iceextenddiff_from_',...
                monthnames(inputs.tozmonth,1,'long'),'_Arcticozoneextremes_over_',seanames{i},...
                num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
            lh = legend([hl,hl2],'relative to TCO lower 20th percentile','relative to TCO upper 20th percentile');
            set(lh,'box','off','fontsize',20,'location','southoutside')
            export_fig(filename,'-png');
            close 1
        end
    end        
    
end    
%close all;
%% plot
createfig('largeportrait','on')
lwidth = 3;
fsize = 18;
colors = cbrewer('qual','Set1',10);
Seanames = {'Barents-Kara sea','Sea of Okhotsk, Bering, Laptet, East Siberian, Chukchi, and Beaufort Seas','Greenland Sea','Sea of Okhotsk and Berling Sea'};
lims = [[-2.5e5 2.5e5];[-2e5 8.5e5];[-2.5e5 1.5e5];[-2e5 2e5]];
yticks = [.5e5,1e5,.5e5,1e5];
% sea ice extent (sum up grid point area for all grid points above inputs.fraction)
%stit  = {'c','d','a','b'};
stit  = {'c','b','a'};
for i = 1:size(Seasplot.lon,1)
    sp(i) = subplot(3,1,i);
    sppos(i,:) = get(sp(i),'position');
    
    %MOVE PLOTS TO MAKE LOOK NICE
    plot([-1 13],[0 0],'--k');
    hold on
%     eh(i) = errorbar([1:size(mons,1)]',seaiceextent.ensemble2.difference(:,i),seaiceextent.ind.std(:,i),'color',colors(1,:),'LineWidth',lwidth-1);    
%     ph(i) = plot([1:size(mons,1)]',seaiceextent.ensemble2.difference(:,i),'LineWidth',lwidth,'Color',...
%         colors(1,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
%     pho(i) = plot([1:size(mons,1)]',seaiceextent.observations.difference(:,i),'LineWidth',lwidth,'Color',...
%         colors(2,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
    
    eh(i) = errorbar(mons,seaiceextent(i).ensdifference(mons),seaiceextent(i).std(mons),'color',colors(1,:),'LineWidth',lwidth-1);    
    ph(i) = plot(mons,seaiceextent(i).ensdifference(mons),'LineWidth',lwidth,'Color',...
        colors(1,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
    pho(i) = plot(mons,obsseaiceextent(i).difference(mons),'LineWidth',lwidth,'Color',...
        colors(2,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
    
    
    if inputs.seasons
        set(gca,'xtick',1:size(mons,2),'xticklabel',monthnames(mons(i,:),1,'single'),...
            'ytick',-8e5:yticks(i):8e5,'yticklabel',-8:yticks(i)/1e5:8,'fontsize',fsize);
    else
        set(gca,'xtick',mons,'xticklabel',monthnames(mons,0,'single'),...
            'ytick',-8e5:yticks(i):8e5,'yticklabel',-8:yticks(i)/1e5:8,'fontsize',fsize);
    end    
    xlim([2.5 12.5])
    ylim(lims(i,:));
    
    ylabel(['Sea ice extent ( ',char(215),'10^5 km^2)'],'fontsize',fsize)
    
    if i == 3 || i == 4
        xlabel('Month','fontsize',fsize);        
    end
    title(Seanames{i},'fontsize',fsize);
    if i == 1
        lh = legend([ph(i),pho(i)],'Model Ensemble','Observations');
        set(lh,'box','off','location','SouthWest','fontsize',fsize);
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

filename = [filedir,'NEW_Seaice_extenddiff_LINEPLOT_from_',...
    monthnames(inputs.tozmonth,1,'long'),'_',inputs.obstouse,'_Arcticozoneextremes_over_',...
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