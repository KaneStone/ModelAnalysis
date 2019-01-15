clear all
% take correlations of singular runs
ClLevel = 'highCl';
tozvar = 'toz';
var = 'TS';

directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
files = dir([directory,'*.nc']);
varfiles = dir([vardirectory,'*.nc']);
lats = [63,90];
% Read in pred run 3d data
count = 1;
toznormanom2 = [];
TSnormanom2 = [];
count2 = 1;
TSnorm2 = zeros(144,96,200);
varmonth = [3,4];
tozmonth = 3;

%%

% for i = 1:length(files)
%     %read in weighted area average data [height,time]
%     [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
%     [~,vardata(i),~] = Read_in_netcdf([vardirectory,varfiles(i).name]);       
%     
%     if i == 1
%         [~,latind(1)] = min(abs(data(i).lat - lats(1)));
%         [~,latind(2)] = min(abs(data(i).lat - lats(2)));
%     end
%     
%     %weighted average    
%     varweighted(i,:) = weightedaverage(squeeze(nanmean(data(i).(tozvar)(:,latind(1):latind(2),:),1)),data(i).lat(latind(1):latind(2)));                            
%     varcombine(i,:,:,:) = vardata(i).TS;
%     
%     % construct year only vector
%     years(i).y = CCMI_years(data(i).date);      
%     
%     dateind = find(years(i).y >= tozdates(1) & years(i).y <= tozdates(2));  
%     
%     toznormanom(i,:) = detrend((varweighted(i,dateind(1)+tozmonth-1:12:dateind(end))'-nanmean(varweighted(i,dateind(1)+tozmonth-1:12:dateind(end))'))./...
%         std(varweighted(i,dateind(1)+tozmonth-1:12:dateind(end))'));
%     toznormanom2 = [toznormanom2,toznormanom(i,:)];
%     
% %     timeperiodhigh = [1995,2016];%[1955,1975]
% 
% %     vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
% %     varfiles = dir([vardirectory,'*.nc']);
% % 
% %     [data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl]...
% %         = predruns_ReadInlayer(vardirectory,varfiles,var,timeperiodhigh,lats,detrend);
%     
%     for k = 1:size(varcombine,2)
%         k
%         %for l = 1:size(varcombine,3)
%             for m = 1:length(varmonth)
%                 TSnormanom(m,i,k,:,:) = detrend((squeeze(varcombine(i,k,:,dateind(1)+varmonth(m)-1:12:dateind(end)))-...                
%                     repmat(nanmean(squeeze(varcombine(i,k,:,dateind(1)+varmonth(m)-1:12:dateind(end))),2),[1,noyears]))./...
%                     repmat(std(squeeze(varcombine(i,k,:,dateind(1)+varmonth(m)-1:12:dateind(end))),1,2),[1,noyears]));                                    
%             end
%             TSnormanom2 = nanmean(TSnormanom,1);
%             for l = 1:size(varcombine,3)
%                 corlatlon(i,k,l) = corr(squeeze(TSnormanom2(1,i,k,l,:)),squeeze(toznormanom(i,:))');
%                 if i == length(files)
%                     if strcmp(ClLevel,'highCl')
%                         tempforcor = [squeeze(TSnormanom2(1,1,k,l,:));squeeze(TSnormanom2(1,2,k,l,:));squeeze(TSnormanom2(1,3,k,l,:));...
%                             squeeze(TSnormanom2(1,4,k,l,:));squeeze(TSnormanom2(1,5,k,l,:));squeeze(TSnormanom2(1,6,k,l,:));...
%                             squeeze(TSnormanom2(1,7,k,l,:));squeeze(TSnormanom2(1,8,k,l,:));squeeze(TSnormanom2(1,9,k,l,:))];
%                     else
%                         tempforcor = [squeeze(TSnormanom2(1,1,k,l,:));squeeze(TSnormanom2(1,2,k,l,:));squeeze(TSnormanom2(1,3,k,l,:));...
%                             squeeze(TSnormanom2(1,4,k,l,:));squeeze(TSnormanom2(1,5,k,l,:));squeeze(TSnormanom2(1,6,k,l,:));...
%                             squeeze(TSnormanom2(1,7,k,l,:));squeeze(TSnormanom2(1,8,k,l,:));squeeze(TSnormanom2(1,9,k,l,:));squeeze(TSnormanom2(1,10,k,l,:))];
%                     end
%                     corlatlon(length(files)+1,k,l) = corr(tempforcor,squeeze(toznormanom2)');
%                 end
%             end
%             
%         %end
%     end    
%     count2 = count2+20;
% end


%% Read in highcl variable

ClLevel = 'highCl';
timeperiodhigh = [1995,2024];%[1955,1975]
timeperiodhighformodel = [1995,2015];%[1955,1975]
noyears = timeperiodhighformodel(2) - timeperiodhighformodel(1)+1;
timeperiodrest = (timeperiodhighformodel(2)+1):timeperiodhigh(2);
noyears2 = timeperiodhigh(2) - timeperiodhighformodel(2);
ifdetrend = 1;
vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[data2.highcl,years2.highcl,composite.highcl,dataMonthArrange.highcl]...
    = predruns_ReadInlayer(vardirectory,varfiles,var,timeperiodhigh,lats,ifdetrend);

monthcombine.(var) = squeeze(nanmean(dataMonthArrange.highcl(:,varmonth,1:noyears,:,:),2));
monthcombinetopredict.(var) = squeeze(nanmean(dataMonthArrange.highcl(:,varmonth,end-(noyears2-1):end,:,:),2));

%% Read in ENSO

[NINO34,IOD,NH] = predruns_calculateENSO(indvarextract,composite,latitude,longitude,indmonth);

%% Read in TOZ highcl and take percentiles
ClLevel = 'highCl';
tozdates = [1995,2024];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
detrend_ozone = 1;
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);

toz_dataMonthArrange.highcltopredict = toz_dataMonthArrange.highcl(:,:,end-(noyears2-1):end);  
toz_dataMonthArrange.highclformodel = toz_dataMonthArrange.highcl(:,:,1:noyears);  

%% latitudes and longitudes
longitude = data2.highcl.lon;
latitude = data2.highcl.lat;

%% taking correlations
for i = 1:size(monthcombine.(var),1)
    for j = 1:size(monthcombine.(var),3)
        [r(i,j,:), p(i,j,:)] = corr(squeeze(monthcombine.(var)(i,:,j,:)),squeeze(toz_dataMonthArrange.highclformodel(i,tozmonth,:)));
    end
end

sig = .05;
ptoplot = p;
ptoplot (ptoplot <= sig) = 0;
ptoplot (ptoplot > sig) = 1;

%% plotting
ind_cp = 0;
if ind_cp
    cbrew = cbrewer('div','RdBu',16);         
    %corlatlon(length(files)+2,:,:) = squeeze(nanmean(corlatlon(1:length(files),:,:),1));    
    %corlatlon(length(files)+3,:,:) = squeeze(nanmean(corlatlon([1,2,3,8,9],:,:),1));    
    if strcmp(ClLevel,'highCl')
        %contourtitle = {'1','2','3','4','5','6','7','8','9','composite','mean','1,2,3,8,9'};
        contourtitle = {'1','2','3','4','5','6','7','8','9'};
    else
        contourtitle = {'1','2','3','4','5','6','7','8','9','10','composite','mean','1,2,3,8,9'};
    end

    rtoplot = permute(r,[1,3,2]);
    ptoplot = permute(ptoplot,[1,3,2]);

    for i = 1:size(r,1);
        subplotmaps(rtoplot(i,:,:),longitude,latitude,{'div','RdBu'},1,ptoplot(i,:,:),12,contourtitle(i),'Longitude','','Correlation','on',...
                [-1,1],22,[-90:10:20],[-90:10:20],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 90],0,'none',1,'Miller Cylindrical');

        latforextract = latitude(49:end-8);

        temp = squeeze(r(i,49:end-8,:));    
        [maxnum(i), maxind(i)] = max(temp(:));
        [minnum(i), minind(i)] = min(temp(:));
        [maxlatind(i), maxlonind(i)] = ind2sub(size(temp),maxind(i));
        [minlatind(i), minlonind(i)] = ind2sub(size(temp),minind(i));
        maxlat(i) = latforextract(maxlatind(i));
        minlat(i) = latforextract(minlatind(i));
        maxlon(i) = longitude(maxlonind(i));
        minlon(i) = longitude(minlonind(i));
        if maxlon(i) > 180
            maxlon(i) = maxlon(i) - 360;
        end
        if minlon(i) > 180
            minlon(i) = minlon(i) - 360;
        end       
        

        m_plot(maxlon(i),maxlat(i),'Marker','p','MarkerSize',30,'MarkerFaceColor',[152,78,163]/255,'MarkerEdgeColor',[152,78,163]/255)
        m_plot(minlon(i),minlat(i),'Marker','p','MarkerSize',30,'MarkerFaceColor',[166,86,40]/255,'MarkerEdgeColor',[166,86,40]/255)
        
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/'...
            ,ClLevel,'_',contourtitle{i},'_individualcorr_novtoz_dects_NH','.png'];
        export_fig(filename,'-png');

    end

    
    
% %% plotting
% cbrew = cbrewer('div','RdBu',16);         
% corlatlon(length(files)+2,:,:) = squeeze(nanmean(corlatlon(1:length(files),:,:),1));    
% corlatlon(length(files)+3,:,:) = squeeze(nanmean(corlatlon([1,2,3,8,9],:,:),1));    
% if strcmp(ClLevel,'highCl')
%     %contourtitle = {'1','2','3','4','5','6','7','8','9','composite','mean','1,2,3,8,9'};
%     contourtitle = {'1','2','3','4','5','6','7','8','9'};
% else
%     contourtitle = {'1','2','3','4','5','6','7','8','9','10','composite','mean','1,2,3,8,9'};
% end
% 
% %rtoplot = permute(r,[1,3,2]);
% 
% for i = 1:size(corlatlon,1);
%     subplotmaps(corlatlon(i,:,:),longitude,latitude,{'div','RdBu'},1,[],12,contourtitle(i),'Longitude','','Correlation','on',...
%             [-1,1],22,[-90:10:20],[-90:10:20],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
%         
%     filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','correlations/maps/',ClLevel,'_',contourtitle{i},'_individualcorr_novtoz_dects_NH2','.png'];
%     export_fig(filename,'-png');
%         
% end
end

%% finding area of largest correlation ouside of 90 
latforextract = latitude(49:end-8);
for i = 1:size(r,1)
    temp = squeeze(r(i,49:end-8,:));    
    [maxnum(i), maxind(i)] = max(temp(:));
    [minnum(i), minind(i)] = min(temp(:));
    [maxlatind(i), maxlonind(i)] = ind2sub(size(temp),maxind(i));
    [minlatind(i), minlonind(i)] = ind2sub(size(temp),minind(i));
    maxlat(i) = latforextract(maxlatind(i));
    minlat(i) = latforextract(minlatind(i));
    maxlon(i) = longitude(maxlonind(i));
    minlon(i) = longitude(minlonind(i));
end

%% take regression

for i = 1:length(files)
    maxlatind(i) = find(latitude == maxlat(i));
    minlatind(i) = find(latitude == minlat(i));
    maxlonind(i) = find(longitude == maxlon(i));
    minlonind(i) = find(longitude == minlon(i));
    predictor(i,:) = squeeze(toz_dataMonthArrange.highclformodel(i,tozmonth,:));    
    topredictrmax(i,:) = squeeze(monthcombine.(var)(i,:,maxlatind(i),maxlonind(i)));
    topredictrmin(i,:) = squeeze(monthcombine.(var)(i,:,minlatind(i),minlonind(i)));
    
    predictorforrest(i,:) = squeeze(toz_dataMonthArrange.highcltopredict(i,tozmonth,:));
    topredictrmaxrest(i,:) = squeeze(monthcombinetopredict.(var)(i,:,maxlatind(i),maxlonind(i)));
    topredictrminrest(i,:) = squeeze(monthcombinetopredict.(var)(i,:,minlatind(i),minlonind(i)));
    
    
    [brmax(i,:),bintrmax(i,:,:),rrmax(i,:),rintrmax(i,:,:),statsrmax(i).s] = regress(topredictrmax(i,:)',[ones(length(predictor),1),...
        predictor(i,:)']);
    [brmin(i,:),bintrmin(i,:,:),rrmin(i,:),rintrmin(i,:,:),statsrmin(i).s] = regress(topredictrmin(i,:)',[ones(length(predictor(i,:)),1),...
        predictor(i,:)']);
    
    regressmodel_rmax(i,:) = brmax(i,1) + (brmax(i,2:end)' * predictor(i,:)');
    regressmodel_rmin(i,:) = brmin(i,1) + (brmin(i,2:end)' * predictor(i,:)');
    
    regressmodelrest_rmax(i,:) = brmax(i,1) + (brmax(i,2:end)' * predictorforrest(i,:))';
    regressmodelrest_rmin(i,:) = brmin(i,1) + (brmin(i,2:end)' * predictorforrest(i,:))';
    
end

%% plotting high correlation areas
cbrew2 = cbrewer('qual','Set1',12);
fsize = 18;
allyears = [timeperiodhighformodel(1):timeperiodhighformodel(2),timeperiodrest];
for i = 1:length(files)    
    createfig('medium','on');
    sp(1) = subplot(2,1,1);
    sp_pos(1,:) = get(sp(1),'position');
    %toplotTS_max = [topredictrmax(i,:),topredictrmaxrest(i,:)];
    %toplot_TCO = [predictor,regressmodelrest_rmax(i,:)];
    [ax,ph1,ph2] = plotyy(timeperiodhighformodel(1):timeperiodhighformodel(2),topredictrmax(i,:),timeperiodhighformodel(1):timeperiodhighformodel(2),predictor(i,:));%,'color',cbrew2(1,:),'LineWidth',3);        
    hold on
    plot(timeperiodrest,topredictrmaxrest(i,:));
    plot(timeperiodrest,regressmodelrest_rmax(i,:));%,'color',cbrew2(1,:),'LineWidth',3);        
    ylimstd1 = std(topredictrmax(i,:))*3;
    ylimmn1 = nanmean(topredictrmax(i,:));
    ylimstd2 = std(predictor(i,:))*3;
    ylimmn2 = nanmean(predictor(i,:));
    set(ax(1),'YLim',[ylimmn1 - ylimstd1, ylimmn1 + ylimstd1])
    set(ax(2),'YLim',[ylimmn2 - ylimstd2, ylimmn2 + ylimstd2])
    set(ph1,'LineWidth',3);
    set(ph2,'LineWidth',3);
    set(ax,'fontsize',fsize);
    title(['maximum correlation: ',num2str(maxlat(i)),'N, ',num2str(maxlon(i)),'E'],'fontsize',fsize+4);
    ylabel(ax(1), 'Surface temperature (K)','fontsize',fsize+2);
    ylabel(ax(2), 'TCO (DU)','fontsize',fsize+2);
    xlabel('Year','fontsize',fsize+2);
    annotation('textbox',sp_pos(1,:),'String',{['r = ', num2str(maxnum(i))]},...
             'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+6,...
            'EdgeColor','none','fontweight','bold')
        
    sp(2) = subplot(2,1,2);
    sp_pos(2,:) = get(sp(2),'position');
    [ax2,ph21,ph22] = plotyy(timeperiodhighformodel(1):timeperiodhighformodel(2),topredictrmin(i,:),timeperiodhighformodel(1):timeperiodhighformodel(2),predictor(i,:));%,'color',cbrew2(1,:),'LineWidth',3);    
    ylimstd21 = std(topredictrmin(i,:))*3;
    ylimmn21 = nanmean(topredictrmin(i,:));
    ylimstd22 = std(predictor(i,:))*3;
    ylimmn22 = nanmean(predictor(i,:));
    set(ax2(1),'YLim',[ylimmn21 - ylimstd21, ylimmn21 + ylimstd21])
    set(ax2(2),'YLim',[ylimmn22 - ylimstd22, ylimmn22 + ylimstd22])
    set(ax2(2),'YDir','Reverse');
    set(ph21,'LineWidth',3);
    set(ph22,'LineWidth',3);
    set(ax2,'fontsize',fsize);
    title(['minimum correlation: ',num2str(minlat(i)),'N, ',num2str(minlon(i)),'E'],'fontsize',fsize+4);
    ylabel(ax2(1), 'Surface temperature (K)','fontsize',fsize+2);
    ylabel(ax2(2), 'TCO (DU)','fontsize',fsize+2);
    xlabel('Year','fontsize',fsize+2);
    % normalise ylimits to 3 standard deviations
    annotation('textbox',sp_pos(2,:),'String',{['r = ', num2str(minnum(i))]},...
             'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+6,...
            'EdgeColor','none','fontweight','bold')    
    
    %ph2 = plotyy(timeperiodhigh(1):timeperiodhigh(2),predictor,'color',cbrew2(2,:),'LineWidth',3);
    %lh = legend([ph,ph2],'Surface temperature','TOZ');        
    
end

%% plotting high correlation areas
cbrew2 = cbrewer('qual','Set1',12);
fsize = 18;
allyears = [timeperiodhighformodel(1):timeperiodhighformodel(2),timeperiodrest];
for i = 1:length(files)    
    createfig('medium','on');
    sp(1) = subplot(2,1,1);
    sp_pos(1,:) = get(sp(1),'position');
    %toplotTS_max = [topredictrmax(i,:),topredictrmaxrest(i,:)];
    %toplot_TCO = [predictor,regressmodelrest_rmax(i,:)];   
    hold on
    plot(timeperiodrest,topredictrmaxrest(i,:));
    hold on
    plot(timeperiodrest,regressmodelrest_rmax(i,:));%,'color',cbrew2(1,:),'LineWidth',3);           
    title(['maximum correlation: ',num2str(maxlat(i)),'N, ',num2str(maxlon(i)),'E'],'fontsize',fsize+4);    
    ylabel('Surface temperature (K)','fontsize',fsize+2);    
    xlabel('Year','fontsize',fsize+2);
    annotation('textbox',sp_pos(1,:),'String',{['r = ', num2str(maxnum(i))]},...
             'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+6,...
            'EdgeColor','none','fontweight','bold')
        
    sp(2) = subplot(2,1,2);
    sp_pos(2,:) = get(sp(2),'position');
    plot(timeperiodrest,topredictrminrest(i,:));
    hold on
    plot(timeperiodrest,regressmodelrest_rmin(i,:));%,'color',cbrew2(1,:),'LineWidth',3);        
    title(['minimum correlation: ',num2str(minlat(i)),'N, ',num2str(minlon(i)),'E'],'fontsize',fsize+4);
    ylabel('Surface temperature (K)','fontsize',fsize+2);    
    xlabel('Year','fontsize',fsize+2);
    % normalise ylimits to 3 standard deviations
    annotation('textbox',sp_pos(2,:),'String',{['r = ', num2str(minnum(i))]},...
             'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+6,...
            'EdgeColor','none','fontweight','bold')    
    
    %ph2 = plotyy(timeperiodhigh(1):timeperiodhigh(2),predictor,'color',cbrew2(2,:),'LineWidth',3);
    %lh = legend([ph,ph2],'Surface temperature','TOZ');        
    
end


