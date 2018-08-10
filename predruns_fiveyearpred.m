function [] = predruns_fiveyearpred(dataVarMonthAve,dataVarMonth,tozdata,inputs,latitude,longitude,pct,differences,correlations)

%% regressio model using data from 1995-2018
modelyears = [1995:2024];
yearsall = [1995,2009];%;2005,2024];
%yearsall = [2010,2024];%;2005,2024];
backwards = 0;
nopredyears = 3;
ifdetrend = 0;
if ifdetrend
    detrendext = 'detrend';
else
    detrendext = 'nodentred';
end

if backwards
    no = yearsall(1)-modelyears(1)-nopredyears;
else
    no = modelyears(end)-yearsall(2)-nopredyears;
end

for l = 1:no
clearvars bdetrend TSdetrended TSleftdetrended bozone corr2use modelfit modelprediction toplot toplotmodel
%years = [1995,2010];
if ~backwards
    years = yearsall+l-1;
else
    years = yearsall-l+1;
end
yearind = modelyears >= years(1) & modelyears <= years(2);
yearindleft = modelyears > years(2) & modelyears <= years(2)+nopredyears;

modelyears2 = years(1):years(end)+nopredyears;
modelyears2 (modelyears2 > 2024) = [];

% if years(1) > modelyears(1)
%     modelyears = years(1):2024;
% else
%     modelyears = [1995:2024];
% end
ifmed = 0;
%% extracing tozdata

tozdataextract = detrend(squeeze(tozdata(:,inputs.tozmonth,:)))';

%% calculate 1995-2018 detrend coefficients

predictors = [(1:sum(yearind))',ones(size(1:sum(yearind)))'];
predictorsleft = sum(yearind)+1:sum(yearind)+sum(yearindleft);

for i = 1:size(dataVarMonthAve,1)
    ozonepredictors = [tozdataextract(yearind,i),ones(size(tozdataextract(yearind,i)))];
    ozoneleft = tozdataextract(yearindleft,i);
    for j = 1:size(dataVarMonthAve,3)
        % detrending        
        if ifdetrend
            bdetrend(i,j,:,:) = predictors\squeeze(dataVarMonthAve(i,yearind,j,:));
            TSdetrended(i,:,j,:) =  squeeze(dataVarMonthAve(i,yearind,j,:)) - squeeze(bdetrend(i,j,1,:))'.*predictors(:,1) - squeeze(bdetrend(i,j,2,:))';
            TSleftdetrended(i,:,j,:) = squeeze(dataVarMonthAve(i,yearindleft,j,:)) - squeeze(bdetrend(i,j,1,:))'.*predictorsleft' - squeeze(bdetrend(i,j,2,:))';
        else
            TSdetrended(i,:,j,:) =  squeeze(dataVarMonthAve(i,yearind,j,:))-squeeze(nanmean(dataVarMonthAve(i,yearind,j,:),2))';
            TSleftdetrended(i,:,j,:) = squeeze(dataVarMonthAve(i,yearindleft,j,:))-squeeze(nanmean(dataVarMonthAve(i,yearind,j,:),2))';
        end
        med(i,:,j) = squeeze(median(TSdetrended(i,:,j,:)));
        if ifmed
            TSdetrended(i,:,j,:) = squeeze(TSdetrended(i,:,j,:)) - med(i,:,j);
            TSleftdetrended(i,:,j,:) = squeeze(TSleftdetrended(i,:,j,:)) - med(i,:,j);
        end
        
        % creating ozone model
        bozone(i,j,:,:) = ozonepredictors\squeeze(TSdetrended(i,:,j,:));
        corr2use(i,j,:) = corr(ozonepredictors(:,1),squeeze(TSdetrended(i,:,j,:)));
        modelfit(i,:,j,:) = squeeze(bozone(i,j,1,:))'.*tozdataextract(yearind,i) + squeeze(bozone(i,j,2,:))';
        modelprediction(i,:,j,:) = squeeze(bozone(i,j,1,:))'.*ozoneleft + squeeze(bozone(i,j,2,:))';
        if ifmed
            modelprediction(i,:,j,:) = squeeze(modelprediction(i,:,j,:)) - med(i,:,j);
        end
        
    end
end

% %% testing
% run = 5;
% lat = 87;
% lon = 22;
% figure;
% 
% hold on
% plot(modelyears,[squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))],'-o','LineWidth',2);
% %plot(predictorsleft,squeeze(TSleftdetrended(run,:,lat,lon)),'LineWidth',2);
% %plot(squeeze(dataVarMonthAve(1,:,87,22))-nanmean(squeeze(dataVarMonthAve(1,:,87,22))))
% plot(modelyears,[squeeze(modelfit(run,:,lat,lon)),squeeze(modelprediction(run,:,lat,lon))],'-o','LineWidth',2);
% plot([years(2),years(2)],[-10,10],'Color','k','LineStyle','--');
% minlim = min([squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))])-...
%     min([squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))])./10;
% maxlim = max([squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))])+...
%     max([squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))])./10;
% ylim([minlim maxlim]);
% lh = legend('Surface temperature anomalies','Model fit and predictions','Fit/Prediction');
% set(lh,'box','off','location','NorthWest');

%plot(predictorsleft,squeeze(modelprediction(run,:,lat,lon)),'LineWidth',2);
%plot(predictorsleft,(tozdataextract(yearindleft,1)-nanmean(tozdataextract(yearindleft,1)))./10);

%% areas of highest correlations

surfacecombine = cat(2,TSdetrended,TSleftdetrended);
modelcombine = cat(2,modelfit,modelprediction);
corr2use2 = permute(correlations.individual,[1,3,2]);

RMSE(l,:,:,:) = sqrt(squeeze(nanmean((TSleftdetrended - modelprediction).^2,2)));

% areas_lons = [290,320;30,120;240,290;40,120];%lons (Greenland,Russia,America,Asia)
% areas_lats = [60,80;50,75;40,60;30,45];%lats

areas_lons = [290,320;30,60;240,290;65,90];%lons (Greenland,Russia,America,Asia)
areas_lats = [60,80;60,75;40,60;30,45];%lats

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([3,4,10,8],:);
mon = 4;
for i = 1:size(areas_lons,1)
    lats = latitude > areas_lats(i,1) & latitude < areas_lats(i,2);
    lons = longitude > areas_lons(i,1) & longitude < areas_lons(i,2);
    latextract = latitude(lats);
    lonextract = longitude(lons);
    [latmesh,lonmesh] = meshgrid(latextract,lonextract);
    lonmesh = lonmesh';
    latmesh = latmesh';
    for k = 1:size(corr2use,1)  
        mult = squeeze(corr2use(k,lats,lons));
         if l == 1           
            mult2 = squeeze(corr2use2(k,lats,lons));
            if i == 2
                [maxval(i,k),maxind(i,k)] = min(mult(:));     
            else
                [maxval(i,k),maxind(i,k)] = max(mult(:));     
            end            
            [maxval2(i,k),maxind2] = max(abs(mult2(:)));        
            lattoplot(l,k,i) = latmesh(maxind(i,k));
            lontoplot(l,k,i) = lonmesh(maxind(i,k)); 
         end
        rval(l,i,k) = mult(maxind(i,k));
        RMSEatloctemp = squeeze(RMSE(l,k,lats,lons));
        RMSEatloc(l,i,k) = RMSEatloctemp(maxind(i,k));
        toplottemp = squeeze(surfacecombine(k,:,lats,lons));        
        toplot(:,k,i) = toplottemp(:,maxind(i,k));
        bvaltemp = squeeze(bozone(k,lats,1,lons));
        bval(i,k) = bvaltemp(maxind(i,k));
        toplottempmodel = squeeze(modelcombine(k,:,lats,lons));        
        toplotmodel(:,k,i) = toplottempmodel(:,maxind(i,k));

    end
end

%% plotting
plotall = 1;
fsize = 18;
if plotall
    for run = 1:9
    fsize = 18;
    titles = {'Greenland','Russia','North America','Asia'};
    cbrew = cbrewer('qual','Set1',10);

    createfig('large','off');
    for i = 1:4
        sp(i) = subplot(3,2,i);
        yyaxis left
        if modelyears2(1) == 2003
            abc = 1;
        end
        plot(modelyears2,toplot(:,run,i),'-o','LineWidth',3,'color',cbrew(1,:))
        hold on
        plot(modelyears2,toplotmodel(:,run,i),'-o','LineWidth',3,'color',cbrew(2,:),'LineStyle','-')    
        set(gca,'Ycolor','k','fontsize',fsize-2);
        ylabel('Temperture anomaly','fontsize',fsize);
        xlabel('Year','fontsize',fsize);
        plot([years(2),years(2)],[-20,10],'Color','k','LineStyle','--');    
        plot([modelyears2(1),modelyears2(end)],[0,0],'Color','k','LineStyle','--');    
        forlim = toplot(:,run,i);    
        ylim([min(forlim(:)) + min(forlim(:)), max(forlim(:)) + max(forlim(:))./10])
        %ylim([min(forlim(:)) - min(forlim(:))./10, max(forlim(:)) + max(forlim(:))./10])

        % bar
        yyaxis right
        toplot2 = sign(toplot(:,run,i));
        toplot2 (sign(toplot(:,run,i)) == sign(toplotmodel(:,run,i))) = 1;
        toplot2 (sign(toplot(:,run,i)) ~= sign(toplotmodel(:,run,i))) = -1;

        bar(modelyears2,toplot2,'LineWidth',3,'Edgecolor','k','FaceColor',cbrew(8,:))
        ylim([-2 8])
        set(gca,'ytick',-1:1:1);
        set(gca,'Ycolor','k');
        title(titles(i),'fontsize',fsize+2);
        ylabel('Anomaly prediction','fontsize',fsize);
        %plot([modelyears(1),modelyears(end)],[-8,-8],'Color','k','LineStyle','--');    


        annotation('textbox',get(sp(i),'position'),'String',['r = ', sprintf('%.2f',rval(l,i,run))],'FitBoxToText','on','LineStyle','none','Fontsize',fsize);
    end

    annotation('textbox',[.01 1 1 0],'String',['Ensemble member ',num2str(run)] ,'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
            'EdgeColor','none','fontweight','bold');    

    cbrew2 = flipud(cbrewer('div','RdBu',21));
    cbrew2 = [cbrew2(1:10,:);cbrew2(12:end,:)];

    toplotcorr = squeeze(corr2use(run,:,:)); 
    toplotcorr = circshift(toplotcorr,size(toplotcorr,2)/2,2);

    sp2 = subplot(3,1,3);  
    m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
    [~,h] = m_contourf(longitude-180,latitude,toplotcorr,-1:.1:1,'LineStyle','none');
    m_coast('color','k','LineWidth',1);
    m_grid('ytick',0:30:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',fsize);            
    set(sp2,'position',[.175 .07 .7 .3]);
    caxis([-1 1]);        
    ch = colorbar;
    set(ch,'YTick',[-1:.2:1],'fontsize',fsize)%,
    colormap(cbrew2);    
    hold on
    for i = 1:numel(lontoplot)
        if lontoplot(i) > 180
            lontoplot(i) = lontoplot(i) - 360;
        end
    end
    for i = 1:4
        m_plot(lontoplot(1,run,i),lattoplot(1,run,i),'LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',4)
    end

    if ifmed
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/No._',num2str(run),'_',num2str(years(1)),'_med_',detrendext,'_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1)];
    else
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/No._',num2str(run),'_',num2str(years(1)),'_',detrendext,'_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1)];
    end
    export_fig(filename,'-png');            

    end
end
end

%%
cbrew = cbrewer('qual','Set1',10);
if ~plotall
    for i = 1:numel(lontoplot)
        if lontoplot(i) > 180
            lontoplot(i) = lontoplot(i) - 360;
        end
    end
end
createfig('large','on');
subplot(2,2,1)
for i = 1:4
    plot(1995:1995+size(RMSEatloc,1)-1,squeeze(nanmean(RMSEatloc(:,i,:),3)),'LineWidth',2,'color',cbrew(i,:))
    hold on
end
title(['RMSE - ',num2str(nopredyears)],'fontsize',20);
set(gca,'fontsize',18);
ylabel('RMSE (K)','fontsize',20);
%xlabel('Year','fontsize',20);
subplot(2,2,2)
for i = 1:4    
    plot(1995:1995+size(RMSEatloc,1)-1,squeeze(nanmean(abs(rval(:,i,:)),3)),'LineWidth',2,'color',cbrew(i,:))
    hold on
end
title(['Correlation - ',num2str(nopredyears),' years'],'fontsize',20);
set(gca,'fontsize',18);
ylabel('r','fontsize',20);
xlabel('Year','fontsize',20);
lh = legend('Greenland','Russia','America','Asia');
set(lh,'box','off');

subplot(2,1,2);
m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
%[~,h] = m_contourf(longitude-180,latitude,toplotcorr,-1:.1:1,'LineStyle','none');
m_coast('color','k','LineWidth',1);
m_grid('ytick',0:30:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',fsize);            
%set(sp2,'position',[.175 .07 .7 .3]);
%caxis([-1 1]);        
%ch = colorbar;
%set(ch,'YTick',[-1:.2:1],'fontsize',fsize)%,
hold on
%colormap(cbrew2);    
for i = 1:4
    m_plot(lontoplot(:,:,i),lattoplot(:,:,i),'LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrew(i,:),'MarkerEdgeColor','k','LineWidth',2)
end

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/RMSE_',num2str(yearsall(1)),'_',num2str(yearsall(2)),'nopredyears-',num2str(nopredyears),'_',monthnames(inputs.varmonthtomean,1,1)];
export_fig(filename,'-png')
end
