function [] = predruns_latsfiveyears_obs(tempin,tozdata,inputs,latitude,longitude)

%% defining model data
dataVarMonthAve = permute(tempin(:,:,12+inputs.varmonthtomean:12:end),[3,2,1]);

%% regressio model using data from 1995-2018
modelyears = [1980:2016];
yearsall = [1980,2010];%;2005,2024];

for l = 1:size(yearsall,1)
clearvars bdetrend TSdetrendedmed TSdetrended TSleftdetrended TSdetrendednomed TSleftdetrendedmed 
clearvars TSleftdetrendednomed bozone modelfit modelprediction modelfitmed modelpredictionmed modelfitnomed modelpredictionnomed
%years = [1995,2010];
years = yearsall(l,:);

if years(1) > modelyears(1)
    modelyears2(l).m = years(1):2016;
else
    modelyears2(l).m = [1980:2016];
end

yearind = modelyears >= years(1) & modelyears <= years(2);
yearindleft = modelyears > years(2);
ifmed = 1;

%% extracting ENSO
clearvars alldata
for i = 1:12    
    alldata(i,:,:,:) = tempin(:,:,12+i:12:end-2);
end
alldata = permute(alldata,[1,4,3,2]);
% read in ENSO 
ENSOyearind = predruns_removeENSOforpred(alldata(:,yearind,:,:),latitude,longitude,1);
ENSOall = predruns_removeENSOforpred(alldata,latitude,longitude,1);
laglength = 3;

%% extracing tozdata

tozdataextract = tozdata;

%% calculate 1995-2018 detrend coefficients
ifdetrend = 1;
predictors = [(1:sum(yearind))',ones(size(1:sum(yearind)))'];
predictorsleft = sum(yearind)+1:sum(yearind)+sum(yearindleft);

for i = 1:1
    
    btozdetrend = predictors\tozdataextract(yearind)';
    tozdetrend = tozdataextract(yearind) - btozdetrend(1)*predictors(:,1)'-btozdetrend(2);
    tozdetrend_left = tozdataextract(yearindleft) - btozdetrend(1)*predictorsleft - btozdetrend(2);
    
    ozonepredictors = [tozdetrend',ones(size(tozdataextract(yearind)))'];    
    for j = 1:size(dataVarMonthAve,2)
        % detrending               
        
        % remove ENSO
        for lag = 1:laglength
            ensopredictors = [squeeze(ENSOyearind(inputs.varmonthtomean-lag+1,:))',ones(size(tozdataextract(yearind)))'];
            benso(j,lag,:,:) = ensopredictors\squeeze(dataVarMonthAve(yearind,j,:));

            ensopredictors_alldata = [squeeze(ENSOall(inputs.varmonthtomean-lag+1,:))',ones(size(ENSOall,2),1)];
            benso_alldata(j,lag,:,:) = ensopredictors_alldata\squeeze(dataVarMonthAve(:,j,:));

        end
        
        %find max lag
        [~,bensomax_ind] = max(abs(benso),[],2);                
        [~,bensomax_ind_alldata] = max(abs(benso_alldata),[],2);   

        for li = 1:size(longitude)
            bensomax(li) = benso(j,bensomax_ind(j,1,1,li),1,li);
            bensomax_alldata(li) = benso_alldata(j,bensomax_ind_alldata(j,1,1,li),1,li);
        end
        if i == 5
            abc = 1
        end
        %removing enso
        for li = 1:size(longitude)
            TSdetrended(:,j,li) = squeeze(dataVarMonthAve(yearind,j,li)) - (squeeze(bensomax(li)).*squeeze(ENSOyearind(inputs.varmonthtomean-bensomax_ind(j,1,1,li)+1,:)))';
            TSleftdetrended(:,j,li) = squeeze(dataVarMonthAve(yearindleft,j,li)) - (squeeze(bensomax_alldata(li)).*squeeze(ENSOall(inputs.varmonthtomean-bensomax_ind_alldata(j,1,1,li)+1,yearindleft))');
        end
                
        if ifdetrend
            bdetrend(j,:,:) = predictors\squeeze(TSdetrended(:,j,:));
            TSdetrended(:,j,:) =  squeeze(TSdetrended(:,j,:)) - squeeze(bdetrend(j,1,:))'.*predictors(:,1) - squeeze(bdetrend(j,2,:))';
            TSleftdetrended(:,j,:) = squeeze(TSleftdetrended(:,j,:)) - squeeze(bdetrend(j,1,:))'.*predictorsleft' - squeeze(bdetrend(j,2,:))';
        else
            TSdetrended(:,j,:) =  squeeze(dataVarMonthAve(yearind,j,:))-squeeze(nanmean(dataVarMonthAve(yearind,j,:),2))';
            TSleftdetrended(:,j,:) = squeeze(dataVarMonthAve(yearindleft,j,:))-squeeze(nanmean(dataVarMonthAve(yearind,j,:),2))';
        end
        
        
%         %find max lag
%         [bensomax(j,:,:),bensomax_ind(j,:,:)] = max(abs(benso(j,lag,:,:)),[],3);
%         [bensomax_alldata(j,:,:),bensomax_ind_alldata(j,:,:)] = max(abs(benso_alldata(j,lag,:,:)),[],3);
% 
%         %removing enso
%         TSdetrended(:,j,:) = squeeze(TSdetrended(:,j,:)) - (squeeze(bensomax(j,1,:)).*squeeze(ENSOyearind(inputs.varmonthtomean-bensomax_ind(j,1,:)+1,:)))';
%         TSleftdetrended(:,j,:) = squeeze(TSleftdetrended(:,j,:)) - (squeeze(bensomax(j,1,:)).*squeeze(ENSOall(inputs.varmonthtomean-bensomax_ind(j,1,:)+1,yearindleft)))';        
        
        % creating ozone model
        bozone(j,:,:) = ozonepredictors\squeeze(TSdetrended(:,j,:));
        corr2use(l,j,:) = corr(ozonepredictors(:,1),squeeze(TSdetrended(:,j,:)));
        modelfit(:,j,:) = squeeze(bozone(j,1,:))'.*tozdetrend' + squeeze(bozone(j,2,:))';
        modelprediction(:,j,:) = squeeze(bozone(j,1,:))'.*tozdetrend_left' + squeeze(bozone(j,2,:))';
        
        if ifmed            
            med(:,j) = squeeze(median(TSdetrended(:,j,:)));        
            TSdetrendedmed(:,j,:) = squeeze(TSdetrended(:,j,:)) - med(:,j)';
            TSleftdetrendedmed(:,j,:) = squeeze(TSleftdetrended(:,j,:)) - med(:,j)';                    
            modelfitmed(:,j,:) = squeeze(modelfit(:,j,:))-med(:,j)';
            modelpredictionmed(:,j,:) = squeeze(modelprediction(:,j,:))-med(:,j)';
            
            TSdetrendednomed(:,j,:) = squeeze(TSdetrended(:,j,:));
            TSleftdetrendednomed(:,j,:) = squeeze(TSleftdetrended(:,j,:));
            modelfitnomed(:,j,:) = squeeze(modelfit(:,j,:));
            modelpredictionnomed(:,j,:) = squeeze(modelprediction(:,j,:));
            
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

surfacecombinemed = cat(1,TSdetrendedmed,TSleftdetrendedmed);
modelcombinemed = cat(1,modelfitmed,modelpredictionmed);

surfacecombinenomed = cat(1,TSdetrendednomed,TSleftdetrendednomed);
modelcombinenomed = cat(1,modelfitnomed,modelpredictionnomed);

% areas_lons = [290,320;30,120;240,290;30,120];%lons (Greenland,Russia,America,Asia)
% areas_lats = [60,80;50,75;30,60;15,45];%lats

% areas_lons = [50,150;230,290;40,90];%lons (Russia,America,Asia)
% areas_lats = [60,75;30,60;30,45];%lats

areas_lons = [50,150;40,90];%lons (Russia,Asia)
areas_lats = [60,75;30,45];%lats

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([4,10,8],:);
mon = 4;

[~,russialatind] = min(abs(latitude - 73));
[~,russialonind] = min(abs(longitude - 120));

for i = 1:size(areas_lons,1)
    lats = latitude > areas_lats(i,1) & latitude < areas_lats(i,2);
    lons = longitude > areas_lons(i,1) & longitude < areas_lons(i,2);
    latextract = latitude(lats);
    lonextract = longitude(lons);
    [latmesh,lonmesh] = meshgrid(latextract,lonextract);
    lonmesh = lonmesh';
    latmesh = latmesh';
    
    mult = squeeze(corr2use(l,lats,lons));
    if l == 1            

        if i == 1
            [maxval(i),maxind(i)] = min(mult(:));     
        else
            [maxval(i),maxind(i)] = max(mult(:));     
        end  
        if i ~= 1
            lattoplot(i) = latmesh(maxind(i));
            lontoplot(i) = lonmesh(maxind(i)); 
        else
            lattoplot(i) = latitude(russialatind(i));
            lontoplot(i) = longitude(russialonind(i));
        end
    end
    if i == 1
        rval(l,i) = corr2use(l,russialatind(i),russialonind(i));        
        toplot(l).t(:,i) = surfacecombinenomed(:,russialatind(i),russialonind(i));
        toplotmed(l).t(:,i) = surfacecombinemed(:,russialatind(i),russialonind(i));

        bvaltemp = squeeze(bozone(lats,1,lons));
        bval(i) = bvaltemp(maxind(i));
        
        toplotmodel(l).t(:,i) = modelcombinenomed(:,russialatind(i),russialonind(i));
        toplotmodelmed(l).t(:,i) = modelcombinemed(:,russialatind(i),russialonind(i));
    else
        rval(l,i) = mult(maxind(i));
        toplottempnomed = squeeze(surfacecombinenomed(:,lats,lons));        
        toplottempmed = squeeze(surfacecombinemed(:,lats,lons));        
        toplot(l).t(:,i) = toplottempnomed(:,maxind(i));
        toplotmed(l).t(:,i) = toplottempmed(:,maxind(i));

        bvaltemp = squeeze(bozone(lats,1,lons));
        bval(i) = bvaltemp(maxind(i));

        toplottempmodelnomed = squeeze(modelcombinenomed(:,lats,lons));        
        toplottempmodelmed = squeeze(modelcombinemed(:,lats,lons));        
        toplotmodel(l).t(:,i) = toplottempmodelnomed(:,maxind(i));
        toplotmodelmed(l).t(:,i) = toplottempmodelmed(:,maxind(i));
    end

end

end

%% saving data
save('/Volumes/MyBook/work/data/predruns/output/obsFigure3/obsFig3.mat','rval','lattoplot','lontoplot','toplot','toplotmodel','modelyears2','latitude','longitude','corr2use');


%% plotting combination
%createfig('large','off');

plot_ind_runs = 0;
plot_specific = 1;

cbrew = cbrewer('qual','Set1',10);
cbrewqual2 = cbrew([4,10,8],:);
cbrewqual3 = cbrewqual2./2;
fsize = 18;

if plot_ind_runs


    titles = {'Russia','North America','Asia'};
    marks = {'s','o','v'};

    for run = 1:9
    fig = figure;
    set(fig,'color','white','position',[100 100 1000 1200],'Visible','off')
    titles = {['Russia - 2005' char(8211) '2020'],['Russia - 1995' char(8211) '2020'],...
        ['Northern America - 2005' char(8211) '2020'],['Northern America - 1995' char(8211) '2020'],...
        ['Asia - 2005' char(8211) '2020'],['Asia - 1995' char(8211) '2020']};
    i = 1;
    for j = 1:6
    %     if  j < 4 
    %         toplotind = 1;
    %     else
    %         toplotind = 2;
    %     end
    if j == 1 || j == 3 || j == 5
        toplotind = 1;
    else
        toplotind = 2;
    end
    %     if j == 1 || j == 4
    %         i = 1;
    %     else
    %     end
        sp(j) = subplot(4,2,j);
        sppos = get(sp(j),'position');

        %yyaxis left
        plot(modelyears2(toplotind).m,toplot(toplotind).t(:,run,i),'-o','LineWidth',3,'color',cbrew(1,:))
        hold on
        plot(modelyears2(toplotind).m,toplotmodel(toplotind).t(:,run,i),'-o','LineWidth',3,'color',cbrew(2,:),'LineStyle','-')    
        set(gca,'Ycolor','k','fontsize',fsize-2);
        if j == 1 || j == 3 || j == 5
            ylabel('Temperture anomaly','fontsize',fsize);
            set(sp(j),'position',[sppos(1)+.025,sppos(2),sppos(3),sppos(4)])
        elseif j == 2 || j == 4 || j == 6
                    set(sp(j),'position',[sppos(1)-.025,sppos(2),sppos(3),sppos(4)])
        end
        if j == 5 || j == 6
            xlabel('Year','fontsize',fsize);

        end
        plot([years(2),years(2)],[-20,10],'Color','k','LineStyle','--');    
        plot([modelyears2(toplotind).m(1),modelyears2(toplotind).m(end)],[0,0],'Color','k','LineStyle','--');    
        forlim = [toplot(toplotind).t(:,run,i);toplotmodel(toplotind).t(:,run,i)]; 
        forlimsort = sort(forlim);
    %     if i == 1 || i == 2
    %         ylim([
    %     elseif i == 3 || i == 4
    %         ylim(
    %     elseif i == 5 || i == 6
    %         ylim(
    %     end
            %ylim([forlimsort(1) + 2*forlimsort(2), forlimsort(end) + forlimsort(end)/2.5])
            ylim([forlimsort(1) + forlimsort(1)./2.5, forlimsort(end) + forlimsort(end)/1.5])

        % bar
    %     yyaxis right
    %     toplot2 = sign(toplotmed(toplotind).t(:,run,i));
    %     toplot2 (sign(toplotmed(toplotind).t(:,run,i)) == sign(toplotmodelmed(toplotind).t(:,run,i))) = 1;
    %     toplot2 (sign(toplotmed(toplotind).t(:,run,i)) ~= sign(toplotmodelmed(toplotind).t(:,run,i))) = -1;
    %     
    %     bar(2021:2024,toplot2(end-3:end),'LineWidth',3,'Edgecolor','k','FaceColor',cbrew(8,:))
    %     ylim([-2 8])
    %     set(gca,'ytick',-1:1:1);
    %     set(gca,'Ycolor','k');

    %     if j == 2 || j == 4 || j == 6
    %         ylabel('Anomaly prediction','fontsize',fsize);
    %     end
    %     %plot([modelyears(1),modelyears(end)],[-8,-8],'Color','k','LineStyle','--');    
         title(titles(j),'fontsize',fsize+2);
        if j == 6
            lh = legend('Surface temperature','Prediction');
            set(lh,'Orientation','horizontal','position',[.025,.16,1,.245],'fontsize',fsize,'box','off');
        end

        annotation('textbox',get(sp(j),'position'),'String',['r = ', sprintf('%.2f',rval(toplotind,i,run))],'FitBoxToText','on','LineStyle','none','Fontsize',fsize);
        if j == 2 
            i = 2;
        elseif j == 4
            i = 3;
        end    
    end

    annotation('textbox',[.01 .98 1 0],'String',['Ensemble member ',num2str(run),', ',monthnames(inputs.tozmonth,0,0),' TCO predicting ',monthnames(inputs.varmonthtomean,0,0),' surface temperatures'],...
        'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
            'EdgeColor','none','fontweight','bold');   

    cbrew2 = flipud(cbrewer('div','RdBu',21));
    cbrew2 = [cbrew2(1:10,:);cbrew2(12:end,:)];

    toplotcorr = squeeze(corr2use(1,run,:,:)); 
    toplotcorr = circshift(toplotcorr,size(toplotcorr,2)/2,2);

    sp2 = subplot(4,1,4);  
    m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
    [~,h] = m_contourf(longitude-180,latitude,toplotcorr,-1:.1:1,'LineStyle','none');
    m_coast('color','k','LineWidth',1);
    m_grid('ytick',0:30:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',fsize);            
    set(sp2,'position',[.19 .03 .7 .25]);
    caxis([-1 1]);
    xlabel('Longitude','fontsize',fsize+2);
    ylab = ylabel('Latitude','fontsize',fsize+2);
    set(ylab,'units','normalized','position',[-.09 .5 0]);
    ch = colorbar;
    set(ch,'YTick',[-1:.2:1],'fontsize',fsize)%,
    set(get(ch,'ylabel'),'string','r','fontsize',fsize+2)
    colormap(cbrew2);    
    hold on
    for i = 1:numel(lontoplot)
        if lontoplot(i) > 180
            lontoplot(i) = lontoplot(i) - 360;
        end
    end
    for i = 1:size(areas_lons,1)
        m_plot(lontoplot(run,i),lattoplot(run,i),'LineStyle','none','Marker',marks{i},'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',4)
    end

    if ifmed
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Model_all_',num2str(run),'_',num2str(years(1)),'_med_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
    else
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Model_all',num2str(run),'_',num2str(years(1)),'_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
    end
    %export_fig(filename,'-pdf');            
    print(filename,'-depsc');     

    end
end



%%
toplotind = 1;
%ru = [5,8];
%runs = [4,7,4,7,4,7];
%runs = [4,7,4,7,4,7];
%runs = [4,7,4,7,4,7];
%runs = repmat(ru,1,3);
if plot_specific

    titles = {'Russia','North America','Asia'};
    marks = {'s','o','v'};

    fig = figure;
    set(fig,'color','white','position',[100 100 1000 1200],'Visible','on')
    titles = {['Observations, ', 'Russia - 1980' char(8211) '2010'],...
        ['Observations, ', 'Northern America - 1980' char(8211) '2010'],...
        ['Observations, ', 'Asia - 1980' char(8211) '2010']};
    for j = 1:3
    %     if  j < 4 
    %         toplotind = 1;
    %     else
    %         toplotind = 2;
    %     end
%     if j == 1 || j == 3 || j == 5
%         toplotind = 1;
%     else
%         toplotind = 1;
%     end
    %     if j == 1 || j == 4
    %         i = 1;
    %     else
    %     end
        sp(j) = subplot(4,1,j);
        sppos(j,:) = get(sp(j),'position');

        %yyaxis left
        plot(modelyears2(toplotind).m,toplot(toplotind).t(:,j),'-o','LineWidth',3,'color',cbrew(1,:))
        hold on
        plot(modelyears2(toplotind).m,toplotmodel(toplotind).t(:,j),'-o','LineWidth',3,'color',cbrew(2,:),'LineStyle','-')    
        set(gca,'Ycolor','k','fontsize',fsize-2);
        if j == 1 || j == 3 || j == 5
            ylabel('Temperture anomaly','fontsize',fsize);
            set(sp(j),'position',[sppos(j,1)+.025,sppos(j,2),sppos(j,3),sppos(j,4)])
        elseif j == 2 || j == 4 || j == 6
                    set(sp(j),'position',[sppos(j,1)-.025,sppos(j,2),sppos(j,3),sppos(j,4)])
        end
        if j == 5 || j == 6
            xlabel('Year','fontsize',fsize);

        end
        
        
        sppos(j,:) = get(sp(j),'position');
        
        plot([years(2),years(2)],[-20,10],'Color','k','LineStyle','--');    
        plot([modelyears2(toplotind).m(1),modelyears2(toplotind).m(end)],[0,0],'Color','k','LineStyle','--');    
        forlim = [toplot(toplotind).t(:,j);toplotmodel(toplotind).t(:,j)]; 
        forlimsort = sort(forlim);

        ylim([forlimsort(1) + forlimsort(1)./2.5, forlimsort(end) + forlimsort(end)/1.5])

         title(titles(j),'fontsize',fsize+2);
        if j == 6
            lh = legend('Surface temperature','Prediction');
            set(lh,'Orientation','horizontal','position',[.025,.16,1,.245],'fontsize',fsize,'box','off');
        end

        annotation('textbox',get(sp(j),'position'),'String',['r = ', sprintf('%.2f',rval(j))],...
            'FitBoxToText','on','LineStyle','none','Fontsize',fsize);
        if j == 2 
            i = 2;
        elseif j == 4
            i = 3;
        end    
    end

    annotation('textbox',[.01 .98 1 0],'String',[monthnames(inputs.tozmonth,0,0),' TCO predicting ',...
        monthnames(inputs.varmonthtomean,0,0),' surface temperatures'],...
        'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
            'EdgeColor','none','fontweight','bold');   

    cbrew2 = flipud(cbrewer('div','RdBu',21));
    cbrew2 = [cbrew2(1:10,:);cbrew2(12:end,:)];

    toplotcorr = squeeze(corr2use(1,:,:)); 
    toplotcorr = circshift(toplotcorr,size(toplotcorr,2)/2,2);

    sp2 = subplot(4,1,4);  
    m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
    [~,h] = m_contourf(longitude-180,latitude,toplotcorr,-.75:.0875:.75,'LineStyle','none');
    m_coast('color','k','LineWidth',1);
    m_grid('ytick',0:30:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',fsize-2);            
    title(['No. ',sprintf('%02d',runs(1)),', correlations'],'fontsize',fsize+2);
    set(sp2,'position',[sppos(1,1) .025 sppos(1,3) .25]);
    caxis([-.75 .75]);
    xlabel('Longitude','fontsize',fsize+2);
    ylab = ylabel('Latitude','fontsize',fsize+2);
    set(ylab,'units','normalized','position',[-.09 .5 0]);
%     ch = colorbar;
%     set(ch,'YTick',[-1:.2:1],'fontsize',fsize)%,
%     set(get(ch,'ylabel'),'string','r','fontsize',fsize+2)
    colormap(cbrew2);    
    hold on
    for i = 1:numel(lontoplot)
        if lontoplot(i) > 180
            lontoplot(i) = lontoplot(i) - 360;
        end
    end
    for i = 1:size(areas_lons,1)
        m_plot(lontoplot(i),lattoplot(i),'LineStyle','none','Marker',marks{i},...
            'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',4)
    end
    
    ch = colorbar;
    %set(ch,'Position',[.89 .15 .0125 .1]);
    set(ch,'YTick',[-75:.175:.75],'fontsize',fsize-2)%,
    set(get(ch,'ylabel'),'string','Correlation','fontsize',fsize+2)
    colormap(cbrew2);    
    hold on
    for i = 1:numel(lontoplot)
        if lontoplot(i) > 180
            lontoplot(i) = lontoplot(i) - 360;
        end
    end    
    
    set(gcf,'Renderer','Painters');
    if ifmed
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Model_specific_',num2str(runs(1)),'_',num2str(runs(2)),'_',num2str(years(1)),'_med_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
    else
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Model_all',num2str(runs(1)),'_',num2str(runs(2)),'_',num2str(years(1)),'_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
    end
    %export_fig(filename,'-pdf');            
    %print(filename,'-depsc');             

end

% for l = 1:size(yearsall,1)
% clearvars bdetrend TSdetrended TSleftdetrended bozone modelfit modelprediction
% clearvars TSdetrendedmed TSleftdetrendedmed modelfitmed modelpredictionmed
% clearvars TSdetrendednomed TSleftdetrendednomed modelfitnomed modelpredictionnomed
% %years = [1995,2010];
% years = yearsall(l,:);
% 
% if years(1) > modelyears(1)
%     modelyears2(l).m = years(1):2016;
% else
%     modelyears2(l).m = [1980:2016];
% end
% yearind = modelyears >= years(1) & modelyears <= years(2);
% yearindleft = modelyears > years(2);
% ifmed = 1;
% %% extracing tozdata
% 
% tozdataextract = tozdata;
% 
% %% calculate 1995-2018 detrend coefficients
% ifdetrend = 1;
% predictors = [(1:sum(yearind))',ones(size(1:sum(yearind)))'];
% predictorsleft = sum(yearind)+1:sum(yearind)+sum(yearindleft);
% 
% 
% 
% ozoneleft = tozdataextract(yearindleft);
% 
% btozdetrend = predictors\tozdataextract(yearind)';
% tozdetrend = tozdataextract(yearind)' - btozdetrend(1)*predictors(:,1)-btozdetrend(2);
% tozdetrend_left = ozoneleft' - btozdetrend(1)*predictorsleft' - btozdetrend(2);
% 
% ozonepredictors = [tozdetrend,ones(size(tozdetrend))];
% 
% 
% for i = 1:1
%     
% 
%     for j = 1:size(dataVarMonthAve,2)
%         % detrending       
%         if ifdetrend
%             bdetrend(i,j,:,:) = predictors\squeeze(dataVarMonthAve(yearind,j,:));
%             TSdetrended(i,:,j,:) =  squeeze(dataVarMonthAve(yearind,j,:)) - squeeze(bdetrend(i,j,1,:))'.*predictors(:,1) - squeeze(bdetrend(i,j,2,:))';
%             TSleftdetrended(i,:,j,:) = squeeze(dataVarMonthAve(yearindleft,j,:)) - squeeze(bdetrend(i,j,1,:))'.*predictorsleft' - squeeze(bdetrend(i,j,2,:))';
%         else
%             TSdetrended(i,:,j,:) =  squeeze(dataVarMonthAve(yearind,j,:))-squeeze(nanmean(dataVarMonthAve(yearind,j,:),1))';
%             TSleftdetrended(i,:,j,:) = squeeze(dataVarMonthAve(yearindleft,j,:))-squeeze(nanmean(dataVarMonthAve(yearind,j,:),1))';
%         end
%         
%         
%         % creating ozone model
%         bozone(i,j,:,:) = ozonepredictors\squeeze(TSdetrended(i,:,j,:));
%         corr2use(l,i,j,:) = corr(ozonepredictors(:,1),squeeze(TSdetrended(i,:,j,:)));
%         modelfit(i,:,j,:) = squeeze(bozone(i,j,1,:))'.*ozonepredictors(:,1) + squeeze(bozone(i,j,2,:))';
%         modelprediction(i,:,j,:) = squeeze(bozone(i,j,1,:))'.*tozdetrend_left + squeeze(bozone(i,j,2,:))';
%                 
%         
%         if ifmed
%             med(i,:,j) = squeeze(median(TSdetrended(i,:,j,:)));        
%             modelfitmed(i,:,j,:) = squeeze(modelfit(i,:,j,:)) - med(i,:,j);
%             modelpredictionmed(i,:,j,:) = squeeze(modelprediction(i,:,j,:)) - med(i,:,j);
%             TSdetrendedmed(i,:,j,:) = squeeze(TSdetrended(i,:,j,:)) - med(i,:,j);
%             TSleftdetrendedmed(i,:,j,:) = squeeze(TSleftdetrended(i,:,j,:)) - med(i,:,j);
%             
%             modelfitnomed(i,:,j,:) = squeeze(modelfit(i,:,j,:));
%             modelpredictionnomed(i,:,j,:) = squeeze(modelprediction(i,:,j,:));
%             TSdetrendednomed(i,:,j,:) = squeeze(TSdetrended(i,:,j,:));% - med(i,:,j);
%             TSleftdetrendednomed(i,:,j,:) = squeeze(TSleftdetrended(i,:,j,:));% - med(i,:,j);
%             
%             
%         end
%         
%     end
% end
% 
% % %% testing
% % run = 5;
% % lat = 87;
% % lon = 22;
% % figure;
% % 
% % hold on
% % plot(modelyears,[squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))],'-o','LineWidth',2);
% % %plot(predictorsleft,squeeze(TSleftdetrended(run,:,lat,lon)),'LineWidth',2);
% % %plot(squeeze(dataVarMonthAve(1,:,87,22))-nanmean(squeeze(dataVarMonthAve(1,:,87,22))))
% % plot(modelyears,[squeeze(modelfit(run,:,lat,lon)),squeeze(modelprediction(run,:,lat,lon))],'-o','LineWidth',2);
% % plot([years(2),years(2)],[-10,10],'Color','k','LineStyle','--');
% % minlim = min([squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))])-...
% %     min([squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))])./10;
% % maxlim = max([squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))])+...
% %     max([squeeze(TSdetrended(run,:,lat,lon)),squeeze(TSleftdetrended(run,:,lat,lon))])./10;
% % ylim([minlim maxlim]);
% % lh = legend('Surface temperature anomalies','Model fit and predictions','Fit/Prediction');
% % set(lh,'box','off','location','NorthWest');
% 
% %plot(predictorsleft,squeeze(modelprediction(run,:,lat,lon)),'LineWidth',2);
% %plot(predictorsleft,(tozdataextract(yearindleft,1)-nanmean(tozdataextract(yearindleft,1)))./10);
% 
% %% areas of highest correlations
% 
% surfacecombinemed = cat(2,TSdetrendedmed,TSleftdetrendedmed);
% modelcombinemed = cat(2,modelfitmed,modelpredictionmed);
% 
% surfacecombinenomed = cat(2,TSdetrendedmed,TSleftdetrendedmed);
% modelcombinenomed = cat(2,modelfitmed,modelpredictionmed);
% 
% % areas_lons = [290,320;30,120;240,290;30,120];%lons (Greenland,Russia,America,Asia)
% % areas_lats = [60,80;50,75;30,60;15,45];%lats
% 
% % areas_lons = [300,320;30,150;240,290;65,90];%lons (Greenland,Russia,America,Asia)
% % areas_lats = [60,80;60,75;30,60;30,45];%lats
% 
% areas_lons = [30,150;240,290;65,90];%lons (Greenland,Russia,America,Asia)
% areas_lats = [60,75;30,60;30,45];%lats
% 
% cbrewqual = cbrewer('qual','Set1',10);
% cbrewqual2 = cbrewqual([3,4,10,8],:);
% mon = 4;
% for i = 1:size(areas_lons,1)
%     lats = latitude > areas_lats(i,1) & latitude < areas_lats(i,2);
%     lons = longitude > areas_lons(i,1) & longitude < areas_lons(i,2);
%     latextract = latitude(lats);
%     lonextract = longitude(lons);
%     [latmesh,lonmesh] = meshgrid(latextract,lonextract);
%     lonmesh = lonmesh';
%     latmesh = latmesh';
%     for k = 1:size(corr2use,2)  
%         mult = squeeze(corr2use(l,k,lats,lons));
%         if l == 1                        
%             if i == 1
%                 [maxval(i,k),maxind(i,k)] = min(mult(:));     
%             else
%                 [maxval(i,k),maxind(i,k)] = max(mult(:));     
%             end                        
%             lattoplot(k,i) = latmesh(maxind(i,k));
%             lontoplot(k,i) = lonmesh(maxind(i,k)); 
%         end
% %         rval(l,i,k) = mult(maxind(i,k));
% %         toplottemp = squeeze(surfacecombine(k,:,lats,lons));        
% %         toplot(l).t(:,k,i) = toplottemp(:,maxind(i,k));
% %         bvaltemp = squeeze(bozone(k,lats,1,lons));
% %         bval(i,k) = bvaltemp(maxind(i,k));
% %         toplottempmodel = squeeze(modelcombine(k,:,lats,lons));        
% %         toplotmodel(l).t(:,k,i) = toplottempmodel(:,maxind(i,k));
%         
%         rval(l,i,k) = mult(maxind(i,k));
%         toplottempnomed = squeeze(surfacecombinenomed(k,:,lats,lons));        
%         toplottempmed = squeeze(surfacecombinemed(k,:,lats,lons));        
%         toplot(l).t(:,k,i) = toplottempnomed(:,maxind(i,k));
%         toplotmed(l).t(:,k,i) = toplottempmed(:,maxind(i,k));
%         
%         bvaltemp = squeeze(bozone(k,lats,1,lons));
%         bval(i,k) = bvaltemp(maxind(i,k));
%         
%         toplottempmodelnomed = squeeze(modelcombinenomed(k,:,lats,lons));        
%         toplottempmodelmed = squeeze(modelcombinemed(k,:,lats,lons));        
%         toplotmodel(l).t(:,k,i) = toplottempmodelnomed(:,maxind(i,k));
%         toplotmodelmed(l).t(:,k,i) = toplottempmodelmed(:,maxind(i,k));
% 
%         
% 
%     end
% end
% 
% %% plotting
% 
% cbrew = cbrewer('qual','Set1',10);
% cbrewqual2 = cbrew([4,10,8],:);
% cbrewqual3 = cbrewqual2./2;
% fsize = 18;
% titles = {['Russia - 1990' char(8211) '2016'],['Russia - 1980' char(8211) '2016'],...
%     ['Northern America - 1990' char(8211) '2016'],['Northern America - 1980' char(8211) '2016'],...
%     ['Asia - 1990' char(8211) '2016'],['Asia - 1980' char(8211) '2016']};
% marks = {'s','o','v'};
% for run = 1:1
% 
% 
% 
% createfig('large','off');
% for i = 1:3
%     sp(i) = subplot(3,2,i);
%     %yyaxis left
%     plot(modelyears2(l).m,toplot(l).t(:,run,i),'-o','LineWidth',3,'color',cbrew(1,:))
%     hold on
%     plot(modelyears2(l).m,toplotmodel(l).t(:,run,i),'-o','LineWidth',3,'color',cbrew(2,:),'LineStyle','-')    
%     set(gca,'Ycolor','k','fontsize',fsize-2);
%     ylabel('Temperture anomaly','fontsize',fsize);
%     xlabel('Year','fontsize',fsize);
%     plot([years(2),years(2)],[-20,10],'Color','k','LineStyle','--');    
%     plot([modelyears2(l).m(1),modelyears2(l).m(end)],[0,0],'Color','k','LineStyle','--');    
%     forlim = [toplot(l).t(:,run,i);toplotmodel(l).t(:,run,i)];    
%     forlimsort = sort(forlim);
%     
%     ylim([forlimsort(2) + forlimsort(2), forlimsort(end) + forlimsort(end)./10])
%     
%     % bar
% %     yyaxis right
% %     toplot2 = sign(toplot(l).t(:,run,i));
% %     toplot2 (sign(toplot(l).t(:,run,i)) == sign(toplotmodel(l).t(:,run,i))) = 1;
% %     toplot2 (sign(toplot(l).t(:,run,i)) ~= sign(toplotmodel(l).t(:,run,i))) = -1;
% %     
% %     bar(modelyears2(l).m,toplot2,'LineWidth',3,'Edgecolor','k','FaceColor',cbrew(8,:))
% %     ylim([-2 8])
% %     set(gca,'ytick',-1:1:1);
% %     set(gca,'Ycolor','k');
%     title(titles(i),'fontsize',fsize+2);
% %     ylabel('Anomaly prediction','fontsize',fsize);
%     %plot([modelyears(1),modelyears(end)],[-8,-8],'Color','k','LineStyle','--');    
%     if i == 4
%         lh = legend('Surface temperature','Prediction');
%         set(lh,'Orientation','horizontal','position',[0,.245,1,.245],'fontsize',fsize,'box','off');
%     end
%     annotation('textbox',get(sp(i),'position'),'String',['r = ', sprintf('%.2f',rval(l,i,run))],'FitBoxToText','on','LineStyle','none','Fontsize',fsize);
% end
% 
% annotation('textbox',[.01 1 1 0],'String',[monthnames(inputs.varmonthtomean,0,0),', Ensemble member ',num2str(run)] ,'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
%         'EdgeColor','none','fontweight','bold');    
% 
% cbrew2 = flipud(cbrewer('div','RdBu',21));
% cbrew2 = [cbrew2(1:10,:);cbrew2(12:end,:)];
% 
% 
% 
% toplotcorr1 = squeeze(corr2use(l,run,:,:)); 
% toplotcorr(l,:,:) = circshift(toplotcorr1,size(toplotcorr1,2)/2,2);
% 
% sp2 = subplot(3,1,3);  
% m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
% [~,h] = m_contourf(longitude-180,latitude,squeeze(toplotcorr(l,:,:)),-1:.1:1,'LineStyle','none');
% m_coast('color','k','LineWidth',1);
% m_grid('ytick',0:30:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',fsize);            
% set(sp2,'position',[.175 .07 .7 .3]);
% caxis([-1 1]);        
% ch = colorbar;
% set(ch,'YTick',[-1:.2:1],'fontsize',fsize)%,
% colormap(cbrew2);    
% hold on
% for i = 1:numel(lontoplot)
%     if lontoplot(i) > 180
%         lontoplot(i) = lontoplot(i) - 360;
%     end
% end
% 
% for i = 1:size(areas_lons,1)
%     m_plot(lontoplot(run,i),lattoplot(run,i),'LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',4)
% end
% 
% if ifmed
%     filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Obs_',num2str(run),'_',num2str(years(1)),'_med_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
% else
%     filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Obs_',num2str(run),'_',num2str(years(1)),'_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
% end
% %export_fig(filename,'-pdf');            
% print(filename,'-depsc');            
% 
% end
% end
% %% plotting combination
% %createfig('large','off');
% fig = figure;
% set(fig,'color','white','position',[100 100 1000 1200],'Visible','off')
% titles = {['Russia - 1990' char(8211) '2016'],['Russia - 1980' char(8211) '2016'],...
%     ['Northern America - 1990' char(8211) '2016'],['Northern America - 1980' char(8211) '2016'],...
%     ['Asia - 1990' char(8211) '2016'],['Asia - 1980' char(8211) '2016']};
% i = 1;
% for j = 1:6
% %     if  j < 4 
% %         toplotind = 1;
% %     else
% %         toplotind = 2;
% %     end
% if j == 1 || j == 3 || j == 5
%     toplotind = 1;
% else
%     toplotind = 2;
% end
% %     if j == 1 || j == 4
% %         i = 1;
% %     else
% %     end
%     sp(j) = subplot(4,2,j);
%     sppos = get(sp(j),'position');
%     
% %     yyaxis left
%     plot(modelyears2(toplotind).m,toplot(toplotind).t(:,run,i),'-o','LineWidth',3,'color',cbrew(1,:))
%     hold on
%     plot(modelyears2(toplotind).m,toplotmodel(toplotind).t(:,run,i),'-o','LineWidth',3,'color',cbrew(2,:),'LineStyle','-')    
%     set(gca,'Ycolor','k','fontsize',fsize-2);
%     if j == 1 || j == 3 || j == 5
%         ylabel('Temperture anomaly','fontsize',fsize);
%         set(sp(j),'position',[sppos(1)+.025,sppos(2),sppos(3),sppos(4)])
%     elseif j == 2 || j == 4 || j == 6
%                 set(sp(j),'position',[sppos(1)-.025,sppos(2),sppos(3),sppos(4)])
%     end
%     if j == 5 || j == 6
%         xlabel('Year','fontsize',fsize);
% 
%     end
%     plot([years(2),years(2)],[-20,10],'Color','k','LineStyle','--');    
%     plot([modelyears2(toplotind).m(1)-1,modelyears2(toplotind).m(end)+1],[0,0],'Color','k','LineStyle','--');    
%     forlim = toplot(toplotind).t(:,run,i);    
%     ylim([min(forlim(:)) + min(forlim(:))./2.5, max(forlim(:)) + max(forlim(:))./1.5])
%     xlim([modelyears2(toplotind).m(1)-1,modelyears2(toplotind).m(end)+1]);
%     % bar
% %     yyaxis right
% %     toplot2 = sign(toplot(toplotind).t(:,run,i));
% %     toplot2 (sign(toplot(toplotind).t(:,run,i)) == sign(toplotmodel(toplotind).t(:,run,i))) = 1;
% %     toplot2 (sign(toplot(toplotind).t(:,run,i)) ~= sign(toplotmodel(toplotind).t(:,run,i))) = -1;
% %     
% %     bar(modelyears2(toplotind).m,toplot2,'LineWidth',3,'Edgecolor','k','FaceColor',cbrew(8,:))
% %     ylim([-2 8])
% %     set(gca,'ytick',-1:1:1);
% %     set(gca,'Ycolor','k');
%     title(titles(j),'fontsize',fsize+2);
% %     if j == 2 || j == 4 || j == 6
% %         ylabel('Anomaly prediction','fontsize',fsize);
% %     end
%     %plot([modelyears(1),modelyears(end)],[-8,-8],'Color','k','LineStyle','--');    
%     if j == 6
%         lh = legend('Surface temperature','Prediction');
%         set(lh,'Orientation','horizontal','position',[.025,.16,1,.245],'fontsize',fsize,'box','off');
%     end
%     
%     annotation('textbox',get(sp(j),'position'),'String',['r = ', sprintf('%.2f',rval(toplotind,i,run))],'FitBoxToText','on','LineStyle','none','Fontsize',fsize);
%     if j == 2 
%         i = 2;
%     elseif j == 4
%         i = 3;
%     end    
% end
% 
% annotation('textbox',[.01 .98 1 0],'String',['Observations of ',monthnames(inputs.tozmonth,0,0),' TCO predicting ',monthnames(inputs.varmonthtomean,0,0),' surface temperatures'],...
%     'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
%         'EdgeColor','none','fontweight','bold');   
% 
% cbrew2 = flipud(cbrewer('div','RdBu',21));
% cbrew2 = [cbrew2(1:10,:);cbrew2(12:end,:)];
%     
% toplotcorr = squeeze(corr2use(1,1,:,:)); 
% toplotcorr = circshift(toplotcorr,size(toplotcorr,2)/2,2);
% 
% sp2 = subplot(4,1,4);  
% m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
% [~,h] = m_contourf(longitude-180,latitude,toplotcorr,-1:.1:1,'LineStyle','none');
% m_coast('color','k','LineWidth',1);
% m_grid('ytick',0:30:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',fsize);            
% set(sp2,'position',[.19 .03 .7 .25]);
% caxis([-1 1]);
% xlabel('Longitude','fontsize',fsize+2);
% ylab = ylabel('Latitude','fontsize',fsize+2);
% set(ylab,'units','normalized','position',[-.09 .5 0]);
% ch = colorbar;
% set(ch,'YTick',[-1:.2:1],'fontsize',fsize)%,
% set(get(ch,'ylabel'),'string','r','fontsize',fsize+2)
% colormap(cbrew2);    
% hold on
% for i = 1:numel(lontoplot)
%     if lontoplot(i) > 180
%         lontoplot(i) = lontoplot(i) - 360;
%     end
% end
% for i = 1:size(areas_lons,1)
%     m_plot(lontoplot(run,i),lattoplot(run,i),'LineStyle','none','Marker',marks{i},'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',4)
% end
% 
% if ifmed
%     filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Obs_all_',num2str(run),'_',num2str(years(1)),'_med_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
% else
%     filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Obs_all',num2str(run),'_',num2str(years(1)),'_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
% end
% %export_fig(filename,'-pdf');            
% print(filename,'-depsc');     



end
