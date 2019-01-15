function [] = predruns_latsfiveyears(alldata,dataVarMonthAve,tozdata,inputs,latitude,longitude)

%% regressio model using data from 1995-2018
modelyears = [1995:2024];
yearsall = [1995,2020;2005,2020];%;2005,2024];
for l = 1:size(yearsall,1)
clearvars bdetrend TSdetrendedmed TSdetrended TSleftdetrended TSdetrendednomed TSleftdetrendedmed 
clearvars TSleftdetrendednomed bozone modelfit modelprediction modelfitmed modelpredictionmed modelfitnomed modelpredictionnomed
%years = [1995,2010];
years = yearsall(l,:);

if years(1) > modelyears(1)
    modelyears2(l).m = years(1):2024;
else
    modelyears2(l).m = [1995:2024];
end
yearind = modelyears >= years(1) & modelyears <= years(2);
yearindleft = modelyears > years(2);
ifmed = 1;
%% extracing tozdata

tozdataextract = squeeze(tozdata(:,inputs.tozmonth,:))';

%% calculate 1995-2018 detrend coefficients
ifdetrend = 1;
predictors = [(1:sum(yearind))',ones(size(1:sum(yearind)))'];
predictorsleft = sum(yearind)+1:sum(yearind)+sum(yearindleft);

% read in ENSO 
ENSOyearind = predruns_removeENSOforpred(alldata(:,:,yearind,:,:),latitude,longitude,0);
ENSOall = predruns_removeENSOforpred(alldata,latitude,longitude,0);
laglength = 3;
% read in alldata

for i = 1:size(dataVarMonthAve,1)
    
    btozdetrend = predictors\tozdataextract(yearind,i);
    tozdetrend = tozdataextract(yearind,i) - btozdetrend(1)*predictors(:,1)-btozdetrend(2);
    tozdetrend_left = tozdataextract(yearindleft,i) - btozdetrend(1)*predictorsleft' - btozdetrend(2);
    
    ozonepredictors = [tozdetrend,ones(size(tozdataextract(yearind,i)))]; 
               
    for j = 1:size(dataVarMonthAve,3)
        % detrending              
        
        % remove ENSO
        for lag = 1:laglength
            ensopredictors = [squeeze(ENSOyearind(i,inputs.varmonthtomean-lag+1,:)),ones(size(tozdataextract(yearind,i)))];
            benso(i,j,lag,:,:) = ensopredictors\squeeze(dataVarMonthAve(i,yearind,j,:));

            ensopredictors_alldata = [squeeze(ENSOall(i,inputs.varmonthtomean-lag+1,:)),ones(size(ENSOall,3),1)];
            benso_alldata(i,j,lag,:,:) = ensopredictors_alldata\squeeze(dataVarMonthAve(i,:,j,:));

        end
        
        %find max lag
        [~,bensomax_ind] = max(abs(benso),[],3);                
        [~,bensomax_ind_alldata] = max(abs(benso_alldata),[],3);   

        for li = 1:size(longitude)
            bensomax(li) = benso(i,j,bensomax_ind(i,j,1,1,li),1,li);
            bensomax_alldata(li) = benso_alldata(i,j,bensomax_ind_alldata(i,j,1,1,li),1,li);
        end
        if i == 5
            abc = 1
        end
        %removing enso
        for li = 1:size(longitude)
            TSdetrended(i,:,j,li) = squeeze(dataVarMonthAve(i,yearind,j,li)) - (squeeze(bensomax(li)).*squeeze(ENSOyearind(i,inputs.varmonthtomean-bensomax_ind(i,j,1,1,li)+1,:)))';
            TSleftdetrended(i,:,j,li) = squeeze(dataVarMonthAve(i,yearindleft,j,li)) - (squeeze(bensomax_alldata(li)).*squeeze(ENSOall(i,inputs.varmonthtomean-bensomax_ind_alldata(i,j,1,1,li)+1,yearindleft))');
        end
             
        if ifdetrend
            bdetrend(i,j,:,:) = predictors\squeeze(TSdetrended(i,:,j,:));
            TSdetrended(i,:,j,:) =  squeeze(TSdetrended(i,:,j,:)) - squeeze(bdetrend(i,j,1,:))'.*predictors(:,1) - squeeze(bdetrend(i,j,2,:))';
            TSleftdetrended(i,:,j,:) = squeeze(TSleftdetrended(i,:,j,:)) - squeeze(bdetrend(i,j,1,:))'.*predictorsleft' - squeeze(bdetrend(i,j,2,:))';
        else
            TSdetrended(i,:,j,:) =  squeeze(TSdetrended(i,:,j,:))-squeeze(nanmean(TSdetrended(i,:,j,:),2))';
            %TSleftdetrended(i,:,j,:) = squeeze(dataVarMonthAve(i,yearindleft,j,:))-squeeze(nanmean(dataVarMonthAve(i,yearind,j,:),2))';
        end
        
%         %find max lag
%         [bensomax(i,j,:,:),bensomax_ind(i,j,:,:)] = max(abs(benso(i,j,lag,:,:)),[],3);
%         [bensomax_alldata(i,j,:,:),bensomax_ind_alldata(i,j,:,:)] = max(abs(benso_alldata(i,j,lag,:,:)),[],3);
% 
%         %removing enso
%         TSdetrended(i,:,j,:) = squeeze(TSdetrended(i,:,j,:)) - (squeeze(bensomax(i,j,1,:)).*squeeze(ENSOyearind(i,inputs.varmonthtomean-bensomax_ind(i,j,1,:)+1,:)))';
%         TSleftdetrended(i,:,j,:) = squeeze(TSleftdetrended(i,:,j,:)) - (squeeze(bensomax(i,j,1,:)).*squeeze(ENSOall(i,inputs.varmonthtomean-bensomax_ind(i,j,1,:)+1,yearindleft)))';

        % creating ozone model
        bozone(i,j,:,:) = ozonepredictors\squeeze(TSdetrended(i,:,j,:));
        corr2use(l,i,j,:) = corr(ozonepredictors(:,1),squeeze(TSdetrended(i,:,j,:)));
        modelfit(i,:,j,:) = squeeze(bozone(i,j,1,:))'.*tozdetrend + squeeze(bozone(i,j,2,:))';
        modelprediction(i,:,j,:) = squeeze(bozone(i,j,1,:))'.*tozdetrend_left + squeeze(bozone(i,j,2,:))';
        
        if ifmed            
            med(i,:,j) = squeeze(median(TSdetrended(i,:,j,:)));        
            TSdetrendedmed(i,:,j,:) = squeeze(TSdetrended(i,:,j,:)) - med(i,:,j);
            TSleftdetrendedmed(i,:,j,:) = squeeze(TSleftdetrended(i,:,j,:)) - med(i,:,j);                    
            modelfitmed(i,:,j,:) = squeeze(modelfit(i,:,j,:))-med(i,:,j);
            modelpredictionmed(i,:,j,:) = squeeze(modelprediction(i,:,j,:))-med(i,:,j);
            
            TSdetrendednomed(i,:,j,:) = squeeze(TSdetrended(i,:,j,:));
            TSleftdetrendednomed(i,:,j,:) = squeeze(TSleftdetrended(i,:,j,:));
            modelfitnomed(i,:,j,:) = squeeze(modelfit(i,:,j,:));
            modelpredictionnomed(i,:,j,:) = squeeze(modelprediction(i,:,j,:));
            
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

surfacecombinemed = cat(2,TSdetrendedmed,TSleftdetrendedmed);
modelcombinemed = cat(2,modelfitmed,modelpredictionmed);

surfacecombinenomed = cat(2,TSdetrendednomed,TSleftdetrendednomed);
modelcombinenomed = cat(2,modelfitnomed,modelpredictionnomed);


% areas_lons = [290,320;30,120;240,290;30,120];%lons (Greenland,Russia,America,Asia)
% areas_lats = [60,80;50,75;30,60;15,45];%lats

% areas_lons = [50,150;230,290;40,90];%lons (Russia,America,Asia)
% areas_lats = [60,75;30,60;30,45];%lats

areas_lons = [50,150;60,120];%lons (Russia,America,Asia)
areas_lats = [60,75;30,45];%lats

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([4,10,8],:);
mon = 4;
for i = 1:size(areas_lons,1)
    lats = latitude > areas_lats(i,1) & latitude < areas_lats(i,2);
    lons = longitude > areas_lons(i,1) & longitude < areas_lons(i,2);
    latextract = latitude(lats);
    lonextract = longitude(lons);
    [latmesh,lonmesh] = meshgrid(latextract,lonextract);
    lonmesh = lonmesh';
    latmesh = latmesh';
    for k = 1:size(corr2use,2)  
        mult = squeeze(corr2use(1,k,lats,lons));
        if l == 1            
            
            if i == 1
                [maxval(i,k),maxind(i,k)] = min(mult(:));     
            else
                [maxval(i,k),maxind(i,k)] = max(mult(:));     
            end                        
            lattoplot(k,i) = latmesh(maxind(i,k));
            lontoplot(k,i) = lonmesh(maxind(i,k)); 
        end
        rval(l,i,k) = mult(maxind(i,k));
        toplottempnomed = squeeze(surfacecombinenomed(k,:,lats,lons));        
        toplottempmed = squeeze(surfacecombinemed(k,:,lats,lons));        
        toplot(l).t(:,k,i) = toplottempnomed(:,maxind(i,k));
        toplotmed(l).t(:,k,i) = toplottempmed(:,maxind(i,k));
        
        bvaltemp = squeeze(bozone(k,lats,1,lons));
        bval(i,k) = bvaltemp(maxind(i,k));
        
        toplottempmodelnomed = squeeze(modelcombinenomed(k,:,lats,lons));        
        toplottempmodelmed = squeeze(modelcombinemed(k,:,lats,lons));        
        toplotmodel(l).t(:,k,i) = toplottempmodelnomed(:,maxind(i,k));
        toplotmodelmed(l).t(:,k,i) = toplottempmodelmed(:,maxind(i,k));

    end
end

end



%% plotting combination
%createfig('large','off');

plot_ind_runs = 0;
plot_specific = 1;

cbrew = cbrewer('qual','Set1',10);
cbrewqual2 = cbrew([4,8],:);
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


%% read in observations

obstoplot = load('/Volumes/ExternalOne/work/data/predruns/output/obsFigure3/obsFig3.mat');

%%
toplotind = 1;
ru = [6,8];
%runs = [4,7,4,7,4,7];
%runs = [4,7,4,7,4,7];
%runs = [4,7,4,7,4,7];
runs = repmat(ru,2,1);
runs = runs(:);
xlims = [1978, 2018;1978, 2018;...
    1993,2026;1993,2026;...
    1993,2026;1993,2026];
labels = {'a','b','c','d','e','f'};
if plot_specific

    titles = {'Russia','North America','Asia'};
    marks = {'s','o','v'};

    fig = figure;
    set(fig,'color','white','position',[100 100 1000 1200],'Visible','on')
    titles = {['Observations, ', ' central Russia - 1980' char(8211) '2010'],...
        ['Observations, ', ' southern Asia - 1980' char(8211) '2010'],...        
        ['No. ',sprintf('%02d',runs(1)),', central Russia - 1995' char(8211) '2020'],...
        ['No. ',sprintf('%02d',runs(2)),', southern Asia - 1995' char(8211) '2020'],...
        ['No. ',sprintf('%02d',runs(3)),', central Russia - 1995' char(8211) '2020'],...
        ['No. ',sprintf('%02d',runs(4)),', southern Asia - 1995' char(8211) '2020']};
    areaind = [1,2,1,2,1,2];
    moveplot = [3,3,2,2,1,1];
    for j = 1:6
    
        sp(j) = subplot(3,2,j);
        sppos(j,:) = get(sp(j),'position');
    
        if j == 1 ||j == 2
            plot(obstoplot.modelyears2(toplotind).m,obstoplot.toplot(toplotind).t(:,areaind(j)),'-o','LineWidth',3,'color',cbrew(1,:))
            hold on
            plot(obstoplot.modelyears2(toplotind).m,obstoplot.toplotmodel(toplotind).t(:,areaind(j)),'-o','LineWidth',3,'color',cbrew(2,:),'LineStyle','-')    
        else
            plot(modelyears2(toplotind).m,toplot(toplotind).t(:,runs(j-2),areaind(j)),'-o','LineWidth',3,'color',cbrew(1,:))
            hold on
            plot(modelyears2(toplotind).m,toplotmodel(toplotind).t(:,runs(j-2),areaind(j)),'-o','LineWidth',3,'color',cbrew(2,:),'LineStyle','-')    
        end
            
        set(gca,'Ycolor','k','fontsize',fsize-2);
        if j == 1 || j == 3 || j == 5
            ylabel('Temperture anomaly','fontsize',fsize);
            set(sp(j),'position',[sppos(j,1)+.025,sppos(j,2),sppos(j,3),sppos(j,4)])
        elseif j == 2 || j == 4 || j == 6
                    set(sp(j),'position',[sppos(j,1)-.025,sppos(j,2),sppos(j,3),sppos(j,4)])
        end
        sppos(j,:) = get(sp(j),'position');
        if j == 1 || j == 2           
            set(sp(j),'position',[sppos(j,1),sppos(j,2)-.075,sppos(j,3),sppos(j,4)])
        elseif j == 3 || j == 4
            set(sp(j),'position',[sppos(j,1),sppos(j,2)-.05,sppos(j,3),sppos(j,4)])
        elseif j == 5 || j == 6
            set(sp(j),'position',[sppos(j,1),sppos(j,2)-.025,sppos(j,3),sppos(j,4)])
        end
        
        if j == 5 || j == 6
            xlabel('Year','fontsize',fsize);

        end
        
        
        sppos(j,:) = get(sp(j),'position');
        
        if j == 1 || j == 2
            plot([2010,2010],[-30,30],'Color','k','LineStyle','--');    
            plot([1980,2020],[0,0],'Color','k','LineStyle','--');            
            forlim = [obstoplot.toplot(toplotind).t(:,areaind(j));obstoplot.toplotmodel(toplotind).t(:,areaind(j))]; 
            forlimsort = sort(forlim);
        else
            plot([years(2),years(2)],[-30,30],'Color','k','LineStyle','--');    
            plot([modelyears2(toplotind).m(1),modelyears2(toplotind).m(end)],[0,0],'Color','k','LineStyle','--');            
            forlim = [toplot(toplotind).t(:,runs(j-2),areaind(j));toplotmodel(toplotind).t(:,runs(j-2),areaind(j))]; 
            forlimsort = sort(forlim);
        end

        ylim([forlimsort(1) + forlimsort(1)./2.5, forlimsort(end) + forlimsort(end)/1.5])

         title(titles(j),'fontsize',fsize);
        if j == 6
            lh = legend('Surface temperature','Prediction');
            set(lh,'Orientation','horizontal','position',[.025,0,1,.05],'fontsize',fsize,'box','off');
        end
        if j == 1 || j == 2
            annotation('textbox',[sppos(j,1),sppos(j,2)-.0175,sppos(j,3:4)],'String',['r = ', sprintf('%.2f',obstoplot.rval(areaind(j)))],...
                'FitBoxToText','on','LineStyle','none','Fontsize',fsize);
        else
            annotation('textbox',[sppos(j,1),sppos(j,2)-.0175,sppos(j,3:4)],'String',['r = ', sprintf('%.2f',rval(toplotind,areaind(j),runs(j-2)))],...
                'FitBoxToText','on','LineStyle','none','Fontsize',fsize);
        end
        if j == 4 
            i = 2;
        %elseif j == 4
        %    i = 3;
        end
        xlim([xlims(j,1) xlims(j,2)]);
        annotation('textbox',[sppos(j,1),sppos(j,2)+.005,sppos(j,3:4)],'String',labels{j},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize,... % 
        'EdgeColor','none','fontweight','bold');    
    end

    
    
    annotation('textbox',[.01 .915 1 0],'String',[monthnames(inputs.tozmonth,0,0),' TCO predicting ',...
        monthnames(inputs.varmonthtomean,0,0),' surface temperatures'],...
        'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
            'EdgeColor','none','fontweight','bold');   

    cbrew2 = flipud(cbrewer('div','RdBu',21));
    cbrew2 = [cbrew2(1:10,:);cbrew2(12:end,:)];

    
    
    set(gcf,'Renderer','Painters');
    if ifmed
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Model_specific_ensoproper',num2str(runs(1)),'_',num2str(runs(3)),'_',num2str(years(1)),'_med_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
    else
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/Model_all',num2str(runs(1)),'_',num2str(runs(2)),'_',num2str(years(1)),'_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
    end
    %export_fig(filename,'-pdf');            
    print(filename,'-depsc');             

    %% plotting contours rvals 
    %toplotcorr = squeeze(corr2use(1,runs(1),:,:)); 
    %toplotcorr = circshift(toplotcorr,size(toplotcorr,2)/2,2);
    
    % interp onto waccm
    for i = 1:length(obstoplot.latitude)
        obscorr2use2(1,i,:) = interp1(obstoplot.longitude,squeeze(obstoplot.corr2use(1,i,:)),longitude);
    end
    for i = 1:length(longitude)
        obscorr2use3(1,1,:,i) = interp1(obstoplot.latitude,squeeze(obscorr2use2(1,:,i)),latitude);
    end
            
    
    toplotcorr = cat(1,obscorr2use3(1,:,:,:),corr2use(1,ru(1),:,:),corr2use(1,ru(2),:,:));
    toplotcorr = squeeze(toplotcorr);
   %toplotcorr = circshift(toplotcorr,size(toplotcorr,2)/2,2);
    
   titles = {['Observations, 1980',char(8211),'2010'],['No. ',sprintf('%02d',runs(1)),', 1995',char(8211),'2020'],['No. ',sprintf('%02d',runs(3)),', 1995',char(8211),'2020']};
   toplotcorr = permute(toplotcorr,[1,3,2]);
    [fig,ch] = subplotmaps(toplotcorr,longitude,latitude,{'div','RdBu'},1,[],16,titles,'Longitude','Latitude','Correlation','on',...
        [-.75,.75],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{'March TCO - April surface temperature correlations'},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
    
    marks = {'s','v'};

    childaxes = get(gcf,'children');
    axesind = [5,4,2];
    for j = 1:3
        for i = 1:size(areas_lons,1)
            fig.CurrentAxes = childaxes(axesind(j));
            hold on
            sppos = get(gca,'position');
            if j == 1
                m_plot(obstoplot.lontoplot(i),obstoplot.lattoplot(i),'LineStyle','none','Marker',marks{i},...
                    'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',4)
            else
                m_plot(lontoplot(runs(j),i),lattoplot(runs(j),i),'LineStyle','none','Marker',marks{i},...
                    'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',4)
            end
            annotation('textbox',[sppos(1)+.11,sppos(2)+.03,sppos(3:4)],'String',labels{j},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize,... % 
                'EdgeColor','none','fontweight','bold');    
        end
    end

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/lines/corrModel_specific_',num2str(runs(1)),'_',num2str(runs(3)),'_',num2str(years(1)),'_med_',num2str(years(2)),monthnames(inputs.varmonthtomean,1,1),'.eps'];
    set(gcf,'Renderer','Painters');
    print(filename,'-depsc');             
    
end
