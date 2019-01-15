% pred runs differences between composite high and low chlorine low and
% high ozone and surface temperature
clear all 

tozvar = 'toz';
var = 'TS';
%area = '6090S';
detrend = 1;
contourplots = 1;
individual_contourplots = 0;
lineplots = 1;
percentile = 10;
lats = [-90,-75];
tozmonth = 10;
varmonth = [12,1,2]; % Can be any number of months in the year (e.g. [12,1,2] for austral summer)
yearplus = 0;
shortnames = 0;
% read in data 
%directqory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
%files = dir([directory,'*.nc']);
%% Read in lowcl variable

ClLevel = 'lowCl';
timeperiodlow = [1955,1976];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
varfilespast = dir([vardirectory,'*.nc']);

[data.lowcl,years.lowcl,composite.lowcl,dataMonthArrange.lowcl]...
    = predruns_ReadInlayer(vardirectory,varfilespast,var,timeperiodlow,lats,1);

%% Read in highcl variable

ClLevel = 'highCl';
timeperiodhigh = [1995,2016];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl]...
    = predruns_ReadInlayer(vardirectory,varfiles,var,timeperiodhigh,lats,1);

%% Read in TOZ highcl and take percentiles
ClLevel = 'highCl';
tozdates = [1995,2015];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats);

[pct_highcl,~] = predruns_varPercentiles(toz_composite.highcl.montharrange,toz_dataMonthArrange.highcl,...
    tozmonth,percentile,length(tozfiles));

%% Read in TOZ lowcl and take percentiles
ClLevel = 'lowCl';
tozpastdates = [1955,1975];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfilespast = dir([directory,'*.nc']);
[toz_data.lowcl,toz_years.lowcl,toz_varweighted.lowcl,toz_composite.lowcl,toz_dataMonthArrange.lowcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfilespast,tozvar,tozpastdates,lats);

[pct_lowcl,~] = predruns_varPercentiles(toz_composite.lowcl.montharrange,toz_dataMonthArrange.lowcl,...
    tozmonth,percentile,length(tozfilespast));

%% Extract high and low variable years, average, and difference

[varextract,varextractmean,vardifference] = predruns_extractpct(dataMonthArrange,pct_highcl,...
    pct_lowcl,varmonth,length(tozfiles),length(tozfilespast));

latitude = data.highcl(1).lat;  
longitude = data.highcl(1).lon;  
longitude = [longitude(end);longitude];

%% two-sampled t-test
for i = 1:length(longitude)-1
    for j = 1:length(latitude)        
        [highcl.h(j,i),highcl.p(j,i),highcl.ci(j,i,:),~] = ttest2(squeeze(varextract.highcl.lowind(:,j,i)),squeeze(varextract.highcl.highind(:,j,i)),'Alpha',.05,'Tail','both');
        [lowcl.h(j,i),lowcl.p(j,i),lowcl.ci(j,i,:),~] = ttest2(squeeze(varextract.lowcl.lowind(:,j,i)),squeeze(varextract.lowcl.highind(:,j,i)),'Alpha',.05,'Tail','both');
    end
end

%% Adding in extra longitude
if contourplots
    TScompositedifference_lowcl = cat(2,vardifference.lowcl(:,end),vardifference.lowcl)';
    TScompositedifference_lowcl = reshape(TScompositedifference_lowcl,[1,size(TScompositedifference_lowcl)]);
    lowcl.h = cat(2,lowcl.h(:,end),lowcl.h)';
    lowcl.h = reshape(lowcl.h,[1,size(lowcl.h)]);
    lowcl.h (lowcl.h == 0) = -1;
    lowcl.h (lowcl.h == 1) = 0;
    
    TScompositedifference_highcl = cat(2,vardifference.highcl(:,end),vardifference.highcl)';
    TScompositedifference_highcl = reshape(TScompositedifference_highcl,[1,size(TScompositedifference_highcl)]);
    highcl.h = cat(2,highcl.h(:,end),highcl.h)';
    highcl.h = reshape(highcl.h,[1,size(highcl.h)]);
    highcl.h (highcl.h == 0) = -1;
    highcl.h (highcl.h == 1) = 0;
    
    if length(varmonth) > 1
        shortnames = 1;            
    end
    
    %% plotting composite average surface temperatures lowcl
    cbrew = cbrewer('div','RdBu',16);         
    
    contourtitle = {[monthnames(varmonth,1,shortnames),' TS difference of ', num2str(100-percentile),'th and ' num2str(percentile),...
        'th ',monthnames(tozmonth,0,0),' TOZ',' percentiles over ',...
        num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S, ',num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2))]};       
    
    subplotmaps(TScompositedifference_lowcl,longitude,latitude,{'div','RdBu'},1,lowcl.h,12,contourtitle,'Longitude','','K','on',...
        [-2,2],18,[-90:10:0],[-90:10:0],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 0],0,'none',1);
    
    if detrend
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','compositedifferences/compositeMaps/',monthnames(tozmonth,0,0),var,'_detrendTS_',...
            num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S_Tperiod-',monthnames(varmonth,1,shortnames),'.png'];
    else
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','compositedifferences/compositeMaps/',var,'_TS_',...
            num2str(timeperiodlow(1)),'-',num2str(timeperiodlow(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S.png'];
    end
    
    export_fig(filename,'-png');


    %% highcl
    cbrew = cbrewer('div','RdBu',16);    
    contourtitle = {[monthnames(varmonth,1,shortnames),' TS difference of ', num2str(100-percentile),'th and ' num2str(percentile),...
        'th ',monthnames(tozmonth,0,0),' TOZ',' percentiles over ',...
        num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S, ',num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2))]};     
    subplotmaps(TScompositedifference_highcl,longitude,latitude,{'div','RdBu'},1,highcl.h,12,contourtitle,'Longitude','latitude','K','on',...
        [-2,2],18,[-90:10:0],[-90:10:0],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 0],0,'none',1);
    
    if detrend
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','compositedifferences/compositeMaps/',monthnames(tozmonth,0,0),var,'_detrendTS_',...
            num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),'_',num2str(100-percentile),...
            'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S_Tperiod-',monthnames(varmonth,1,shortnames),'.png'];
    else
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','compositedifferences/compositeMaps/',var,'_TS_',...
            num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),'_',num2str(100-percentile),'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S.png'];
    end
    
    export_fig(filename,'-png');
end
%%
TScompositedifference_individual_highcl = permute(cat(3,TScompositedifference_individual_highcl,TScompositedifference_individual_highcl(:,:,1)),[1,3,2]);
%% here
indtitles = {'No. 1','No. 2','No. 3','No. 4','No. 5','No. 6','No. 7','No. 8','No. 9'};
if individual_contourplots
    cbrew = cbrewer('div','RdBu',16);
    subplotmaps(TScompositedifference_individual_highcl,longitude,latitude,{'div','RdBu'},1,[],12,indtitles,'Longitude','','K','on',...
        [-4,4],18,[-90:10:0],[-90:10:0],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 0],0,'none',1);
    
    subplotmaps(TScompositedifference_individual_lowcl,longitude,latitude,{'div','RdBu'},1,[],12,indtitles,'Longitude','','K','on',...
        [-4,4],18,[-90:10:0],[-90:10:0],[longitude(1:10:end)],[longitude(1:10:end)],'',1,[0 360],[-90 0],0,'none',1);
    
end

if lineplots
    
%% plotting line plots of detrended temperature
    cbrew = cbrewer('qual','Paired',10);
    month = 12;
    fsize = 14;

%% longitude surface temperature line plots
createfig('large','on');
cllevel = 'TSdata_highcl';
if strcmp(cllevel,'TSdata_lowcl')
    timeperiod = timeperiodlow;
    toplot = TSdata_lowcl;        
    indtousehigh = highind_lowcl;
    indtouselow = lowind_lowcl;
    differencetoplot = TScompositedifference_lowcl;
    limits = [1:timeperiod(2)-timeperiod(1)+1:size(TSdata_lowcl,3),size(TSdata_lowcl,3)];
else strcmp(cllevel,'TSdata_highcl')
    timeperiod = timeperiodhigh;
    toplot = TSdata_highcl;
    indtousehigh = highind_highcl;
    indtouselow = lowind_highcl;
    differencetoplot = TScompositedifference_highcl;
    limits = [1:timeperiod(2)-timeperiod(1)+1:size(TSdata_highcl,3),size(TSdata_highcl,3)];
end
    
lattoplot = -30;
lontoplot = 120;
[~,latind] = min(abs(lattoplot-latitude));
[~,lonind] = min(abs(lontoplot-longitude));
fsize = 18;

titles = {[num2str(lattoplot),'{\circ}','S'],'Low ozone years-high ozone years',...
    ['Detrended temperature time series at ',num2str(lattoplot),'N and ',num2str(lontoplot),'E']};

for i = 1:3
    
    if i == 1
        subplot(3,1,i)
        ph = plot(longitude(2:end),squeeze(toplot(:,latind,:)),'color',[.7 .7 .7],'LineWidth',2);
        hold on
        phhighTOZ = plot(longitude(2:end),squeeze(toplot(:,latind,indtousehigh(ozonemonth).a)),'color',cbrew(9,:),'LineWidth',2);        
        phlowTOZ = plot(longitude(2:end),squeeze(toplot(:,latind,indtouselow(ozonemonth).a)),'color',cbrew(7,:),'LineWidth',2);
        %phhighTOZmean = plot(longitude(2:end),nanmean(squeeze(toplot(:,latind,indtousehigh(ozonemonth).a)),2),'color',cbrew(10,:),'LineWidth',3);
        %phlowTOZmean = plot(longitude(2:end),nanmean(squeeze(toplot(:,latind,indtouselow(ozonemonth).a)),2),'color',cbrew(8,:),'LineWidth',3);
        phhighTOZmean = plot(longitude(2:end),TScomposite_high_highcl(:,latind),'color',cbrew(10,:),'LineWidth',3);
        phlowTOZmean = plot(longitude(2:end),TScomposite_low_highcl(:,latind),'color',cbrew(8,:),'LineWidth',3);
    elseif i == 2
        subplot(3,1,i)
        differencetoplot = squeeze(differencetoplot);
        if size(differencetoplot,1) == length(longitude)
            differencetoplot(1,:) = [];
        end
        plot(longitude(2:end),differencetoplot(:,latind),'LineWidth',2)
    else
        count = 1;
        subplot(3,1,3)
        for j = 1:length(datahighcl)
            TSph(j) = plot(squeeze(toplot(lonind,latind,count(j):count(j)+timeperiodhigh(2)-timeperiodhigh(1))),'LineWidth',2,'color',[.7 .7 .7],'Marker','*');
            uistack(TSph(j),'bottom');
            hold on
            
            %indseparated = indtouselow(ozonemonth).a(>= limts(j) & <= limits(j+1))
            
            if j == 1
                TSphlow = plot(indtouselow(ozonemonth).a(indtouselow(ozonemonth).a<(count(j)+timeperiod(2)-timeperiod(1))),...
                    squeeze(toplot(lonind,latind,indtouselow(ozonemonth).a(indtouselow(ozonemonth).a<(count+timeperiod(2)-timeperiod(1))))),...
                    'LineStyle','none','Marker','o','markerFaceColor',cbrew(8,:),'MarkerEdgeColor',cbrew(8,:));
                TSphhigh = plot(indtousehigh(ozonemonth).a(indtousehigh(ozonemonth).a<(count+timeperiod(2)-timeperiod(1)))-(j-1)*(timeperiod(2)-timeperiod(1)+1),...
                    squeeze(toplot(lonind,latind,indtousehigh(ozonemonth).a(indtousehigh(ozonemonth).a<(count+timeperiod(2)-timeperiod(1))))),...
                    'LineStyle','none','Marker','o','markerFaceColor',cbrew(10,:),'MarkerEdgeColor',cbrew(10,:));
            else
                plot(indtouselow(ozonemonth).a(indtouselow(ozonemonth).a<(count(j)+timeperiod(2)-timeperiod(1)+1) &...
                    indtouselow(ozonemonth).a>(count(j-1)+timeperiod(2)-timeperiod(1)+1))-...
                    (j-1)*(timeperiod(2)-timeperiod(1)+1),...
                    squeeze(toplot(lonind,latind,indtouselow(ozonemonth).a(indtouselow(ozonemonth).a<(count(j)+timeperiod(2)-timeperiod(1)+1) &...
                    indtouselow(ozonemonth).a>(count(j-1)+timeperiod(2)-timeperiod(1)+1)))),...
                    'LineStyle','none','Marker','o','markerFaceColor',cbrew(8,:),'MarkerEdgeColor',cbrew(8,:));
                plot(indtousehigh(ozonemonth).a(indtousehigh(ozonemonth).a<(count(j)+timeperiod(2)-timeperiod(1)+1) &...
                    indtousehigh(ozonemonth).a>(count(j-1)+timeperiod(2)-timeperiod(1)+1))-...
                    (j-1)*(timeperiod(2)-timeperiod(1)+1),...
                    squeeze(toplot(lonind,latind,indtousehigh(ozonemonth).a(indtousehigh(ozonemonth).a<(count(j)+timeperiod(2)-timeperiod(1)+1) &...
                    indtousehigh(ozonemonth).a>(count(j-1)+timeperiod(2)-timeperiod(1)+1)))),...
                    'LineStyle','none','Marker','o','markerFaceColor',cbrew(10,:),'MarkerEdgeColor',cbrew(10,:));
            end
                
            count(j+1) = count(j)+timeperiod(2)-timeperiod(1)+1;
        end                
        
    end
        
    if i == 1 || i == 2
        set(gca,'xtick',0:40:360,'xticklabel',0:40:360,'fontsize',fsize-2);
        xlabel('Longitude (degrees)','fontsize',fsize);
        xlim([-5 365])    
    else
        set(gca,'xtick',1:5:timeperiod(2)-timeperiod(1)+1,'xticklabel',timeperiod(1):5:timeperiod(2),'fontsize',fsize-2);
        xlabel('Year','fontsize',fsize);
        xlim([0 timeperiod(2)-timeperiod(1)+2])    
        %ylim([min(squeeze(nanmean(toplot(12,:,latind,:),2)))-min(squeeze(nanmean(toplot(12,:,latind,:),2)))/1000,...
        %    max(squeeze(nanmean(toplot(12,:,latind,:),2)))+max(squeeze(nanmean(toplot(12,:,latind,:),2)))/1000]);
    end
    ylabel('Temperature (K)','fontsize',fsize);    
    title(titles{i},'fontsize',fsize+2)                    
    if i == 1
        lh = legend([phhighTOZmean, phlowTOZmean],[num2str(100-percentile),' percentile'],[num2str(percentile),' percentile']);
        set(lh,'fontsize',fsize-2,'box','off','location','SouthEast');
    end
end
annotation('textbox','string',...
    [num2str(timeperiod(1)),'-',num2str(timeperiod(2)),',{ }',num2str(100-percentile),' and ',num2str(percentile),' TOZ percentiles over ',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S, ',num2str(abs(lattoplot)),'S']...
    ,'fontsize',fsize+4,'fontweight','bold','LineStyle','none','position',[.225,.9,.1,.1])%,'fontsize',fsize+4,'color','k','xoff',0,'yoff',.04);
    
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','compositedifferences/temperatureLinePlots/Lineplot','_detrendedTS_',...
     num2str(timeperiod(1)),'-',num2str(timeperiod(2)),'_',num2str(100-percentile),'and',...
     num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S',...
     '_',num2str(abs(lattoplot)),'S_',num2str(abs(lontoplot)),'E_',monthnames(ozonemonth),'TOZ.pdf'];
export_fig(filename,'-pdf')
        
end
