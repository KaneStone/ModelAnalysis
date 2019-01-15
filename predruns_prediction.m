%% Create linear regression model and estiamte temperature based on toz

% prediction runs correlations to surface temperature

clear all 

tozvar = 'toz';
var = 'TS';
pslvar = 'PSL';
%area = '6090S';
detrend = 1;
detrend_ozone = 1;
contourplots = 1;
individual_contourplots = 0;
lineplots = 1;
percentile = 50;
lats = [-90,-75];
lats2 = [-30,-10];
tozmonth = 11;
varmonth = [12,1,2]; % Can be any number of months in the year (e.g. [12,1,2] for austral summer)
indmonth = 11; % Can be any number of months in the year (e.g. [12,1,2] for austral summer)
pslvarmonth = 11;
if length(varmonth) > 1
    shortnames = 1;            
end
yearplus = 0;
shortnames = 0;
stringtogether = 0;
sig = .05;
strapme = 0;

% read in data 
%directqory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
%files = dir([directory,'*.nc']);
%% Read in lowcl variable

ClLevel = 'lowCl';
timeperiodlow = [1955,1976];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
varfilespast = dir([vardirectory,'*.nc']);

[data.lowcl,years.lowcl,composite.lowcl,dataMonthArrange.lowcl]...
    = predruns_ReadInlayer(vardirectory,varfilespast,var,timeperiodlow,lats,detrend);

%% Read in highcl variable

ClLevel = 'highCl';
timeperiodhigh = [1995,2016];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl]...
    = predruns_ReadInlayer(vardirectory,varfiles,var,timeperiodhigh,lats,detrend);

%% Read in lowcl PSL

ClLevel = 'lowCl';
timeperiodlow = [1955,1976];%[1955,1975]

pslvardirectory = ['/Volumes/ExternalOne/work/data/predruns/',pslvar,'/',ClLevel,'/'];
pslvarfilespast = dir([pslvardirectory,'*.nc']);

[psldata.lowcl,pslyears.lowcl,pslcomposite.lowcl,psldataMonthArrange.lowcl]...
    = predruns_ReadInlayer(pslvardirectory,pslvarfilespast,pslvar,timeperiodlow,lats,detrend);

%% Read in highcl PSL

ClLevel = 'highCl';
timeperiodhigh = [1995,2016];%[1955,1975]

pslvardirectory = ['/Volumes/ExternalOne/work/data/predruns/',pslvar,'/',ClLevel,'/'];
pslvarfiles = dir([pslvardirectory,'*.nc']);

[psldata.highcl,pslyears.highcl,pslcomposite.highcl,psldataMonthArrange.highcl]...
    = predruns_ReadInlayer(pslvardirectory,pslvarfiles,pslvar,timeperiodhigh,lats,detrend);

%% Read in TOZ highcl and take percentiles
ClLevel = 'highCl';
tozdates = [1995,2015];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);

[pct_highcl,tozextract.highcl] = predruns_varPercentiles(toz_composite.highcl.montharrange,toz_dataMonthArrange.highcl,...
    tozmonth,percentile,length(tozfiles));

%% Read in TOZ lowcl and take percentiles
ClLevel = 'lowCl';
tozpastdates = [1955,1975];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfilespast = dir([directory,'*.nc']);
[toz_data.lowcl,toz_years.lowcl,toz_varweighted.lowcl,toz_composite.lowcl,toz_dataMonthArrange.lowcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfilespast,tozvar,tozpastdates,lats,detrend_ozone);

[pct_lowcl,tozextract.lowcl] = predruns_varPercentiles(toz_composite.lowcl.montharrange,toz_dataMonthArrange.lowcl...
    ,tozmonth,percentile,length(tozfilespast));


%% Read in TOZ highcl and take percentiles
ClLevel = 'highCl';
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
[toz_data2.highcl,toz_years2.highcl,toz_varweighted2.highcl,toz_composite2.highcl,toz_dataMonthArrange2.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats2,detrend_ozone);

%% Read in TOZ lowcl and take percentiles
ClLevel = 'lowCl';
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
[toz_data2.lowcl,toz_years2.lowcl,toz_varweighted2.lowcl,toz_composite2.lowcl,toz_dataMonthArrange2.lowcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfilespast,tozvar,tozpastdates,lats2,detrend_ozone);

%% extracting percentiles

% toz
tozextract.lowcl.combine = [toz_composite.lowcl.montharrange(tozmonth,pct_lowcl.highind.a),...
    toz_composite.lowcl.montharrange(tozmonth,pct_lowcl.lowind.a)]';
tozextract.highcl.combine = [toz_composite.highcl.montharrange(tozmonth,pct_highcl.highind.a),...
    toz_composite.highcl.montharrange(tozmonth,pct_highcl.lowind.a)]';

tozextract2.lowcl.combine = [toz_composite2.lowcl.montharrange(tozmonth,pct_lowcl.highind.a),...
    toz_composite2.lowcl.montharrange(tozmonth,pct_lowcl.lowind.a)]';
tozextract2.highcl.combine = [toz_composite2.highcl.montharrange(tozmonth,pct_highcl.highind.a),...
    toz_composite2.highcl.montharrange(tozmonth,pct_highcl.lowind.a)]';

lowcl_tozcorr = corr(tozextract.lowcl.combine,tozextract2.lowcl.combine);
highcl_tozcorr = corr(tozextract.highcl.combine,tozextract2.highcl.combine);

% variable
[varextract,varextractmean,vardifference] = predruns_extractpct(dataMonthArrange,pct_highcl,...
    pct_lowcl,varmonth,length(tozfiles),length(tozfilespast));

[indvarextract,indvarextractmean,indvardifference] = predruns_extractpct(dataMonthArrange,pct_highcl,...
    pct_lowcl,indmonth,length(tozfiles),length(tozfilespast));

[pslvarextract,pslvarextractmean,pslvardifference] = predruns_extractpct(psldataMonthArrange,pct_highcl,...
    pct_lowcl,pslvarmonth,length(tozfiles),length(tozfilespast));

longitude = data(1).lowcl.lon;
latitude = data(1).lowcl.lat;

%% Calculate southern Oscillation index

[SOI] = predruns_calculateSOI(pslvarextract,pslcomposite,latitude,longitude,indmonth);

%% Calculate NINO

[NINO34,IOD,NH] = predruns_calculateENSO(indvarextract,composite,latitude,longitude,indmonth);

%% create regression model using both ozone latitude averages (currently testing)

amt = 10;
lattp = -45; %-62 -36 -30 -65
lontp = 110; % 20 150 245 240å
[~,lattpind] = min(abs(lattp-latitude));
[~,lontpind] = min(abs(lontp-longitude));
userandom = 0;
extractbestfits = 1;
cltouse = 'highcl';
fsize = 10;
predictors = 'ozoneonly'; %SOIonly,ozoneonly,all,polarozone,midlatozone
findneutralSOI = 1;
%% use random 10

if userandom
    createfig('large','on');
    nosubplots = 6;
    for i = 1:nosubplots;

        randind = ceil(rand(10,1)*length(tozextract.(cltouse).combine));
        restind = 1:length(tozextract.(cltouse).combine);
        restind (restind(randind)) = [];
        var_combine = [varextract.(cltouse).highind(:,lattpind,lontpind);varextract.(cltouse).lowind(:,lattpind,lontpind)];
        varrandextract = var_combine(restind);
        tozrandextract = tozextract.(cltouse).combine(restind);
        tozrandextract2 = tozextract2.(cltouse).combine(restind);

        varrest = var_combine(randind);
        tozrest = tozextract.(cltouse).combine(randind);
        tozrest2 = tozextract2.(cltouse).combine(randind);

        [b,bint,r,rint,stats] = regress(varrandextract,[ones(length(tozrandextract),1),...
            tozrandextract,tozrandextract2]);

        % find outliers 
        outlierind = rint(:,1) < 0;
        varrandextract = varrandextract(outlierind);
        tozrandextract = tozrandextract(outlierind);
        tozrandextract2 = tozrandextract2(outlierind);

        % do regression with outliers removed
        [b2,bint2,r2,rint2,stats2] = regress(varrandextract,[ones(length(tozrandextract),1),...
            tozrandextract,tozrandextract2]);

        test = b(2)*tozrest+b(3)*tozrest2;

        test1 = (b(2)*tozrest+b(3)*tozrest2...
            -nanmean(b(2)*tozrest+b(3)*tozrest2))*4 ...
            +nanmean(b(2)*tozrest+b(3)*tozrest2);

        test2 = (b2(2)*tozrest+b2(3)*tozrest2...
            -nanmean(b2(2)*tozrest+b2(3)*tozrest2)) ...
            +nanmean(b2(2)*tozrest+b2(3)*tozrest2);

        t_rest = b(1)+test;
        t_rest1 = b(1)+test1;
        t_rest2 = b2(1)+test2;

        %plot(t_rest,'-o');
        sp(i) = subplot(nosubplots/2,2,i);
        sp_pos(i,:) = get(sp(i),'position');
        box on
        hold on
        plot(varrest,'k-o','LineWidth',2);
        plot(t_rest,'-o','LineWidth',2);
        %plot(t_rest1,'-o','LineWidth',2);
        plot(t_rest2,'-o','LineWidth',2);

        cor = corr(varrest,b(1)+test);
        %cor1 = corr(varrest,b(1)+test1);
        cor2 = corr(varrest,b2(1)+test2);
        if i == 1
            lh = legend('WACCM surface temperatures','Regression model prediction','Regression model with removed outliers');
            set(lh,'position',[.5,.02,.05,.05],'fontsize',fsize,'box','off');
        end
        title(['Random sample ',num2str(i)],'fontsize',fsize+2); 
        set(gca,'fontsize',fsize);
        xlabel('Year','fontsize',fsize+2);
        ylabel('Surface temperature (K)','fontsize',fsize+2);

         annotation('textbox',[sp_pos(i,:)],'String',{['r = ', num2str(cor)];['r = ', num2str(cor2)]},...
             'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize,...
            'EdgeColor','none','fontweight','bold')
    end
        annotation('textbox',[0 .90 1 .1],'String',{['Regression model prediction of surface temperatures, ',num2str(lattp),'S, ',num2str(lontp),'E']},...
            'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,...
            'EdgeColor','none','fontweight','bold')
        
        randfortit = rand(1)*1000;

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','regressionPredictions/',monthnames(tozmonth,0,0),tozvar,'_detrendTS_',...
        num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2)),'_',num2str(100-percentile),...
        'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),...
        'and',num2str(abs(lats2(1))),'-',num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,shortnames),'_',num2str(randfortit),'.pdf'];
    export_fig(filename,'-pdf');
        
end

%%

if strcmp(cltouse,'highcl')
    tptouse = timeperiodhigh;
else
    tptouse = timeperiodlow;
end    
    

toz1 = (tozextract.(cltouse).combine - nanmean(tozextract.(cltouse).combine))./std(tozextract.(cltouse).combine,1);
toz2 = tozextract2.(cltouse).combine; 

divby = 3;
divby2 = 2;
if extractbestfits
    predictorOrder = {'polarOzone','midLatOzone','SOI','NINO34','IOD'};    
    predictorSwitch = logical([1,0,0,0,0]);
    allPredictors = [toz1,toz2,SOI.(cltouse),NINO34.(cltouse),IOD.(cltouse)];
    allPredictors = allPredictors(:,predictorSwitch);
    predictorOrder = predictorOrder(predictorSwitch);
    
    var_combine = [varextract.(cltouse).highind(:,lattpind,lontpind);varextract.(cltouse).lowind(:,lattpind,lontpind)];
    
    [b,bint,r,rint,stats] = regress(var_combine,[ones(length(toz2),1),...
        allPredictors]);
    regressmodel = b(1) + b(2:end)' * allPredictors';
    
    createfig('largeportrait','on')
    sp(1) = subplot(2,2,1);
    sp_pos(1,:) = get(sp(1),'position');
    hold on
    ph(1) = plot(var_combine,'-o','LineWidth',3);
    ph(2) = plot(regressmodel,'-o','LineWidth',3);
    
    set(gca,'fontsize',fsize+4);
    xlabel('','fontsize',fsize+6);
    ylabel('Temperature (K)','fontsize',fsize+6);
    title(['Upper and lower ', num2str(percentile),'th percentiles and regression model fit'],'fontsize',fsize+5);
    box on
    
    pearsonsr = corr(regressmodel',var_combine);
    spearmansr = corr(regressmodel',var_combine,'type','spearman');
    kendallsr = corr(regressmodel',var_combine,'type','kendall');
    
    
    annotation('textbox',sp_pos(1,:),'String',{['r = ', num2str(pearsonsr)];['{\rho} = ', num2str(spearmansr)];['{\tau} = ', num2str(kendallsr)]},...
             'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+6,...
            'EdgeColor','none','fontweight','bold')
    
    lh = legend(ph,'Model surface temperatures','Regression model fit');
    set(lh,'fontsize',fsize+4,'location','NorthEast','box','off');
        
    %% plot residuals
    [sorted_r, sort_rindex] = sort(r,'ascend');
    [sorted_rabs, sort_rindexabs] = sort(abs(r),'ascend');
    sp(2) = subplot(2,2,2);
    plot(r,'LineWidth',3)
    hold on
    plot(sorted_r,'LineWidth',3);
    set(gca,'fontsize',fsize+4);
    xlabel('','fontsize',fsize+6);
    ylabel('Residual (Regression Model - Data (degrees))','fontsize',fsize+6);
    title('Sorted residuals','fontsize',fsize+6);
    
    %filename
    %%
    sp(3) = subplot(2,2,3);
    box on
    residual_min = min(r(sort_rindexabs(1:length(sort_rindexabs)/divby)));
    residual_min2 = min(r(sort_rindexabs(1:length(sort_rindexabs)/(divby/divby2))));
    residual_max = max(r(sort_rindexabs(1:length(sort_rindexabs)/divby)));
    residual_max2 = max(r(sort_rindexabs(1:length(sort_rindexabs)/(divby/divby2))));
    residual1 = find(r >= residual_min & r <= residual_max);
    residual2 = find(r <= residual_min2);
    residual3 = find(r >= residual_max2);
    residual4 = find(r <= residual_min2 | r >= residual_max2);
    
    
    regressmodel_residual1 = b(1)+b(2:end)' * allPredictors(residual1,:)';
        
    sp_pos(3,:) = get(sp(3),'position');
    hold on
%     plot(var_combine(residual1),'-o','LineWidth',3);
%     plot(b(1)+b(2)*toz1(residual1)+b(3)*toz2(residual1),'-o','LineWidth',3);
    plot(var_combine(residual1),'-o','LineWidth',3);
    plot(regressmodel_residual1,'-o','LineWidth',3);
    
    set(gca,'fontsize',fsize+4);
    xlabel('','fontsize',fsize+6);
    ylabel('Residual (Regression Model - Data (degrees))','fontsize',fsize+6);
    title(['Data points with residuals between ', sprintf('%.3f',(residual_min)),' and ',sprintf('%.3f',residual_max),' degrees'],'fontsize',fsize+6);
    pearsonsr_rest = corr(regressmodel_residual1',var_combine(residual1));
    spearmansr_rest = corr(regressmodel_residual1',var_combine(residual1),'type','spearman');
    kendallsr_rest = corr(regressmodel_residual1',var_combine(residual1),'type','kendall');
    
    annotation('textbox',sp_pos(3,:),'String',{['r = ', num2str(pearsonsr_rest)];['{\rho} = ', num2str(spearmansr_rest)];['{\tau} = ', num2str(kendallsr_rest)]},...
             'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+6,...
            'EdgeColor','none','fontweight','bold')
    
        %%
    
    sp(4) = subplot(2,2,4);
    sp_pos(4,:) = get(sp(4),'position');
    box on
    hold on
    neutralSOIindex = find(abs(SOI.(cltouse)) < .25);    
    
    regressmodel_residual4 = b(1) + b(2:end)' * allPredictors(residual4,:)';
    
    plot(var_combine(residual4),'-o','LineWidth',3);
    plot(regressmodel_residual4,'-o','LineWidth',3);
    
    set(gca,'fontsize',fsize+4);
    xlabel('','fontsize',fsize+6);
    ylabel('Residual (Regression Model - Data (degrees))','fontsize',fsize+6);
    title(['Data points with residuals greater than ', sprintf('%.3f',residual_max2),' and less than ',sprintf('%.3f',residual_min2),' degrees'],'fontsize',fsize+6);
    pearsonsr_outliers = corr(regressmodel_residual4',var_combine(residual4));
    spearmansr_outliers = corr(regressmodel_residual4',var_combine(residual4),'type','spearman');
    kendallsr_outliers = corr(regressmodel_residual4',var_combine(residual4),'type','kendall');
    
    annotation('textbox',sp_pos(4,:),'String',{['r = ', num2str(pearsonsr_outliers)];['{\rho} = ', num2str(spearmansr_outliers)];['{\tau} = ', num2str(kendallsr_outliers)]},...
             'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+6,...
            'EdgeColor','none','fontweight','bold')

    if length(varmonth) > 1
        stringtogether = 1;
        shortnames = 1;
    end
    predForTitle = strcat(predictorOrder,',{ }');
    annotation('textbox',[0 .90 1 .1],'String',{['Regression model (',monthnames(tozmonth,0,0),'{ }',[predForTitle{:}],') prediction of ',...
        monthnames(varmonth,stringtogether,shortnames),' surface temperatures, ',num2str(lattp),'S, ',num2str(lontp),'E, ',num2str(tptouse(1)),'-',num2str(tptouse(2))]},...
                'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,...
                'EdgeColor','none','fontweight','bold')

    %%
    if ~exist(['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','regression/testcases/',...
            num2str(abs(lattp)),'S',num2str(abs(lontp)),'E/'],'dir')
        mkdir(['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','regression/testcases/',...
            num2str(abs(lattp)),'S',num2str(abs(lontp)),'E/']);
         if ~exist(['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','regression/testcases/',...
                num2str(abs(lattp)),'S',num2str(abs(lontp)),'E/Figures/'],'dir')
            mkdir(['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','regression/testcases/',...
                num2str(abs(lattp)),'S',num2str(abs(lontp)),'E/Figures/']);
         end
    end
    predForFilename = strcat(predictorOrder,'_');
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','regression/testcases/',...
        num2str(abs(lattp)),'S',num2str(abs(lontp)),'E','/Figures/',monthnames(tozmonth,0,0),tozvar,'_detrendTS_',...
        num2str(tptouse(1)),'-',num2str(tptouse(2)),'_',num2str(100-percentile),...
        'and',num2str(percentile),'percentile_',monthnames(varmonth,1,shortnames),'_',[predForFilename{:}],'.pdf'];

    export_fig(filename,'-pdf');

%     %%
%     filename2 = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','regression/testcases/',...
%         num2str(abs(lattp)),'S',num2str(abs(lontp)),'E','/indexoutput/',monthnames(tozmonth,0,0),tozvar,'_detrendTS_',...
%         num2str(tptouse(1)),'-',num2str(tptouse(2)),'_',num2str(100-percentile),...
%         'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),...
%         'and',num2str(abs(lats2(1))),'-',num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,shortnames),'_',predictors,'.mat'];
% 
%     save(filename2,'residual1','residual2','residual3');
end
%% plotting TS maps of 

res.(['r',num2str(abs(lattp)),'S',num2str(abs(lontp)),'E']).good = residual1;
res.(['r',num2str(abs(lattp)),'S',num2str(abs(lontp)),'E']).low = residual2;
res.(['r',num2str(abs(lattp)),'S',num2str(abs(lontp)),'E']).high = residual3;

var_combineformap = [varextract.(cltouse).highind;varextract.(cltouse).lowind];
difffromlow = squeeze(nanmean(var_combineformap(residual1,:,:))) - squeeze(nanmean(var_combineformap(residual2,:,:)));
difffromhigh = squeeze(nanmean(var_combineformap(residual1,:,:))) - squeeze(nanmean(var_combineformap(residual3,:,:)));
highfromlow = squeeze(nanmean(var_combineformap(residual3,:,:))) - squeeze(nanmean(var_combineformap(residual2,:,:)));


toplot = permute(cat(3,difffromlow,difffromhigh),[3,2,1]);

longitudetoplot = [longitude(end);longitude];
toplot = cat(2,toplot(:,end,:),toplot);

cbrew = cbrewer('div','RdBu',16);         
    
subplottitles = {[monthnames(varmonth,1,shortnames),'{ }',var,' dif. from negative residuals'],...
    [monthnames(varmonth,1,shortnames),'{ }',var,' dif. from positive residuals']};

mtitle = ['Regressors: ',monthnames(tozmonth,0,0),'{ }', [predForTitle{:}],' ',num2str(abs(lattp)),'S, ',num2str(abs(lontp)),'E, ',num2str(tptouse(1)),'-',num2str(tptouse(2))];

subplotmaps(toplot,longitudetoplot,latitude,{'div','RdBu'},1,[],10,subplottitles,'Longitude','','Difference','on',...
    [-1.5,1.5],18,[-90:15:15],[-90:15:15],[longitude(1:10:end)],[longitude(1:10:end)],mtitle,1,[0 360],[-90 20],0,'none',1);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/','regression/testcases/',...
    num2str(abs(lattp)),'S',num2str(abs(lontp)),'E','/Figures/Differences_',monthnames(tozmonth,0,0),tozvar,'_detrendTS_',...
    num2str(tptouse(1)),'-',num2str(tptouse(2)),'_',num2str(100-percentile),...
    'and',num2str(percentile),'percentile_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),...
    'and',num2str(abs(lats2(1))),'-',num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,shortnames),'_',[predForFilename{:}],'.png'];

export_fig(filename,'-png');
%%
