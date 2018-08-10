% Read in all data and take regression
clear variables
clc
close all
%% Read in All HighCl data

%% import SWOOSH
Stimeperiod = [2000 2014];
%[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201611-latpress-10deg-L31.nc');
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedeqfillo3q_swoosh-v02.6-198401-201712-latpress-2.5deg-L31.nc');
sfields = fieldnames(SWOOSH);
SWOOSH.(sfields{1}) = cat(3,SWOOSH.(sfields{1}),zeros(size(SWOOSH.(sfields{1}),1),size(SWOOSH.(sfields{1}),2),1));
SWOOSH.(sfields{1}) (SWOOSH.(sfields{1}) == 0) = NaN;
%read in SWOOSH regression functions
[SWOOSHregfun] = ozoneRegression_SWOOSHregfun(Stimeperiod);

SWOOSHyears = repmat(1984:2017,[12,1]);
SWOOSHextract = permute(SWOOSH.(sfields{1})(:,:,SWOOSHyears >= Stimeperiod(1) & SWOOSHyears <= Stimeperiod(2)),[2,1,3]);
%% import highCl
highCllevel = [SWOOSH.level;[.9,.8,.7,.6,.5,.4,.3,.2,.1]'];
highClTimperiod = [2000 2009];
var = 'T';
if strcmp(var,'T')
    vartitle = 'temperature';
elseif strcmp(var,'Z3')
    vartitle = 'geopotential Height';
elseif strcmp(var,'U')
    vartitle = 'zonal Wind';
elseif strcmp(var,'O3')
    vartitle = 'ozone';
end
[highClData,WAClat] = ReadInHighClRegression(highCllevel,highClTimperiod,'highCl',var);

%% take regression of highcl
clearvars O3Anomaly predictors b
for i = 1:length(highClData)
    tic;
   [b(i),predictors(i,:,:),O3Anomaly(i)] = ...
       ozoneRegressionTrends(highClData(i).(var),highClData(i).Uregfun,0,0,highClData(i).NINO34,...
       highClData(i).HFsouth,highClData(i).HFnorth,0,0,highClTimperiod,0,0,0,var); % 
   
   toc;
end

%% plot residuals
plotres = 0;
if plotres
    lats = [-30 30];
    lev = [25,150];
    latind = WAClat >= lats(1) & WAClat <= lats(2);   
    levind = highCllevel >= lev(1) & highCllevel <= lev(2);    
    levind2 = find(levind);    
    
    clearvars highclO3 highclO3_Vertave highclO3vert
    for i = 1:10      
        for j = 1:sum(levind)        
            highclO3_anomaly(j,i,:) = weightedaverage(squeeze(O3Anomaly(i).percent(:,levind2(j),latind))',WAClat(latind));    
            highclO3_anomaly_res(j,i,:) = weightedaverage(squeeze(O3Anomaly(i).percent_residuals_months(:,levind2(j),latind))',WAClat(latind));    
        end
        highclO3_anomaly_vertave(i,:) = squeeze(nanmean(highclO3_anomaly(:,i,:),1));
        highclO3_anomaly_res_vertave(i,:) = squeeze(nanmean(highclO3_anomaly_res(:,i,:),1));            
    end
    
    createfig('medium','on')
    for i = 1:10
        if i == 10
            plot(highclO3_anomaly_vertave(i,:),'r','LineWidth',2);
        else
            plot(highclO3_anomaly_vertave(i,:),'b','LineWidth',2);
        end
        hold on
    end
    for i = 1:10
        if i ~= 10            
            plot(highclO3_anomaly_res_vertave(i,:),'color',[.6 .6 .6],'LineWidth',2);
        else
            plot(highclO3_anomaly_res_vertave(i,:),'color','k','LineWidth',2);
        end
        hold on
    end


    %% constructing correlation matrix
    clearvars corr_raw corr_res
    for i = 1:9
        for j = 1:9        
            corr_raw(i,j) = corr(highclO3_anomaly_vertave(i,:)',highclO3_anomaly_vertave(j,:)');
            corr_res(i,j) = corr(highclO3_anomaly_res_vertave(i,:)',highclO3_anomaly_res_vertave(j,:)');
        end
    end

    highclO3_anomaly_vertave(10,:) = [];
    highclO3_anomaly_res_vertave(10,:) = [];
    %% plot histogram
    fig = figure;
    fsize = 7;
    set(fig,'color','white','position',[100 100 500 350]);
    histogram(corr_raw,-.2:.1:.5,'Normalization','probability')
    hold on
    histogram(corr_res,-.2:.1:.5,'Normalization','probability')
    set(gca,'fontsize',fsize+1)
    ylabel('Percent','fontsize',fsize+2);
    xlabel('Correlation between ensemble members','fontsize',fsize+2);
    title('Correlation matrix values between 30{\circ}S and 30{\circ}N, and 150 and 25 hPa','fontsize',fsize+4);
    lh = legend('Ozone anomalies','Ozone anomaly residuals after regression');
    set(lh,'box','off','fontsize',fsize+2);

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/SolarandTrendsPaper/Draft/Review/MyReply/NewFigures/Randomness_pdf'];

    print(gcf, '-dpdf', '/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/SolarandTrendsPaper/Draft/Review/MyReply/NewFigures/Randomness_pdf.pdf'); 

    %%

    fig = figure;
    fsize = 7;
    set(fig,'color','white','position',[100 100 500 350]);
    histogram(highclO3_anomaly_vertave(:),-17:1:17,'Normalization','probability')
    hold on
    histogram(highclO3_anomaly_res_vertave(:),-17:1:17,'Normalization','probability')
    set(gca,'fontsize',fsize+1)
    ylabel('Percent','fontsize',fsize+2);
    xlabel('Correlation between ensemble members','fontsize',fsize+2);
    title('Correlation matrix values between 30{\circ}S and 30{\circ}N, and 150 and 25 hPa','fontsize',fsize+4);
    lh = legend('Ozone anomalies','Ozone anomaly residuals after regression');
    set(lh,'box','off','fontsize',fsize+2);
    ylim([0 .20]);
    print(gcf, '-dpdf', '/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/SolarandTrendsPaper/Draft/Review/MyReply/NewFigures/Randomness_pdf_ofactualvalues.pdf');

%%

    corr_raw (corr_raw >= .9) = NaN;
    corr_res (corr_res >= .9) = NaN;

    corr_raw = [corr_raw,corr_raw(:,1)];
    corr_raw = [corr_raw;corr_raw(1,:)];

    corr_res = [corr_res,corr_res(:,1)];
    corr_res = [corr_res;corr_res(1,:)];

    fig = figure;
    set(fig,'position',[100 100 500 1000],'color','white')
    cbrew = cbrewer('seq','Reds',26);
    cbrew = cbrew(4:end-3,:);
    cbrew2 = cbrewer('seq','Blues',10);
    cbrew2  = flipud(cbrew2(4:7,:));
    cbrewtouse = [cbrew2;cbrew];
    %cbrew = flipud(cbrew);

    sh(1) = subplot(2,1,1);
    pcolor(1:10,1:10,corr_raw);
    colormap(cbrewtouse);
    caxis([-.1 .5])
    subpos(1,:) = get(sh(1),'position');
    set(gca,'fontsize',16);
    %xlabel('No.','fontsize',20);
    ylabel('No.','fontsize',18);
    set(gca,'color',[.6 .6 .6],'xtick',1.5:1:9.5,'xticklabel',1:1:9,'ytick',1.5:1:9.5,'yticklabel',1:1:9);
    set(sh(1),'position',[subpos(1,1)-.05,subpos(1,2)-.1,subpos(1,3:4)-.05]);
    title('Ozone anomalies','fontsize',18);

    sh(2) = subplot(2,1,2);
    pcolor(1:10,1:10,corr_res);
    colormap(cbrewtouse);
    caxis([-.1 .5])
    subpos(2,:) = get(gca,'position');


    ylabel('No.','fontsize',18);
    xlabel('No.','fontsize',18);
    set(gca,'fontsize',16);
    set(gca,'color',[.6 .6 .6],'xtick',1.5:1:9.5,'xticklabel',1:1:9,'ytick',1.5:1:9.5,'yticklabel',1:1:9);
    set(sh(2),'position',[subpos(2,1)-.05,subpos(2,2),subpos(2,3:4)-.05]);
    ch = colorbar;
    set(ch,'position',[subpos(2,1)+subpos(2,3)-.08,subpos(2,2),.05,.665]);
    set(get(ch,'ylabel'),'string','Correlation','fontsize',18)
    set(gca,'color',[.6 .6 .6]);

    title('Ozone anomaly residuals after regression','fontsize',18);

    annotation('textbox',[.01 .87 .95 0],'String','Correlation matrices','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',22,... % 
        'EdgeColor','none','fontweight','bold');    
    annotation('textbox',[.01 .84 .95 0],'String','(between 30{\circ}S and 30{\circ}N, and 150 and 25 hPa)','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',20,... % 
        'EdgeColor','none','fontweight','bold');    

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/SolarandTrendsPaper/Draft/Review/MyReply/NewFigures/Randomness'];
    export_fig(filename,'-pdf');

    %% Correlation matrix of ozone anomalies between 30S and 30N, and 150 and 25 hPa

    rawxcorr = xcorr(highclO3_anomaly_vertave(10,:));
    rawxcorr2 = xcorr(highclO3_anomaly_res_vertave(10,:));

    rawxcorr = rawxcorr(ceil(length(rawxcorr)./2):end);
    rawxcorr2 = rawxcorr2(ceil(length(rawxcorr2)./2):end);

    figure
    plot(rawxcorr)
    hold on
    plot(rawxcorr2)

end
%%

O3Anomaly(11).ppmv_residuals_months = nanmean(cat(4,O3Anomaly(1:9).ppmv_residuals_months),4);
O3Anomaly(11).percent_residuals_months = nanmean(cat(4,O3Anomaly(1:9).percent_residuals_months),4);
%O3Anomaly(11).percent_residuals_months_AR1 = nanmean(cat(4,O3Anomaly(1:9).percent_residuals_months_AR1),4);
O3Anomaly(11).percent_residuals = nanmean(cat(4,O3Anomaly(1:9).percent_residuals),4);
%O3Anomaly(11).regmodel_months = nanmean(cat(4,O3Anomaly(1:9).regmodel_months),4);
% O3Anomaly(11).regmodel = nanmean(cat(4,O3Anomaly(1:9).regmodel),4);

%%

%ensavestd(O3Anomaly);

%%

[btest,balltest,ball_pvalue] = ozoneRegressionEnsAve(O3Anomaly(11).percent_residuals_months,highClTimperiod(2) - highClTimperiod(1));
[btest2,balltest2] = ozoneRegressionEnsAve(O3Anomaly(11).percent_residuals_months,highClTimperiod(2) - highClTimperiod(1));

sig = .01;
ball_pvalue_toplot = ball_pvalue;
ball_pvalue_toplot (ball_pvalue_toplot <= sig) = 0;
ball_pvalue_toplot (ball_pvalue_toplot > sig) = 1;

%% individual global trends
for i = 1:9
    if strcmp(var,'T') || strcmp(var,'Z3') || strcmp(var,'U')    
        %[bhclind(i).b,ballhclind(i).b] = ozoneRegressionEnsAve(O3Anomaly(i).ppmv_residuals_months,highClTimperiod(2) - highClTimperiod(1));
        [bhclind(i).b,ballhclind(i).b] = ozoneRegressionEnsAve(O3Anomaly(i).ppmv,highClTimperiod(2) - highClTimperiod(1));
    else
         %[bhclind(i).b,ballhclind(i).b] = ozoneRegressionEnsAve(O3Anomaly(i).percent_residuals_months,highClTimperiod(2) - highClTimperiod(1));
         [bhclind(i).b,ballhclind(i).b] = ozoneRegressionEnsAve(O3Anomaly(i).percent,highClTimperiod(2) - highClTimperiod(1));
    end
end

%%
bcat = cat(4,ballhclind(:).b);
bindtoplot = permute(squeeze(bcat(:,:,2,:)*120),[3,2,1]);
%btoplot(2,:,:) = balltest2(:,:,2)*120;
%btoplot = permute(btoplot,[1,3,2]);

bstd = std(bindtoplot,1,1);

prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

titles = {'No. 1','No. 2','No. 3','No. 4','No. 5','No. 6','No. 7','No. 8','No. 9'};

mtit = {['Individual ensemble member global ',vartitle,' linear trends over ', num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))]};

if strcmp(var,'T') 
    clim = [-2 2];
    ctitle = 'K/decade';
elseif strcmp(var,'Z3') 
    clim = [-100 100];
    ctitle = 'm/decade';
elseif strcmp(var,'U')   
    clim = [-2 2];
    ctitle = 'm/s/decade';
else
    clim = [-12 12];
    ctitle = 'Percent/decade';
end


%bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
% bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
%bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
%bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
subplotmaps(bindtoplot,WAClat,log(highCllevel),{'div','RdBu'},1,[],16,titles,'Latitude','Pressure (hPa)',ctitle,'on',...
    clim,22,-90:30:90,-90:30:90,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(10) log(300)],1,'-',0,'');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','Ensave_indmem_',var,'_',num2str(highClTimperiod(1)),'-',num2str(highClTimperiod(2))];
set(gcf,'Renderer','Painters')

print(filename,'-depsc');
export_fig(filename,'-png');

%%
subplotmaps(bstd,WAClat,log(highCllevel),{'div','RdBu'},1,[],20,{['Ens-ave standard deviation in linear trends ',num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))]},'Latitude','Pressure (hPa)','Percent/decade','on',...
    [-4 4],22,-90:30:90,-90:30:90,...
    fliplr(logprestick),fliplr(presticklabel),{''},1,[-90 90],[log(1) log(200)],1,'-',0,'');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','Ensave_indmem_std',num2str(highClTimperiod(1)),'-',num2str(highClTimperiod(2))];

set(gcf,'position',[100 100 900 700])

%export_fig(filename,'-png');
export_fig(filename,'-pdf');



%%
clearvars watemp bmonth
%lats = [20 60];
lats = [-80 -50;50 80];
%lats = [80 90];


%bonth = [month,pressure,lat,beta]
for i = 1:12
    for j = 1:size(O3Anomaly(11).percent_residuals_months,2)
        for k = 1:2
            latind = WAClat >= lats(k,1) & WAClat <= lats(k,2);
            watemp(:,i,j,k) = weightedaverage(permute(squeeze(O3Anomaly(11).percent_residuals_months(i:12:end,j,latind)),[2,1]),WAClat(latind));
            [bmonth(i,j,k,:),~,~,~,bmonthstats] = regress(watemp(:,i,j,k),[ones(size(watemp,1),1),[1:size(watemp,1)]']);
            bmonth_pvalue(1,i,j,k) = bmonthstats(3);
        end
    end
    
end

plotHIGHCL = 0;
if plotHIGHCL
    %% testing with contour plot
    bmonthtoplot = permute(bmonth(:,:,2),[3,1,2])*10;
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {['Monthly ozone trends between ',num2str(lats(1)),' and ',num2str(lats(2)),'{\circ}','N']};

    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(bmonthtoplot,1:12,log(highCllevel),{'div','RdBu'},1,[],22,titles,'Month','Pressure (hPa)','Percent/decade','on',...
        [-4 4],22,1:12,1:12,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[1 12],[log(.1) log(30)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/',...
        'Highclmonthregression_',num2str(highClTimperiod(1)),'-',num2str(highClTimperiod(2)),num2str(lats(1)),'and',num2str(lats(2)),'N'];

    export_fig(filename,'-png');

    %%
    lat = 10;
    lev = 27;
    mon = 9;
    figure
    for i = 1:11
        hold on    
        if i == 11
            plot(O3Anomaly(i).percent_residuals_months(:,lev,lat),'LineWidth',3)
        else
            plot(O3Anomaly(i).percent_residuals_months(:,lev,lat))
        end
    end

    figure
    for i = 1:11
        hold on
        if i == 11
            plot(O3Anomaly(i).percent_residuals(:,lev,lat),'LineWidth',3)
        else
            plot(O3Anomaly(i).percent_residuals(:,lev,lat))
        end
    end

    figure
    plot(O3Anomaly(11).percent_residuals_months(mon:12:end,lev,lat),'LineWidth',3)
    hold on
    plot(O3Anomaly(11).percent_residuals(mon:12:end,lev,lat),'--','LineWidth',3)

    %plot(O3Anomaly(1).residuals(:,27,80))

    %% testing with contour plot
    clearvars btoplot
    btoplot(1,:,:) = balltest(:,:,2)*120;
    btoplot(2,:,:) = balltest2(:,:,2)*120;
    btoplot = permute(btoplot,[1,3,2]);

    prestick = [300,200,100,90:-10:10,9:-1:1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1};
    logprestick = log(prestick);

    titles = {'Ensemble average 1998-2024 yearly ozone trends','2'};

    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(btoplot,WAClat,log(highCllevel),{'div','RdBu'},1,[],22,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
        [-4 4],18,-90:10:90,-90:10:90,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[-90 90],[log(.1) log(300)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','Highclmonthregression_',num2str(highClTimperiod(1)),'-',num2str(highClTimperiod(2))];

    export_fig(filename,'-png');
end


%% import Chem-only and MAM
SDtimeperiod = [2000 2014];
[SDWaccmData] = ReadInSDWACCM(highCllevel,SDtimeperiod);

%% take regression of SD WACCM (chem-only and MAM)

[bchemonly,predictorschemonly,O3Anomalychemonly] = ...
    ozoneRegressionTrends(SDWaccmData(2).O3,SDWaccmData(2).U,SDWaccmData(2).solar,SDWaccmData(2).SPE,...
    SDWaccmData(2).NINO34,SDWaccmData(2).HFS,SDWaccmData(2).HFN,SDWaccmData(2).NO2,0,SDtimeperiod,1,0,0);%SDWaccmData(2).SPEintpres 



%% Take regression of residuals
[btest,ballchemonly,ballchem_pvalue] = ozoneRegressionEnsAve(O3Anomalychemonly.percent_residuals_months,SDtimeperiod(2) - SDtimeperiod(1));
[btest2,ballchemonly2] = ozoneRegressionEnsAve(O3Anomalychemonly.percent_residuals,SDtimeperiod(2) - SDtimeperiod(1));

sig = .01;
ballchem_pvalue_toplot = ballchem_pvalue;
ballchem_pvalue_toplot (ballchem_pvalue_toplot <= sig) = 0;
ballchem_pvalue_toplot (ballchem_pvalue_toplot > sig) = 1;

 %%
    clearvars watemp_chemonly bmonth_chemonly watemp_chemonly2 bmonth_chemonly2
    %lats = [70 90];
    lats = [-80 -50; 50 80];
latsforpresmonth = lats;


    for i = 1:12
        for j = 1:size(O3Anomalychemonly.percent_residuals_months,2)
            for k = 1:2
                latind = WAClat >= lats(k,1) & WAClat <= lats(k,2);
            
            watemp_chemonly(:,i,j,k) = weightedaverage(permute(squeeze(O3Anomalychemonly.percent_residuals_months(i:12:end,j,latind)),[2,1]),WAClat(latind));
            watemp_chemonly2(:,i,j,k) = weightedaverage(permute(squeeze(O3Anomalychemonly.percent(i:12:end,j,latind)),[2,1]),WAClat(latind));
            [bmonth_chemonly(i,j,k,:),~,~,~,bmonthstatsc1] = regress(watemp_chemonly(:,i,j,k),[ones(size(watemp_chemonly,1),1),[1:size(watemp_chemonly,1)]']);
            [bmonth_chemonly2(i,j,k,:),~,~,~,bmonthstatsc2] = regress(watemp_chemonly2(:,i,j,k),[ones(size(watemp_chemonly2,1),1),[1:size(watemp_chemonly2,1)]']);
            
            bmonth_pvalue(2,i,j,k) = bmonthstatsc1(3);
            bmonth_pvalue(3,i,j,k) = bmonthstatsc2(3);
            end
        end

    end




%%

% figure
% plot(O3Anomalychemonly.percent_regmodel(:,27,95))
% hold on
% plot(O3Anomalychemonly.percent(:,27,95));
% 
% figure
% plot(O3Anomalychemonly.percent_regmodel_months(:,27,95))
% hold on
% plot(O3Anomalychemonly.percent(:,27,95));
% 
% figure
% plot(O3Anomalychemonly.percent_residuals_months(12:12:end,20,92));
% hold on
% plot(O3Anomalychemonly.percent_regmodel_months(12:12:end,20,92))
% hold on
% plot(O3Anomalychemonly.percent(12:12:end,20,92));

plotCHEMONLY = 0;
if plotCHEMONLY
    %%
    mon = 10;
    lev = 27;
    latt = 85;
    figure
    plot(O3Anomalychemonly.percent_residuals_months(mon:12:end,lev,latt));
    hold on
    %plot(O3Anomalychemonly.percent_residuals(mon:12:end,lev,latt));

    %plot(O3Anomalychemonly.percent_residuals(mon:12:end,lev,latt));
    plot(O3Anomalychemonly.percent_regmodel_months(mon:12:end,lev,latt))
    plot(O3Anomalychemonly.percent_regmodel(mon:12:end,lev,latt))
    hold on
    plot(O3Anomalychemonly.percent(mon:12:end,lev,latt));

    lh = legend('res','modelmonths','model','data');
    set(lh,'box','off');


   
    %% testing with contour plot
    clearvars bmonthtoplot

    bmonthtoplot(1,:,:) = permute(bmonth_chemonly(:,:,1,2),[3,1,2])*10;
    %bmonthtoplot(2,:,:) = permute(bmonth2(:,:,2),[3,1,2])*10;
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {['Monthly ozone trends between ', num2str(lats(1)),' and ',num2str(lats(2)),'{\circ}','N']};

    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(bmonthtoplot,1:12,log(highCllevel),{'div','RdBu'},1,[],22,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
        [-4 4],22,1:12,1:12,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[1 12],[log(.1) log(30)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/',...
        'Chemonly_monthregression_nosolarSPE_5080N_2_',num2str(SDtimeperiod(1)),'-',num2str(SDtimeperiod(2)),num2str(lats(1)),'and',num2str(lats(2)),'N'];

    %export_fig(filename,'-png');

    %% testing with contour plot
    clearvars btoplot
    btoplot(1,:,:) = ballchemonly(:,:,2)*120;
    btoplot(2,:,:) = ballchemonly2(:,:,2)*120;
    btoplot = permute(btoplot,[1,3,2]);

    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {'Chem-only 2000-2014 yearly ozone trends','2'};

    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(btoplot,WAClat,log(highCllevel),{'div','RdBu'},1,[],22,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
        [-4 4],22,-90:10:90,-90:10:90,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[-90 90],[log(.1) log(30)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','Chemonlymonthregression2_',num2str(SDtimeperiod(1)),'-',num2str(SDtimeperiod(2))];

   % export_fig(filename,'-png');
end

%% MAM

[bMAM,predictorsMAM,O3AnomalyMAM] = ...
    ozoneRegressionTrends(SDWaccmData(1).O3,SDWaccmData(1).U,SDWaccmData(2).solar,...
    0,SDWaccmData(1).NINO34, SDWaccmData(1).HFS,SDWaccmData(1).HFN,SDWaccmData(1).NO2*1e9,0,SDtimeperiod,0,0,0);%SDWaccmData(2).SPEintpres, SDWaccmData(1).HFS,SDWaccmData(1).HFN,SDtimeperiod  SDWaccmData(1).NINO34

% [bMAM,predictorsMAM,O3AnomalyMAM] = ...
%     ozoneRegressionTrends(SDWaccmData(1).O3,0,SDWaccmData(2).solar,...%SDWaccmData(2).solar
%     0,0,0,0,SDWaccmData(1).NO2*1e9,SDtimeperiod,0,0);%SDWaccmData(2).SPEintpres, SDWaccmData(1).HFS,SDWaccmData(1).HFN,SDtimeperiod 


[btest,ballMAM] = ozoneRegressionEnsAve(O3AnomalyMAM.percent_residuals_months,SDtimeperiod(2) - SDtimeperiod(1));
[btest2,ballMAM2] = ozoneRegressionEnsAve(O3AnomalyMAM.percent,SDtimeperiod(2) - SDtimeperiod(1));

%%

clearvars watemp bmonth
%lats = [50 80];
lats = [-90 -80];
%lats = [80 90];

latind = WAClat >= lats(1) & WAClat <= lats(2);

for i = 1:12
    for j = 1:size(O3AnomalyMAM.percent_residuals_months,2)
        watemp(:,i,j) = weightedaverage(permute(squeeze(O3AnomalyMAM.percent_residuals_months(i:12:end,j,latind)),[2,1]),WAClat(latind));
        watemp2(:,i,j) = weightedaverage(permute(squeeze(O3AnomalyMAM.percent(i:12:end,j,latind)),[2,1]),WAClat(latind));
        bmonth(i,j,:) = regress(watemp(:,i,j),[ones(size(watemp,1),1),[1:size(watemp,1)]']);
        bmonth2(i,j,:) = regress(watemp2(:,i,j),[ones(size(watemp2,1),1),[1:size(watemp2,1)]']);
    end
    
end

plotMAM = 0;
if plotMAM
    %% 
    mon = 6;
    lev = 25;
    latt = 2;
    figure
    plot(O3AnomalyMAM.percent_residuals_months(mon:12:end,lev,latt),'-o'); 
    hold on
    %plot(O3AnomalyMAM.percent_residuals(mon:12:end,lev,latt));
    plot(O3AnomalyMAM.percent_regmodel_months(mon:12:end,lev,latt))    
    %plot(O3AnomalyMAM.percent_regmodel(mon:12:end,lev,latt))    
    plot(O3AnomalyMAM.percent(mon:12:end,lev,latt));

    %legend('res months','res','model months','model','data')
    legend('res months','model months','data')

    figure
    plot(squeeze(SDWaccmData(1).NO2(lev,latt,mon:12:end)),'-o');
    hold on
    plot(squeeze(SDWaccmData(2).NO2(lev,latt,mon:12:end)));
    
    %% testing with contour plot
    clearvars btoplot
    btoplot(1,:,:) = ballMAM(:,:,2)*120;
    %btoplot(2,:,:) = ballMAM2(:,:,2)*120;
    btoplot = permute(btoplot,[1,3,2]);

    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {'Ensemble average 2000-2014 yearly ozone trends','2'};

    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(btoplot,WAClat,log(highCllevel),{'div','RdBu'},1,[],22,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
        [-4 4],22,-90:10:90,-90:10:90,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[-90 90],[log(.1) log(300)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','MAMmonthregression2_',num2str(SDtimeperiod(1)),'-',num2str(SDtimeperiod(2))];

    export_fig(filename,'-png');


    %%

    bmonthtoplot(1,:,:) = permute(bmonth(:,:,2),[3,1,2])*10;
    %bmonthtoplot(2,:,:) = permute(bmonth2(:,:,2),[3,1,2])*10;
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {['Monthly ozone trends between ', num2str(lats(1)),' and ',num2str(lats(2)),'{\circ}','N']};

    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(bmonthtoplot,1:12,log(highCllevel),{'div','RdBu'},1,[],22,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
        [-4 4],22,1:12,1:12,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[1 12],[log(.1) log(30)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/',...
        'MAMmonthregression_',num2str(SDtimeperiod(1)),'-',num2str(SDtimeperiod(2)),num2str(lats(1)),'and',num2str(lats(2)),'N'];

    export_fig(filename,'-png');
end


%% SWOOSH

% [bSWOOSH,predictorsSWOOSH,O3AnomalySWOOSH] = ...
%     ozoneRegressionTrends(SWOOSHextract,SDWaccmData(1).U,SDWaccmData(2).solar,...
%     0,SDWaccmData(1).NINO34,SDWaccmData(1).HFS,SDWaccmData(1).HFN,SDWaccmData(1).NO2*1e9,SDtimeperiod,0,0,0);%SDWaccmData(2).SPEintpres,, SDWaccmData(1).NO2*1e9

[bSWOOSH,predictorsSWOOSH,O3AnomalySWOOSH] = ...
    ozoneRegressionTrends(SWOOSHextract,SWOOSHregfun.singaporedata,...
    SWOOSHregfun.solardata,0,SWOOSHregfun.MEIdata,0,0,0,0,Stimeperiod,0,0,0);%SDWaccmData(2).SPEintpres,, SDWaccmData(1).NO2*1e9

% [bSWOOSH,predictorsSWOOSH,O3AnomalySWOOSH] = ...
%     ozoneRegressionTrends(SWOOSHextract,SDWaccmData(1).U,SDWaccmData(2).solar,...
%     0,SDWaccmData(1).NINO34,SDWaccmData(1).HFS,SDWaccmData(1).HFN,0,SDtimeperiod,0,0,0);%SDWaccmData(2).SPEintpres,, SDWaccmData(1).NO2*1e9

% [bSWOOSH,predictorsSWOOSH,O3AnomalySWOOSH] = ...
%     ozoneRegressionTrends(SWOOSHextract,0,0,...
%     0,0,0,0,0,Stimeperiod,0,1,0);%SDWaccmData(2).SPEintpres,, SDWaccmData(1).NO2*1e9

%%

[btest,ballSWOOSH] = ozoneRegressionEnsAve(O3AnomalySWOOSH.percent_residuals_months,Stimeperiod(2) - Stimeperiod(1));
[btest2,ball2SWOOSH] = ozoneRegressionEnsAve(O3AnomalySWOOSH.percent,Stimeperiod(2) - Stimeperiod(1));
[btest3,ball3SWOOSH] = ozoneRegressionEnsAve(O3AnomalySWOOSH.percent_residuals,Stimeperiod(2) - Stimeperiod(1));

clearvars watemp bmonth
lats = [10 40];
%lats = [80 90];

latind = SWOOSH.lat >= lats(1) & SWOOSH.lat <= lats(2);

for i = 1:12
    for j = 1:size(O3AnomalySWOOSH.percent_residuals_months,2)
        watemp(:,i,j) = weightedaverage(permute(squeeze(O3AnomalySWOOSH.percent_residuals_months(i:12:end,j,latind)),[2,1]),WAClat(latind));
        watemp2(:,i,j) = weightedaverage(permute(squeeze(O3AnomalySWOOSH.percent_residuals(i:12:end,j,latind)),[2,1]),WAClat(latind));
        bmonth(i,j,:) = regress(watemp(:,i,j),[ones(size(watemp,1),1),[1:size(watemp,1)]']);
        bmonth2(i,j,:) = regress(watemp2(:,i,j),[ones(size(watemp2,1),1),[1:size(watemp2,1)]']);
    end
    
end

%% 
% mon = 12;
% lev = 27;
% latt = 83;
% figure
% plot(O3AnomalyMAM.percent_residuals_months(mon:12:end,lev,latt));
% hold on
% plot(O3AnomalyMAM.percent_regmodel_months(mon:12:end,lev,latt))
% hold on
% plot(O3AnomalyMAM.percent(mon:12:end,lev,latt));
% 
% legend('res','model','data')

plotSWOOSH = 0;
if plotSWOOSH

    %% testing with contour plot
    clearvars btoplot
    btoplot(1,:,:) = ballSWOOSH(:,:,2)*120;
    %btoplot(2,:,:) = ball2SWOOSH(:,:,2)*120;
    btoplot(2,:,:) = ball3SWOOSH(:,:,2)*120;
    btoplot = permute(btoplot,[1,3,2]);

    btoplot(:,[1,2,end-1,end],:) = NaN;
    
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {['SWOOSH ', num2str(Stimeperiod(1)), char(8211), num2str(Stimeperiod(2)), ' yearly ozone trends (full regression)'],...
        ['SWOOSH ', num2str(Stimeperiod(1)), char(8211), num2str(Stimeperiod(2)), ' yearly ozone trends (direct linear trend)']};%'SWOOSH 1984-1999 yearly ozone trends (full regression)',

    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(btoplot,SWOOSH.lat,log(SWOOSH.level),{'div','RdBu'},1,[],20,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
        [-8 8],22,-90:10:90,-90:10:90,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[-90 90],[log(.1) log(200)],1,'-',0,'');
    %set(gcf,'position',[100 100 900 700])
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','SWOOSHmonthregression2_',num2str(Stimeperiod(1)),'-',num2str(Stimeperiod(2))];

    export_fig(filename,'-pdf');



    %%
    clearvars bmonthtoplot
    bmonthtoplot(1,:,:) = permute(bmonth(:,:,2),[3,1,2])*10;
    %bmonthtoplot(2,:,:) = permute(bmonth2(:,:,2),[3,1,2])*10;
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {['Monthly ozone trends between ', num2str(lats(1)),' and ',num2str(lats(2)),'{\circ}','N']};

    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(bmonthtoplot,1:12,log(SWOOSH.level),{'div','RdBu'},1,[],22,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
        [-4 4],22,1:12,1:12,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[1 12],[log(.1) log(30)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/',...
        'SWOOSHmonthregression_',num2str(SDtimeperiod(1)),'-',num2str(SDtimeperiod(2)),num2str(lats(1)),'and',num2str(lats(2)),'N'];

    export_fig(filename,'-png');
end


%% create line plot
createlineplot = 1;
if createlineplot
    lats = [-69 -63;63 67];

    lev = [22,27];
    mon = [10,12];
    
    titles = {[num2str(abs(lats(1,2))),char(8211),num2str(abs(lats(1,1))),'{\circ}S, ',sprintf('%.1f',SWOOSH.level(lev(1))),' hPa, ',monthnames(mon(1),0,0)],...
        [num2str(abs(lats(2,1))),char(8211),num2str(abs(lats(2,2))),'{\circ}N, ',sprintf('%.1f',SWOOSH.level(lev(2))),' hPa, ',monthnames(mon(2),0,0)]};
    fsize = 18;
    lwidth = 3;
    colors = cbrewer('qual','Set1',10);
    plot_bothhem = 0;
    
    figure;
    set(gcf,'position',[100 100 600 900],'color','white');
    if plot_bothhem    
        for k = 1:2

            sp(i) = subplot(2,1,k);
            latind = WAClat >= lats(k,1) & WAClat <= lats(k,2);

            if k == 2
                sppos = get(sp(i),'position');
                set(sp(i),'position',[sppos(1),sppos(2)+.05,sppos(3:4)]);
            end

            NO2temp = weightedaverage(squeeze(SDWaccmData(2).NO2(lev(k),latind,mon(k):12:end)),WAClat(latind))+1e-8;

            NO2anompercent = (NO2temp - nanmean(NO2temp))./nanmean(NO2temp)*100/5;

            latindSWOOSH = SWOOSH.lat >= lats(k,1) & SWOOSH.lat <= lats(k,2);
            SWOOSHlatmean = weightedaverage(squeeze(O3AnomalySWOOSH.percent(mon(k):12:end,lev(k),latindSWOOSH))',SWOOSH.lat(latindSWOOSH));
            SWOOSHlatmean (SWOOSHlatmean == 0) = NaN;
            Chemonlylatmean = weightedaverage(squeeze(O3Anomalychemonly.percent(mon(k):12:end,lev(k),latind))',WAClat(latind));
            Chemonlylatmean_regress = weightedaverage(squeeze(O3Anomalychemonly.percent_residuals_months(mon(k):12:end,lev(k),latind))',WAClat(latind));
            MAMlatmean = weightedaverage(squeeze(O3AnomalyMAM.percent(mon(k):12:end,lev(k),latind))',WAClat(latind));


            plot(Chemonlylatmean,'--','LineWidth',lwidth,'color',colors(3,:))
            hold on
            plot(Chemonlylatmean_regress,'-','LineWidth',lwidth,'color','k')                     
            plot(NO2anompercent,':','LineWidth',lwidth,'color',[.5 .5 .5]);
            plot(MAMlatmean,'LineWidth',lwidth,'color',colors(2,:))        
            plot(SWOOSHlatmean,'LineWidth',lwidth,'color',colors(1,:))        
            set(gca,'fontsize',fsize+2,'xtick',1:2:15,'xticklabel',2000:2:2015)
            if k == 2
                xlabel('Year','fontsize',fsize+2);
            else
                set(gca,'ytick',-10:2:6,'yticklabel',-10:2:6);
            end
            ylabel('Anomaly (%)','fontsize',fsize+2);
            title(titles{k},'fontsize',fsize+4)
            %title('Simulated and observed solar proton effects on ozone')
            xlim([0 16])
            if k == 1
                lh = legend('O_3 Chem-only (WSE)','O_3 Chem-only','NO_2 Chem-only (%/20)','O_3 Chem-Dyn-Vol','O_3 SWOOSH');
                set(lh,'box','off','fontsize',fsize-2,'location','southeast');
            end

        end

        annotation('textbox',[0 .905 1 .1],'String','Simulated and observed solar proton effects','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+4,...
            'EdgeColor','none','fontweight','bold')
    

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/BothHem_LinePlots_select'];
    else
    %%
    sp(1) = subplot(2,1,1);
    latind = WAClat >= lats(2,1) & WAClat <= lats(2,2);

    NO2temp = weightedaverage(squeeze(SDWaccmData(2).NO2(lev(k),latind,mon(k):12:end)),WAClat(latind))+1e-8;
    NO2anompercent = (NO2temp - nanmean(NO2temp))./nanmean(NO2temp)*100/8;
    latindSWOOSH = SWOOSH.lat >= lats(k,1) & SWOOSH.lat <= lats(k,2);
    SWOOSHlatmean = weightedaverage(squeeze(O3AnomalySWOOSH.percent(mon(k):12:end,lev(k),latindSWOOSH))',SWOOSH.lat(latindSWOOSH));
    SWOOSHlatmean (SWOOSHlatmean == 0) = NaN;
    Chemonlylatmean = weightedaverage(squeeze(O3Anomalychemonly.percent(mon(k):12:end,lev(k),latind))',WAClat(latind));
    Chemonlylatmean_regress = weightedaverage(squeeze(O3Anomalychemonly.percent_residuals_months(mon(k):12:end,lev(k),latind))',WAClat(latind));
    MAMlatmean = weightedaverage(squeeze(O3AnomalyMAM.percent(mon(k):12:end,lev(k),latind))',WAClat(latind));

    plot(Chemonlylatmean,'--','LineWidth',lwidth,'color',colors(3,:))
    hold on
    plot(MAMlatmean,'LineWidth',lwidth,'color',colors(2,:))        
    plot(SWOOSHlatmean,'LineWidth',lwidth,'color',colors(1,:))        
    set(gca,'fontsize',fsize+2,'xtick',1:2:15,'xticklabel',2000:2:2015)
    ylabel('Ozone Anomaly (%)','fontsize',fsize+2);
    title(titles{2},'fontsize',fsize+4)
    xlim([0 16]);
    ylim([-20 15]);
    lh = legend('O_3 Chem-only (WSE)','O_3 Chem-Dyn-Vol','O_3 SWOOSH');
    set(lh,'box','off','fontsize',fsize-2,'location','southeast');
    sp(2) = subplot(2,1,2);    
    
    plot(Chemonlylatmean,'--','LineWidth',lwidth,'color',colors(3,:))
    hold on
    plot(Chemonlylatmean_regress,'-','LineWidth',lwidth,'color','k')                     
    plot(NO2anompercent,':','LineWidth',lwidth,'color',[.5 .5 .5]);
    set(gca,'fontsize',fsize+2,'xtick',1:2:15,'xticklabel',2000:2:2015)
    xlabel('Year','fontsize',fsize+2);
    xlim([0 16]);
    ylim([-15 15]);
    ylabel('Ozone Anomaly (%)','fontsize',fsize+2);
    %title(titles{2},'fontsize',fsize+4)
    sppos = get(sp(2),'position');
    set(sp(2),'position',[sppos(1),sppos(2)+.07,sppos(3:4)]);
    lh = legend('O_3 Chem-only (WSE)','O_3 Chem-only','NO_2 Chem-only (%/8)');
    set(lh,'box','off','fontsize',fsize-2,'location','southeast');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/LinePlots_select'];
    
    annotation('textbox',[0 .905 1 .1],'String','Simulated and observed solar proton effects','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+4,...
            'EdgeColor','none','fontweight','bold')
    
    end
end

export_fig(filename,'-png');
export_fig(filename,'-pdf');

%%

createPaperFigures = 1;
if createPaperFigures
    clearvars Globaltrends
    Globaltrends(1,:,:) = balltest2(:,:,2)*120;
    Globaltrends(2,:,:) = ballchemonly(:,:,2)*120;
    Globaltrends = permute(Globaltrends,[1,3,2]);
    clearvars btoplot            
    
    ptoplot = permute(cat(1,ball_pvalue_toplot,ballchem_pvalue_toplot),[1,3,2]);
    
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);
    
    mtit = {'Global ozone linear trends'};
    
    titles = {['Ensemble average 1998' char(8211) '2024 ozone trends'],['Chem-only 2000' char(8211) '2014 ozone trends']};
    
    subplotmaps(Globaltrends,WAClat,log(highCllevel),{'div','RdBu'},1,[],22,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
        [-4 4],22,-90:10:90,-90:10:90,...
        fliplr(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.1) log(30)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','Globaltrends_',num2str(highClTimperiod(1)),'-',num2str(highClTimperiod(2))];

    export_fig(filename,'-png');
    export_fig(filename,'-pdf');

    % plot month/pressure near poles    
    
    %% testing with contour plot
    clearvars Polesig_toplot
    
    Poletrends(1,:,:) = permute(bmonth(:,:,1,2),[3,1,2])*10;
    Poletrends(2,:,:) = circshift(permute(bmonth(:,:,2,2),[3,1,2])*10,[0,6,0]);
    Poletrends(3,:,:) = permute(bmonth_chemonly(:,:,1,2),[3,1,2])*10;
    Poletrends(4,:,:) = circshift(permute(bmonth_chemonly(:,:,2,2),[3,1,2])*10,[0,6,0]);
    Poletrends(5,:,:) = permute(bmonth_chemonly2(:,:,1,2),[3,1,2])*10;
    Poletrends(6,:,:) = circshift(permute(bmonth_chemonly2(:,:,2,2),[3,1,2])*10,[0,6,0]);
    
    Polesig(1,:,:) = bmonth_pvalue(1,:,:,1);
    Polesig(2,:,:) = circshift(bmonth_pvalue(1,:,:,2),[0,6,0]);
    Polesig(3,:,:) = bmonth_pvalue(2,:,:,1);
    Polesig(4,:,:) = circshift(bmonth_pvalue(2,:,:,2),[0,6,0]);
    Polesig(5,:,:) = bmonth_pvalue(3,:,:,1);
    Polesig(6,:,:) = circshift(bmonth_pvalue(3,:,:,2),[0,6,0]);
    
    sig = .01;
    Polesig_toplot = Polesig;
    Polesig_toplot (Polesig_toplot <= sig) = 0;
    Polesig_toplot (Polesig_toplot > sig) = 1;
    
    bmonthtoplot = permute(bmonth(:,:,2),[3,1,2])*10;
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    %titles = {['Monthly ozone trends between ',num2str(lats(1)),' and ',num2str(lats(2)),'{\circ}','N']};
    titles = {['Ens-ave, 1998',char(8211),'2024, ',num2str(abs(latsforpresmonth(1,2))),char(8211),num2str(abs(latsforpresmonth(1,1))),'{\circ}','S']...
        ,['Ens-ave, 1998',char(8211),'2024, ',num2str(abs(latsforpresmonth(2,1))),char(8211),num2str(abs(latsforpresmonth(2,2))),'{\circ}','N']...
        ,['Chem-only, 2000',char(8211),'2014, ',num2str(abs(latsforpresmonth(1,2))),char(8211),num2str(abs(latsforpresmonth(1,1))),'{\circ}','S']...
        ,['Chem-only, 2000',char(8211),'2014, ',num2str(abs(latsforpresmonth(2,1))),char(8211),num2str(abs(latsforpresmonth(2,2))),'{\circ}','N']...
        ,['Chem-only (WSE), 2000',char(8211),'2014, ',num2str(abs(latsforpresmonth(1,2))),char(8211),num2str(abs(latsforpresmonth(1,1))),'{\circ}','S']...
        ,['Chem-only (WSE), 2000',char(8211),'2014, ',num2str(abs(latsforpresmonth(2,1))),char(8211),num2str(abs(latsforpresmonth(2,2))),'{\circ}','N']};

    xticks = {'J','F','M','A','M','J','J','A','S','O','N','D';...
        'J','A','S','O','N','D','J','F','M','A','M','J';...
        'J','F','M','A','M','J','J','A','S','O','N','D';...
        'J','A','S','O','N','D','J','F','M','A','M','J';...
        'J','F','M','A','M','J','J','A','S','O','N','D';...
        'J','A','S','O','N','D','J','F','M','A','M','J'};
    
    mtit = {'Monthly ozone linear trends'}; 
    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(Poletrends,1:12,log(highCllevel),{'div','RdBu'},1,[],17,titles,'Month','Pressure (hPa)','Percent/decade','on',...
        [-4 4],22,1:12,xticks,...
        fliplr(logprestick),fliplr(presticklabel),mtit ,1,[1 12],[log(.1) log(30)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/',...
        'MonthlyTrends_',num2str(highClTimperiod(1)),'-',num2str(highClTimperiod(2)),num2str(latsforpresmonth(1)),'and',num2str(latsforpresmonth(2)),'N'];

    export_fig(filename,'-png');
    export_fig(filename,'-pdf');
    print(filename,'-depsc');
    
end