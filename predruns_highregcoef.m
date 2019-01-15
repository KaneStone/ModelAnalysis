% plot high correlation areas and take predictions
clear all
remove_ENSO = 0;
SHreg = load('/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefsNovembertoz_TS_detrend1995-2016_90-75S_Tperiod-DecJanFeb_.mat');
NHreg = load('/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefsMarchtoz_TS_detrend1995-2016_75-90S_Tperiod-MarApr_rmENSO.mat');

%% Finding maxmimum and minimum correlation spots
% highcl case only
% lonland = [110, 150; 165, 180; 15, 35; 280,320]; % Australia, NZ, Africa, South America
% latland = [-45, -32; -47, -33; -38, -28; -55,0];

lonland = [110, 150; 15, 35]; % Australia, NZ, Africa, South America
latland = [-45,-32; -38, -28];

noval = 1;
for i = 1:size(SHreg.r.ind.highcl)
    for j = 1:size(lonland,1)
        lonind = SHreg.lons >= lonland(j,1) & SHreg.lons <= lonland(j,2);        
        latind = SHreg.lats >= latland(j,1) & SHreg.lats <= latland(j,2);        
        highcl_reg = squeeze(SHreg.r.ind.highcl(i,lonind,latind));        
        [srt,srtind] = sort(highcl_reg(:));        
        maxval(i,j) = srt(end);
        minval(i,j) = srt(1);
        [maxcoor(i,j,:,1),maxcoor(i,j,:,2)] = ind2sub(size(highcl_reg),srtind(end-noval+1:end));
        [mincoor(i,j,:,1),mincoor(i,j,:,2)] = ind2sub(size(highcl_reg),srtind(1:noval));
        
        latextract = SHreg.lats(latind);
        lonextract = SHreg.lons(lonind);        
        latsofmin(i,j) = latextract(squeeze(mincoor(i,j,:,2)));
        lonsofmin(i,j) = lonextract(squeeze(mincoor(i,j,:,1)));        
        latsofmax(i,j) = latextract(squeeze(maxcoor(i,j,:,2)));
        lonsofmax(i,j) = lonextract(squeeze(maxcoor(i,j,:,1)));
    end
end
lonsofmax2 = lonsofmax;
lonsofmin2 = lonsofmin;
lonsofmin2 (lonsofmin2 > 180) = (lonsofmin2 (lonsofmin2 > 180)) - 360;
lonsofmax2 (lonsofmax2 > 180) = (lonsofmax2 (lonsofmax2 > 180)) - 360;

%% Read in highcl Surface temperature
ClLevel = 'highCl';
timeperiodhigh = [1995,2016];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/','TS','/',ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl]...
    = predruns_ReadInlayer(vardirectory,varfiles,'TS',timeperiodhigh,[-90,-75],1);

%% Read in TOZ highcl and take percentiles
% years to estimate = 2016:2023

ClLevel = 'highCl';
tozdates = [1995,2015];
directory = ['/Volumes/ExternalOne/work/data/predruns/','TOZ','/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,'toz',tozdates,[-95,-75],0);

%% Read in highcl Surface temperature
ClLevel = 'highCl';
timeperiodhigh = [2016,2024];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/','TS','/',ClLevel,'/'];
varfiles = dir([vardirectory,'*.nc']);

[data_forpred.highcl,years_forpred.highcl,composite_forpred.highcl,dataMonthArrange_forpred.highcl]...
    = predruns_ReadInlayer(vardirectory,varfiles,'TS',timeperiodhigh,[-90,-75],1);

temp_forpred = nanmean(cat(5,squeeze(dataMonthArrange_forpred.highcl(:,12,1:end-1,:,:)),squeeze(dataMonthArrange_forpred.highcl(:,1,2:end,:,:)),...
    squeeze(dataMonthArrange_forpred.highcl(:,2,2:end,:,:))),5);

%% Read in TOZ highcl and take percentiles
% years to estimate = 2016:2023

ClLevel = 'highCl';
tozdates = [2016,2023];
directory = ['/Volumes/ExternalOne/work/data/predruns/','TOZ','/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data_forpred.highcl,toz_years_forpred.highcl,toz_varweighted_forpred.highcl,toz_composite_forpred.highcl,toz_dataMonthArrange_forpred.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,'toz',tozdates,[-95,-75],1);

%%
    varmonth = [12,1,2];
for j = 1:length(varmonth)
    %varmonths
    if varmonth(j) <= 2
        dataMonthArrange_varmonth(:,j,:,:,:) = squeeze(dataMonthArrange.highcl(:,varmonth(j),2:end,:,:));            
        dataMonthArrange_varmonth_forpred(:,j,:,:,:) = squeeze(dataMonthArrange_forpred.highcl(:,varmonth(j),2:end,:,:));            
    else
        dataMonthArrange_varmonth(:,j,:,:,:) = squeeze(dataMonthArrange.highcl(:,varmonth(j),1:end-1,:,:));            
        dataMonthArrange_varmonth_forpred(:,j,:,:,:) = squeeze(dataMonthArrange_forpred.highcl(:,varmonth(j),1:end-1,:,:));
    end

end

%% ENSO
if remove_ENSO

    varmonth = [12,1,2];
    latlimits = [-5 5];
    lonlimits = [190 240];

    latindex = SHreg.lats >= latlimits(1) & SHreg.lats <= latlimits(2);
    lonindex = SHreg.lons >= lonlimits(1) & SHreg.lons <= lonlimits(2);

    for j = 1:length(varmonth)
        count = 0;
        NINO_mn = squeeze(nanmean(dataMonthArrange.highcl(:,:,:,latindex,lonindex),5));
        NINO_mn = squeeze(nanmean(NINO_mn,4));
        NINOallmean = nanmean(NINO_mn,3);
        NINOstd = std(NINO_mn,1,3);

        NINO_mn_forpred = squeeze(nanmean(dataMonthArrange_forpred.highcl(:,:,:,latindex,lonindex),5));
        NINO_mn_forpred = squeeze(nanmean(NINO_mn_forpred,4));
        NINOallmean_forpred = nanmean(NINO_mn_forpred,3);
        NINOstd_forpred = std(NINO_mn_forpred,1,3);

        for k = 1:4
            if varmonth(j) <=2
                if varmonth(j) - count <= 0
                    varmontemp = varmonth(j)-count+12;
                    NINO34(j,k,:,:) = detrend(((squeeze(NINO_mn(:,varmontemp,1:end-1)) - ...
                        NINOallmean(:,varmontemp))./NINOstd(:,varmontemp))');    
                    NINO34_forpred(j,k,:,:) = detrend(((squeeze(NINO_mn_forpred(:,varmontemp,1:end-1)) - ...
                        NINOallmean_forpred(:,varmontemp))./NINOstd_forpred(:,varmontemp))');    
                else                    
                    NINO34_forpred(j,k,:,:) = detrend(((squeeze(NINO_mn_forpred(:,varmonth(j)-count,2:end)) - ...
                        NINOallmean_forpred(:,varmonth(j)-count))./NINOstd_forpred(:,varmonth(j)-count))');    
                    NINO34(j,k,:,:) = detrend(((squeeze(NINO_mn(:,varmonth(j)-count,2:end)) - ...
                        NINOallmean(:,varmonth(j)-count))./NINOstd(:,varmonth(j)-count))');    
                end               
            else                
                NINO34_forpred(j,k,:,:) = detrend(((squeeze(NINO_mn_forpred(:,varmonth(j)-count,1:end-1)) - ...
                    NINOallmean_forpred(:,varmonth(j)-count))./NINOstd_forpred(:,varmonth(j)-count))');    
                NINO34(j,k,:,:) = detrend(((squeeze(NINO_mn(:,varmonth(j)-count,1:end-1)) - ...
                    NINOallmean(:,varmonth(j)-count))./NINOstd(:,varmonth(j)-count))');    
            end
             count = count+1;
        end
    end

%%
% NINO_mn = squeeze(nanmean(dataMonthArrange.highcl(:,11,:,latindex,lonindex),5));
% NINO_mn = squeeze(nanmean(NINO_mn,3));
% NINOallmean = nanmean(NINO_mn,2);
% NINOstd = std(NINO_mn,1,2);
% NINO34 = detrend(((NINO_mn - NINOallmean)./NINOstd)');    
% temp = nanmean(cat(5,squeeze(dataMonthArrange.highcl(:,12,1:end-1,:,:)),squeeze(dataMonthArrange.highcl(:,1,2:end,:,:)),...
%     squeeze(dataMonthArrange.highcl(:,2,2:end,:,:))),5);
% 
% % nino for pred
% NINO_mn = squeeze(nanmean(dataMonthArrange_forpred.highcl(:,11,:,latindex,lonindex),5));
% NINO_mn = squeeze(nanmean(NINO_mn,3));
% NINOallmean = nanmean(NINO_mn,2);
% NINOstd = std(NINO_mn,1,2);
% NINO34_forpred = detrend(((NINO_mn - NINOallmean)./NINOstd)');    

%% removing ENSO
% finding maximum regression lag time

    

    %%

    for j = 1:size(dataMonthArrange_varmonth,1) % members
        for k = 1:size(dataMonthArrange_varmonth,2) % months
            for l = 1:size(dataMonthArrange_varmonth,4) % latitudes
                for m = 1:size(dataMonthArrange_varmonth,5) % longitudes
                    for lag = 1:size(NINO34,2) % lag
                        [b(j,k,l,m,lag,:)] = regress(squeeze(dataMonthArrange_varmonth(j,k,:,l,m)),...
                            [ones(size(NINO34,3),1),squeeze(NINO34(k,lag,:,j))]);                        
                        [b_forpred(j,k,l,m,lag,:)] = regress(squeeze(dataMonthArrange_varmonth_forpred(j,k,:,l,m)),...
                            [ones(size(NINO34_forpred,3),1),squeeze(NINO34_forpred(k,lag,:,j))]);                        
                        % finding largest lag correlation
                    end
                    [~,llc(j,k,l,m)] = max(abs(squeeze(b(j,k,l,m,:,2))));
                    [~,llc_forpred(j,k,l,m)] = max(abs(squeeze(b_forpred(j,k,l,m,:,2))));
                    blag(j,k,l,m,:) = b(j,k,l,m,llc(j,k,l,m),:);
                    blag_forpred(j,k,l,m,:) = b_forpred(j,k,l,m,llc_forpred(j,k,l,m),:);
                    dataVarMonth(j,:,k,l,m) = squeeze(dataMonthArrange_varmonth(j,k,:,l,m)) - ...
                        squeeze(b(j,k,l,m,llc(j,k,l,m),2))*squeeze(NINO34(k,llc(j,k,l,m),:,j));                                    
                    dataVarMonth_forpred(j,:,k,l,m) = squeeze(dataMonthArrange_varmonth_forpred(j,k,:,l,m)) - ...
                        squeeze(b_forpred(j,k,l,m,llc_forpred(j,k,l,m),2))*squeeze(NINO34_forpred(k,llc_forpred(j,k,l,m),:,j));                                    
                end
            end
        end        
    end

    dataVarMonthAve = squeeze(nanmean(dataVarMonth,3));
    dataVarMonthAve_forpred = squeeze(nanmean(dataVarMonth_forpred,3));
else
    
    dataVarMonthAve = squeeze(nanmean(dataMonthArrange_varmonth,2));
    dataVarMonthAve_forpred = squeeze(nanmean(dataMonthArrange_varmonth_forpred,2));
end

% %%
% load('NINO.highcl.mat')
% for j = 1:size(temp,1)
%     for k = 1:size(temp,3)
%         for l = 1:size(temp,4)
%             [b(j,k,l,:)] = regress(squeeze(temp(j,:,k,l))',[ones(length(NINO34)-1,1),NINO34(1:end-1,j)]);                        
%             temp(j,:,k,l) = squeeze(temp(j,:,k,l))' - squeeze(b(j,k,l,2))*NINO34(1:end-1,j);
%             temp_forpred_rmENSO(j,:,k,l) = squeeze(temp_forpred(j,:,k,l))' - squeeze(b(j,k,l,2))*NINO34_forpred(1:end-1,j);            
%         end
%     end
% end      

%%
testing_enso_correlation = 0;
if testing_enso_correlation
    toznov = squeeze(toz_dataMonthArrange.highcl(:,11,:));
    for i = 1:9
        NINO_toz_corr(i) = corr(toznov(i,:)',NINO34(:,i));
    end
    figure;
    plot(NINO_toz_corr);
    figure
    plot((toznov(i,:)' - nanmean(toznov(i,:),2))./std(toznov(i,:),1,2))
    hold on
    plot(NINO34(:,i) - nanmean(NINO34(:,i))./std(NINO34(:,i)))
end

%% Use regression beta values to estimate surface temperature
regyears = [2016,2023];
yearind = find(toz_years_forpred.highcl(1).y(1:12:end) >= regyears(1) & toz_years_forpred.highcl(1).y(1:12:end) <= regyears(2)); 
for i = 1:size(toz_dataMonthArrange.highcl,1)
    for j = 1:size(latsofmax,2)
        latind = SHreg.lats == latsofmax(i,j);
        lonind = SHreg.lons == lonsofmax(i,j);
        regpred(:,i,j) = squeeze(toz_dataMonthArrange_forpred.highcl(i,11,:)).*...
            SHreg.bpred.ind.highcl(i,lonind,latind,2)+SHreg.bpred.ind.highcl(i,lonind,latind,1);
%         temperature(:,i,j) = nanmean(cat(2,squeeze(dataMonthArrange.highcl(i,12,1:end-1,latind,lonind)),...
%             squeeze(dataMonthArrange.highcl(i,1,2:end,latind,lonind)),...
%             squeeze(dataMonthArrange.highcl(i,2,2:end,latind,lonind))),2);
%         temperature_rmENSO(:,i,j) = temp(i,:,latind,lonind);
        
        %SHreg.bpred.ind.highcl(i,lonind,latind,2)
        corrtemp(i,j) = SHreg.r.ind.highcl(i,lonind,latind);
        
        %[benso(j,k,l,:)] = regress(squeeze(temp.(fields{i})(j,:,k,l))',[ones(length(NINO34.(fields{i}))-1,1),NINO34.(fields{i})(1:end-1,j)]);                        
        %            dataVarMonthAve.(fields{i})(j,:,k,l) = squeeze(temp.(fields{i})(j,:,k,l))' - squeeze(b(j,k,l,2))*NINO34.(fields{i})(1:end-1,j);  
        
    end
end

%%
plotregpred = 1;
%titles = {'Australia', 'New Zealand', 'Africa', 'South America'};
titles = {'Australia','Africa'};
titles2 = titles;
for i = 1:length(titles)
    %titles2{i} (strfind(titles,' ')) = '';
    titles2{i} (titles2{i} == ' ') = '';
end

%%
if remove_ENSO
    ENSOext = '_rmESNO';
else
    ENSOext = '';
end
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

if plotregpred
    for area = 1:2
        for i = 1:9
            latind = SHreg.lats == latsofmax(i,area);
            lonind = SHreg.lons == lonsofmax(i,area);
            fig = createfig('medium','on');
            %regpredplot = (regpred(:,i,area) - nanmean(regpred(:,i,area)))./std(regpred(:,i,area));
            %temp_forpredplot = (squeeze(temp_forpred(i,:,latind,lonind)) - nanmean(squeeze(temp_forpred(i,:,latind,lonind))))./std(squeeze(temp_forpred(i,:,latind,lonind)));
            temp_forpredplot = squeeze(dataVarMonthAve_forpred(i,:,latind,lonind));
            tempmean = nanmean(squeeze(dataVarMonthAve_forpred(i,:,latind,lonind)));
            tempstd = std(squeeze(dataVarMonthAve_forpred(i,:,latind,lonind)));
            regpredplot = (regpred(:,i,area) - nanmean(regpred(:,i,area)))./std(regpred(:,i,area))*tempstd+tempmean;
            plot(regyears(1):regyears(2),regpredplot,'LineWidth',3);
            hold on
            plot(regyears(1):regyears(2),temp_forpredplot,'LineWidth',3);
            if remove_ENSO
                title([titles{area},' No. ',num2str(i),' - ENSO removed'],'fontsize',20);
            else
                title([titles{area},' No. ',num2str(i)],'fontsize',20);
            end
            ylabel('Temperature (K)','fontsize',18);
            xlabel('Year','fontsize',18);
            set(gca,'fontsize',16);
            legend('90-75S toz regression prediction','Surface temperature');            
            filename = [titles2{area},'_',sprintf('%02d',i),ENSOext]; 
            filedir = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/InitialPredictions/',filename];
            export_fig(filedir,'-png');                 
%             S = dbstack('-completenames');
%             S(1).file
%             figLog(mfilename('fullpath'),filename);
        end
    end
end

%%
plotozonetemp = 0;
if plotozonetemp
    for i = 1:9
        fig = createfig('medium','on');
        toztemp = squeeze((toz_dataMonthArrange.highcl(i,11,1:end))-...
            nanmean(toz_dataMonthArrange.highcl(i,11,1:end),3))./std(toz_dataMonthArrange.highcl(i,11,1:end),1,3);
        temptemp = (temperature(1:21,i,1) - nanmean(temperature(1:21,i,1)))./std(temperature(1:21,i,1));
        temptemp_rmENSO = (temperature_rmENSO(1:21,i,1) - nanmean(temperature_rmENSO(1:21,i,1)))./std(temperature_rmENSO(1:21,i,1));
        plot(tozdates(1):tozdates(2),toztemp,'r','LineWidth',3)
        hold on
        plot(tozdates(1):tozdates(2),temptemp,'b','LineWidth',3);
        plot(tozdates(1):tozdates(2),temptemp_rmENSO,'g','LineWidth',3);
        corrofpoint(i,1) = corr(toztemp,temptemp);
        corrofpoint(i,2) = corr(toztemp,temptemp_rmENSO);
        corrofpoint(i,3) = corrtemp(i,1);
        legend('toz','TS','TS_rmENSO');
    end
end

%% testing
% %% histogram
% load('/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefs_Justins');
% toplot = polar_correlations(:,61:end);
% % createfig('medium','on')
% % h = histogram(sort(toplot(:)),50,'Normalization','Probability');
% yERA = -1:.01:1;
% mu = nanmean(toplot(:));
% sigma = nanstd(toplot(:));    
% fERA = exp(-(yERA-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% % createfig('medium','on')
% %plot(y,f,'LineWidth',2);
% 
% cbrew = cbrewer('qual','Set1',10);
% createfig('medium','on')
% for i = 1:9
%     subplot(3,3,i)
%     toplot = SHreg.r.ind.lowGHG(i,:,1:48);
%     h(i) = histogram(sort(toplot(:)),50,'Normalization','Probability');
%     y(i,:) = -1:.01:1;
%     mu(i,:) = nanmean(toplot(:));
%     sigma(i,:) = nanstd(toplot(:));    
%     f(i,:) = exp(-(y(i,:)-mu(i,:)).^2./(2*sigma(i,:)^2))./(sigma(i,:)*sqrt(2*pi));
%     title(num2str(i))
% end
% 
% toplot = SHreg.r.composite(3,:,1:48);
% yc = -1:.01:1;
% mu = nanmean(toplot(:));
% sigma = nanstd(toplot(:));    
% fc = exp(-(yc-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% 
% toplot = SHreg.r.ensave(3,:,1:48);
% mu = nanmean(toplot(:));
% sigma = nanstd(toplot(:));    
% fe = exp(-(yc-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
% 
% createfig('medium','on')
% for i = 1:9
%    ph(i) = plot(y(i,:),f(i,:),'LineWidth',2,'color',cbrew(i,:));
%    hold on
% end
% ph(i+1) = plot(yERA,fERA,'LineWidth',2);
% ph(i+2) = plot(yc,fc,'LineWidth',2);
% ph(i+3) = plot(yc,fe,'LineWidth',2);
% legend(ph)

%%
cbrew = cbrewer('qual','Set1',10);
createfig('medium','on');
m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[-90 0]);


% for j = 1:size(lonland,1)
%     for k = 1:[3,5,9]%size(SHreg.r.ind.highcl,1)
%         h(j,k) = m_plot(lonsofmin2(k,j),latsofmin(k,j),'o');
%         hold on
%         set(h(j,k),'MarkerFaceColor',cbrew(k,:),'MarkerEdgeColor',cbrew(k,:),'MarkerSize',10);
%     end
% end
% 
% m_coast('color','k','LineWidth',3);
% m_grid('ytick',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90],'XaxisLocation','bottom','fontsize',18);
% title('Min');
% 
% filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/InitialPredictions/','Minlocations',ENSOext];
% export_fig(filename,'-pdf');
% 
% %
createfig('medium','on');
for j = 1:size(lonland,1)
    for k = 1:size(SHreg.r.ind.highcl,1)
        h(j,k) = m_plot(lonsofmax2(k,j),latsofmax(k,j),'o');
        hold on
        set(h(j,k),'MarkerFaceColor',cbrew(k,:),'MarkerEdgeColor',cbrew(k,:),'MarkerSize',10,'LineStyle','none');
    end
end
%lh = legend(h(1,:),'1','2','3','4','5','6','7','8','9');
lh = legend('3','5','9');
set(lh,'location','EastOutside','box','off');
m_coast('color','k','LineWidth',3);
m_grid('ytick',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90],'XaxisLocation','bottom','fontsize',18);
title('Max');

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/InitialPredictions/','Maxlocations',ENSOext];
export_fig(filename,'-pdf');
%%
% rtoplot = nanmean(SHreg.r.ind.highcl([1,6,3,9],:,:),1);
% ptoplot = nanmean(SHreg.p.ind.highcl([1,6,3,9],:,:),1);
% subplotmaps(rtoplot,SHreg.lons,SHreg.lats,{'div','RdBu'},1,ptoplot,16,{['No.',num2str(i)]},'Longitude','latitude','Correlation','on',...
%             [-.5 .5],18,[-90:10:90],[-90:10:90],[SHreg.lons(1:10:end)],[SHreg.lons(1:10:end)],'',1,[0 360],[-90,0],0,'none',1,'Miller Cylindrical');
