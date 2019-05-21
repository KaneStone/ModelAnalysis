% Compare daily total column ozone and 10 hPa 60N zonal wind.
clear variables

%% inputs
tozmonth = 3;
lats = [63,90];
detrend = 1;
days = [day(datetime('1-Jan-2018'),'dayofyear'),day(datetime('1-Feb-2018'),'dayofyear'),...
    day(datetime('1-Mar-2018'),'dayofyear'),day(datetime('1-Apr-2018'),'dayofyear'),...
    day(datetime('1-May-2018'),'dayofyear'),day(datetime('1-Jun-2018'),'dayofyear'),...
    day(datetime('1-Jul-2018'),'dayofyear'),day(datetime('1-Aug-2018'),'dayofyear'),...
    day(datetime('1-Sep-2018'),'dayofyear'),day(datetime('1-Oct-2018'),'dayofyear'),...
    day(datetime('1-Nov-2018'),'dayofyear'),day(datetime('1-Dec-2018'),'dayofyear'),366];

days_shift = [1,days(8)-days(7)+1,days(9)-days(7)+1,days(10)-days(7)+1,days(11)-days(7)+1,...
    days(12)-days(7)+1,366-days(7)+1,366-days(7)+31+1,366-days(7)+31+28+1,366-days(7)+31+28+31+1,...
    366-days(7)+31+28+31+30+1,366-days(7)+31+28+31+30+31+1,366];

cbrew = cbrewer('qual','Paired',10);
plotozoneextremes = 0;

%%  Read in toz
tozdirectory = '/Volumes/ExternalOne/work/data/predruns/toz/highCl/daily/6090N/';
files = dir([tozdirectory,'*.nc']);
tic;
for i = 1:length(files)
    
    [~,tozdata(i),~] = Read_in_netcdf([tozdirectory,files(i).name]);
    if i == 1
        [~,latind(1)] = min(abs(tozdata(i).lat - lats(1)));
        [~,latind(2)] = min(abs(tozdata(i).lat - lats(2)));
        years = CCMI_years(tozdata(i).date,1);      
        yearsUnique = unique(years);
    end
    varweighted(i,:) = weightedaverage(squeeze(nanmean(tozdata(i).toz(:,latind(1):latind(2),:),1)),tozdata(i).lat(latind(1):latind(2)));      
    %shift to start of July
    varweighted_shift(i,:) = circshift(varweighted(i,:),365-days(7)-1);
    varweighted_shift(i,end-365:end) = NaN;
    
    %take monthly average
    count = 1;
    for k = 1:length(yearsUnique)
        tozdaily_years_shift(i,k,:) = varweighted_shift(i,count:count+364);
        tozdaily_years(i,k,:) = varweighted(i,count:count+364);
        count = count+365;
        for j = 1:12
            tozdaily_monthlymean(i,j,k) = nanmean(tozdaily_years(i,k,days_shift(j):days_shift(j+1)-1),3);
        end
    end
end
toc;

%% Read in wind at 10 hPa and 60N

%% Read in data

winddirectory = '/Volumes/ExternalOne/work/data/predruns/U/highCl/daily/';

windfiles = dir([winddirectory,'*.nc']);

for i = 1:length(files)
    [~,winddata,~] = Read_in_netcdf([winddirectory,windfiles(i).name]);
    
    %combine into 4d array
    
    wind_dataall(i,:,:,:,:) = winddata.U;
    if i == 1
        pressure = winddata.lev;
        latitude = winddata.lat;
    end
end

[~,windlatind] = min(abs(60 - latitude));
[~,windpresind] = min(abs(10 - pressure));

winddaily = squeeze(nanmean(wind_dataall(:,:,windlatind,windpresind,:),2));

winddaily_shift = circshift(winddaily,[0,365 - days(7)-1]);
winddaily_shift(:,end-365:end) = NaN;

count = 1;
for k = 1:length(yearsUnique)
    winddaily_years_shift(:,k,:) = winddaily_shift(:,count:count+364);
    winddaily_years(:,k,:) = winddaily(:,count:count+364);
    count = count+365;
    for j = 1:12
        winddaily_monthlymean(:,j,k) = nanmean(winddaily_years(:,k,days_shift(j):days_shift(j+1)-1),3);
    end
end

%% Find SSWs and SPV indexes

% restructure wind array so that the months go from July-June
%winddaily_years_shift = circshift(winddaily_years,[0,0,365 - days(7)+1]); % this is wrong
%tozdaily_years_shift = circshift(tozdaily_years,[0,0,365 - days(7)+1]);
yearsUnique_shift = [yearsUnique(2:end),yearsUnique(end)+1];

SSWwind = -2;
SPVwind = 40;

for j = 1:size(winddaily_years_shift,1)
    count1 = 1;
    count = 1;
    for i = 1:size(winddaily_years_shift,2)
        SSWind(j,i).w = find(squeeze(winddaily_years_shift(j,i,days_shift(6):days_shift(10)-1)) < SSWwind);
        SSWind(j,i).w2 = SSWind(j,i).w - 31;
        decind = SSWind(j,i).w2 <= 0;
        SSWind(j,i).w2(decind) = SSWind(j,i).w2(decind)+365; 
        % find difference
        diffdates = diff(SSWind(j,i).w);    
        warmings = find(diffdates > 20,1);        

        if ~isempty(SSWind(j,i).w)
            diffdates = [1;diffdates];
        end
        %find central date (first date that winds become easterly);

        for k = 1:length(diffdates)
            if k == 1 || diffdates(k) >= 20 
                Centraldates(j).SSW(count,:,:,:,:,:,:) = datevec(datenum(yearsUnique_shift(i),1,SSWind(j,i).w2(k)));
                count = count+1;
            end
        end

        SSWSPVyears.SSW(j,i) = ~isempty(SSWind(j,i).w);       
        %-------------------------------------------------------------------------------------------
        
        SPVind(j,i).w = find(squeeze(winddaily_years_shift(j,i,days_shift(6):days_shift(10)-1)) > SPVwind);
        
        SPVind(j,i).w2 = SPVind(j,i).w - 31;
        decind = SPVind(j,i).w2 <= 0;
        SPVind(j,i).w2(decind) = SPVind(j,i).w2(decind)+365; 
        
        %find difference
        SPVdiffdates = diff(SPVind(j,i).w);    
        

        if ~isempty(SPVind(j,i).w)
            SPVdiffdates = [1;SPVdiffdates];
        end
                
        for k = 1:length(SPVdiffdates)
            if k == 1 || SPVdiffdates(k) >= 20 
                Centraldates(j).SPV(count1,:,:,:,:,:,:) = datevec(datenum(yearsUnique_shift(i),1,SPVind(j,i).w2(k)));
                count1 = count1+1;
            end
        end
                
        SSWSPVyears.SPV(j,i) = ~isempty(SPVind(j,i).w);
    end
    
    % find upper and lower March 20th percentiles
    upperpct(j) = prctile(squeeze(tozdaily_monthlymean(j,3,:)),80);
    upperind(j,:) = find(squeeze(tozdaily_monthlymean(j,3,:)) >= upperpct(j));
    lowerpct(j) = prctile(squeeze(tozdaily_monthlymean(j,3,:)),20);    
    lowerind(j,:) = find(squeeze(tozdaily_monthlymean(j,3,:)) <= lowerpct(j));
end

SSWSPVyears.SSW;
SSWSPVyears.SPV;

SSWSPVyears.SSWonly = SSWSPVyears.SSW; 
SSWSPVyears.SPVonly = SSWSPVyears.SPV; 

SSWSPVyears.SSWonly (SSWSPVyears.SSW == SSWSPVyears.SPV) = 0;
SSWSPVyears.SPVonly (SSWSPVyears.SPV == SSWSPVyears.SSW) = 0;


%% Wind first years is actually second year
for i = 1:10
    additionyear = 0;
    SSWyearsforcom = find(SSWSPVyears.SSW(i,:));
    SPVyearsforcom = find(SSWSPVyears.SPV(i,:));
    ozoneind.upper.noSSW(i).a = setdiff(upperind(i,:)+additionyear,SSWyearsforcom);
    ozoneind.upper.noSPV(i).a = setdiff(upperind(i,:)+additionyear,SPVyearsforcom);
    ozoneind.lower.noSSW(i).a = setdiff(lowerind(i,:)+additionyear,SSWyearsforcom);
    ozoneind.lower.noSPV(i).a = setdiff(lowerind(i,:)+additionyear,SPVyearsforcom);
    ozoneind.upper.SSW(i).a = intersect(upperind(i,:)+additionyear,SSWyearsforcom);
    ozoneind.upper.SPV(i).a = intersect(upperind(i,:)+additionyear,SPVyearsforcom);
    ozoneind.lower.SSW(i).a = intersect(lowerind(i,:)+additionyear,SSWyearsforcom);
    ozoneind.lower.SPV(i).a = intersect(lowerind(i,:)+additionyear,SPVyearsforcom);
    
    ozoneind.upper.noother(i).a = setdiff(ozoneind.upper.SSW(i).a,ozoneind.upper.SPV(i).a);
    ozoneind.lower.noother(i).a = setdiff(ozoneind.lower.SPV(i).a,ozoneind.lower.SSW(i).a);
end

save(['/Volumes/ExternalOne/work/data/predruns/output/TCOwindextremes/coincidence_wind',num2str(SPVwind),'and',num2str(SSWwind),'.mat'],'ozoneind');

%%
no.lowernoSPV = sum(logical([ozoneind.lower.noSPV(:).a]));
no.lowerSPV = sum(logical([ozoneind.lower.SPV(:).a]));
no.uppernoSSW = sum(logical([ozoneind.upper.noSSW(:).a]));
no.upperSSW = sum(logical([ozoneind.upper.SSW(:).a]));
no.uppernoother = sum(logical([ozoneind.upper.noother(:).a]));
no.lowernoother = sum(logical([ozoneind.lower.noother(:).a]));
%no.upperSPV = sum(logical([ozoneind.upper.SPV(:).a]));
%% plot ozone and wind during ozone extremes

if plotozoneextremes
    fsize = 18;

    for ensmem = 1:10
        createfig('large','on');
    for j = 1:2
         subplot(2,1,j)
         if j == 1
             toplot = upperind;
             titleext = 'Upper';
         else
             toplot = lowerind;
             titleext = 'Lower';
         end
        for i = 1:size(toplot(ensmem,:),2)

            yeartoplot = yearsUnique(toplot(ensmem,i));
            yearind(i) = find(yearsUnique_shift == yeartoplot+1);
            fsize = 18;
            monslabel = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

            yyaxis left
            hold on    
            plot(squeeze(winddaily_years_shift(ensmem,yearind(i),:)),'-','LineWidth',2,'color',cbrew(1,:));
            if i == size(toplot(ensmem,:),2)
               phWind = plot(nanmean(squeeze(winddaily_years_shift(ensmem,yearind,:))),'-','LineWidth',3,'color',cbrew(2,:));
            end
            if i == 1
                plot([0,400],[0,0],'--k');
                plot([185,185],[-40,60],'--k');
                ylim([-25,55]);
                set(gca,'xtick',days_shift,'xticklabel',monslabel,'fontsize',fsize,'color','k');
                title([titleext,' ozone extremes'],'fontsize',fsize+2);
                ylabel(['Zonal wind at 10 hPa and 60',char(176),'N (m/s)'],'fontsize',fsize+2);
                if j == 2
                    xlabel('Month','fontsize',fsize+2);
                end
                box on
            end


            yyaxis right
            hold on
            plot(squeeze(tozdaily_years_shift(ensmem,yearind(i),:)),'-','LineWidth',2,'color',cbrew(5,:));    
            if i == size(toplot(ensmem,:),2)
                phTCO = plot(nanmean(squeeze(tozdaily_years_shift(ensmem,yearind,:))),'-','LineWidth',3,'color',cbrew(6,:));
                if j == 1
                    lh = legend([phWind,phTCO],'Zonal wind','TCO');
                    set(lh,'fontsize',fsize+2,'box','off','location','NorthWest');            
                end
            end

            if i == 1
                %plot upper and lower pct
                plot([0,400],[upperpct(ensmem),upperpct(ensmem)],'--','color',[.7 .7 .7],'LineWidth',2)
                plot([0,400],[lowerpct(ensmem),lowerpct(ensmem)],'--','color',[.7 .7 .7],'LineWidth',2)
                ylabel('Total column ozone (DU)','fontsize',fsize+2);            
            end

            ylim([280,520]);
        end
    end

    annotation('textbox',[0 1 1 0],'String',['Ensemble member ', num2str(ensmem)],'FitBoxToText',...
        'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');   

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/EnsMem',sprintf('%02d',ensmem),'_zonalwind_TCO_duringOzoneExtremes.pdf'];

    export_fig(filename,'-pdf');
    end
end

%% plot individual ozone and wind during wind extremes
plotwindextremes = 0;
if plotwindextremes
    fsize = 18;

    for ensmem = 1:10
        createfig('large','on');
    for j = 1:2
         subplot(2,1,j)
         if j == 1
             toplot = find(SSWSPVyears.SSWonly(ensmem,:));
             titleext = 'SSW';
         else
             toplot = find(SSWSPVyears.SPVonly(ensmem,:));
             titleext = 'SPV';
         end
        for i = 1:length(toplot)

            yeartoplot = yearsUnique(toplot(i));
            yearind(i) = find(yearsUnique_shift == yeartoplot+1);
            fsize = 18;
            monslabel = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

            yyaxis left
            hold on    
            plot(squeeze(winddaily_years_shift(ensmem,yearind(i),:)),'-','LineWidth',2,'color',cbrew(1,:));
            if i == length(toplot)
               phWind = plot(nanmean(squeeze(winddaily_years_shift(ensmem,yearind,:))),'-','LineWidth',3,'color',cbrew(2,:));
            end
            if i == 1
                plot([0,400],[0,0],'--k');
                plot([185,185],[-40,60],'--k');
                ylim([-25,55]);
                set(gca,'xtick',days_shift,'xticklabel',monslabel,'fontsize',fsize);
                title([titleext,' years'],'fontsize',fsize+2);
                ylabel(['Zonal wind at 10 hPa and 60',char(176),'N (m/s)'],'fontsize',fsize+2);
                if j == 2
                    xlabel('Month','fontsize',fsize+2);
                end
                box on
            end


            yyaxis right
            hold on
            plot(squeeze(tozdaily_years_shift(ensmem,yearind(i),:)),'-','LineWidth',2,'color',cbrew(5,:));    
            if i == length(toplot)
                phTCO = plot(nanmean(squeeze(tozdaily_years_shift(ensmem,yearind,:))),'-','LineWidth',3,'color',cbrew(6,:));
                if j == 1
                    lh = legend([phWind,phTCO],'Zonal wind','TCO');
                    set(lh,'fontsize',fsize+2,'box','off','location','NorthWest');            
                end
            end

            if i == 1
                %plot upper and lower pct
                plot([0,400],[upperpct(ensmem),upperpct(ensmem)],'--','color',[.7 .7 .7],'LineWidth',2)
                plot([0,400],[lowerpct(ensmem),lowerpct(ensmem)],'--','color',[.7 .7 .7],'LineWidth',2)
                ylabel('Total column ozone (DU)','fontsize',fsize+2);            
            end

            ylim([280,520]);
        end
    end

    annotation('textbox',[0 1 1 0],'String',['Ensemble member ', num2str(ensmem)],'FitBoxToText',...
        'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');   

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/EnsMem',sprintf('%02d',ensmem),'_zonalwind_TCO_duringWindExtremes.pdf'];

    export_fig(filename,'-pdf');
    end
end

%% plot ensemble ozone and wind during wind extremes
fsize = 18;

createfig('large','on');
%condition = 'O3andWind';
condition = 'O3noWind';
%condition = 'O3andOpposite'; 
%condition = 'O3noOpposite';
for p = 1:2
    createfig('large','on');
    if p == 1
        condition = 'O3andWind';
    else
        condition = 'O3noWind';
    end
    isemp = [];
    for j = 1:2
        meantoplot = [];
        meantoplotozone = [];
        subplot(2,1,j)
        %indtomean = 10;
    for ensmem = 1:10

        if j == 1
            if strcmp(condition,'O3andWind')
                %toplot = ozoneind.upper.SSW(ensmem).a(:);         
                toplot = ozoneind.upper.noother(ensmem).a(:);         
                titleext = 'Years with upper March TCO extremes and sudden stratospheric warmings ';
                filext = 'zonalwind_TCOextremes_and_Windextremes';

                for l = 1:10
                    %isemp(l) = isempty(ozoneind.upper.SSW(l).a);
                    isemp(l) = isempty(ozoneind.upper.noother(l).a);
                end
                indtomean = find(~isemp,1,'last');

            elseif strcmp(condition,'O3noWind')
                toplot = ozoneind.upper.noSSW(ensmem).a(:);         
                titleext = 'Years with upper March TCO extremes but no sudden stratospheric warmings';
                filext = 'zonalwind_TCOextremes_no_Windextremes';

                for l = 1:10
                    %isemp(l) = isempty(ozoneind.upper.noSSW(l).a);
                    isemp(l) = isempty(ozoneind.upper.noSSW(l).a);
                end
                indtomean = find(~isemp,1,'last');


            elseif strcmp(condition,'O3andOpposite')
                toplot = ozoneind.upper.SPV(ensmem).a(:);         
                titleext = 'SPV and upper ozone';
                filext = 'zonalwind_TCOextremes_opposite_Windextremes';
            elseif strcmp(condition,'O3noOpposite')
                toplot = ozoneind.upper.noSPV(ensmem).a(:);         
                titleext = 'no SPV and upper ozone';
                filext = 'zonalwind_TCOextremes_noopposite_windextremes';
            end
        else
            if strcmp(condition,'O3andWind')
                %toplot = ozoneind.lower.SPV(ensmem).a(:);         
                toplot = ozoneind.lower.noother(ensmem).a(:);         
                titleext = 'Years with lower March TCO extremes and strong polar vortex events';

                for l = 1:10
                    %isemp(l) = isempty(ozoneind.lower.SPV(l).a);
                    isemp(l) = isempty(ozoneind.lower.noother(l).a);
                end
                indtomean = find(~isemp,1,'last');                        

            elseif strcmp(condition,'O3noWind')
                toplot = ozoneind.lower.noSPV(ensmem).a(:);         
                titleext = 'Years with lower March TCO extremes but no strong polar vortex events';

                for l = 1:10
                    isemp(l) = isempty(ozoneind.lower.noSPV(l).a);
                end
                indtomean = find(~isemp,1,'last');            

            elseif strcmp(condition,'O3andOpposite')
                toplot = ozoneind.lower.SSW(ensmem).a(:);         
                titleext = 'Sudden stratospheric warmings and upper March TCO extremes';
            elseif strcmp(condition,'O3noOpposite')
                toplot = ozoneind.lower.noSSW(ensmem).a(:);         
                titleext = 'Strong polar vortex and lower March TCO extremes';
            end
        end
        for i = 1:length(toplot)

            yeartoplot = yearsUnique(toplot(i));
            yearind(i) = find(yearsUnique_shift == yeartoplot+1);
            fsize = 18;
            monslabel = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

            yyaxis left
            hold on    
            plot(squeeze(winddaily_years_shift(ensmem,yearind(i),:)),'-','LineWidth',2,'color',cbrew(1,:));
            if i == length(toplot)
                if length(yearind) == 1
                    meantoplot = [meantoplot;squeeze(winddaily_years_shift(ensmem,yearind,:))'];
                else
                    meantoplot = [meantoplot;squeeze(winddaily_years_shift(ensmem,yearind,:))];
                end
               %phWind = plot(nanmean(squeeze(winddaily_years_shift(ensmem,yearind,:))),'-','LineWidth',3,'color',cbrew(2,:));
                if ensmem == indtomean
                    meantoplotcombine(j,p,:) = nanmean(meantoplot);
                    phWind = plot(nanmean(meantoplot),'-','LineWidth',3,'color',cbrew(2,:));
                end
            end
            if i == 1
    %             plot([0,400],[0,0],'--k');
    %             plot([185,185],[-40,60],'--k');
                ylim([-25,60]);
                set(gca,'xtick',days_shift(1:end-1),'xticklabel',monslabel,'fontsize',fsize,'Ycolor','k');
                xlim([-5 370]);
                title([titleext],'fontsize',fsize+2);
                ylabel(['Zonal wind at 10 hPa and 60',char(176),'N (m/s)'],'fontsize',fsize+2);
                if j == 2
                    xlabel('Month','fontsize',fsize+2);
                end
                box on
            end


            yyaxis right
            hold on
            plot(squeeze(tozdaily_years_shift(ensmem,yearind(i),:)),'-','LineWidth',2,'color',cbrew(5,:));    
            if i == length(toplot)
                if length(yearind) == 1
                    meantoplotozone = [meantoplotozone;squeeze(tozdaily_years_shift(ensmem,yearind,:))'];
                else
                    meantoplotozone = [meantoplotozone;squeeze(tozdaily_years_shift(ensmem,yearind,:))];
                end
                if ensmem == indtomean
                    meantoplot_oz_combine(j,p,:) = nanmean(meantoplotozone);
                    phTCO = plot(nanmean(meantoplotozone),'-','LineWidth',3,'color',cbrew(6,:));
                end
                if j == 1 && ensmem == indtomean
                    lh = legend([phWind,phTCO],'Zonal wind','TCO');
                    set(lh,'fontsize',fsize+2,'box','off','location','NorthWest');            
                end
            end

            if i == 1
                %plot upper and lower pct
    %             plot([0,400],[upperpct(ensmem),upperpct(ensmem)],'--','color',[.7 .7 .7],'LineWidth',2)
    %             plot([0,400],[lowerpct(ensmem),lowerpct(ensmem)],'--','color',[.7 .7 .7],'LineWidth',2)
                ylabel('Total column ozone (DU)','fontsize',fsize+2);
                set(gca,'Ycolor','k');
                ylim([260,530]);
            end

            
        end
        clearvars yearind
    end

    annotation('textbox',[0 1 1 0],'String',['Ensemble composite over', ' 1995',char(8211),'2024'],'FitBoxToText',...
        'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');   

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/EnsComposite','_',filext,'.pdf'];

    export_fig(filename,'-pdf');
    end
end

% plot means only
createfig('large','on');
hold on
count = 1;
lstyle = {'-','--','-','--'};
titextcombine = {'Upper ozone','Lower ozone'};
for j = 1:size(meantoplotcombine,1)
    subplot(2,1,j)
    for p = 1:size(meantoplot_oz_combine)
        
        yyaxis left
        ph(p) = plot(squeeze(meantoplotcombine(j,p,:)),'LineStyle',lstyle{count},'LineWidth',3,'color',cbrew(2,:));
        
        if p == 1
            
%             plot([0,400],[0,0],'--k');
%             plot([185,185],[-40,60],'--k');
            ylim([-10,40]);
            set(gca,'xtick',days_shift(1:end-1),'xticklabel',monslabel,'fontsize',fsize,'Ycolor','k');
            xlim([-5 370]);
            title([titextcombine{j}],'fontsize',fsize+2);
            ylabel(['Zonal wind at 10 hPa and 60',char(176),'N (m/s)'],'fontsize',fsize+2);
            if j == 2
                xlabel('Month','fontsize',fsize+2);
            end
            box on
            
        end
                        
        hold on
        yyaxis right
        phoz(p) = plot(squeeze(meantoplot_oz_combine(j,p,:)),'LineStyle',lstyle{count},'LineWidth',3,'color',cbrew(6,:));
        
        if p == 1
            ylabel('Total column ozone (DU)','fontsize',fsize+2);            
            ylim([260,520]);
            set(gca,'Ycolor','k');
        end
                                
        count = count+1;
        
        if p == 2 && j == 1
            lh = legend([ph,phoz],'Zonal wind, SSWs','Zonal wind, no SSWs','Upper ozone, SSWs','Upper ozone, no SSWs');
            set(lh,'fontsize',fsize,'box','off','location','NorthWest');            
        elseif p == 2 && j == 2
            lh = legend([ph,phoz],'Zonal wind, SPVs','Zonal wind, no SPVs','Upper ozone, SPVs','Upper ozone, no SPVs');
            set(lh,'fontsize',fsize,'box','off','location','NorthWest');            
        end                
    end
end

annotation('textbox',[0 1 1 0],'String',['Coincidences of ozone and wind extremes:', ' 1995',char(8211),'2024'],'FitBoxToText',...
        'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');   

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/EnsCompositeCoincidence','.pdf'];

export_fig(filename,'-pdf');

%plot(meantoplotcombine(j,l,:)

%% plot individual ozone and wind during both ozone and wind extremes
fsize = 18;
clearvars yearind
for ensmem = 1:10
    createfig('large','on');
for j = 1:2
     subplot(2,1,j)
     if j == 1
         toplot = ozoneind.upper.SSW(ensmem).a(:);
         %toplot = SSWyearsforcom;
         titleext = 'SSW and upper ozone';
     else
         toplot = ozoneind.lower.SPV(ensmem).a(:);
         titleext = 'SPV and lower ozone';
     end
    for i = 1:length(toplot)
        yeartoplot = yearsUnique(toplot(i));
        yearind(i) = find(yearsUnique_shift == yeartoplot+1);
        
        fsize = 18;
        monslabel = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

        yyaxis left
        hold on    
        plot(squeeze(winddaily_years_shift(ensmem,yearind(i),:)),'-','LineWidth',2,'color',cbrew(1,:));
        if i == length(toplot) && i ~= 1
           phWind = plot(nanmean(squeeze(winddaily_years_shift(ensmem,yearind,:))),'-','LineWidth',3,'color',cbrew(2,:));
        end
        if i == 1
            plot([0,400],[0,0],'--k');
            plot([185,185],[-40,60],'--k');
            ylim([-25,55]);
            set(gca,'xtick',days_shift,'xticklabel',monslabel,'fontsize',fsize);
            title([titleext],'fontsize',fsize+2);
            ylabel(['Zonal wind at 10 hPa and 60',char(176),'N (m/s)'],'fontsize',fsize+2);
            if j == 2
                xlabel('Month','fontsize',fsize+2);
            end
            box on
        end


        yyaxis right
        hold on
        plot(squeeze(tozdaily_years_shift(ensmem,yearind(i),:)),'-','LineWidth',2,'color',cbrew(5,:));    
        if i == length(toplot) && i ~= 1
            phTCO = plot(nanmean(squeeze(tozdaily_years_shift(ensmem,yearind,:))),'-','LineWidth',3,'color',cbrew(6,:));
            if j == 1
                lh = legend([phWind,phTCO],'Zonal wind','TCO');
                set(lh,'fontsize',fsize+2,'box','off','location','NorthWest');            
            end
        end

        if i == 1
            %plot upper and lower pct
            plot([0,400],[upperpct(ensmem),upperpct(ensmem)],'--','color',[.7 .7 .7],'LineWidth',2)
            plot([0,400],[lowerpct(ensmem),lowerpct(ensmem)],'--','color',[.7 .7 .7],'LineWidth',2)
            ylabel('Total column ozone (DU)','fontsize',fsize+2);            
        end
        
        ylim([280,520]);
    end
    clearvars yearind
end

annotation('textbox',[0 1 1 0],'String',['Ensemble member ', num2str(ensmem)],'FitBoxToText',...
    'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
    'EdgeColor','none','fontweight','bold');   

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/EnsMem',sprintf('%02d',ensmem),'_zonalwind_TCO_individual.pdf'];

export_fig(filename,'-pdf');
end
