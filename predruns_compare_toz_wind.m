% Compare daily total column ozone and 10 hPa 60N zonal wind.
clear variables

%% inputs

SSWwind = 0;
SPVwind = 40;

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
%     varweighted_shift(i,end-365:end) = NaN;
    winddaily_shift(:,1:365) = NaN;
    
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

%%
testingtoz = 0;
if testingtoz
    [~,test,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/predruns/toz/highCl/daily/temp/out.nc');
    [~,test2,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/predruns/toz/highCl/TOZ_b.e11.BWTREFC2.f19_g16.ccmi34.HighCl.002.nc');
    [~,latind2(1)] = min(abs(test.lat - lats(1)));
    [~,latind2(2)] = min(abs(test.lat - lats(2)));

    testwa = weightedaverage(squeeze(nanmean(test.toz(:,latind2(1):latind2(2),:),1)),test.lat(latind2(1):latind2(2)));  
    test2wa = weightedaverage(squeeze(nanmean(test2.toz(:,latind2(1):latind2(2),:),1)),test2.lat(latind2(1):latind2(2)));  

    for j = 1:12
        testma(j,:) = testwa(j:12:end-1);
        test2ma(j,:) = test2wa(j:12:end);
    end
    testma2 = squeeze(tozdaily_monthlymean(1,:,:));
    test2ma2 = squeeze(tozdaily_monthlymean(2,:,:));
    figure; plot(testma(12,:)); hold on; plot(testma2(12,:),'--')
    figure; plot(test2ma(12,:)); hold on; plot(test2ma2(12,:),'--')
end
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
% winddaily_shift(:,end-365:end) = NaN;
winddaily_shift(:,1:365) = NaN;

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
Centraldates = [];
SSWSPVmonth.SSW = zeros(10,30);
SSWSPVmonth.SPV = zeros(10,30);
SSWSPVyears.SPV = zeros(10,30);
SSWSPVyears.SSW = zeros(10,30);
for j = 1:size(winddaily_years_shift,1)
    count1 = 1;
    count = 1;
    for i = 1:size(winddaily_years_shift,2)
        SSWind(j,i).w = find(squeeze(winddaily_years_shift(j,i,days_shift(6):days_shift(13)-1)) < SSWwind); % days_shift(6) = December 1st, days_shift(10) = Apr 1st
        SSWind(j,i).w2 = SSWind(j,i).w - 31;
        decind = SSWind(j,i).w2 <= 0;
        SSWind(j,i).w2(decind) = SSWind(j,i).w2(decind)+365; 
        % find difference
        diffdates = diff(SSWind(j,i).w);    
        warmings = [];

        if ~isempty(SSWind(j,i).w)
            diffdates = [1;diffdates];
            warmings = [1;find(diffdates > 20,1)];       
            finalwarm = find(diffdates >= 20,1,'last');
            diffdates(finalwarm:end) = 1;
        end
        %find central date (first date that winds become easterly);
        
        
        
        for k = 1:length(diffdates)
            if SSWind(j,i).w(k) <= 150
            if k == 1 || diffdates(k) >= 20                
                
                diffdates2 = diffdates(k+1:end);
                if j == 1 && i == 15
                    abc = 1
                end
                if sum(diffdates2(diffdates2>=20)) == 0
                    % does it still go positive in the future?
                    
                    if sum(diffdates2(diffdates2>=2)) >=1 && SSWind(j,i).w(k) < 150 % stating if direction changes again before April 1st then not the final warming!
                        % not a final warming, just a shorter return to
                        % positive period
                    else
                        % final warming
                        break;
                    end
                end
                if i == 26
                    abc = 1;
                end
                Centraldates(j).SSW(count,:,:,:,:,:,:) = datevec(datenum(yearsUnique_shift(i),1,SSWind(j,i).w2(k)));
                Centraldates(j).SSWday(count) = SSWind(j,i).w(k);
                if SSWSPVmonth.SSW(j,i) == 0 
                    SSWSPVmonth.SSW(j,i) = ~isempty(SSWind(j,i).w).*Centraldates(j).SSW(count,2);
                end
                SSWSPVyears.SSW(j,i) = ~isempty(SSWind(j,i).w);  % check here!!!!
                count = count+1;
            end 
            end
        end

        %SSWSPVyears.SSW(j,i) = ~isempty(SSWind(j,i).w);  
              
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
                Centraldates(j).SPVday(count1) = SPVind(j,i).w(k);
                SSWSPVmonth.SPV(j,i) = ~isempty(SPVind(j,i).w).*Centraldates(j).SPV(count1,2);
                count1 = count1+1;
            end
        end
                
        SSWSPVyears.SPV(j,i) = ~isempty(SPVind(j,i).w);
    end
    
    % find upper and lower March 20th percentiles
    upperpct(j) = prctile(squeeze(tozdaily_monthlymean(j,3,:)),80);
    upperind(j,:) = squeeze(tozdaily_monthlymean(j,3,:)) >= upperpct(j);
    upperind2(j,:) = find(squeeze(tozdaily_monthlymean(j,3,:)) >= upperpct(j));
    lowerpct(j) = prctile(squeeze(tozdaily_monthlymean(j,3,:)),20);    
    lowerind(j,:) = squeeze(tozdaily_monthlymean(j,3,:)) <= lowerpct(j);
    lowerind2(j,:) = find(squeeze(tozdaily_monthlymean(j,3,:)) <= lowerpct(j));
end

SSWSPVyears.SSW;
SSWSPVyears.SPV;

SSWSPVyears.SSWonly = SSWSPVyears.SSW; 
SSWSPVyears.SPVonly = SSWSPVyears.SPV; 

% finding SSW and SPV months during ozone extremes.
SSWSPVmonth.SSWduringozone = SSWSPVmonth.SSW;
SSWSPVmonth.SPVduringozone = SSWSPVmonth.SPV;

SSWSPVmonth.SSWduringozone (SSWSPVmonth.SSWduringozone ~= upperind.*SSWSPVmonth.SSW) = 0;
SSWSPVmonth.SPVduringozone (SSWSPVmonth.SPVduringozone ~= lowerind.*SSWSPVmonth.SPV) = 0;

% 
SSWSPVyears.SSWonly (SSWSPVyears.SSW == SSWSPVyears.SPV) = 0;
SSWSPVyears.SPVonly (SSWSPVyears.SPV == SSWSPVyears.SSW) = 0;

SSWSPVyears.both = SSWSPVyears.SSW;
SSWSPVyears.both (SSWSPVyears.both ~= SSWSPVyears.SPV) = 0;

SSWfreq = sum(SSWSPVyears.SSWonly(:))./length(SSWSPVyears.SSWonly(:));
SPVfreq = sum(SSWSPVyears.SPVonly(:))./length(SSWSPVyears.SSWonly(:));

%%
SSWlogicalday = zeros(10,30);
SPVlogicalday = zeros(10,30);
ys = 1995:2024;
for i = 1:10
    for j = 1:length(Centraldates(i).SSW(:,1))
        for k = 1:30
            if Centraldates(i).SSW(j,1) == ys(k)
                SSWlogicalday(i,k) = Centraldates(i).SSWday(j);                  
            end            
        end
    end
    
    for j = 1:length(Centraldates(i).SPV(:,1))
        for k = 1:30
            if Centraldates(i).SPV(j,1) == ys(k)
                SPVlogicalday(i,k) = Centraldates(i).SPVday(j);                  
            end
        end
    end
    
end
SSWlogicalday = [SSWlogicalday(:,end),SSWlogicalday(:,1:end-1)];
SPVlogicalday = [SPVlogicalday(:,end),SPVlogicalday(:,1:end-1)];

SSWlogicalday_SSWonly = SSWlogicalday;
SPVlogicalday_SPVonly = SPVlogicalday;
SSWlogicalday_SSWonly (SSWSPVyears.SSW == SSWSPVyears.SPV) = 0;
SPVlogicalday_SPVonly (SSWSPVyears.SSW == SSWSPVyears.SPV) = 0;

%convert days from December to logical array
%% histogram

% all SSWs
SSWSPVmonth.SSW2 = SSWSPVmonth.SSW;
SSWSPVmonth.SPV2 = SSWSPVmonth.SPV;
SSWSPVmonth.SSW2 (SSWSPVmonth.SSW2 == 0) = [];
SSWSPVmonth.SPV2 (SSWSPVmonth.SPV2 == 0) = [];

SSWSPVmonth.SSW2 (SSWSPVmonth.SSW2 == 12) = 0;
SSWSPVmonth.SPV2 (SSWSPVmonth.SPV2 == 12) = 0;

SSWcentraldateaverage = nanmean(cat(2,Centraldates(:).SSWday));
SPVcentraldateaverage = nanmean(cat(2,Centraldates(:).SPVday));

% SSWs during ozone extremes
SSWSPVmonth.SSWduringozone2 = SSWSPVmonth.SSWduringozone;
SSWSPVmonth.SPVduringozone2 = SSWSPVmonth.SPVduringozone;
SSWSPVmonth.SSWduringozone2 (SSWSPVmonth.SSWduringozone2 == 0) = [];
SSWSPVmonth.SPVduringozone2 (SSWSPVmonth.SPVduringozone2 == 0) = [];

SSWSPVmonth.SSWduringozone2 (SSWSPVmonth.SSWduringozone2 == 12) = 0;
SSWSPVmonth.SPVduringozone2 (SSWSPVmonth.SPVduringozone2 == 12) = 0;

SSWSPVmonth.SSWduringozone2 (SSWSPVmonth.SSWduringozone2 == 12) = 0;
SSWSPVmonth.SPVduringozone2 (SSWSPVmonth.SPVduringozone2 == 12) = 0;

% days during ozone extremes
SSWlogicalday_do = SSWlogicalday;
SPVlogicalday_do = SPVlogicalday;
SSWlogicalday_do (SSWlogicalday_do ~= upperind.*SSWlogicalday) = 0;
SPVlogicalday_do (SPVlogicalday_do ~= lowerind.*SPVlogicalday) = 0;

SSWlogicalday_do2 = SSWlogicalday_do;
SPVlogicalday_do2 = SPVlogicalday_do;

SSWlogicalday_do2 (SSWlogicalday_do2 == 0) = [];
SPVlogicalday_do2 (SPVlogicalday_do2 == 0) = [];

% no coincident days during ozone extremes
SSWlogicalday_do3 = SSWlogicalday_SSWonly;
SPVlogicalday_do3 = SPVlogicalday_SPVonly;
SSWlogicalday_do3 (SSWlogicalday_do3 ~= upperind.*SSWlogicalday_SSWonly) = 0;
SPVlogicalday_do3 (SPVlogicalday_do3 ~= lowerind.*SPVlogicalday_SPVonly) = 0;

SSWlogicalday_do4 = SSWlogicalday_do3;
SPVlogicalday_do4 = SPVlogicalday_do3;

SSWlogicalday_do4 (SSWlogicalday_do4 == 0) = [];
SPVlogicalday_do4 (SPVlogicalday_do4 == 0) = [];

% all events days

SSWlogicalday2 = SSWlogicalday;
SPVlogicalday2 = SPVlogicalday;

SSWlogicalday2 (SSWlogicalday2 == 0) = [];
SPVlogicalday2 (SPVlogicalday2 == 0) = [];

% no coincidence days
SSWlogicalday_SSWonly2 = SSWlogicalday_SSWonly;
SPVlogicalday_SPVonly2 = SPVlogicalday_SPVonly;

SSWlogicalday_SSWonly2 (SSWlogicalday_SSWonly2 == 0) = [];
SPVlogicalday_SPVonly2 (SPVlogicalday_SPVonly2 == 0) = [];

%%

% figure; histogram(SSWlogicalday2,10); hold on; histogram(SPVlogicalday2,10); title('All events');
% figure; histogram(SSWlogicalday_SSWonly2,10); hold on; histogram(SPVlogicalday_SPVonly2,10); title('No coincidence');
 figure; histogram(SSWlogicalday_do2,10); hold on; histogram(SPVlogicalday_do2,10); title('All events during ozone extremes');
ylim([0 10]);
xlim([-10 150]);
set(gca,'xtick',days_shift(6:10)-days_shift(6),'xticklabel',{'December','January','February','March','April'});

%% 
figure; histogram(SSWlogicalday2,10); hold on; histogram(SPVlogicalday2,10); title('No coincidence during ozone extremes');
ylim([0 30]);
xlim([-10 150]);
set(gca,'xtick',days_shift(6:10)-days_shift(6),'xticklabel',{'December','January','February','March','April'});
%%
%figure; histogram(SSWlogicalday_do4,10); hold on; histogram(SPVlogicalday_do4,10); title('No coincidence during ozone extremes');
%% find indives which came last in years that have coincident SSW and SPV
plotboth = 0;

if plotboth
count = 1;

for i = 1:size(SSWSPVyears.both,1)
    for j = 1:size(SSWSPVyears.both,2)
        if SSWSPVyears.both(i,j) == 1
            co.SSW = find(squeeze(winddaily_years_shift(i,j,days_shift(6):days_shift(10)-1)) < SSWwind); % days_shift(6) = December 1st, days_shift(10) = Apr 1st
            co.SPV = find(squeeze(winddaily_years_shift(i,j,days_shift(6):days_shift(10)-1)) > SPVwind); % days_shift(6) = December 1st, days_shift(10) = Apr 1st
            co.SSWmax = max(co.SSW);
            co.SPVmax = max(co.SPV);
            if co.SSWmax > co.SPVmax
                co.lastSSW(i,j) = 1;
            else
                co.lastSSW(i,j) = 0;
            end
            if co.SSWmax < co.SPVmax
                co.lastSPV(i,j) = 1;
            else
                co.lastSPV(i,j) = 0;
            end
                
%             createfig('medium','on')
%             yyaxis left
%             plot(days_shift(6):days_shift(11)-1,squeeze(winddaily_years_shift(i,j,days_shift(6):days_shift(11)-1)),'linewidth',2)
%             ylabel('zonal wind (m/s)','fontsize',20);
%             yyaxis right
%             hold on    
%             plot(days_shift(6):days_shift(11)-1,squeeze(tozdaily_years_shift(i,j,days_shift(6):days_shift(11)-1)),'linewidth',2)        
%             xlim([days_shift(6)-10,days_shift(11)+10]);
%             plot([0,400],[upperpct(i),upperpct(i)],'--','color',[.7 .7 .7],'LineWidth',2)
%             plot([0,400],[lowerpct(i),lowerpct(i)],'--','color',[.7 .7 .7],'LineWidth',2)           
%             ylabel('toz','fontsize',20);
%             set(gca,'xtick',days_shift(6:10),'xticklabel',['Dec';'Jan';'Feb';'Mar';'Apr'],'fontsize',20);            
%             filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/Both/EnsMemBoth',sprintf('%02d',i),'_',sprintf('%02d',j),'_zonalwind_TCO_duringOzoneExtremes.pdf'];
%                 
%             export_fig(filename,'-pdf');
            
            % close 1
            %count = count+1;
        end       
    end
end
end
%% Wind first years is actually second year
for i = 1:10
    additionyear = 0;
    SSWyearsforcom = find(SSWSPVyears.SSW(i,:));
    SPVyearsforcom = find(SSWSPVyears.SPV(i,:));
    ozoneind.upper.noSSW(i).a = setdiff(upperind2(i,:)+additionyear,SSWyearsforcom);
    ozoneind.upper.noSPV(i).a = setdiff(upperind2(i,:)+additionyear,SPVyearsforcom);
    ozoneind.lower.noSSW(i).a = setdiff(lowerind2(i,:)+additionyear,SSWyearsforcom);
    ozoneind.lower.noSPV(i).a = setdiff(lowerind2(i,:)+additionyear,SPVyearsforcom);
    ozoneind.upper.SSW(i).a = intersect(upperind2(i,:)+additionyear,SSWyearsforcom);
    ozoneind.upper.SPV(i).a = intersect(upperind2(i,:)+additionyear,SPVyearsforcom);
    ozoneind.lower.SSW(i).a = intersect(lowerind2(i,:)+additionyear,SSWyearsforcom);
    ozoneind.lower.SPV(i).a = intersect(lowerind2(i,:)+additionyear,SPVyearsforcom);
    
    ozoneind.upper.noother(i).a = setdiff(ozoneind.upper.SSW(i).a,ozoneind.upper.SPV(i).a);
    ozoneind.lower.noother(i).a = setdiff(ozoneind.lower.SPV(i).a,ozoneind.lower.SSW(i).a);
end

%find central dates

%%
no.lowernoSPV = sum(logical([ozoneind.lower.noSPV(:).a]));
no.lowerSPV = sum(logical([ozoneind.lower.SPV(:).a]));
no.uppernoSSW = sum(logical([ozoneind.upper.noSSW(:).a]));
no.upperSSW = sum(logical([ozoneind.upper.SSW(:).a]));
no.uppernoother = sum(logical([ozoneind.upper.noother(:).a]));
no.lowernoother = sum(logical([ozoneind.lower.noother(:).a]));


save(['/Volumes/ExternalOne/work/data/predruns/output/TCOwindextremes/coincidence_wind',num2str(SPVwind),'and',num2str(SSWwind),'.mat'],'ozoneind','no');

%no.upperSPV = sum(logical([ozoneind.upper.SPV(:).a]));
%% plot ozone and wind during ozone extremes

if plotozoneextremes
    fsize = 18;

    for ensmem = 1:10
        createfig('large','on');
    for j = 1:2
         subplot(2,1,j)
         if j == 1
             toplot = upperind2;
             titleext = 'Upper';
         else
             toplot = lowerind2;
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
                set(gca,'xtick',days_shift,'xticklabel',monslabel,'fontsize',fsize);
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
        clearvars yearind
    end

    annotation('textbox',[0 1 1 0],'String',['Ensemble member ', num2str(ensmem)],'FitBoxToText',...
        'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');   

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/EnsMem',sprintf('%02d',ensmem),'_zonalwind_TCO_duringOzoneExtremes.pdf'];

    export_fig(filename,'-pdf');
    end
end

%% plot individual ozone and wind during wind extremes
plotwindextremes = 1;
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
        clearvars yearind
    end

    annotation('textbox',[0 1 1 0],'String',['Ensemble member ', num2str(ensmem)],'FitBoxToText',...
        'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');   

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/EnsMem',sprintf('%02d',ensmem),'_zonalwind_TCO_duringWindExtremes.pdf'];

    export_fig(filename,'-pdf');
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
                filext = ['zonalwind_TCOextremes_and_Windextremes_',num2str(SSWwind),'SSW','_',num2str(SPVwind),'SPV'];

                for l = 1:10
                    %isemp(l) = isempty(ozoneind.upper.SSW(l).a);
                    isemp(l) = isempty(ozoneind.upper.noother(l).a);
                end
                indtomean = find(~isemp,1,'last');

            elseif strcmp(condition,'O3noWind')
                toplot = ozoneind.upper.noSSW(ensmem).a(:);         
                titleext = 'Years with upper March TCO extremes but no sudden stratospheric warmings';
                filext = ['zonalwind_TCOextremes_no_Windextremes_',num2str(SSWwind),'SSW','_',num2str(SPVwind),'SPV'];

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
            lh = legend([ph,phoz],'Zonal wind, SPVs','Zonal wind, no SPVs','Lower ozone, SPVs','Lower ozone, no SPVs');
            set(lh,'fontsize',fsize,'box','off','location','NorthWest');            
        end                
    end
end

annotation('textbox',[0 1 1 0],'String',['Coincidences of ozone and wind extremes:', ' 1995',char(8211),'2024'],'FitBoxToText',...
        'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');   

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/EnsCompositeCoincidence','_',num2str(SSWwind),'SSW','_',num2str(SPVwind),'SPV','.pdf'];

export_fig(filename,'-pdf');

%plot(meantoplotcombine(j,l,:)

%% plot individual ozone and wind during both ozone extremes
plotind = 1;
meantoplotwindupper = [];
meantoplotwindlower = [];
meantoplottozlower = [];
meantoplottozupper = [];
meantoplotwind = [];
meantoplottoz = [];
if plotind
    fsize = 18;
    clearvars yearind
    figure;
    set(gcf,'position',[100 100 700 800],'color','white');
    for ensmem = 1:10
        
        for j = 1:2


             toplot_upper = upperind2(ensmem,:);
             toplotens = upperind;
             %toplot = SSWyearsforcom;
             titleext = 'upper ozone';

             toplot_lower = lowerind2(ensmem,:);
             toplotens = lowerind;
             titleext = 'lower ozone';

            for i = 1:length(toplot_upper)
                yeartoplotlower = yearsUnique(toplot_lower(i));
                yeartoplotupper = yearsUnique(toplot_upper(i));
                yearind_upper(i) = find(yearsUnique_shift == yeartoplotupper+1);
                yearind_lower(i) = find(yearsUnique_shift == yeartoplotlower+1);

                fsize = 18;
                monslabel = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

                %yyaxis left
                hold on    
                if j == 1

                    meantoplotwindupper(:,:,ensmem) = [meantoplotwind;squeeze(winddaily_years_shift(ensmem,upperind2(ensmem,:),:))'];
                    meantoplotwindlower(:,:,ensmem) = [meantoplotwind;squeeze(winddaily_years_shift(ensmem,lowerind2(ensmem,:),:))'];
                else                
                    meantoplottozupper(:,:,ensmem) = [meantoplottoz;squeeze(tozdaily_years_shift(ensmem,upperind2(ensmem,:),:))'];
                    meantoplottozlower(:,:,ensmem) = [meantoplottoz;squeeze(tozdaily_years_shift(ensmem,lowerind2(ensmem,:),:))'];
                end                                   
            end
        end
    end
    
    
% plotting
for j = 1:2
    subplot(2,1,j)
   
    if j == 1
        plot(meantoplottozlower(:,:),'-','LineWidth',2,'color',cbrew(1,:));
        hold on    
        plot(meantoplottozupper(:,:),'-','LineWidth',2,'color',cbrew(5,:));    
        phl = plot(nanmean(meantoplottozlower(:,:),2),'-','LineWidth',2,'color',cbrew(2,:));
        phh = plot(nanmean(meantoplottozupper(:,:),2),'-','LineWidth',2,'color',cbrew(6,:));
        titleext = 'Arctic average TCO';
    else
        plot(meantoplotwindlower(:,:),'-','LineWidth',2,'color',cbrew(1,:));
        hold on
        plot(meantoplotwindupper(:,:),'-','LineWidth',2,'color',cbrew(5,:));    
        plot(nanmean(meantoplotwindlower(:,:),2),'-','LineWidth',2,'color',cbrew(2,:));
        plot(nanmean(meantoplotwindupper(:,:),2),'-','LineWidth',2,'color',cbrew(6,:));        
        titleext = ['Zonal wind at 60',char(176),'N and 10hPa'];
    end
    plot([0,400],[0,0],'--k','LineWidth',2);
    
    if j == 1               
        ylim([260,530]);
    else                    
        ylim([-30,80]); 
    end
    xlim([-5 370]);
    set(gca,'xtick',days_shift,'xticklabel',monslabel,'fontsize',fsize);
    title([titleext],'fontsize',fsize+2);
    ylabel('m/s','fontsize',fsize+2);
    if j == 2
        xlabel('Month','fontsize',fsize+2);
    end
    box on
    if j == 1
        %plot upper and lower pct
        plot([0,400],[nanmean(upperpct),nanmean(upperpct)],'--','color','k','LineWidth',2)
        plot([0,400],[nanmean(lowerpct),nanmean(lowerpct)],'--','color','k','LineWidth',2)
        ylabel('Dobson Units','fontsize',fsize+2);            
        lh = legend([phh,phl],'Ozone upper 20th percentile years','Ozone lower 20th percentile years');
        set(lh,'box','off','location','NorthWest','fontsize',fsize);
    else
        plot([0,400],[40 40],'--','color','k','LineWidth',2)
    end
    
end
    
annotation('textbox',[0 1 1 0],'String',['Ozone extreme years'],'FitBoxToText',...
    'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
    'EdgeColor','none','fontweight','bold');   

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/EnsMem',sprintf('%02d',ensmem),'_TCOonly_all.pdf'];

export_fig(filename,'-pdf');
end

%% plotting standard deviations
ozoneupperstandarddeviations = nanstd(meantoplottozupper(:,:),0,2);
ozonelowerstandarddeviations = nanstd(meantoplottozlower(:,:),0,2);

windupperstandarddeviations = nanstd(meantoplotwindupper(:,:),0,2);
windlowerstandarddeviations = nanstd(meantoplotwindlower(:,:),0,2);

tozall = permute(tozdaily_years_shift,[3,1,2]);
tozallstd = nanstd(tozall(:,:),0,2);

windall = permute(winddaily_years_shift,[3,1,2]);
windallstd = nanstd(windall(:,:),0,2);

cbrew2 = cbrewer('qual','Set1',10);

figure;
set(gcf,'position',[100 100 700 800],'color','white');
subplot(2,1,1);
pi1 = plot(ozonelowerstandarddeviations,'color',cbrew(2,:),'LineWidth',4,'LineStyle','-');
hold on
pi2 = plot(ozoneupperstandarddeviations,'color',cbrew(6,:),'LineWidth',4,'LineStyle','-');
pi3 = plot(tozallstd,'color','k','LineWidth',4,'LineStyle','-.');
ylim([0,30]);     
xlim([-5 370]);

set(gca,'fontsize',fsize-2)
set(gca,'xtick',days_shift(1:end-1),'xticklabel',monslabel);
ylabel('Standard deviation (Dobson units)','fontsize',fsize);

title('Model ensemble Arctic average TCO','fontsize',fsize+2)


lh = legend([pi2,pi1,pi3],'Ozone upper 20th percentile years','Ozone lower 20th percentile years','All years');
set(lh,'box','off','location','NorthWest','fontsize',fsize);

%wind
subplot(2,1,2)
box on
hold on
plot(windlowerstandarddeviations,'color',cbrew(2,:),'LineWidth',4,'LineStyle','-');
plot(windupperstandarddeviations,'color',cbrew(6,:),'LineWidth',3,'LineStyle','-');
plot(windallstd,'color','k','LineWidth',4,'LineStyle','-.');

ylim([0,15]); 
xlim([-5 370]);

set(gca,'fontsize',fsize)
set(gca,'xtick',days_shift(1:end-1),'xticklabel',monslabel);    
ylabel('standard deviation (m/s)','fontsize',fsize);
xlabel('Month','fontsize',fsize);
title(['Model ensemble zonal wind at 60',char(176),'N and 10hPa'],'fontsize',fsize+2)

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/Model_Observedtozandwind.pdf'];

export_fig(filename,'-pdf');



%% plot ozone wind for extremes that don't correspond have coincident SSW and SPV
plotindother = 1;
meantoplotwindlower = [];
meantoplottozlower = [];
meantoplotwindupper = [];
meantoplottozupper = [];
fsize = 18;
monslabel = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};
if plotindother
    fsize = 18;
    clearvars yearind
    createfig('large','on');
    for ensmem = 1:10
        
       


             toplot_upper = ozoneind.upper.noSSW(ensmem).a;
             
             %toplot = SSWyearsforcom;
             titleext = 'upper ozone';

             toplot_lower = ozoneind.lower.noSPV(ensmem).a;
             
             titleext = 'lower ozone';            
            for i = 1:length(toplot_upper)
                
                yeartoplotupper = yearsUnique(toplot_upper(i));
                yearind_upper(i) = find(yearsUnique_shift == yeartoplotupper+1);                                                            

                meantoplotwindupper = cat(1,meantoplotwindupper,squeeze(winddaily_years_shift(ensmem,ozoneind.upper.noSSW(ensmem).a(i),:))');
                meantoplottozupper = cat(1,meantoplottozupper,squeeze(tozdaily_years_shift(ensmem,ozoneind.upper.noSSW(ensmem).a(i),:))');
                                                 
            end
            
            for i = 1:length(toplot_lower)
                
                yeartoplotlower = yearsUnique(toplot_lower(i));
                yearind_lower(i) = find(yearsUnique_shift == yeartoplotlower+1);
                
                meantoplotwindlower = cat(1,meantoplotwindlower,squeeze(winddaily_years_shift(ensmem,ozoneind.lower.noSPV(ensmem).a(i),:))');
                meantoplottozlower = cat(1,meantoplottozlower,squeeze(tozdaily_years_shift(ensmem,ozoneind.lower.noSPV(ensmem).a(i),:))');
                
            end
        
    end
    % plotting
for j = 1:2
    subplot(2,1,j)
   
    if j == 1
        plot(meantoplottozlower','-','LineWidth',2,'color',cbrew(1,:));
        hold on    
        plot(meantoplottozupper','-','LineWidth',2,'color',cbrew(5,:));    
        phl = plot(nanmean(meantoplottozlower,1),'-','LineWidth',2,'color',cbrew(2,:));
        phh = plot(nanmean(meantoplottozupper,1),'-','LineWidth',2,'color',cbrew(6,:));
        
        titleext = 'Arctic average TCO';
    else
        plot(meantoplotwindlower','-','LineWidth',2,'color',cbrew(1,:));
        hold on
        plot(meantoplotwindupper','-','LineWidth',2,'color',cbrew(5,:));    
        plot(nanmean(meantoplotwindlower,1),'-','LineWidth',2,'color',cbrew(2,:));
        plot(nanmean(meantoplotwindupper,1),'-','LineWidth',2,'color',cbrew(6,:));
        titleext = ['Zonal wind at 60',char(176),'N and 10hPa'];
    end
    
    if j == 1
        ylim([260,530]);     
        plot([0,400],[nanmean(upperpct),nanmean(upperpct)],'--','color','k','LineWidth',2)
        plot([0,400],[nanmean(lowerpct),nanmean(lowerpct)],'--','color','k','LineWidth',2)
        %plot([185,185],[-40,600],'--k','LineWidth',2);
        lh = legend([phh,phl],'Ozone upper 20th percentile years','Ozone lower 20th percentile years');
        set(lh,'box','off','location','NorthWest','fontsize',fsize);
        ylabel('Dobson Units','fontsize',fsize+2);            
    else
    
        plot([0,400],[0,0],'--k','LineWidth',2);
        plot([0,400],[40,40],'--k','LineWidth',2);
        %plot([185,185],[-40,600],'--k','LineWidth',2);

        ylim([-30,60]);
        
    end
    xlim([-5 370]);
    set(gca,'xtick',days_shift,'xticklabel',monslabel,'fontsize',fsize);
    title([titleext],'fontsize',fsize+2);
    
    if j == 2
        xlabel('Month','fontsize',fsize+2);
        ylabel('m/s','fontsize',fsize+2);
    end
    box on
    
end
    
annotation('textbox',[0 1 1 0],'String','Ozone extreme years that do not correspond to SSW or SPV years','FitBoxToText',...
    'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
    'EdgeColor','none','fontweight','bold');   

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/TCOonly_nocoincidentwithwind.pdf'];

export_fig(filename,'-pdf');
end

%%
plotindother = 1;
meantoplotwindlower = [];
meantoplottozlower = [];
meantoplotwindupper = [];
meantoplottozupper = [];
fsize = 18;
monslabel = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};
if plotindother
    fsize = 18;
    clearvars yearind
    createfig('large','on');
    for ensmem = 1:10
        
       


             toplot_upper = ozoneind.upper.SSW(ensmem).a;
             
             %toplot = SSWyearsforcom;
             titleext = 'upper ozone';

             toplot_lower = ozoneind.lower.SPV(ensmem).a;
             
             titleext = 'lower ozone';            
            for i = 1:length(toplot_upper)
                
                yeartoplotupper = yearsUnique(toplot_upper(i));
                yearind_upper(i) = find(yearsUnique_shift == yeartoplotupper+1);                                                            

                meantoplotwindupper = cat(1,meantoplotwindupper,squeeze(winddaily_years_shift(ensmem,ozoneind.upper.SSW(ensmem).a(i),:))');
                meantoplottozupper = cat(1,meantoplottozupper,squeeze(tozdaily_years_shift(ensmem,ozoneind.upper.SSW(ensmem).a(i),:))');
                                                 
            end
            
            for i = 1:length(toplot_lower)
                
                yeartoplotlower = yearsUnique(toplot_lower(i));
                yearind_lower(i) = find(yearsUnique_shift == yeartoplotlower+1);
                
                meantoplotwindlower = cat(1,meantoplotwindlower,squeeze(winddaily_years_shift(ensmem,ozoneind.lower.SPV(ensmem).a(i),:))');
                meantoplottozlower = cat(1,meantoplottozlower,squeeze(tozdaily_years_shift(ensmem,ozoneind.lower.SPV(ensmem).a(i),:))');
                
            end
        
    end
    % plotting
for j = 1:2
    subplot(2,1,j)
   
    if j == 1
        plot(meantoplottozlower','-','LineWidth',2,'color',cbrew(1,:));
        hold on    
        plot(meantoplottozupper','-','LineWidth',2,'color',cbrew(5,:));    
        phl = plot(nanmean(meantoplottozlower,1),'-','LineWidth',2,'color',cbrew(2,:));
        phh = plot(nanmean(meantoplottozupper,1),'-','LineWidth',2,'color',cbrew(6,:));
        
        titleext = 'Arctic average TCO';
    else
        plot(meantoplotwindlower','-','LineWidth',2,'color',cbrew(1,:));
        hold on
        plot(meantoplotwindupper','-','LineWidth',2,'color',cbrew(5,:));    
        plot(nanmean(meantoplotwindlower,1),'-','LineWidth',2,'color',cbrew(2,:));
        plot(nanmean(meantoplotwindupper,1),'-','LineWidth',2,'color',cbrew(6,:));
        titleext = ['Zonal wind at 60',char(176),'N and 10hPa'];
    end
    
    if j == 1
        ylim([260,530]);     
        plot([0,400],[nanmean(upperpct),nanmean(upperpct)],'--','color','k','LineWidth',2)
        plot([0,400],[nanmean(lowerpct),nanmean(lowerpct)],'--','color','k','LineWidth',2)
        plot([185,185],[-40,600],'--k','LineWidth',2);
        lh = legend([phh,phl],'Ozone upper 20th percentile years','Ozone lower 20th percentile years');
        set(lh,'box','off','location','NorthWest','fontsize',fsize);
        ylabel('Dobson Units','fontsize',fsize+2);            
    else
    
        plot([0,400],[0,0],'--k','LineWidth',2);
        plot([0,400],[40,40],'--k','LineWidth',2);
        plot([185,185],[-40,600],'--k','LineWidth',2);

        ylim([-30,60]);
        
    end
    xlim([-5 370]);
    set(gca,'xtick',days_shift,'xticklabel',monslabel,'fontsize',fsize);
    title([titleext],'fontsize',fsize+2);
    
    if j == 2
        xlabel('Month','fontsize',fsize+2);
        ylabel('m/s','fontsize',fsize+2);
    end
    box on
    
end
    
annotation('textbox',[0 1 1 0],'String','Ozone extreme years that correspond to SSW or SPV years','FitBoxToText',...
    'on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
    'EdgeColor','none','fontweight','bold');   

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/TCOzonalwind/TCO_coincidentwithwind.pdf'];

export_fig(filename,'-pdf');
end
