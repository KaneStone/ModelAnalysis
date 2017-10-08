% Tropopause adjusted height ozonesonde time series and WACCM volcanic
% influence at 50N

% Stations

% - Valentia
% - De Bilt
% - Lindenberg
% - Legionowo
% - Bratt's Lake
% - Kelowna
% - Edmonton
% - Goose Bay

addpath('/Users/kanestone/work/projects/ozonesonde/code');

%stations = {'Churchill'};
stations = {'GooseBay','Lindenberg','Legionowo','DeBilt','Valentia','BrattsLake',...
    'Kelowna','Alert','EurekaLab','NyAlesund','Resolute','Lerwick','Churchill'};

directory = '/Users/kanestone/work/projects/ozonesonde/';

years = 2000:2014;
numyears = 15;

%altitudes = 0:200:30000;
%altsmid = 100:200:30000;

altitudes = 0:1000:30000;
altsmid = 500:1000:30000;


for i = 1:length(stations)
    year = '2000';
    datecount = 1;
    for k = 1:numyears
        [ozonesonde_data.(stations{i}),date, next_year] = read_ozonesonde_raw(stations{i},'ECC',...
            [directory,stations{i},'/ecc/'],year);    
        if next_year
            year = num2str(str2double(year)+1);
            year_missing(k) = 1;              
        else
            year_missing(k) = 0;
            for j = 1:length(ozonesonde_data.(stations{i}))
                vmr(j).profile = ozonesonde_data.(stations{i})(j).profiles(:,2)./...
                    (ozonesonde_data.(stations{i})(j).profiles(:,1));                
                month(i,datecount) = str2double(date(j,5:6));       
                year1(i,datecount) = str2double(year);
                %lapse(j).l = diff(ozonesonde_data.(stations{i})(j).profiles(:,3));
                Altitude(j).profile = -(log(ozonesonde_data.(stations{i})(j).profiles(:,1)./...
                    1013.25)*8.31.*273.15)./(9.806*.0289);
                
                %calculating regular altitudes
                for p = 1:length(altitudes)-1
                    %[altind,altvalue] = find(Altitude(j).profile
                    vmr_altaverage(j,p) = nanmean(vmr(j).profile(Altitude(j).profile > altitudes(p) ...
                        & Altitude(j).profile < altitudes(p+1)));
                    temp_altaverage(j,p) = nanmean(ozonesonde_data.(stations{i})(j).profiles(Altitude(j).profile > altitudes(p) ...
                        & Altitude(j).profile < altitudes(p+1),3));
                end
                
                %calculating tropopause
                %lapse(j,:) = diff(temp_altaverage(j,40:end))*5;
                lapse(j,:) = diff(temp_altaverage(j,3:end));
                count = 1;
                if j == 24
                    test = 1;
                end
                for p = 1:size(lapse,2)-3;                    
                    if lapse(j,p) >= -2
                        lapse_candidate(count) = p; 
                        aveabove(count) = nanmean(lapse(j,p+1:p+2));
                        count = count+1;
%                     else
%                         lapse_candidate = [];
%                         aveabove = [];
                    end
                end               
                if ~exist('aveabove','var')
                    troptemp = [];
                else troptemp = lapse_candidate(aveabove >= -2);
                end
                if isempty(troptemp)
                    tropindex(datecount) = NaN;                                
                else
                    %tropindex(datecount) = troptemp(1)+39;                                
                    tropindex(datecount) = troptemp(1)+2;                                
                end
                   
                clear lapse_candidate aveabove
                
                if isnan(tropindex(datecount))
                    vmr2to5above(i,datecount) = NaN;
                else
                    vmr2to5above(i,datecount) = vmr_altaverage(j,tropindex(datecount)+2);
                    %vmr2to5above(i,datecount) = nanmean(vmr_altaverage(j,tropindex(datecount)+1:tropindex(datecount)+2))*10;
                end
                
                %%
%                 figure
%                 plot(temp_altaverage*-.005,1:length(altitudes)-1)
%                 hold on
%                 plot([0 .45],[tropindex tropindex],'--k')
%                 plot(vmr_altaverage,1:length(altitudes)-1)
                
                datecount = datecount+1;
            end
            year = num2str(str2double(year)+1); 
            clearvars vmr Altitude vmr_altaverage temp_altaverage lapse
        end
    end
end

%% Plotting all monthly values

matstand = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

monthstitle = {'January','February','March','April','May','June','July','August',...
    'September','October','November','December'};
fsize = 20;
mon2plot = 8:10;
yeartoplot = 2000:2015;
for j = 1:length(stations)
    figure;
    set(gcf,'color','white','position',[100 100 1000 700]);
    for i = 1:numyears
        for k = 1:length(mon2plot)
            s1(k) = scatter(repmat(i,1,length(vmr2to5above(j,year1(j,:) == yeartoplot(i) & month(j,:) == mon2plot(k)))),...
                vmr2to5above(j,year1(j,:) == yeartoplot(i) & month(j,:) == mon2plot(k)),80,matstand(k,:),'filled');
            hold on
        end
    end
    set(gca,'xtick',1:1:15,'xticklabel',2000:1:2015,'fontsize',fsize)
    xlabel('Year','fontsize',fsize+2);
    ylabel('ppm','fontsize',fsize+2);
    xlim([0 16])
    %ylim([.05 .5]);
    title(stations{j},'fontsize',fsize+4);
    lh = legend(s1,monthstitle(mon2plot));
        set(lh,'Orientation','horizontal','position',[.51 .01 .01 .01],'box','off','fontsize',fsize+2);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/Ozonesondes/',stations{j},'_',...
        monthstitle{mon2plot},'tropadjusted.pdf'];
    export_fig(filename,'-pdf');
end
