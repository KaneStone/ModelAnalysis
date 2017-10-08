
%% plot TS surface temperatures
clear all
ClLevel = 'highcl';
TSdirectory = ['/Volumes/MyBook/work/data/predruns/','TS','/',ClLevel,'/'];

TSdirectory = ['/Volumes/MyBook/work/data/predruns/','TS','/',ClLevel,'/'];
TSfiles = dir([TSdirectory,'*.nc']);

latitudestoplot = [-80];
longitudestoplot = 0:30:330;

timeperiod = [1995,2015];%[1955,1975]
for i = 1:length(TSfiles)
    %read in TS data
    [~,TSdata(i),~] = Read_in_netcdf([TSdirectory,TSfiles(i).name]);  
    
    years(i).y = CCMI_years(TSdata(i).date);
    dateindfirst = find(years(i).y == timeperiod(1),1);
    dateindlast = find(years(i).y == timeperiod(2),1,'last');    

    for j = 1:12
        for k = 1:length(TSdata(i).lat)
            TSmean1 = nanmean(TSdata(i).TS(:,k,dateindfirst+j-1:12:dateindlast),3);
            TSdata_rearrange_detrend(i,j,k,:,:) = detrend(squeeze(TSdata(i).TS(:,k,dateindfirst+j-1:12:dateindlast)),'linear')...
                + repmat(TSmean1,[1,size(TSdata(i).TS(:,k,dateindfirst+j-1:12:dateindlast),3)]);
        end
    end
end

%% plotting
fsize = 18;
createfig('medium','on')
count = 1;
montoplot = 12;
for i = 1:length(latitudestoplot)
    [~,latind] = min(abs(latitudestoplot(i)-TSdata(1).lat));
    for j = 1:length(longitudestoplot)
        [~,lonind] = min(abs(longitudestoplot(j)-TSdata(1).lon));
        ph(count) = plot(squeeze(nanmean(TSdata_rearrange_detrend(:,montoplot,latind,lonind,:),1)));
        forleg{count} = [num2str(abs(latitudestoplot(i))),'S',', ',num2str(longitudestoplot(j)),'E'];
        count = count+1;  
        hold on
        
    end
end

set(gca,'xtick',1:5:30,'xticklabel',timeperiod(1):5:timeperiod(1)+30,'fontsize',fsize);
xlabel('Year','fontsize',fsize+2);
ylabel('Temperature (K)','fontsize',fsize+2);
lh = legend(ph,forleg);
set(lh,'location','SouthEast');