clear all
ClLevel = 'highCl';
tozvar = 'TOZ';
var = 'TS';
directory = ['/Volumes/MyBook/work/data/predruns/',tozvar,'/',ClLevel,'/'];
TSdirectory = ['/Volumes/MyBook/work/data/predruns/','TS','/',ClLevel,'/'];
files = dir([directory,'*.nc']);
TSfiles = dir([TSdirectory,'*.nc']);
lat = -30;
countlow = 1;
counthigh = 1;
month = 10;
timeperiod = [1995,2015];
for i = 1:length(files)
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    [~,TSdata(i),~] = Read_in_netcdf([TSdirectory,TSfiles(i).name]);
    
    years(i).y = CCMI_years(TSdata(i).date);
    dateindfirst = find(years(i).y == timeperiod(1),1);
    dateindlast = find(years(i).y == timeperiod(2),1,'last');
    
    for j = 1:12
        TSdata_rearrange(i,:,:,j,:) = TSdata(i).TS(:,:,dateindfirst+j-1:12:dateindlast);
        TOZ_rearrange(i,:,:,j,:) = data(i).toz(:,:,dateindfirst+j-1:12:dateindlast);
    end    
    
    %for j = 1:length(TSdata(i).lon)
        for k = 1:length(TSdata(i).lat)
            for l = 10
                r(i,:,k,l) = corr(squeeze(TSdata_rearrange(i,:,k,l+2,:))',squeeze(TOZ_rearrange(i,:,k,l,:))');
            end
        end
    %end
end   

%% 
for i = 1:9
    figure;
    contourf(squeeze(r(i,:,:,10)));
end