function [data,years,data_composite,dataMonthArrange] = predruns_ReadInlayer(directory,files,...
    var,dates,lats,ifdetrend)

% Read in pred run 3d data
count = 1;
for i = 1:length(files)
    %read in weighted area average data [height,time]
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    if i == 1
        [~,latind(1)] = min(abs(data(i).lat - lats(1)));
        [~,latind(2)] = min(abs(data(i).lat - lats(2)));
    end    
    
    % construct year only vector
    years(i).y = CCMI_years(data(i).date);      
    
    %constructing composite
    if i == 1
        dateindfirst = find(years(1).y == dates(1),1);
        dateindlast = find(years(1).y == dates(2),1,'last');
    end
    data_composite.data(:,:,count:count+dateindlast-dateindfirst) = data(i).(var)(:,:,dateindfirst:dateindlast);
    count = count+size(data(i).(var)(:,:,dateindfirst:dateindlast),3);   
    
    dateind = find(years(i).y >= dates(1) & years(i).y <= dates(2));  
    for j = 1:12                
        if ifdetrend
            for k = 1:length(data(1).lat)
                dataMonthArrange(i,j,:,k,:) = detrend(squeeze(data(i).TS(:,k,dateind(1)+j-1:12:dateind(end)))')...
                    + repmat(nanmean(squeeze(data(i).TS(:,k,dateind(1)+j-1:12:dateind(end)))),[length(data(i).lon),1])';
            end
        else
            dataMonthArrange(i,j,:,:,:) = data(i).TS(:,:,dateind(1)+j-1:12:dateind(end));
        end        
    end
end

for j = 1:12
    if ifdetrend
         bp = 1:size(data_composite.data,3)/length(files):length(files)*size(data_composite.data,3)/length(files);
        for k = 1:length(data(1).lat)
            data_composite.montharrange(j,:,k,1:length(data_composite.data)/12) = ...
                (detrend(squeeze(data_composite.data(:,k,j:12:end))','linear',bp)...
                +repmat(squeeze(nanmean(data_composite.data(:,k,j:12:end))),[1,length(data(i).lon)]))';
        end
    else
        data_composite.montharrange(j,:,:,1:length(data_composite.data)/12) = data_composite.data(:,:,j:12:end);
    end
end

end
    