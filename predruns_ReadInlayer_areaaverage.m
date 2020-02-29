function [data,years,varweighted,data_composite,dataMonthArrange] = ...
    predruns_ReadInlayer_areaaverage(directory,files,var,dates,lats,ifdetrend,varmonth)

% definging end point
ifmore = varmonth > 12;

if sum(ifmore) > 0
    ext = 12;
else
    ext = 0;
end

% Read in pred run 3d data
count = 1;
for i = 1:length(files)
    %read in weighted area average data [height,time]
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    if i == 1
        [~,latind(1)] = min(abs(data(i).lat - lats(1)));
        [~,latind(2)] = min(abs(data(i).lat - lats(2)));
    end
    %weighted average    
    varweighted(i,:) = weightedaverage(squeeze(nanmean(data(i).(var)(:,latind(1):latind(2),:),1)),data(i).lat(latind(1):latind(2)));                            
    
    % construct year only vector
    if i == 1
        years(i).y = CCMI_years(data(i).date,0);      
    else
        years(i).y = CCMI_years(data(i).date,1);      
    end
    
    %constructing composite
    dateindfirst = find(years(i).y == dates(1),1);
    dateindlast = find(years(i).y == dates(2),1,'last')-ext;
    data_composite.data(count:count+dateindlast-dateindfirst) = varweighted(i,dateindfirst:dateindlast);
    count = count+size(data(i).(var)(:,:,dateindfirst:dateindlast),3);        
    
end
    
for i = 1:size(varweighted,1)    
    dateind = find(years(i).y >= dates(1) & years(i).y <= dates(2));  
    dateind = dateind(1:end-ext);
    for j = 1:12    
        if ifdetrend            
            dataMonthArrange(i,j,:) = detrend(varweighted(i,dateind(1)+j-1:12:dateind(end))) + ...
                squeeze(nanmean(varweighted(i,dateind(1)+j-1:12:dateind(end))));                     
        else
             dataMonthArrange(i,j,:) = varweighted(i,dateind(1)+j-1:12:dateind(end));
        end
    end
end


for j = 1:12
    if ifdetrend
        bp = 1:length(data_composite.data)/length(files):length(files)*length(data_composite.data)/length(files);
        data_composite.montharrange(j,1:length(data_composite.data)/12) = ...
            detrend(data_composite.data(j:12:end),'linear',bp) + nanmean(data_composite.data(j:12:end));
    else
        data_composite.montharrange(j,1:length(data_composite.data)/12) = data_composite.data(j:12:end);
    end
end

end
    
