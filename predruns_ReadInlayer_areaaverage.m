function [data,years,varweighted,data_composite,dataMonthArrange] = predruns_ReadInlayer(directory,files,var,dates,lats)

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
    years(i).y = CCMI_years(data(i).date);      
    
    %constructing composite
    dateindfirst = find(years(i).y == dates(1),1);
    dateindlast = find(years(i).y == dates(2),1,'last');
    data_composite.data(count:count+dateindlast-dateindfirst) = varweighted(i,dateindfirst:dateindlast);
    count = count+size(data(i).(var)(:,:,dateindfirst:dateindlast),3);        
    
end
    
for i = 1:size(varweighted,1)    
    dateind = find(years(i).y >= dates(1) & years(i).y <= dates(2));  
    for j = 1:12        
        dataMonthArrange(i,j,:) = varweighted(i,dateind(1)+j-1:12:dateind(end));
    end
end

for j = 1:12
   data_composite.montharrange(j,1:length(data_composite.data)/12) = data_composite.data(j:12:end);
end

end
    