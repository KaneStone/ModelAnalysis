function [data,pressure,years,dataRegPres,regPres,dataMonthArrange,dataRegPres_composite]...
    = predruns_ReadIn(directory,files,var,dates)

% Read in pred run 4d data
count = 1;
for i = 1:length(files)
    %read in weighted area average data [height,time]
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    
    %calculate pressure from hybrid_height
    pressure(i).p = permute(repmat(data(i).hyam.*data(i).P0,1,size(data(i).(var),1),size(data(i).(var),3)) + ...
       repmat(data(i).hybm,1,size(data(i).(var),1),size(data(i).(var),3)) ...
       .* permute(repmat(data(i).PS,1,1,length(data(i).lev)),[3,1,2]),[2,1,3])./100; 
    
    %weighted average
    for j = 1:size(data(i).(var),2)
        varweighted(i).wa(j,:) = weightedaverage(squeeze(data(i).(var)(:,j,:)),data(i).lat);        
        presweighted(i).wa(j,:) = weightedaverage(squeeze(pressure(i).p(:,j,:)),data(i).lat);        
    end     
        
    % construct year only vector
    years(i).y = CCMI_years(data(i).date);
    
    % interpolate onto regular pressure
    [dataRegPres(i,:,:),regPres(i,:)] = intRegPres(varweighted(i).wa,presweighted(i).wa);
    
    % composite
    dateindfirst = find(years(i).y == dates(1),1);
    dateindlast = find(years(i).y == dates(2),1,'last');
    dataRegPres_composite.data(1:size(dataRegPres,2),count:count+dateindlast-dateindfirst) = squeeze(dataRegPres(i,:,dateindfirst:dateindlast));
    count = count+size(data(i).(var)(:,:,dateindfirst:dateindlast),3);        
end

dataRegPres(length(files)+1,:,:) = nanmean(dataRegPres,1);
regPres(length(files)+1,:) = nanmean(regPres,1);
years(length(files)+1).y = years(1).y;

for i = 1:size(dataRegPres,1)    
    dateind = find(years(i).y >= dates(1) & years(i).y <= dates(2));  
    for j = 1:12        
        dataMonthArrange(i,j,:,:) = dataRegPres(i,:,dateind(1)+j-1:12:dateind(end));
    end
end

for i = 1:12
    dataRegPres_composite.MonthArrange(i,:,:) = dataRegPres_composite.data(:,i:12:end);
end

end