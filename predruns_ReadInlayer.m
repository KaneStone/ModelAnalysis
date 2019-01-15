function [data,years,data_composite,dataMonthArrange,dataMonthArrangeMean] = predruns_ReadInlayer(directory,files,...
    var,dates,ifdetrend)

% Read in pred run 3d data
count = 1;
for i = 1:length(files)
    %read in weighted area average data [height,time]
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);    
    
    % construct year only vector
    if isfield(data(i),'date')
        years(i).y = CCMI_years(data(i).date,1);      
    else        
        temp = repmat(dates(1):dates(2),12,1);
        years(i).y = temp(:);
    end
    
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
                if i == 7
                    a = 1;
                end
                dataMonthArrange(i,j,:,k,:) = detrend(squeeze(data(i).(var)(:,k,dateind(1)+j-1:12:dateind(end)))')...
                    + repmat(nanmean(squeeze(data(i).(var)(:,k,dateind(1)+j-1:12:dateind(end))),2),[1,dates(2) - dates(1)+1])';
                
%                 dataMonthArrange(i,j,:,k,:) = detrend(squeeze(data(i).(var)(:,k,dateind(1)+j-1:12:dateind(end)))') + ...
%                     nanmean(squeeze(data(i).(var)(:,k,dateind(1)+j-1:12:dateind(end))'));
                dataMonthArrangeMean(i,j,:,k,:) = repmat(nanmean(squeeze(data(i).(var)(:,k,dateind(1)+j-1:12:dateind(end))),2),[1,dates(2) - dates(1)+1])';
                
                
            end
        else
            dataMonthArrange(i,j,:,:,:) = data(i).(var)(:,:,dateind(1)+j-1:12:dateind(end));
            dataMonthArrangeMean = [];
        end        
    end
end

for j = 1:12
    if ifdetrend
         bp = 1:size(data_composite.data,3)/length(files):length(files)*size(data_composite.data,3)/length(files);
        for k = 1:length(data(1).lat)
            data_composite.montharrange(j,:,k,1:length(data_composite.data)/12) = ...
                (detrend(squeeze(data_composite.data(:,k,j:12:end))','linear',bp)...
                +repmat(squeeze(nanmean(data_composite.data(:,k,j:12:end),3)),[1,length(data_composite.data)/12])')';
        end
    else
        data_composite.montharrange(j,:,:,1:length(data_composite.data)/12) = data_composite.data(:,:,j:12:end);
    end
end

if ~ifdetrend
    dataMonthArrange = permute(dataMonthArrange,[1,2,5,4,3]);

end
    
