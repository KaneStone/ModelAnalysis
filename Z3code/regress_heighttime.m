function [b,bint] = regress_heighttime(data,years,year_range,remove2002)

% data is a structure of similar runs
% data is in the form of [x,time]
% year_range is the range over the regression takes place
% remove2002 is a switch

names = fieldnames(data); 
yearnames = fieldnames(years);
for i = 1:length(data)
    
    if remove2002
        dateind_temp = find(years(i).y >= year_range(1) & years(i).y <= year_range(2));    
        dateind_temp (years(i).y(dateind_temp) >= 2002 & years(i).y(dateind_temp) < 2003) = [];
        dateind(i,:) = dateind_temp;
    else
        dateind(i,:) = find(years(i).(yearname{1}) >= year_range(i) & years(i).(yearnames{1}) <= year_range(2));         
    end    
    dataatdate(i,:,:) = data(i).(names{1})(:,dateind(i,:));    
    for j = 1:12
        dataatdate_rearrange(i,j,:,:) = dataatdate(i,:,j:12:end);        
        for k = 1:size(dataatdate_rearrange,3)
            [b(i,j,k,:),bint(i,j,k,:,:)] = ... %(run,month,pressure,b)
                regress(squeeze(dataatdate_rearrange(i,j,k,:)),...
                [ones(1,size(dataatdate_rearrange(i,j,k,:),4));...
                1:size(dataatdate_rearrange(i,j,:,:),4)]');                            
        end        
    end       
end
end
