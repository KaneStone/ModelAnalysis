function [years,PSdata,datefiles,PSfiles] = ReadinCESMDateandPS(date,PS)

datedirectory = '/Volumes/ExternalOne/work/data/CESM-CCMI/date/';
PSdirectory = '/Volumes/ExternalOne/work/data/CESM-CCMI/PS/';
if date
    datefiles = dir([datedirectory,'*.nc*']);
    for i = 1:length(datefiles);    
        [~,datedata(i).d,~] = Read_in_netcdf([datedirectory,datefiles(i).name]);    
        temp(i).y = num2str(datedata(i).d.date);
        for j = 1:length(temp(i).y)
            years(i).y(j,:) = str2double(temp(i).y(j,1:4));
        end
        years(i).y = circshift(years(i).y,1);
        years(i).y(1) = years(i).y(2); 
    end
else
    years = [];
    datefiles = [];
end
if PS
    PSfiles = dir([PSdirectory,'*.nc*']);
    for i = 1:length(PSfiles);    
        [~,PSdata(i).p,~] = Read_in_netcdf([PSdirectory,PSfiles(i).name]);    
    end
else
    PSdata = [];
    PSfiles = [];
end

end
