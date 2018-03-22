
function [ta,years,ta_ensall,ensallyears,lon,lat] = Read_in_CCMIACCESS_fortrends(directory,files,var)

for i = 1:length(files)

    st = regexp(files(i).name,'ta')+3;
    fn = regexp(files(i).name,'19')-2;
    nameprint{i,1} = files(i).name(st:fn);        
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    i
    if i == 14
        a = 1;
    end
    for j = 1:length(data(i).lon)
        ta(i,j).w = squeeze(data(i).(var)(j,:,:,:));
        ta(i,j).w = [ta(i,j).w, zeros(size(ta(i,j).w,1),1740-size(ta(i,j).w,2))];
        ta(i,j).w (ta(i,j).w == 0) = NaN;
    end
%         if size(data(i).time,1) < 900
%             years_temp = repmat([1950:2020],[12,1]);
%             years(i).y = years_temp(:);
% 
%         else
    
    years(i).y = data(i).ccmi_year;
%         end
end

for j = 1:length(data(1).lon)
    ta_ensall(1,j).w = nanmean(cat(3,ta(1,j).w),3);
    ta_ensall(2,j).w = nanmean(cat(3,ta(2:3,j).w),3);
    ta_ensall(3,j).w = nanmean(cat(3,ta(4:5,j).w),3);
    ta_ensall(4,j).w = nanmean(cat(3,ta(6:7,j).w),3);    
end
ensallyears(1).y = years(1).y;
ensallyears(2).y = years(2).y;
ensallyears(3).y = years(2).y;
ensallyears(4).y = years(2).y;

lon = data(1).lon;
lat = data(1).lat;

end
