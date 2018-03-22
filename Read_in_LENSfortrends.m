function [ta,years,ta_ensall,ensallyears,lon,lat] = Read_in_LENSfortrends(directory,files,var)

for i = 1:length(files)

    st = regexp(files(i).name,'ta')+3;
    fn = regexp(files(i).name,'19')-2;
    nameprint{i,1} = files(i).name(st:fn);        
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    i
    for j = 1:length(data(i).lon)
        ta(i,j).w = squeeze(data(i).(var)(j,:,:,:));
        %ta(i,j).w = [ta(i,j).w, zeros(size(ta(i,j).w,1),1812-size(ta(i,j).w,2))];
        %ta(i,j).w (ta(i,j).w == 0) = NaN;
    end
%         if size(data(i).time,1) < 900
%             years_temp = repmat([1950:2020],[12,1]);
%             years(i).y = years_temp(:);
% 
%         else
    
    years(i).y = CCMI_years(data(i).date,1);
%         end
end

% take ensemble averages
for j = 1:length(data(1).lon)
    ta_ensall(1,j).w = nanmean(cat(3,ta(:,j).w),3);    
    ensallyears(1,j).y = years(1).y;    
end

lon = data(1).lon;
lat = data(1).lat;

end
