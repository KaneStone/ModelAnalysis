% Read in and quickly plot temperature
clear all

directory = '/Volumes/My Book for Mac/work/data/CESM-CCMI/T/zonalmean/6090S_T_b.e11.BWTSENC2fODS1960.f19_g16.ccmi34.001.nc';
[data.attributes, data.data,data.info] = Read_in_netcdf(directory);
data.data.date = squeeze(data.data.date(1,1,:));
years = num2str(data.data.date);
for j = 1:length(years)
        yearsfinal(j,:) = str2double(years(j,1:4));
end
yearsfinal = circshift(yearsfinal,1);
yearsfinal(1) = yearsfinal(2); 

lat = -75;
pres = 2;
pres2 = 100;
month = 7;

[~,latind] = min(abs(data.data.lat-lat));
[~,presind] = min(abs(data.data.lev-pres));
[~,presind2] = min(abs(data.data.lev-pres2));

extract = squeeze(data.data.T(1,latind,presind,month:12:end));
extract2 = squeeze(data.data.T(1,latind,presind2,month:12:end));

% plotting
subplot(2,1,1);
plot(yearsfinal(month:12:end),extract);
subplot(2,1,2);
plot(yearsfinal(month:12:end),extract2);

