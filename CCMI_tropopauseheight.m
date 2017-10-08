% Calculate CCMI tropopause pressure;
includeERA = 1;
variable = 'T';
area = '7590S';
% trenddates = [19790201 19990101; 19990201 20150101];
% finding all files
directory = ['/Volumes/MyBook/work/data/CESM-CCMI/',variable,'/zonalmean/',area,'/'];
Z3directory = ['/Volumes/MyBook/work/data/CESM-CCMI/','Z3','/zonalmean/',area,'/'];
datedirectory = ['/Volumes/MyBook/work/data/CESM-CCMI/date/'];
files = dir([directory,'*.nc*']);
Z3files = dir([Z3directory,'*.nc*']);
datefiles = dir([datedirectory,'*.nc*']);

%% Reading in ERA-Interim
ERApressure = ncread('/Volumes/MyBook/work/data/ERA-Interim/Z3/Z3_ERA-Interim.nc','level');
if includeERA
    if strcmp(variable,'Z3')
        eravar = 'z';
        scale = 9.86;
    elseif strcmp(variable,'T')
        eravar = 't';
        scale = 1;
    end
    ERAdata = ncread(['/Volumes/MyBook/work/data/ERA-Interim/',variable,'/',variable,'_ERA-Interim.nc'],eravar)./scale;
    
    ERAlatitude = ncread('/Volumes/MyBook/work/data/ERA-Interim/Z3/Z3_ERA-Interim.nc','latitude');
    if strcmp(area,'6090S')
        ERAlatindex = find(ERAlatititude >= -90 && ERAlatititude <= -60);
    elseif strcmp(area,'6090N')
        ERAlatindex = find(ERAlatitude >= 60 & ERAlatitude <= 90);        
    elseif strcmp(area,'7590S')
        ERAlatindex = find(ERAlatitude >= -90 & ERAlatitude <= -75);        
    end
    
    ERAzonalmean = squeeze(nanmean(ERAdata(:,ERAlatindex,:,:),1));

    for j = 1:size(ERAzonalmean,2)
        Z3weighted_ERA(j,:) = weightedaverage(squeeze(ERAzonalmean(:,j,:)),ERAlatitude(ERAlatindex));        
    end
end


%%

for i = 1:length(files);    
    [~, data(i),~] = Read_in_netcdf([directory,files(i).name]);
    [~, Z3data(i),~] = Read_in_netcdf([Z3directory,Z3files(i).name]);
    [~, date(i),~] = Read_in_netcdf([datedirectory,datefiles(i).name]);    
    data(i).(variable) = squeeze(data(i).(variable));
    pressure(i).wa = (repmat(data(i).ap.*100000,1,size(data(i).(variable),2)) + ...
        repmat(data(i).b,1,size(data(i).(variable),2)) ...
        .* double(repmat(data(i).PS,1,length(data(i).lev))'))./100;        

    years = num2str(date(i).date);
    for j = 1:length(years)
        yearsfinal(i).y(j,:) = str2double(years(j,1:4));
    end
    yearsfinal(i).y = circshift(yearsfinal(i).y,1);
    yearsfinal(i).y(1) = yearsfinal(i).y(2); 
    
    % converting from geopotential height to geometric height
    GMH(i).h = (Z3data(i).Z3*6356000)./(6356000-Z3data(i).Z3);
    
end

%% interpolating onto regular km grid

ht = 3000:1000:100000;

for i = 1:length(files); 
    for j = 1:size(GMH(i).h,2)
        temp_interp(i).t(:,j) = interp1(double(GMH(i).h(:,j)),double(data(i).T(:,j)),ht);
    end
end
%% find turning points
for i = 1:length(files);     
    tp(i).t = diff(temp_interp(i).t(5:end,:)) > -2;
    for j = 1:size(GMH(i).h,2)
        lowest(i).w(j,:) = find(tp(i).t(:,j),10);        
        %difflow = diff(lowest_to_third(i).w(j,:));
        %if 
        tropopauseheight(i).w(j) = ht(lowest(i).w(j,1)+5);
    end    
end
