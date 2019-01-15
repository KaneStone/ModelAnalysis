% Read in area weighted average data, put onto regular pressure levels 
clear all
% include ERA-Interim data

includeERA = 1;
variable = 'Z3';
area = '7590S';
% trenddates = [19790201 19990101; 19990201 20150101];
% finding all files
directory = ['/Volumes/ExternalOne/work/data/CESM-CCMI/',variable,'/zonalmean/',area,'/'];
datedirectory = ['/Volumes/ExternalOne/work/data/CESM-CCMI/date/'];
files = dir([directory,'*.nc*']);
datefiles = dir([datedirectory,'*.nc*']);

%% Reading in ERA-Interim
ERApressure = ncread('/Volumes/ExternalOne/work/data/ERA-Interim/Z3/Z3_ERA-Interim.nc','level');
if includeERA
    if strcmp(variable,'Z3')
        eravar = 'z';
        scale = 9.86;
    elseif strcmp(variable,'T')
        eravar = 't';
        scale = 1;
    end
    ERAdata = ncread(['/Volumes/ExternalOne/work/data/ERA-Interim/',variable,'/',variable,'_ERA-Interim.nc'],eravar)./scale;
    
    ERAlatitude = ncread('/Volumes/ExternalOne/work/data/ERA-Interim/Z3/Z3_ERA-Interim.nc','latitude');
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

    for j = 1:12
        Z3_ERA_rearrange(j,:,:) = Z3weighted_ERA(:,j:12:240);
        Z3_ERA_rearrange2_1(j,:,:) = Z3weighted_ERA(:,j+240:12:276);
        Z3_ERA_rearrange2_2(j,:,:) = Z3weighted_ERA(:,j+288:12:456);
        Z3_ERA_rearrange2 = cat(3,Z3_ERA_rearrange2_1,Z3_ERA_rearrange2_2);
        %Z3_ERA_rearrange = circshift(Z3_ERA_rearrange,[6,0,0]);
        for k = 1:size(Z3_ERA_rearrange,2)
            [bERA(j,k,:),bintERA(j,k,:,:)] = ... %(month,pressure,b)
                regress(squeeze(Z3_ERA_rearrange(j,k,:)),[ones(1,size(Z3_ERA_rearrange(j,k,:),3));1:size(Z3_ERA_rearrange(j,:,:),3)]');                            
            [bERA2(j,k,:),bintERA2(j,k,:,:)] = ... %(month,pressure,b)
                regress(squeeze(Z3_ERA_rearrange2(j,k,:)),[ones(1,size(Z3_ERA_rearrange2(j,k,:),3));1:size(Z3_ERA_rearrange2(j,:,:),3)]');                            
        end        
    end
end


%%

for i = 1:length(files);    
    [~, data(i),~] = Read_in_netcdf([directory,files(i).name]);
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
end

%% take mirror trends and plot

pres = [1000:-50:150 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.05:.15 ... 
    .1:-.005:.015 ...
    .01:-.0005:.0015 ...
    .001:-.00005:.00015];

% first put onto regular grid by interpolating pressure

for i = 1:length(data)
    for j = 1:size(data(i).(variable),2)        
        regp(i).wa(:,j)  = interp1(log(pressure(i).wa(:,j)),data(i).(variable)(:,j),log(double(ERApressure)),'linear','extrap');
    end   
end

%%
for i = 1:length(files)
    filenames{i,:} = files(i).name;
end
save(['/Volumes/ExternalOne/work/data/CESM-CCMI/',variable,'/output/',variable,'_',area,'.mat'],'regp','ERApressure','yearsfinal','filenames');

if includeERA
    ERAInterim = Z3weighted_ERA;
    count = 1;
    yearcount = 1979;
    for i = 1:floor(size(Z3weighted_ERA,2)/12)
        yearsERA(count:count+11) = repmat(yearcount,12,1);
        yearcount = yearcount+1;
        count = count+12;
    end
    yearsERA = [yearsERA,yearsERA(end)+1];
    save(['/Volumes/ExternalOne/work/data/CESM-CCMI/',variable,'/output/',variable,'_ERA-Interim_',area,'.mat'],'ERAInterim','ERApressure','yearsERA');
end


