%function [] = ReadOMI()
clear all
directory = '/Users/kanestone/work/projects/OMI/';
year = '2008';
files = dir([directory,'Y',year,'/','*.txt']);    

latitudes = [-89.5:89.5];
longitudes = [-179.5:179.5];
data_daily = zeros(length(longitudes),length(latitudes),366);
mon3 = zeros(1,366);
days = [31,28,31,30,31,30,31,31,30,31,30,31];

% reading in files
for i = 1:length(files)        
    mon(i,:) = files(i).name(18:19);
    mon2(i) = str2double(mon(i,:));
    day(i,:) = files(i).name(20:21);
    day2(i) = str2double(day(i,:));
    doy(i) = datevec2doy([str2double(year),mon2(i),day2(i),0,0,0]);
    mon3(doy(i)) = mon2(i);
    fid = fopen([directory,'Y',year,'/',files(i).name]);
    eof = 0;
    start = 0;
    count = 1;
    count1 = 1;
    count2 = 1;
    while eof == 0
        if start
            line = fgetl(fid);  
            if strcmp(line(35:37),'lat')
                numofmeas = 10;               
                for j = 1:numofmeas
                    data_daily(count1,count2,doy(i)) = str2double(line(1+count:1+count+2));
                    count = count+3;
                    count1 = count1+1;
                end
                count = 1;            
                count1 = 1;
                count2 = count2+1;
            else
                numofmeas = (length(line)-1)/3;
                for j = 1:numofmeas
                    data_daily(count1,count2,doy(i)) = str2double(line(1+count:1+count+2));
                    count = count+3;
                    count1 = count1+1;
                end
                count = 1;
            end                                                                    
        else
            line = fgetl(fid);            
        end
        
        if strcmp(line(1:10),' Latitudes');
            start = 1;
        elseif strcmp(line(42:46),' 89.5')
            eof = 1;
        end
    end
    i
    fclose(fid);
end

data_daily (data_daily == 0) = NaN;

%monthly averaging
data_monthly = zeros(size(data_daily,1),size(data_daily,2),12);

for i = 1:12
    data_monthly(:,:,i) = nanmean(data_daily(:,:,mon3 == i),3);
end

% for i = 1:length(days)
%     if i == 1
%         data_daily_monthly(:,:,i) = nanmean(data_daily(:,:,1:days(i)),3);
%     else
%         data_daily_monthly(:,:,i) = nanmean(data_daily(:,:,sum(days(1:i-1))+1:sum(days(1:i))),3);
%     end
% end
    
save(['/Users/kanestone/work/projects/OMI/output/','OMI_',year,'.mat'],'data_monthly','data_daily','latitudes','longitudes');

%end
%Read in omi