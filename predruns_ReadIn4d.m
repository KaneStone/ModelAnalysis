function [data,years,data_at_level,data_composite,dataMonthArrange]...
    = predruns_ReadIn4d(directory,files,var,dates,level,ifdetrend)

data_at_level = zeros(length(files),144,96,(dates(2)-dates(1)+1)*12);
% Read in pred run 4d data
count = 1;
filename = ['/Volumes/MyBook/work/data/predruns/output/',var,'_',num2str(level),'hPa',files(1).name(end-19:end-13),'.mat'];
for i = 1:length(files)
    %read in weighted area average data [height,time]
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    
    %calculate pressure from hybrid_height    
    pressure(i).p = permute(repmat(data(i).hyam.*data(i).P0,1,size(data(i).(var),1),size(data(i).(var),2),size(data(i).(var),4)) + ...
       repmat(data(i).hybm,1,size(data(i).(var),1),size(data(i).(var),2),size(data(i).(var),4)) ...
       .* permute(repmat(data(i).PS,1,1,1,length(data(i).lev)),[4,1,2,3]),[2,3,1,4])./100;     
        
    % construct year only vector
    years(i).y = CCMI_years(data(i).date,1);
    
    %test = interpn(pressure(i).p,1:144,1:96,500,1:300);
    
    % interpolate onto regular pressure
    
%     %testing
%     for l = 1:size(pressure(i).p,4)
%         scatteredInterpolant(data(i).lon,data(i).lat,data(i).lev,log(squeeze(pressure(i).p(:,:,:,l))));
%     end
    
    if ~exist(filename,'file')
        
        for j = 1:size(pressure(i).p,1)
            j
            for k = 1:size(pressure(i).p,2)
                for l = 1:size(pressure(i).p,4)
                    data_at_level(i,j,k,l) = interp1qr(log(squeeze(pressure(i).p(j,k,:,l))),squeeze(data(i).(var)(j,k,:,l)),log(level));
                end
            end
        end            
    end
end
if ~exist(filename,'file')
    save(filename,'data_at_level')
else
    load(filename);
end

for i = 1:size(data_at_level,1)    
    
    %constructing composite
    if i == 1
        dateindfirst = find(years(1).y == dates(1),1);
        dateindlast = find(years(1).y == dates(2),1,'last');
    end
    
    data_composite.data(:,:,count:count+dateindlast-dateindfirst) = data_at_level(i,:,:,dateindfirst:dateindlast);
    
    count = count+size(data_at_level(:,:,:,dateindfirst:dateindlast),4);   
    
    dateind = find(years(i).y >= dates(1) & years(i).y <= dates(2));  
    for j = 1:12       
        if ifdetrend
            for k = 1:length(data(1).lat)
                dataMonthArrange(i,j,:,k,:) = detrend(squeeze(data_at_level(i,:,k,dateind(1)+j-1:12:dateind(end)))') + ...
                    repmat(nanmean(squeeze(data_at_level(i,:,k,dateind(1)+j-1:12:dateind(end))),2),[1,dates(2) - dates(1)+1])';
            end
        else
            dataMonthArrange(i,j,:,:,:) = squeeze(data_at_level(i,:,:,dateind(1)+j-1:12:dateind(end)));
        end
    end
end

for i = 1:12
    if ifdetrend
        bp = 1:size(data_composite.data,3)/length(files):length(files)*size(data_composite.data,3)/length(files);
        for k = 1:length(data(1).lat)
            data_composite.montharrange(i,:,k,1:length(data_composite.data)/12) = (detrend(squeeze(data_composite.data(:,k,j:12:end))','linear',bp)...
                +repmat(squeeze(nanmean(data_composite.data(:,k,j:12:end),3)),[1,length(data_composite.data)/12])')';
        end
    else
        data_composite.montharrange(j,:,:,1:length(data_composite.data)/12) = data_composite.data(:,:,j:12:end);
    end 
end
if ~ifdetrend
    dataMonthArrange = permute(dataMonthArrange,[1,2,5,4,3]);
end
end

