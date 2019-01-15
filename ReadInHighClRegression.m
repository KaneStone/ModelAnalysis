function [highClData,WAClat] = ReadInHighClRegression(intpres,timeperiod,type,variable)

%% Read in highcl runs
directory = ['/Volumes/ExternalOne/work/data/predruns/',variable,'/',type,'/zonalmean/'];
files = dir([directory,'*.nc']);
for i = 1:length(files)
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);    
    if i == 1
        highclyears = 1995:2024;
        highclyears = repmat(highclyears,[12,1]);
        highclyears = highclyears(:);        
    end
    data(i).(variable) = data(i).(variable)(:,:,highclyears >= timeperiod(1) & highclyears <= timeperiod(2));
    data(i).PS = data(i).PS(:,highclyears >= timeperiod(1) & highclyears <= timeperiod(2));
    
    pressure(i).p =  permute(repmat(data(i).hyam*100000,[1,size(data(i).PS)]),[2,1,3]) + ...
        permute(repmat(data(i).hybm,[1,size(data(i).PS)]),[2,1,3]) .* ...
        double(permute(repmat(data(i).PS,[1,1,length(data(i).lev)]),[1,3,2]));           
    for j = 1:size(data(i).(variable),1)
        [highClData(i).(variable)(:,j,:),~] = intRegPres(squeeze(data(i).(variable)(j,:,:)),...
            squeeze(pressure(i).p(j,:,:))./100,intpres);
    end
end

WAClat = data(1).lat;

%% Read in Uwind
levs = [30,10];
directory = ['/Volumes/ExternalOne/work/data/predruns/U/',type,'/'];
filesU = dir([directory,'*.nc']);
for i = 1:length(filesU)    
    [~,highclU(i),~] = Read_in_netcdf([directory,filesU(i).name]);
    
    if i == 1        
        [~,levind] = min(abs(levs-highclU(i).lev));
    end
    highClData(i).Uregfun = highclU(i).U(levind,:);
    highClData(i).Uregfun = highClData(i).Uregfun(:,highclyears >= timeperiod(1) & highclyears <= timeperiod(2));
    
end

%% Read in heat flux
levs = [100];
directoryS = '/Volumes/ExternalOne/work/data/predruns/VTH3d/highCl/south/';
directoryN = '/Volumes/ExternalOne/work/data/predruns/VTH3d/highCl/north/';
filesHF_south = dir([directoryS,'*.nc']);
filesHF_north = dir([directoryN,'*.nc']);
for i = 1:length(filesHF_south)
    if i == 1
        [~,levind] = min(abs(levs-highclU(i).lev));
    end
    [~,highclHFsouth(i),~] = Read_in_netcdf([directoryS,filesHF_south(i).name]);    
    [~,highclHFnorth(i),~] = Read_in_netcdf([directoryN,filesHF_north(i).name]);    
    
    pressurei_s(i).p =  repmat(highclHFsouth(i).hyai*100000,[1,size(highclHFsouth(i).PS)]) + ...
        repmat(highclHFsouth(i).hybi,[1,size(highclHFsouth(i).PS)]) .* ...
        double(permute(repmat(highclHFsouth(i).PS,[1,length(highclHFsouth(i).ilev)]),[2,1]));                       
    
    pressurei_n(i).p =  repmat(highclHFnorth(i).hyai*100000,[1,size(highclHFnorth(i).PS)]) + ...
        repmat(highclHFnorth(i).hybi,[1,size(highclHFnorth(i).PS)]) .* ...
        double(permute(repmat(highclHFnorth(i).PS,[1,length(highclHFnorth(i).ilev)]),[2,1]));                       
    
    latsouth = data(1).lat > -80 & data(1).lat < -40;
    latnorth = data(1).lat > 40 & data(1).lat < 80;
    
%     for j = 1:size(pressure(i).p,2)
%         PHFS(i).s(j,:) = weightedaverage(squeeze(pressurei_s(i).p(latsouth,j,:)),data(1).lat(latsouth));
%         PHFN(i).n(j,:) = weightedaverage(squeeze(pressurei_p(i).p(latnorth,j,:)),data(1).lat(latnorth));
%     end
    % interpolate to pressure level
    for j = 1:size(highclHFsouth(i).VTH3d,2)        
        HFS(i).h(j) = interp1(log(pressurei_s(i).p(30:60,j)./100),highclHFsouth(i).VTH3d(30:60,j),log(levs));
        HFN(i).h(j) = interp1(log(pressurei_n(i).p(30:60,j)./100),highclHFnorth(i).VTH3d(30:60,j),log(levs));
    end
    
    highClData(i).HFsouth = highclHFsouth(i).VTH3d(levind,highclyears >= timeperiod(1) & highclyears <= timeperiod(2));
    highClData(i).HFnorth = highclHFnorth(i).VTH3d(levind,highclyears >= timeperiod(1) & highclyears <= timeperiod(2));
    
    HFSextract(i).h = HFS(i).h(highclyears >= timeperiod(1) & highclyears <= timeperiod(2));
    HFNextract(i).h = HFN(i).h(highclyears >= timeperiod(1) & highclyears <= timeperiod(2));
end

%% Read in Surface Pressure and construct ENSO

directory = ['/Volumes/ExternalOne/work/data/predruns/TS/',type,'/'];
filesTS = dir([directory,'*.nc']);
for i = 1:length(filesTS)
    [~,highclTS(i),~] = Read_in_netcdf([directory,filesTS(i).name]);    
    highClData(i).TS = highclTS(i).TS(:,:,highclyears >= timeperiod(1) & highclyears <= timeperiod(2));
end

highClData(length(filesTS)+1).TS = nanmean(cat(4,highClData(:).TS),4);

latlimits = [-5 5];
lonlimits = [190 240];

latindex = highclTS(1).lat >= latlimits(1) & highclTS(1).lat <= latlimits(2);
lonindex = highclTS(1).lon >= lonlimits(1) & highclTS(1).lon <= lonlimits(2);

for i = 1:length(filesTS)+1
    for j = 1:12
        NINO_mn(:,j,:) = squeeze(nanmean(highClData(i).TS(lonindex,latindex,j:12:end),1));
        NINO_mn2(j,:) = squeeze(nanmean(NINO_mn(:,j,:),1));        
    end    
    NINO_mn3 = NINO_mn2(:);
    NINOallmean = nanmean(NINO_mn3);
    NINOstd = std(NINO_mn3,1);
    highClData(i).NINO34 = detrend((NINO_mn3 - NINOallmean)./NINOstd);
end

%% taking ensemble averages
highClData(10).(variable) = nanmean(cat(4,highClData(:).(variable)),4);
highClData(10).Uregfun = nanmean(cat(3,highClData(:).Uregfun),3);
highClData(10).HFsouth = nanmean(cat(3,highClData(:).HFsouth),3);
highClData(10).HFnorth = nanmean(cat(3,highClData(:).HFnorth),3);

end
