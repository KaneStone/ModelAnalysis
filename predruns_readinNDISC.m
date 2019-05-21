function [data] = predruns_readinNDISC(datatype)

if strcmp(datatype,'seaiceindex')
    directory = '/Volumes/ExternalOne/work/data/NSIDC_seaice/seaiceindex/';
    files = dir([directory,'*.csv']);
    for i = 1:length(files)
        datatemp = importdata([directory,files(i).name]);
        data(:,:,i) = datatemp.data;
    end
    data (data == -9999) = NaN;
elseif strcmp(datatype,'Goddard')
    directory = '/Volumes/ExternalOne/work/data/NSIDC_seaice/';
    files = dir([directory,'*.nc']);
    for i = 1:length(files)
        [data(i).seaiceconcentration] = ncread([directory,files(i).name],'goddard_bt_seaice_conc_monthly'); %goddard_merged_seaice_conc_monthly
        [data(i).time] = ncread([directory,files(i).name],'time');
        [data(i).latitude] = ncread([directory,files(i).name],'latitude');
        [data(i).longitude] = ncread([directory,files(i).name],'longitude');
        data(i).seaiceconcentration (data(i).seaiceconcentration < 0) = NaN;
    end    
elseif strcmp(datatype,'ERA')
    directory = '/Volumes/ExternalOne/work/data/ERA-Interim/SEAICE/';
    files = dir([directory,'*.nc']);
    for i = 1:length(files)
        [data(i).seaiceconcentration] = ncread([directory,files(i).name],'siconc');
        [data(i).time] = ncread([directory,files(i).name],'time');
        [data(i).latitude] = ncread([directory,files(i).name],'latitude');
        [data(i).longitude] = ncread([directory,files(i).name],'longitude');
        data(i).seaiceconcentration (data(i).seaiceconcentration < 0) = NaN;
    end    
else    
    
    directory = '/Volumes/ExternalOne/work/data/NSIDC_seaice/G10010_SIBT1850_v1.1/';
    files = dir([directory,'*.nc']);
    for i = 1:length(files)
        [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    end    
end
end
