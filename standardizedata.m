function [standardout,regpres,latitudes,longitudes,years] = standardizedata(filename,var)

%% Read in 4d data and output standrd data

[~,data,~] = Read_in_netcdf(filename);

years = CCMI_years(data.date);
latitudes = data.lat;
longitudes = data.lon;

pressure = permute(repmat(data.hyam.*data.P0,1,size(data.(var),1),size(data.(var),2),...
    size(data.(var),4)) + repmat(data.hybm,1,size(data.(var),1),size(data.(var),2),...
    size(data.(var),4)) .* permute(repmat(data.PS,1,1,1,length(data.lev)),[4,1,2,3]),[2,3,1,4])./100;

%% integrate to regular pressure

for i = 1:size(data.(var),1)
    for j = 1:size(data.(var),2)
        [standardout(i,j,:,:),regpres] = intRegPres(squeeze(data.(var)(i,j,:,:)),squeeze(pressure(i,j,:,:)));
    end
end

end
