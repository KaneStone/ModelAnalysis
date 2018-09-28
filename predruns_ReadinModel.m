function [surfacedata,tozdata] = predruns_ReadinModel(inputs,seaice)

% Reads in the model surface variable and total ozone

%% Read in data (will read in multiple simulation groups (highC,lowCl,lowGHG)
for i = 1:length(inputs.ClLevel)
    %% Read in surface temperature or other similar variable
    vardirectory = ['/Volumes/MyBook/work/data/predruns/',inputs.var,'/',inputs.ClLevel{i},'/'];
    varfiles = dir([vardirectory,'*.nc']);
    if seaice
        varfiles(1) = []
    end

    [surfacedata.(inputs.ClLevel{i}).data,surfacedata.(inputs.ClLevel{i}).years,surfacedata.(inputs.ClLevel{i}).composite,...
        surfacedata.(inputs.ClLevel{i}).dataMonthArrange,surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean]...
        = predruns_ReadInlayer(vardirectory,varfiles,inputs.var,inputs.timeperiodvar,inputs.detrend);

    %% Read in TOZ or other similar variable

    tozdirectory = ['/Volumes/MyBook/work/data/predruns/',inputs.tozvar,'/',inputs.ClLevel{i},'/'];
    tozfiles = dir([tozdirectory,'*.nc']);

    [tozdata.(inputs.ClLevel{i}).data,tozdata.(inputs.ClLevel{i}).years,tozdata.(inputs.ClLevel{i}).varweighted,...
        tozdata.(inputs.ClLevel{i}).toz_composite,tozdata.(inputs.ClLevel{i}).dataMonthArrange] = ...
        predruns_ReadInlayer_areaaverage(tozdirectory,tozfiles,inputs.tozvar,inputs.timeperiodtoz,inputs.lats,...
        inputs.detrend_ozone);
end

end
