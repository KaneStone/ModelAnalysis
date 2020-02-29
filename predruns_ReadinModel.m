function [surfacedata,tozdata] = predruns_ReadinModel(inputs)

% Reads in the model surface variable and total ozone

%% Read in data (will read in multiple simulation groups (highC,lowCl,lowGHG)
for i = 1:length(inputs.ClLevel)
    %% Read in surface temperature or other similar variable
    %vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',inputs.var,'/',inputs.ClLevel{i},'/'];
    vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',inputs.var,'/',inputs.ClLevel{i},'/'];
    
    varfiles = dir([vardirectory,'*.nc']); % making sure I am using the same data as the ozone temperature paper
    if strcmp(inputs.var,'ICEFRAC') 
        varfiles(1) = [];
        [surfacedata.(inputs.ClLevel{i}).data,surfacedata.(inputs.ClLevel{i}).years,surfacedata.(inputs.ClLevel{i}).composite,...
            surfacedata.(inputs.ClLevel{i}).dataMonthArrange,surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean]...
            = predruns_ReadInlayer(vardirectory,varfiles,inputs.var,inputs.timeperiodvar,inputs.detrend);
    elseif strcmp(inputs.var,'SNOWHICE')
        varfiles(1) = []; 
        [surfacedata.(inputs.ClLevel{i}).data,surfacedata.(inputs.ClLevel{i}).years,surfacedata.(inputs.ClLevel{i}).composite,...
            surfacedata.(inputs.ClLevel{i}).dataMonthArrange,surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean]...
            = predruns_ReadInlayer(vardirectory,varfiles,inputs.var,inputs.timeperiodvar,inputs.detrend);                    
        
        inputs.var = 'SNOWHLND';        
        vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',inputs.var,'/',inputs.ClLevel{i},'/'];    
        varfiles = dir([vardirectory,'*.nc']);
        
        varfiles(1) = []; % making sure I am using the same data as the ozone temperature paper
        [data2,~,composite2,...
            dataMonthArrange2,dataMonthArrangeMean2]...
            = predruns_ReadInlayer(vardirectory,varfiles,inputs.var,inputs.timeperiodvar,inputs.detrend);
        for j = 1:size(varfiles,1)
            surfacedata.(inputs.ClLevel{i}).data(j).SNOWHICE = surfacedata.(inputs.ClLevel{i}).data(j).SNOWHICE+data2(j).SNOWHLND;                        
        end
%         surfacedata.(inputs.ClLevel{i}).composite.data = surfacedata.(inputs.ClLevel{i}).composite.data + composite2.data;
%         surfacedata.(inputs.ClLevel{i}).composite.montharrange = surfacedata.(inputs.ClLevel{i}).composite.data + composite2.data;
%         surfacedata.(inputs.ClLevel{i}).dataMonthArrange = surfacedata.(inputs.ClLevel{i}).dataMonthArrange + dataMonthArrange2;
%         surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean = surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean + dataMonthArrangeMean2;
        
        surfacedata.(inputs.ClLevel{i}).composite.data = composite2.data;
        surfacedata.(inputs.ClLevel{i}).composite.montharrange = composite2.montharrange;
        surfacedata.(inputs.ClLevel{i}).dataMonthArrange = dataMonthArrange2;
        surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean = dataMonthArrangeMean2;
        
        datasize = size(surfacedata.(inputs.ClLevel{i}).dataMonthArrange);
        
        % read in landfrac
        [~,lf,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/predruns/landfrac/LANDFRAC_b.e11.BWTREFC2.f19_g16.ccmi32.001.cam.h0.2097-08.nc');
        lfrep = permute(repmat(lf.LANDFRAC,[1,1,datasize(1:3)]),[3,4,5,2,1]);
        
        % remove data from ocean areas
        surfacedata.(inputs.ClLevel{i}).dataMonthArrange (lfrep <= .2) = NaN;        
        surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean (lfrep <= .2) = NaN;
    elseif strcmp(inputs.var,'PRECSL')
        varfiles(1) = []; 
        [surfacedata.(inputs.ClLevel{i}).data,surfacedata.(inputs.ClLevel{i}).years,surfacedata.(inputs.ClLevel{i}).composite,...
            surfacedata.(inputs.ClLevel{i}).dataMonthArrange,surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean]...
            = predruns_ReadInlayer(vardirectory,varfiles,inputs.var,inputs.timeperiodvar,inputs.detrend);              
        
        inputs.var = 'PRECSC';        
        vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',inputs.var,'/',inputs.ClLevel{i},'/'];    
        varfiles = dir([vardirectory,'*.nc']);
        varfiles(1) = [];
        [data2,~,composite2,...
            dataMonthArrange2,dataMonthArrangeMean2]...
            = predruns_ReadInlayer(vardirectory,varfiles,inputs.var,inputs.timeperiodvar,inputs.detrend);
        for j = 1:size(varfiles,1)
            surfacedata.(inputs.ClLevel{i}).data(j).PRECSL = surfacedata.(inputs.ClLevel{i}).data(j).PRECSL+data2(j).PRECSC;                        
        end
        
        surfacedata.(inputs.ClLevel{i}).composite.data = surfacedata.(inputs.ClLevel{i}).composite.data + composite2.data;
        surfacedata.(inputs.ClLevel{i}).composite.montharrange = surfacedata.(inputs.ClLevel{i}).composite.data + composite2.data;
        surfacedata.(inputs.ClLevel{i}).dataMonthArrange = surfacedata.(inputs.ClLevel{i}).dataMonthArrange + dataMonthArrange2;
        surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean = surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean + dataMonthArrangeMean2;
        
        
    elseif strcmp(inputs.var,'TS')
        varfiles(1) = [];
        [surfacedata.(inputs.ClLevel{i}).data,surfacedata.(inputs.ClLevel{i}).years,surfacedata.(inputs.ClLevel{i}).composite,...
            surfacedata.(inputs.ClLevel{i}).dataMonthArrange,surfacedata.(inputs.ClLevel{i}).dataMonthArrangeMean]...
            = predruns_ReadInlayer(vardirectory,varfiles,inputs.var,inputs.timeperiodvar,inputs.detrend);      
    else
        surfacedata = [];
    end

    %% Read in TOZ or other similar variable
    
    %tozdirectory = ['/Volumes/ExternalOne/work/data/predruns/',inputs.tozvar,'/',inputs.ClLevel{i},'/'];
    tozdirectory = ['/Volumes/ExternalOne/work/data/predruns/',inputs.tozvar,'/',inputs.ClLevel{i},'/'];
    tozfiles = dir([tozdirectory,'*.nc']);
    tozfiles(1) = []; % making sure I am using the same data as the ozone temperature paper
    [tozdata.(inputs.ClLevel{i}).data,tozdata.(inputs.ClLevel{i}).years,tozdata.(inputs.ClLevel{i}).varweighted,...
        tozdata.(inputs.ClLevel{i}).toz_composite,tozdata.(inputs.ClLevel{i}).dataMonthArrange] = ...
        predruns_ReadInlayer_areaaverage(tozdirectory,tozfiles,inputs.tozvar,inputs.timeperiodtoz,inputs.lats,...
        inputs.detrend_ozone,inputs.varmonth);
end

end
