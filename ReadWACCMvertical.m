function [NumberDensity, MolConc, TAdata, Pressure, Altitude, GMH, Latitudes, Longitudes, datedata] = ...
    ReadWACCMvertical(variable,temporal_period,commondir,prescoordinfile,ND)
% Read in WACCM 4D data and calculate associated pressure values
% temporal_period = 'monthly','daily'
% commondir = common location of files'
% prescoordinfile = specidy whether hybrid height coordinates are in
    % variable file.

Longitudes = [];

if strcmp(temporal_period,'daily')
    dailydir = 'daily/';
else
    dailydir = [];
end

if strcmp(variable,'Extinction')
    names = {'MAM','VCMAM'};
%elseif strcmp(variable,'SAD_SULFC') || strcmp(variable,'Temperature')
%    names = {'CCMI','MAM1990','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
else
    names = {'CCMI','MAM1990','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
    %names = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
end

%file locations
directory = [commondir,variable,'/',dailydir];
files = dir([directory,variable,'*']);

TAdirectory = [commondir,'Temperature/',dailydir];
TAfiles = dir([TAdirectory,'TA*']);

if prescoordinfile
    
    Z3directory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/Z3/',dailydir];
    Z3files = dir([Z3directory,'Z3*']);

    datedir = ['/Users/kanestone/work/projects/WACCM/netcdffiles/date/',dailydir];
    datefiles = dir([datedir,'date*']);
    
else
    HCdirectory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/hybrid_ht_conversion/',dailydir];
    HCfiles = dir([HCdirectory,'hybrid*']);

    PSdirectory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/PS/',dailydir];
    PSfiles = dir([PSdirectory,'PS*']);

    Z3directory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/Z3/',dailydir];
    Z3files = dir([Z3directory,'Z3*']);

    datedir = ['/Users/kanestone/work/projects/WACCM/netcdffiles/date/',dailydir];
    datefiles = dir([datedir,'date*']);
end

for i = 1:length(files)
    %read in O3
    [info.(names{i}), data.(names{i}), attributes.(names{i})] = ...
        Read_in_netcdf([directory,files(i).name]);
    
    if strcmp(variable,'NONO2')
        %data.(names{i}).NONO2 = permute(data.(names{i}).NONO2,[3,2,1]);
    end
    
    %read in temperature
    [TAinfo.(names{i}), TAdata.(names{i}), TAattributes.(names{i})] = ...
        Read_in_netcdf([TAdirectory,TAfiles(i).name]);
    
    if prescoordinfile
        HCdata.(names{i}) = data.(names{i});
        PSdata.(names{i}).PS = data.(names{i}).PS;
        HCdata.(names{i}).hybm = HCdata.(names{i}).b;
        %read in geopotenial height        
    else
        %read in hybrid_ht conversion parameters
        [HCinfo.(names{i}), HCdata.(names{i}), HCattributes.(names{i})] = ...
            Read_in_netcdf([HCdirectory,HCfiles(i).name]);

        %read in surface pressure
        [PSinfo.(names{i}), PSdata.(names{i}), PSattributes.(names{i})] = ...
            Read_in_netcdf([PSdirectory,PSfiles(i).name]);
    end
    %read in geopotenial height
    [Z3info.(names{i}), Z3data.(names{i}), Z3attributes.(names{i})] = ...
        Read_in_netcdf([Z3directory,Z3files(i).name]);

    %read in dates
    [dateinfo.(names{i}), datedata.(names{i}), dateattributes.(names{i})] = ...
        Read_in_netcdf([datedir,datefiles(i).name]);
        
    %calculating pressure
    if strcmp(variable,'Extinction')
        %Pressure.(names{i}) = zeros(size(data.(names{i}).EXTINCTdn));    
        Pressure.(names{i}) = zeros(size(data.(names{i}).EXTxASYMdn));    
    else
        Pressure.(names{i}) = zeros(size(data.(names{i}).(variable)));    
    end
    
    %Broken here
    if strcmp(variable,'Extinction')
        %numberdims = ndims(data.(names{i}).EXTINCTdn);
        numberdims = ndims(data.(names{i}).EXTxASYMdn);
    else
        numberdims = ndims(data.(names{i}).(variable));
    end
    
    if numberdims > 3
        permute_coords = [2,3,1,4];
        PSrepcoords = [1,1,1];
        PScoords = [1,2,4,3];
    else
        permute_coords = [2,1,3];
        PSrepcoords = [1,1];
        PScoords = [1,3,2];
    end
    if prescoordinfile
        hyamp = permute(repmat(HCdata.(names{i}).ap * 100000,[1,size(PSdata.(names{i}).PS)]),permute_coords);
    else
        hyamp = permute(repmat(HCdata.(names{i}).hyam .* HCdata.(names{i}).P0,[1,size(PSdata.(names{i}).PS)]),permute_coords);
    end
    hybm = permute(repmat(HCdata.(names{i}).hybm,[1,size(PSdata.(names{i}).PS)]),permute_coords);
    %P0 = permute(repmat(HCdata.(names{i}).P0,[size(HCdata.(names{i}).hyam,1),size(PSdata.(names{i}).PS)]),permute_coords);
    PS = double(permute(repmat(PSdata.(names{i}).PS,[PSrepcoords,size(HCdata.(names{i}).hybm,1)]),PScoords));
    Pressure.(names{i}) = hyamp + hybm .* PS;        
    
    Altitude.(names{i}) = -(log(Pressure.(names{i})./101325)*8.31.*273.15)./(9.806*.0289);        

    % converting from geopotential height to geometric height
    GMH.(names{i}) = (Z3data.(names{i}).Z3*6356000)./(6356000-Z3data.(names{i}).Z3);
    
    %converting to number density
    if strcmp(variable,'SAD_SULFC') || strcmp(variable,'HNO3_STS') || strcmp(variable,'HNO3_GAS') ...
            || strcmp(variable,'HNO3_NAT') || strcmp(variable,'HNO3_TOTAL') || strcmp(variable,'Extinction') ...
            || strcmp(variable,'SAD_ICE')
        NumberDensity.(names{i}) = [];
    else
        if ND
            NumberDensity.(names{i}) = vmr2conc(data.(names{i}).(variable),TAdata.(names{i}).T,Pressure.(names{i}),variable,'conc');
        else
            NumberDensity.(names{i}) = []
        end
    end
    
    if strcmp(variable,'Extinction')
        %MolConc.(names{i}) = data.(names{i}).EXTINCTdn;
        MolConc.(names{i}) = data.(names{i}).EXTxASYMdn;
    else
        MolConc.(names{i}) = data.(names{i}).(variable);
    end
    
    Latitudes = data.(names{i}).lat;  
    if numberdims > 3
        Longitudes = data.(names{i}).lon;    
    end
end

end