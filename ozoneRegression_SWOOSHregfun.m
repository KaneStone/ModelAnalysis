function SWOOSHRegfun = ozoneRegression_SWOOSHregfun(timeperiod)

       %% Read in Multi-variate ENSO index

       MEI = importdata('/Volumes/MyBook/work/data/SWOOSH/SWOOSHregfunctions/MEI.dat');
       MEI_dateindex = MEI.data(:,1) >= timeperiod(1) & MEI.data(:,1) <= timeperiod(2);
       MEIdata = MEI.data(MEI_dateindex,2:end);
       MEIdata = permute(MEIdata,[2,1]);
       SWOOSHRegfun.MEIdata = MEIdata(:);
   
    %% Read in QBO proxies (30 and 10 hPa winds)
    
    syears = [1984:2017]; 
    singapore10 = importdata('/Volumes/MyBook/work/data/SWOOSH/SWOOSHregfunctions/qbo10.dat');
    singapore10 = reshape(singapore10,[12,length(singapore10)./12])';    
    singapore30 = importdata('/Volumes/MyBook/work/data/SWOOSH/SWOOSHregfunctions/qbo30.dat');
    singapore30 = reshape(singapore30,[12,length(singapore30)./12])';
    s_dateindex = syears >= timeperiod(1) & syears <= timeperiod(2);
    singapore10data = singapore10(s_dateindex,:);
    singapore30data = singapore30(s_dateindex,:);
    singapore10data = permute(singapore10data,[2,1]);
    singapore30data = permute(singapore30data,[2,1]);
    SWOOSHRegfun.singaporedata(1,:) = singapore10data(:);
    SWOOSHRegfun.singaporedata(2,:) = singapore30data(:);

    %% import solar cycle
    
    solar = importdata('/Volumes/MyBook/work/data/SWOOSH/SWOOSHregfunctions/solarcycle.dat');
    solar_dateindex = solar.data(:,1) >= timeperiod(1) & solar.data(:,1) <= timeperiod(2);
    SWOOSHRegfun.solardata = solar.data(solar_dateindex,4);
    
    %% Read in Aerosol
    aerosol = ncread('/Volumes/MyBook/work/data/SWOOSH/SWOOSHregfunctions/tau_map_2012-12.nc','tau');
    %[~,aerosol,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/SWOOSHregfunctions/tau_map_2012-12.nc');
    aerdate = repmat(1850:2012,[12,1]);
    aerdate = aerdate(:);
    aerdateindex = aerdate >= timeperiod(1) & aerdate <= timeperiod(2);
    SWOOSHRegfun.aerosols = nanmean(aerosol(:,aerdateindex),1);
    SWOOSHRegfun.aerosols = [SWOOSHRegfun.aerosols,ones(1,length(SWOOSHRegfun.solardata) - length(SWOOSHRegfun.aerosols))*SWOOSHRegfun.aerosols(1)]';
    SWOOSHRegfun.aerosols (SWOOSHRegfun.aerosols == 0) = SWOOSHRegfun.aerosols(1); 
        
end


