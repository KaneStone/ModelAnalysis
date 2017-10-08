%% Solar cycle analysis

%% Read in ozone data from specified dynamics runs

commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';

[NumberDensity, MolConc, TAdata, Pressure, Altitude, GMH, Latitudes, Longitudes, datedata] = ...
    ReadWACCMvertical('O3','monthly',commondir,0);

%% Put chemonly onto regular pressure

[cofissts_regpres,pres] = presinterpWACCM(Pressure.Chemonlynoleap,NumberDensity.Chemonlynoleap);


%%

plot(squeeze(cofissts_regpres(10,39,7:12:end)));
%figure;
%plot(log(pres),squeeze(cofissts_regpres(10,:,u)));

