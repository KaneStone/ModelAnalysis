% WACCM doy and daily time series plots mimicking ozonesonde locations
clear all
%% Reading in WACCM;
year_extract = '2015';
variable = 'O3';
%waclat = -90;
%pres = 150;
wacyear = str2double(year_extract);
wactime = [str2double([num2str(wacyear),'0101']) str2double([num2str(wacyear),'1231'])];
wactime2 = [str2double([num2str(wacyear),'0201']) str2double([num2str(wacyear+1),'0101'])];
pres = 130;
wacpres = pres*100;
%import WACCM sulfc
commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[WACCM_ND, WACCM_MolConc, ~, WACCM_Pressure, ~, ~, WACCM_Latitudes, WACCMdatedata] = ReadWACCMvertical('O3','monthly',commondir,0);
if strcmp(variable,'O3');
    [~,WACCM_daily_data,~] = Read_in_netcdf(...
        '/Volumes/My Book for Mac/work/data/WACCM/O3/daily/O3_f.e11.FWTREFC1SD.f19.f19.ccmi30.mam.004f.cam.h5_43_226_merged.nc');
    [~,WACCM_daily_data_VC,~] = Read_in_netcdf(...
        '/Volumes/My Book for Mac/work/data/WACCM/O3/daily/O3_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.004f.cam.h5_43_226_merged.nc');
    factor = 1e6;
elseif strcmp(variable,'T');    
    [~,WACCM_daily_data,~] = Read_in_netcdf(...
        '/Volumes/My Book for Mac/work/data/WACCM/Temperature/daily/TA_f.e11.FWTREFC1SD.f19.f19.ccmi30.mam.004f.cam.h5_merged_43_226hPa.nc');   
    factor = 1;
end
WACCM_daily_data_VC = WACCM_daily_data;
[~,wacpresindex] = min(abs(wacpres-squeeze(WACCM_Pressure.MAM(1,:,1))));

wacdays = [365,366,365,365,365,366,365,365,365,366,365,365,365,366,365,365,365];

waclat = [-90,-70.7,-69,-68.6,-64.2]; % South Pole, Neumayer, Syowa, Davis, Marambio
waclon = [158,351.7,39.6,78,303.4];

%waclat = [-73,-73,-73,-73,-73,-73]; % South Pole, Neumayer, Syowa, Davis, Marambio
%waclon = [0,60,120,180,240,300];

for i = 1:length(waclat)
    [~,waclatindex(i)] = min(abs(waclat(i)-WACCM_daily_data.lat));
    [~,waclonindex(i)] = min(abs(waclon(i)-WACCM_daily_data.lon));
    daily(i,:,:) = squeeze(WACCM_daily_data.(variable)(waclonindex(i),waclatindex(i),:,:));
    daily_VC(i,:,:) = squeeze(WACCM_daily_data_VC.(variable)(waclonindex(i),waclatindex(i),:,:));
end

daily_zm = squeeze(nanmean(WACCM_daily_data.(variable)));
daily2011 = daily_zm(waclatindex,WACCM_daily_data.date >= wactime(1) & WACCM_daily_data.date <= wactime(2));
%monthly2011 = squeeze(WACCM_MolConc.MAM(waclatindex,wacpresindex, WACCMdatedata.MAM.date >= wactime2(1) & WACCMdatedata.MAM.date <= wactime2(2)));

%dailyzm_at_lat = squeeze(daily(waclatindex,:,:));

%%
for j = 1:length(waclat)
    
    for i = 1:size(WACCM_daily_data.(variable),4)    
        inddate = num2str(WACCM_daily_data.date(i));
        indyear(i) = str2double(inddate(1:4));    
        indmonth(i) = str2double(inddate(5:6));
        dailyzm_at_pres(j,i) = interp1(log(WACCM_daily_data.lev),squeeze(daily(j,:,i)),log(pres))*factor;
        %dailyzm_at_pres(j,i) = interp1(log(WACCM_daily_data.lev),squeeze(daily(j,:,i)),log(pres))*factor;
        dailyzm_at_pres_VC(j,i) = interp1(log(WACCM_daily_data_VC.lev),squeeze(daily_VC(j,:,i)),log(pres))*factor;
    end
    yearin = min(indyear);
    doyWACCM(j).a = zeros(max(indyear)-min(indyear)+1,366);
    dailyzm_at_pres_out(j).a = zeros(max(indyear)-min(indyear)+1,366);
    dailyzm_at_pres_out_VC(j).a = zeros(max(indyear)-min(indyear)+1,366);
    for i = 1:max(indyear)-min(indyear)+1;
        dailyzm_at_pres_out(j).a(i,1:length(indyear(indyear == yearin))) = dailyzm_at_pres(j,indyear == yearin);
        dailyzm_at_pres_out_VC(j).a(i,1:length(indyear(indyear == yearin))) = dailyzm_at_pres_VC(j,indyear == yearin);
        doyWACCM(j).a(i,1:length(indyear(indyear == yearin))) = 1:length(indyear(indyear == yearin));
        yearin = yearin+1;
    end        
end

for j = 1:length(waclat)
    for i = 1:12
        xyear(j).a(i).m = indyear(indmonth == i);
        ozonemonth(j).a(i).m = dailyzm_at_pres(j,indmonth == i);
        ozonemonth_VC(j).a(i).m = dailyzm_at_pres_VC(j,indmonth == i);
    end
end

% plotting separate time series and interpolation methods

figure;
plot(daily(j,:,i))

%stations = {'South Pole','Syowa','Davis'};
stations = {'South Pole','Neumayer','Syowa','Davis','Marambio'};
pressures = [150,155];
startyear = 1999;
numyears = 17;
yearstohighlight = 2015;
fsize = 22;
wac = 1;
months = 11;
combine = 1;

if strcmp(variable,'O3')
    save(['/Users/kanestone/work/projects/WACCM/output/','WACCMdaily_forozonesondes_',num2str(pres),'.mat'],'xyear','ozonemonth',...
        'ozonemonth_VC','doyWACCM','dailyzm_at_pres_out','dailyzm_at_pres_out_VC');
    %save(['/Users/kanestone/work/projects/WACCM/output/','WACCMdaily_forstandardlat_',num2str(pres),'.mat'],'xyear','ozonemonth',...
    %    'ozonemonth_VC','doyWACCM','dailyzm_at_pres_out','dailyzm_at_pres_out_VC');
else strcmp(variable,'T')
    save(['/Users/kanestone/work/projects/WACCM/output/','WACCMdaily_forozonesondes_Temperature_',num2str(pres),'.mat'],'xyear','ozonemonth',...
        'ozonemonth_VC','doyWACCM','dailyzm_at_pres_out','dailyzm_at_pres_out_VC');
end

%plot_ozonesonde_timeseries(months,stations,pressures,xyear,...
%    ozonemonth,startyear,waclat,yearstohighlight,numyears,fsize,'MAM',combine)
%plot_ozonesonde_timeseries(months,stations,pressures,xyear,...
%    ozonemonth_VC,startyear,waclat,yearstohighlight,numyears,fsize,'VCMAM',combine)
