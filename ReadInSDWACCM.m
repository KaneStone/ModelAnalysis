function [SDWaccmData,Latitudes,Longitudes,intpres,Altitude,ND,Pressure] = ReadInSDWACCM(intpres,timeperiod)

%% Read in WACCM O3 data
commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[ND, MolConc, ~, Pressure, Altitude, ~, Latitudes, Longitudes, datedata] = ReadWACCMvertical('O3','monthly',commondir,0,1);

MolConc = rmfield(MolConc, 'MAM1990');
datedata = rmfield(datedata, 'MAM1990');
Pressure = rmfield(Pressure, 'MAM1990');
runs = fields(MolConc);
runs = {runs{2},runs{6}};

for i = 1:2
    CCMIy = CCMI_years(datedata.(runs{i}).date,1);
    dateind = CCMIy >= timeperiod(1) & CCMIy <= timeperiod(2);
    dateindHF = CCMIy >= timeperiod(1) & CCMIy <= timeperiod(2);
    for j = 1:size(MolConc.Chemonlynoleap,1)
        [SDWaccmData(i).O3(:,j,:),~] = intRegPres(squeeze(MolConc.(runs{i})(j,:,dateind)),...
            squeeze(Pressure.(runs{i})(j,:,dateind))./100,intpres);
    end
    SDWaccmData(i).O3nointerp = squeeze(MolConc.(runs{i})(:,:,dateind));
end

%% Read in U wind if needed

[~,U,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/U/wa_-5_5_U_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.mam.004f.cam.h0zm_merged_c160105.nc');
lev = [10,50];
[~,levind] = min(abs(U.lev - lev));
SDWaccmData(1).U = squeeze(U.U(levind,dateind));
SDWaccmData(2).U = repmat(squeeze(U.U(levind,1:12)),[1,timeperiod(2) - timeperiod(1)+1]);

%% Read in U wind if needed

[~,HFS,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/VTH3d/wa_-80_-40_VTH3d_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.mam.004f.cam.h0zm_merged_c160105.nc');
[~,HFN,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/VTH3d/wa_40_80_VTH3d_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.mam.004f.cam.h0zm_merged_c160105.nc');
lev = 50;
[~,levind] = min(abs(HFS.ilev - lev));
SDWaccmData(1).HFS = squeeze(HFS.VTH3d(levind,dateindHF));
SDWaccmData(1).HFN = squeeze(HFN.VTH3d(levind,dateindHF));
SDWaccmData(2).HFS = repmat(SDWaccmData(1).HFS(1,1:12),[1,size(SDWaccmData(1).HFS,2)./12]); % up to here.
SDWaccmData(2).HFN = repmat(SDWaccmData(1).HFN(1,1:12),[1,size(SDWaccmData(1).HFN,2)./12]); % up to here.

%% Read in andconstruct ENSO if needed
[~,TS,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/TS/merged_TS_f.e11.FWTREFC1SD.f19.f19.ccmi30.mam.004f.cam.h0.1999-2014.nc');

SDWaccmData(1).TS = squeeze(TS.TS(:,:,dateind));

latlimits = [-5 5];
lonlimits = [190 240];

latindex = TS.lat >= latlimits(1) & TS.lat <= latlimits(2);
lonindex = TS.lon >= lonlimits(1) & TS.lon <= lonlimits(2);


for j = 1:12
    NINO_mn(:,j,:) = squeeze(nanmean(SDWaccmData(1).TS(lonindex,latindex,j:12:end),1));
    NINO_mn2(j,:) = squeeze(nanmean(NINO_mn(:,j,:),1));        
end    
NINO_mn3 = NINO_mn2(:);
NINOallmean = nanmean(NINO_mn3);
NINOstd = std(NINO_mn3,1);
SDWaccmData(1).NINO34 = detrend((NINO_mn3 - NINOallmean)./NINOstd);

%% Read in Solar

[~,solar,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/solar/f107a/f107a_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.1999.004f.cam.h0zm_merged_c160105.nc');
SDWaccmData(2).solar = solar.f107a(dateind);

%% Read in SPE

[~,SPE,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/SPE/spes_1963-2014_c150717.nc');

SPEde = SPE.Prod(:,SPE.date >= 19990101 & SPE.date <= 20160101);
SPEde_date = SPE.date(SPE.date >= 19990101 & SPE.date <= 20160101);

temp = num2str(SPEde_date);
for j = 1:length(temp)
    mons(j) = str2double(temp(j,5:6));
    years(j) = str2double(temp(j,1:4));
end
clearvars temp

count = 1;
for i = timeperiod(1)-1:timeperiod(2)
    for j = 1:12
        SPE_monthmean(count,j,:) = nanmean(SPEde(:,years == i & mons == j),2);        
    end    
    count = count+1;
end
SDWaccmData(2).SPE = permute(SPE_monthmean,[3,2,1]);
SDWaccmData(2).SPE = SDWaccmData(2).SPE(:,:);
clearvars count

[SDWaccmData(2).SPEintpres,~] = intRegPres(SDWaccmData(2).SPE,...
            repmat(SPE.pressure,[1,size(SDWaccmData(2).SPE,2)]),intpres);

%% Read in NO2 regression function

load('/Volumes/MyBook/work/data/regressionOutput/NO2forregression.mat');

SDWaccmData(1).NO2 = squeeze(SDanomalllats(1,:,:,:));
SDWaccmData(2).NO2 = squeeze(SDanomalllats(2,:,:,:));
SDWaccmData(1).NO2nointerp = squeeze(SDanomalllats_nointerp(1,:,:,:));
SDWaccmData(2).NO2nointerp = squeeze(SDanomalllats_nointerp(2,:,:,:));
        
end