% chem only gas phase trends, while removing SPE and solar cycle.

%trends go from 2000-2014

%% import WACCM data
commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[NumberDensity, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes, Longitudes,datedata] = ReadWACCMvertical('O3','monthly',commondir,0);
[~,TCO,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/TCO/TOZ_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.1999.fsst.nl.004f_merged_c160407.nc');

MolConc = rmfield(MolConc, 'MAM1990');
datedata = rmfield(datedata, 'MAM1990');
Pressure = rmfield(Pressure, 'MAM1990');
runs = fields(MolConc);

%% import SWOOSH
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201611-latpress-10deg-L31.nc');
%SWOOSHlatindex = find(SWOOSH.lat > lats(1) & SWOOSH.lat < lats(2));
SWOOSH.combinedo3q = SWOOSH.combinedo3q(:,:,1:end-11);
for i = 1:12
    SWOOSH.montharrange(:,:,:,i) = SWOOSH.combinedo3q(:,:,i:12:end);
end
SWOOSHyears = 1984:2015;
%SWOOSHtimeindex = find(SWOOSHyears >= compareyears(1) & SWOOSHyears <= compareyears(2));

%% Read in Solar

[~,solara,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/solar/f107a/f107a_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.1999.004f.cam.h0zm_merged_c160105.nc');
[~,solar,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/solar/f107/f107_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.1999.004f.cam.h0zm_merged_c160105.nc');

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

%%
count = 1;
for i = 1999:2014
    for j = 1:12
        SPE_monthmean(count,j,:) = nanmean(SPEde(:,years == i & mons == j),2);
        SPE_monthmedian(count,j,:) = nanmedian(SPEde(:,years == i & mons == j),2);
        SPE_monthsum(count,j,:) = sum(SPEde(:,years == i & mons == j),2);
    end
    SPE_yearmean(count,:) = nanmean(SPEde(:,years == i),2);
    SPE_yearsum(count,:) = sum(SPEde(:,years == i),2);
    count = count+1;
end
clearvars count

%% remove solar and SPE
chemonly = MolConc.Chemonlynoleap;
MAM = MolConc.MAM;
% for i = 1:size(chemonly
% [b,bint,~,~,stats] = regress(

% plot(squeeze(chemonly(81,34,12:12:end)))
% hold on
% plot(squeeze(MAM(81,34,12:12:end)))
%chemonly = cat(3,chemonly,NaN);
solara.f107a = [solara.f107a;NaN];
solar.f107 = [solar.f107;NaN];
%% solarmonths
for i = 1:12
    solarmonths(:,i) = solar.f107(i:12:end-12);
    solarmonthsa(:,i) = solara.f107a(i:12:end-12);
    chemonlymontharrange(:,:,:,i) = chemonly(:,:,i:12:end-11);
end

%% plot WACCM Latitude Pressure trends

latPresTrends(MolConc.Chemonlynoleap,solarmonths,SPE_monthmean,Pressure.Chemonlynoleap,datedata.Chemonlynoleap.date,SWOOSH.level,SPE.pressure);


%% testing
%totest = squeeze(chemonlymontharrange(81,34,4:end,12));
totest = squeeze(chemonlymontharrange(10,34,4:end,5));
totestSWOOSH = squeeze(SWOOSH.montharrange(2,27,SWOOSHyears >= 2002 & SWOOSHyears <= 2014,5));
% normalizing 

totest_n = (totest - nanmean(totest))./std(totest);
totestSWOOSH_n = (totestSWOOSH - nanmean(totestSWOOSH))./std(totestSWOOSH);
%SPE_totest_n = (SPE_monthsum(:,11,22) - nanmean(SPE_monthsum(:,11,22)))./std(SPE_monthsum(:,11,22));
%SPE_totest_n = (nanmean(SPE_monthmedian(:,9:11,23),2) - nanmean(nanmean(SPE_monthmedian(:,9:11,23),2)))./std(nanmean(SPE_monthmedian(:,9:11,23),2));
%SPE_totest_n = (nanmean(SPE_monthsum(:,9:11,23),2) - nanmean(nanmean(SPE_monthsum(:,9:11,23),2)))./std(nanmean(SPE_monthsum(:,9:11,23),2));
%SPE_totest_n = (nanmean(SPE_monthmean(4:end,9:11,22),2) - nanmean(nanmean(SPE_monthmean(4:end,9:11,22),2)))./std(nanmean(SPE_monthmean(4:end,9:11,22),2));
SPE_totest_n = (nanmean(SPE_monthmean(4:end,2:4,22),2) - nanmean(nanmean(SPE_monthmean(4:end,2:4,22),2)))./std(nanmean(SPE_monthmean(4:end,2:4,22),2));
%SPE_totest_n2 = (SPE_monthsum(:,10,22) - nanmean(SPE_monthsum(:,10,22)))./std(SPE_monthsum(:,10,22));
%SPE_totest_n3 = (SPE_monthsum(:,9,22) - nanmean(SPE_monthsum(:,9,22)))./std(SPE_monthsum(:,9,22));
%SPE_totest_n = (SPE_yearsum(:,22) - nanmean(SPE_yearsum(:,22)))./std(SPE_yearsum(:,22));
solar_totest_n = (nanmean(solarmonths(4:end,10),2) - nanmean(nanmean(solarmonths(4:end,10),2)))./std(nanmean(solarmonths(4:end,10),2));

%b = regress(totest_n,[ones(length(totest_n),1),SPE_totest_n,SPE_totest_n2,SPE_totest_n3,solar_totest_n,]); % solar.f107(12:12:end)
b = regress(totest_n,[ones(length(totest_n),1),SPE_totest_n,solar_totest_n,[1:length(totest_n)]']); % solar.f107(12:12:end)
%b = regress(totest_n,[ones(length(totest_n),1),SPE_totest_n]); % solar.f107(12:12:end)
bSWOOSH = regress(totestSWOOSH_n,[ones(length(totestSWOOSH_n),1),SPE_totest_n,solar_totest_n,[1:length(totestSWOOSH_n)]']); % solar.f107(12:12:end)
%totest_srm = totest_n - b(2)*SPE_totest_n - b(3)*SPE_totest_n2 - b(4)*SPE_totest_n3 - b(5)*solar_totest_n;
totest_srm = totest_n - b(2)*SPE_totest_n - b(3)*solar_totest_n;
totestSWOOSH_srm = totestSWOOSH_n - bSWOOSH(2)*SPE_totest_n - bSWOOSH(3)*solar_totest_n;

%% plotting

cbrew = cbrewer('qual','Paired',12);
createfig('medium','on');
hold on
plot(SPE_totest_n,'color',cbrew(6,:),'LineWidth',2)
plot(solar_totest_n,'color',cbrew(8,:),'LineWidth',2)

plot(totest_n,'color',cbrew(1,:),'LineWidth',2);
plot(totestSWOOSH_n,'color',cbrew(3,:),'LineWidth',2);

plot(totest_srm,'color',cbrew(2,:),'LineWidth',4);
plot(totestSWOOSH_srm,'color',cbrew(4,:),'LineWidth',4);
legend

