% looking at NOX
clear all
[~,data,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/NOX/NOX_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.1999.fsst.nl.004f.cam.h0zm_merged_c160407.nc');
[~,SPE,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/SPE/spes_1963-2014_c150717.nc');
[~,temp,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/reactionrates/OddOx_lossrates_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.1999.fsst.nl.004f.cam.h0zm_merged_c160407.nc');
HOx = temp.OddOx_HOx_Loss;
SPEde = SPE.Prod(:,SPE.date >= 20000101 & SPE.date <= 20141201);
SPEde_date = SPE.date(SPE.date >= 20000101 & SPE.date <= 20141201);

temp = num2str(SPEde_date);
for j = 1:length(temp)
    mons(j) = str2double(temp(j,5:6));
    years(j) = str2double(temp(j,1:4));
end
clearvars temp


%% import NO2

commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[~, MolConc, ~, Pressure, ~, ~, Latitudes, Longitudes, datedata] = ReadWACCMvertical('NO2','monthly',commondir,0,1);

MolConc = rmfield(MolConc, 'MAM1990');
datedata = rmfield(datedata, 'MAM1990');
Pressure = rmfield(Pressure, 'MAM1990');
runs = fields(MolConc);
runs = {runs{2},runs{6}};

NO2chemonly = MolConc.(runs{2});

load('/Volumes/MyBook/work/data/regressionOutput/NO2forregression.mat');

NO2anom = squeeze(SDanomalllats_nointerp(2,:,:,:));
NO2anom = permute(NO2anom,[2,1,3]);

%% deseasonalise NO2 based on 2006-2010
for i = 1:12
    NO2anom2(:,:,i,:) = NO2chemonly(:,:,12+i:12:end-11) - nanmean(NO2chemonly(:,:,108+i:12:156),3);
end
NO2anom_reshape = NO2anom2(:,:,:);
NO2datedata = datedata.Chemonlynoleap.date(13:end-11);
%%
count = 1;
for i = 2000:2014
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

SPEextract = SPE_monthsum(:,:,22);
SPEextract = SPEextract';
SPEextract = SPEextract(:);

%% average over polar regions.
lats = [-90 -60;60 90];
pres = 2;
[~,presind] = min(abs(data.lev - pres));
for i = 1:size(lats,1)
    latind = data.lat >= lats(i,1) & data.lat <= lats(i,2);    
    data_wa_NOx(i,:) = weightedaverage(squeeze(data.NOX(latind,presind,:)),data.lat(latind)); %weighted average lat mean, pressure, time    
    data_wa_NO2(i,:) = weightedaverage(squeeze(NO2anom(latind,presind,:)),data.lat(latind)); %weighted average lat mean, pressure, time    
    data_wa_NO2_2(i,:) = weightedaverage(squeeze(NO2anom_reshape(latind,presind,:)),data.lat(latind)); %weighted average lat mean, pressure, time    
end

%%
cbrew = cbrewer('qual','Set1',10);
fsize = 18;
SPEdedateind = 1:179/5449:180;
SPEdedateind = SPEdedateind(1:end-1);
%SPEde (SPEde <= 22) = NaN;
%SPEextract (SPEextract <= 10) = NaN;
%datetoplot = [20000101,200
lwidth = 3;
createfig('medium','on')
yyaxis right
%plot(SPEextract,'LineWidth',3);
%plot(SPEextract,'LineWidth',2,'color','k','LineStyle','--');
plot(SPEdedateind,SPEde(22,:),'LineWidth',2,'color','k','LineStyle','--');
set(gca, 'YScale', 'log')
set(gca,'fontsize',fsize);
ylim([-1000 4200]);
ylabel('Ion pair production rate (/cm3/sec)','fontsize',fsize+2);
set(gca,'YColor','k')
hold on
yyaxis left
ylim([-5e-9 21e-9])
plot(data_wa_NO2_2(1,:)','LineWidth',3,'color',cbrew(1,:),'LineStyle','-')
hold on
plot(data_wa_NO2_2(2,:)','LineWidth',3,'color',cbrew(2,:),'LineStyle','-')
ylabel('NO2 anomaly from 2008-2011 (mol/mol)','fontsize',fsize+2);
set(gca,'xtick',1:24:180,'xticklabel',2000:2:2014);
set(gca,'YColor','k')
title('Comparison of solar proton events and polar NO_2 at 2 hPa','fontsize',fsize+4);
lh = legend('South Polar NO_2 anomaly (60-90{\circ}S)','North Polar NO_2 anomaly (60-90{\circ}N)','Ion pair production rate');
set(lh,'location','NorthEast','box','off');
xlabel('Year','fontsize',fsize+2);

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','SPENO2comparisonlog'];
 
export_fig(filename,'-pdf');

%%

lat = 10;
mon = 7;
SPElagmean = nanmean(squeeze(SPE_monthmean(:,mon-5:mon,22)),2);
%SPElagmean2 = squeeze(SPE_monthmean(:,mon,25));

figure
plot((squeeze(data.NOX(lat,32,mon:12:end))-nanmean(squeeze(data.NOX(lat,32,mon:12:end))))./std(squeeze(data.NOX(lat,32,mon:12:end))),'k');
hold on
plot((SPElagmean - nanmean(SPElagmean))./std(SPElagmean));
%plot((SPElagmean2 - nanmean(SPElagmean2))./std(SPElagmean2));
plot((squeeze(HOx(lat,34,mon:12:end))-nanmean(squeeze(HOx(lat,34,mon:12:end))))./std(squeeze(HOx(lat,34,mon:12:end))),'--');
