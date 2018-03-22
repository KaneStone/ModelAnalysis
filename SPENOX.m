% looking at NOX
clear all
[~,data,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/NOX/NOX_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.1999.fsst.nl.004f.cam.h0zm_merged_c160407.nc');
[~,SPE,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/SPE/spes_1963-2014_c150717.nc');
[~,temp,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/reactionrates/OddOx_lossrates_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.1999.fsst.nl.004f.cam.h0zm_merged_c160407.nc');
HOx = temp.OddOx_HOx_Loss;
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
