% Read in total Cl

clear variables

[~,highcl,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/predruns/CLOY/highcl/ensave_zm_CLOY_b.e11.BWTREFC2.f19_g16.ccmi34.HighCl.nc');
[~,lowcl,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/predruns/CLOY/lowcl/ensave_zm_CLOY_b.e11.BWTREFC2.f19_g16.ccmi34.LowCl.nc');

highcl_mean = nanmean(highcl.CLOY,3)*1e9;
lowcl_mean = nanmean(lowcl.CLOY,3)*1e9;

%% plotting

lwidth = 3;
fsize = 18;
createfig('medium','on');

plot(highcl.lat,highcl_mean(:,27),'LineWidth',lwidth);
set(gca,'fontsize',fsize,'xtick',[-90:15:90])
ylabel('ppb (highCl)','fontsize',fsize+2);
xlabel('latitude','fontsize',fsize+2);
hold on
yyaxis right
plot(highcl.lat,lowcl_mean(:,27),'LineWidth',lwidth);
xlim([-95 95]);
lh = legend('1995-2024 zonal mean (HighCl)','1955-1979 zonal mean (LowCl)');
set(lh,'fontsize',fsize','box','off')
ylabel('ppb (lowCl)','fontsize',fsize+2);
title('Cly at 1 hPa','fontsize',fsize+2)

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/Cly_highcl_and_lowcl','.pdf',];
export_fig(filename,'-pdf');
