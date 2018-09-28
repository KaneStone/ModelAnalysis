% Compare SWOOSH anomalies

Stimeperiod = [2000 2016];

[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedanomfillo3q_swoosh-v02.6-198401-201712-latpress-2.5deg-L31.nc');
sfields = fieldnames(SWOOSH);

[~,SWOOSH2,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedeqfillo3q_swoosh-v02.6-198401-201712-latpress-2.5deg-L31.nc');

sfields2 = fieldnames(SWOOSH2);

[~,SWOOSH3,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201712-latpress-2.5deg-L31.nc');

sfields3 = fieldnames(SWOOSH3);

%% create anomalies
lat = -75;
[~,latind] = min(abs(SWOOSH.lat - lat));
pres = 30;
[~,presind] = min(abs(SWOOSH.level - pres));
S1anom = (squeeze(SWOOSH.combinedanomfillo3q(latind,presind,:)) - nanmean(squeeze(SWOOSH.combinedanomfillo3q(latind,presind,:))))./...
    nanmean(squeeze(SWOOSH.combinedanomfillo3q(latind,presind,:)))*100;

S2anom = (squeeze(SWOOSH2.combinedeqfillo3q(latind,presind,:)) - nanmean(squeeze(SWOOSH2.combinedeqfillo3q(latind,presind,:))))./...
    nanmean(squeeze(SWOOSH2.combinedeqfillo3q(latind,presind,:)))*100;

S3anom = (squeeze(SWOOSH3.combinedo3q(latind,presind,:)) - nanmean(squeeze(SWOOSH3.combinedo3q(latind,presind,:))))./...
    nanmean(squeeze(SWOOSH3.combinedo3q(latind,presind,:)))*100;

%% plotting
lw = 3;
fsize = 18;
createfig('medium','on');

ph = plot(1:length(S1anom),S1anom,'LineWidth',lw);
hold on
ph2 = plot(1:length(S2anom),S2anom,'LineWidth',lw);
ph3 = plot(1:length(S3anom),S3anom,'LineWidth',lw);
set(gca,'fontsize',fsize);
xlabel('Year','fontsize',fsize+2);
ylabel('O3 anomaly (%)','fontsize',fsize+2);
title(['SWOOSH dataset comparison at ',num2str(SWOOSH.lat(latind)),'{\circ}','N and ',num2str(SWOOSH.level(presind)),' hPa'],'fontsize',fsize+4);
set(gca,'xtick',1:60:500,'xticklabel',1984:5:2020);
lh = legend('SWOOSH anomaly fill','SWOOSH eq. lat fill','SWOOSH raw');
set(lh,'location','SouthEast','box','off');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/SWOOSHdatasetcomparison_',num2str(SWOOSH.lat(latind)),'N','_',num2str(SWOOSH.level(presind)),'hPa'];
export_fig(filename,'-pdf');