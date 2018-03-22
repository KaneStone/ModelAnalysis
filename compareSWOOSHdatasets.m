% Compare SWOOSH anomalies

Stimeperiod = [2000 2016];

[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedanomfillo3q_swoosh-v02.6-198401-201712-latpress-2.5deg-L31.nc');
sfields = fieldnames(SWOOSH);

[~,SWOOSH2,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedeqfillo3q_swoosh-v02.6-198401-201712-latpress-2.5deg-L31.nc');
sfields2 = fieldnames(SWOOSH2);

%% create anomalies
lat = 62;
S1anom = (squeeze(SWOOSH.combinedanomfillo3q(lat,27,12:12:end)) - nanmean(squeeze(SWOOSH.combinedanomfillo3q(lat,27,12:12:end))))./...
    nanmean(squeeze(SWOOSH.combinedanomfillo3q(lat,27,12:12:end)))*100;

S2anom = (squeeze(SWOOSH2.combinedeqfillo3q(lat,27,12:12:end)) - nanmean(squeeze(SWOOSH2.combinedeqfillo3q(lat,27,12:12:end))))./...
    nanmean(squeeze(SWOOSH2.combinedeqfillo3q(lat,27,12:12:end)))*100;

%% plotting
lw = 3;
fsize = 18;
createfig('medium','on');

ph = plot(1984:2017,S1anom,'LineWidth',lw);
hold on
ph2 = plot(1984:2017,S2anom,'LineWidth',lw);
set(gca,'fontsize',fsize);
xlabel('Year','fontsize',fsize+2);
ylabel('O3 anomaly (%)','fontsize',fsize+2);
title(['SWOOSH dataset comparison at ',num2str(SWOOSH.lat(lat)),'{\circ}','N and ',num2str(SWOOSH.level(27)),' hPa'],'fontsize',fsize+4);

lh = legend('SWOOSH anomaly fill','SWOOSH eq. lat fill');
set(lh,'location','SouthEast','box','off');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/SWOOSHdatasetcomparison'];
export_fig(filename,'-pdf');