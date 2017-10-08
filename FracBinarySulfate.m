% determining fraction of binary sulfate


[STS_ND, STS_MolConc, ~, STS_Pressure, ~, ~, STS_Latitudes, STS_datedata] = ...
    ReadWACCMvertical('HNO3_STS','monthly');
[HNO3GAS_ND, HNO3GAS_MolConc, ~, HNO3GAS_Pressure, ~, ~, HNO3GAS_Latitudes, HNO3GAS_datedata] = ...
    ReadWACCMvertical('HNO3_GAS','monthly');
[HNO3TOTAL_ND, HNO3TOTAL_MolConc, ~, HNO3TOTAL_Pressure, ~, ~, HNO3TOTAL_Latitudes, HNO3TOTAL_datedata] = ...
    ReadWACCMvertical('HNO3_TOTAL','monthly');
[HNO3NAT_ND, HNO3NAT_MolConc, ~, HNO3NAT_Pressure, ~, ~, HNO3NAT_Latitudes, HNO3NAT_datedata] = ...
    ReadWACCMvertical('HNO3_NAT','monthly');

lats = [-90 -30];
time = [20150201 20160101];
yearprint = num2str(time(1));
yearprint = yearprint(1:4);
latind = find(STS_Latitudes > lats(1) & STS_Latitudes < lats(2));
timeind = find(STS_datedata.MAM.date >= time(1) & STS_datedata.MAM.date <= time(2));

STS_date_lat = STS_MolConc.MAM(latind,:,timeind);
STS_date_lat_pressure = STS_Pressure.MAM(latind,:,timeind);
GAS_date_lat = HNO3GAS_MolConc.MAM(latind,:,timeind);
GAS_date_lat_pressure = HNO3GAS_Pressure.MAM(latind,:,timeind);
TOTAL_date_lat = HNO3TOTAL_MolConc.MAM(latind,:,timeind);
TOTAL_date_lat_pressure = HNO3TOTAL_Pressure.MAM(latind,:,timeind);
NAT_date_lat = HNO3NAT_MolConc.MAM(latind,:,timeind);
NAT_date_lat_pressure = HNO3NAT_Pressure.MAM(latind,:,timeind);

STS_date_lat_VC = STS_MolConc.VCMAM(latind,:,timeind);
STS_date_lat_pressure_VC = STS_Pressure.VCMAM(latind,:,timeind);
GAS_date_lat_VC = HNO3GAS_MolConc.VCMAM(latind,:,timeind);
GAS_date_lat_pressure_VC = HNO3GAS_Pressure.VCMAM(latind,:,timeind);
TOTAL_date_lat_VC = HNO3TOTAL_MolConc.VCMAM(latind,:,timeind);
TOTAL_date_lat_pressure_VC = HNO3TOTAL_Pressure.VCMAM(latind,:,timeind);
NAT_date_lat_VC = HNO3NAT_MolConc.VCMAM(latind,:,timeind);
NAT_date_lat_pressure_VC = HNO3NAT_Pressure.VCMAM(latind,:,timeind);

%% interp onto regular pressure levels

pres = [1000:-10:110 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.1:.1];

preslog = log(pres);

prestick = [1000 500 250 100 50 25 10 5 2.5 1 .5 .25 .1 .05 .025 .01 .005 .0025 ...
    .001 .0005 .00025 .0001];

presticklog = log(prestick);

for i = 1:size(STS_date_lat_pressure,1)
    for j = 1:size(STS_date_lat_pressure,3)
        STS_date_lat_interp(i,:,j) = interp1(log(squeeze(STS_date_lat_pressure(i,:,j))./100),...
            squeeze(STS_date_lat(i,:,j)),log(pres));
        STS_date_lat_interp_VC(i,:,j) = interp1(log(squeeze(STS_date_lat_pressure_VC(i,:,j))./100),...
            squeeze(STS_date_lat_VC(i,:,j)),log(pres));
        GAS_date_lat_interp(i,:,j) = interp1(log(squeeze(GAS_date_lat_pressure(i,:,j))./100),...
            squeeze(GAS_date_lat(i,:,j)),log(pres));
        GAS_date_lat_interp_VC(i,:,j) = interp1(log(squeeze(GAS_date_lat_pressure_VC(i,:,j))./100),...
            squeeze(GAS_date_lat_VC(i,:,j)),log(pres));
        NAT_date_lat_interp(i,:,j) = interp1(log(squeeze(NAT_date_lat_pressure(i,:,j))./100),...
            squeeze(NAT_date_lat(i,:,j)),log(pres));
        NAT_date_lat_interp_VC(i,:,j) = interp1(log(squeeze(NAT_date_lat_pressure_VC(i,:,j))./100),...
            squeeze(NAT_date_lat_VC(i,:,j)),log(pres));
        TOTAL_date_lat_interp(i,:,j) = interp1(log(squeeze(TOTAL_date_lat_pressure(i,:,j))./100),...
            squeeze(TOTAL_date_lat(i,:,j)),log(pres));
        TOTAL_date_lat_interp_VC(i,:,j) = interp1(log(squeeze(TOTAL_date_lat_pressure_VC(i,:,j))./100),...
            squeeze(TOTAL_date_lat_VC(i,:,j)),log(pres));
    end
end

GASminSTS = GAS_date_lat_interp - STS_date_lat_interp;
GASminSTS_VC = GAS_date_lat_interp_VC - STS_date_lat_interp_VC;

STSGASratio = STS_date_lat_interp./GAS_date_lat_interp;
STSGASratio_VC = STS_date_lat_interp_VC./GAS_date_lat_interp_VC;

monthstitle = {'January','February','March','April','May','June','July','August',...
    'September','October','November','December'};

%% plotting vertical profiles
fsize = 18;
lattoplot = [1,3,5,7,9,11];
month = 6;
lwidth = 2;
createfig('large');
for i = 1:length(lattoplot)
    sp(i) = subplot(3,2,i);
    MAM = plot(GAS_date_lat_interp(lattoplot(i),:,month),preslog,'LineWidth',lwidth,'color','k');
    hold on
    GAS = plot(GAS_date_lat_interp(lattoplot(i),:,month),preslog,'LineWidth',lwidth);
    STS = plot(STS_date_lat_interp(lattoplot(i),:,month),preslog,'LineWidth',lwidth);
    NAT = plot(NAT_date_lat_interp(lattoplot(i),:,month),preslog,'LineWidth',lwidth);
    TOTAL = plot(TOTAL_date_lat_interp(lattoplot(i),:,month),preslog,'LineWidth',lwidth);
    ax = gca;
    ax.ColorOrderIndex = 1;
    VCMAM = plot(GAS_date_lat_interp_VC(lattoplot(i),:,month),preslog,'--','LineWidth',lwidth,'color','k');
    plot(GAS_date_lat_interp_VC(lattoplot(i),:,month),preslog,'--','LineWidth',lwidth)
    plot(STS_date_lat_interp_VC(lattoplot(i),:,month),preslog,'--','LineWidth',lwidth)
    plot(NAT_date_lat_interp_VC(lattoplot(i),:,month),preslog,'--','LineWidth',lwidth)
    plot(TOTAL_date_lat_interp_VC(lattoplot(i),:,month),preslog,'--','LineWidth',lwidth)

    set(gca,'Ydir','reverse');
    set(gca,'ytick',fliplr(presticklog),'YtickLabel',fliplr(prestick),'fontsize',fsize)
    if mod(i,2)
        ylabel('Pressure (hPa)','fontsize',fsize+2);
    end
    if i >= length(lattoplot)-1
        xlabel('Mol/mol','fontsize',fsize+2);
    end
    title([sprintf('%.1f',abs(HNO3GAS_Latitudes(lattoplot(i)))),'S'],'fontsize',fsize+4)
    
    ylim([log(1) log(1000)])
    if i == length(lattoplot)
        lh = legend([MAM, VCMAM, GAS, STS, NAT, TOTAL],...
            'MAM','VC-MAM','Gas phase','STS condensed','NAT condensed','Total');
        set(lh,'fontsize',fsize+2,'box','off','orientation','horizontal','position',[.5 .01 .01 .01]);
    end
    
end
pmtit = mtit([monthstitle{month},' 2015 HNO3'],'xoff',-.00,'yoff',.04,'fontsize',fsize+10);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMHNO3/',yearprint,'_',monthstitle{month},'_profile.pdf'];
export_fig(filename,'-pdf')

%% colorbar

cbrew = cbrewer('seq','YlOrBr',10);

cmin = 0;
cmax = 10e-9;
interval = 1e-9;
STScmin = 0;
STScmax = 2.2e-9;
STSinterval = .2e-9;
ratmin = 0;
ratmax = 1;
ratint = .1;

%% plotting MAM

fsize = 18;
month = 11;

createfig('large');
sp = subplot(3,2,1);
%contourf(STS_date_lat(:,:,month)
contourf(STS_Latitudes(latind),preslog,STS_date_lat_interp(:,:,month)',STScmin:STSinterval:STScmax);
set(gca,'Ydir','reverse');
set(gca,'ytick',fliplr(presticklog),'YtickLabel',fliplr(prestick),'fontsize',fsize)
ylabel('Pressure (hPa)','fontsize',fsize+2);
%xlabel('Latitude','fontsize',fsize+2);
title('MAM STS condensed HNO3','fontsize',fsize+4);
ch = colorbar;
cbaxloc = get(ch,'Position');
suboutpos = get(sp,'OuterPosition');            
colormap(cbrew);  
caxis([STScmin STScmax]);

subplot(3,2,3);
contourf(HNO3GAS_Latitudes(latind),preslog,GAS_date_lat_interp(:,:,month)',[cmin:interval:cmax]);
set(gca,'Ydir','reverse');
set(gca,'ytick',fliplr(presticklog),'YtickLabel',fliplr(prestick),'fontsize',fsize)
ylabel('Pressure (hPa)','fontsize',fsize+2);
%xlabel('Latitude','fontsize',fsize+2);
title('MAM Gas phase HNO3','fontsize',fsize+4)
ch = colorbar;
cbaxloc = get(ch,'Position');
suboutpos = get(sp,'OuterPosition');            
colormap(cbrew);  
caxis([cmin cmax]);

subplot(3,2,5);
%contourf(HNO3GAS_Latitudes(latind),preslog,GASminSTS(:,:,month)',[cmin:interval:cmax]);
contourf(HNO3GAS_Latitudes(latind),preslog,STSGASratio(:,:,month)',[ratmin:ratint:ratmax]);
set(gca,'Ydir','reverse');
set(gca,'ytick',fliplr(presticklog),'YtickLabel',fliplr(prestick),'fontsize',fsize)
ylabel('Pressure (hPa)','fontsize',fsize+2);
xlabel('Latitude','fontsize',fsize+2);
title('MAM STS/Gas','fontsize',fsize+4)
ch = colorbar;
cbaxloc = get(ch,'Position');
suboutpos = get(sp,'OuterPosition');            
colormap(cbrew);  
caxis([ratmin ratmax]);

% Begin VC MAM--------

subplot(3,2,2);
contourf(STS_Latitudes(latind),preslog,STS_date_lat_interp_VC(:,:,month)',[STScmin:STSinterval:STScmax])
set(gca,'Ydir','reverse');
set(gca,'ytick',fliplr(presticklog),'YtickLabel',fliplr(prestick),'fontsize',fsize)
%ylabel('Pressure (hPa)','fontsize',fsize+2);
%xlabel('Latitude','fontsize',fsize+2);
title('VC-MAM STS condensed HNO3','fontsize',fsize+4)
ch = colorbar;
cbaxloc = get(ch,'Position');
suboutpos = get(sp,'OuterPosition');            
colormap(cbrew);  
caxis([STScmin STScmax]);

subplot(3,2,4);
contourf(HNO3GAS_Latitudes(latind),preslog,GAS_date_lat_interp_VC(:,:,month)',[cmin:interval:cmax])
set(gca,'Ydir','reverse');
set(gca,'ytick',fliplr(presticklog),'YtickLabel',fliplr(prestick),'fontsize',fsize)
%ylabel('Pressure (hPa)','fontsize',fsize+2);
%xlabel('Latitude','fontsize',fsize+2);
title('VC-MAM Gas phase HNO3','fontsize',fsize+4)
ch = colorbar;
cbaxloc = get(ch,'Position');
suboutpos = get(sp,'OuterPosition');            
colormap(cbrew);  
caxis([cmin cmax]);

subplot(3,2,6);
%contourf(HNO3GAS_Latitudes(latind),preslog,GASminSTS_VC(:,:,month)')
contourf(HNO3GAS_Latitudes(latind),preslog,STSGASratio_VC(:,:,month)',[ratmin:ratint:ratmax])
set(gca,'Ydir','reverse');
set(gca,'ytick',fliplr(presticklog),'YtickLabel',fliplr(prestick),'fontsize',fsize)
%ylabel('Pressure (hPa)','fontsize',fsize+2);
xlabel('Latitude','fontsize',fsize+2);
title('VC-MAM STS/Gas','fontsize',fsize+4)
ch = colorbar;
cbaxloc = get(ch,'Position');
suboutpos = get(sp,'OuterPosition');            
colormap(cbrew);  
caxis([ratmin ratmax]);

mtit([monthstitle{month}, ' 2015 HNO3'],'xoff',-.00,'yoff',.04,'fontsize',fsize+10);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMHNO3/',yearprint,'_',monthstitle{month},'_LatPres.png'];
export_fig(filename,'-png')
