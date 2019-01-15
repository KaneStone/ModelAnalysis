%% Looking at Methyl Chloroform

%Calculate EESC
clear all
%constants
timeseries = [1765:2500];
alpha = 65; %For polar latitudes
air = 28.97; %(g mol-1)
en = 0;
plot_EESC = 0;
plot_regression = 1;
cbrew = cbrewer('qual','Set1',5);
%read in WMO 2011 file
fid = fopen('/Volumes/ExternalOne/work/data/rcp6.0_table_wmo2011'); 
while ~en
    line = fgetl(fid);
    if strcmp(line(1:10),'      Year')
        en = 1;
    end
end
data = fscanf(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',[21 Inf])';

ODS_data = data(1:end,5:21);
%ODS_data_EESC = ODS_data;

%calculating EESC using Newman (2007) 5.5 mean age of air fractional values

%The equation for each species is as follows:
%(no. atoms)*MMR*(air/species)*(fractional halogen release)*(Br scaling factor if needed)

ODS_data_EESC(:,1) = 3*ODS_data(:,1)*(air./137.37)*.99;         %CFC-11                 (CCl3F)
ODS_data_EESC(:,2) = 2*ODS_data(:,2)*(air./120.91)*.86;         %CFC-12                 (CCl2F2)
ODS_data_EESC(:,3) = 3*ODS_data(:,3)*(air./187.37)*.90;         %CFC-113                (CCl2FCClF2) 
ODS_data_EESC(:,4) = 2*ODS_data(:,4)*(air./170.92)*.40;         %CFC-114                (CClF2CClF2)
ODS_data_EESC(:,5) = 1*ODS_data(:,5)*(air./154.47)*.15;         %CFC-115                (CClF2CF3)
ODS_data_EESC(:,6) = 4*ODS_data(:,6)*(air./153.82)*1;           %Carbon tetrachlorine   (CCl4)
ODS_data_EESC(:,7) = 3*ODS_data(:,7)*(air./133.4)*.99;          %Methyl chloroform      (CH3CCl3)
ODS_data_EESC(:,8) = 1*ODS_data(:,8)*(air./86.47)*.41;          %HCFC-22                (CHClF3)
ODS_data_EESC(:,9) = 2*ODS_data(:,9)*(air./116.95)*.90;         %HCFC-141b              (CH3CCl2F)
ODS_data_EESC(:,10) = 1*ODS_data(:,10)*(air./100.5)*.29;        %HCFC-142b              (CH3CClF2)
ODS_data_EESC(:,11) = 1*ODS_data(:,11)*(air./165.36)*1;         %Halon-1211             (CBrClF2)
ODS_data_EESC(:,12) = 1*ODS_data(:,11)*(air./165.36)*1*alpha;   %Halon-1211             (CBrClF2)
ODS_data_EESC(:,13) = 2*ODS_data(:,12)*(air./209.82)*1*alpha;   %Halon-1202             (CBr2F2)
ODS_data_EESC(:,14) = 1*ODS_data(:,13)*(air./148.91)*.8*alpha;  %Halon-1301             (CBrF3)
ODS_data_EESC(:,15) = 2*ODS_data(:,14)*(air./259.8)*1*alpha;    %Halon-2402             (CBrF2CBrF2)
ODS_data_EESC(:,16) = 1*ODS_data(:,15)*(air./94.94)*.99*alpha;  %Methyl bromide         (CH3Br)
ODS_data_EESC(:,17) = 1*ODS_data(:,16)*(air./50.49)*.91;        %Methyl chloride        (CH3Cl)
%ODS_data_EESC(:,18) = 1*ODS_data(:,17)*.1;                     %HFC134A                (CH2FCF3)

EESC = sum(ODS_data_EESC,2);%+((.39e-9*(air./165.36))*.99); %(VSL Bromine);
EESC_wmcf = sum(ODS_data_EESC(:,[1:6,8:17]),2);%+((.39e-9*(air./165.36))*.99); %(VSL Bromine);
EESC_wmcf2 = sum(ODS_data_EESC(:,[1:6,11:17]),2);%+((.39e-9*(air./165.36))*.99); %(VSL Bromine);

EESC_recent = EESC(timeseries >= 1955 & timeseries <= 2100);
EESC_wmcf_recent = EESC_wmcf(timeseries >= 1955 & timeseries <= 2100); 
EESC_wmcf_recent2 = EESC_wmcf2(timeseries >= 1955 & timeseries <= 2100); 

gradient = diff(EESC_recent);
gradient2 = diff(EESC_wmcf_recent);
gradient3 = diff(EESC_wmcf_recent2);

%% find turnaround periods

period = 10;

[test,maxlocEESC] = max(EESC_recent);
[test,maxlocEESC2] = max(EESC_wmcf_recent);
[test,maxlocEESC3] = max(EESC_wmcf_recent2);

diffallpast = (EESC_recent(maxlocEESC-period)-EESC_recent(maxlocEESC))*1e9;
diffallfuture = (EESC_recent(maxlocEESC+period)-EESC_recent(maxlocEESC))*1e9;

diffallpast2 = (EESC_wmcf_recent(maxlocEESC-period)-EESC_wmcf_recent(maxlocEESC))*1e9;
diffallfuture2 = (EESC_wmcf_recent(maxlocEESC+period)-EESC_wmcf_recent(maxlocEESC))*1e9;

diffallpast3 = (EESC_wmcf_recent2(maxlocEESC-period)-EESC_wmcf_recent2(maxlocEESC))*1e9;
diffallfuture3 = (EESC_wmcf_recent2(maxlocEESC+period)-EESC_wmcf_recent2(maxlocEESC))*1e9;

%% difference in every 10 year period
period2 = 10;
for i = 1:14
    EESC_recent_change(i) = (EESC_recent(period2*i+1) - EESC_recent(period2*(i-1)+1))*1e9;
end

lbl = {'70-60','80-70','90-80','00-90','10-00','20-10','30-20','40-30','50-40','60-50','70-60','80-70','90-80','00-90'};

plot(EESC_recent_change)
set(gca,'xtick',1:1:14,'xticklabel',lbl);

%% plotting
cbrew = cbrewer('qual','Set1',10);

createfig('large','on')
subplot(2,1,1)
ph = plot(1954.5:2099.5,EESC_recent,'LineWidth',3,'color',cbrew(1,:));
hold on
phw = plot(1954.5:2099.5,EESC_wmcf_recent,'LineWidth',3,'color',cbrew(2,:));
phw2 = plot(1954.5:2099.5,EESC_wmcf_recent2,'LineWidth',3,'color',cbrew(3,:));
plot(repmat(1954.5+maxlocEESC-1,1,2),[1e-9,5e-9],'LineStyle','--','LineWidth',2,'color',cbrew(1,:))
plot(repmat(1954.5+maxlocEESC2-1,1,2),[1e-9,5e-9],'LineStyle','--','LineWidth',2,'color',cbrew(2,:))
plot(repmat(1954.5+maxlocEESC3-1,1,2),[1e-9,5e-9],'LineStyle','--','LineWidth',2,'color',cbrew(3,:))

set(gca,'xtick',1954.5:10:2099.5,'xticklabel',1960:10:2100,'fontsize',18);
xlabel('Year','fontsize',20);
ylabel('ppbv','fontsize',20);
title('EESC','fontsize',22);
lh = legend('EESC','EESC without CH3CCl3','EESC without CH3CCl3 and HCFCs');
set(lh,'box','off','fontsize',18);

subplot(2,1,2)
ph = plot(1954.5:2099.5,EESC_recent,'LineWidth',3,'color',cbrew(1,:));
hold on
phw = plot(1954.5:2099.5,EESC_wmcf_recent,'LineWidth',3,'color',cbrew(2,:));
phw2 = plot(1954.5:2099.5,EESC_wmcf_recent2,'LineWidth',3,'color',cbrew(3,:));
plot(repmat(1954.5+maxlocEESC-1,1,2),[1e-9,5e-9],'LineStyle','--','LineWidth',2,'color',cbrew(1,:))
plot(repmat(1954.5+maxlocEESC2-1,1,2),[1e-9,5e-9],'LineStyle','--','LineWidth',2,'color',cbrew(2,:))
plot(repmat(1954.5+maxlocEESC3-1,1,2),[1e-9,5e-9],'LineStyle','--','LineWidth',2,'color',cbrew(3,:))

plot(repmat(1954.5+maxlocEESC-1-period,1,2),[1e-9,5e-9],'LineStyle','--','LineWidth',2,'color','k')
plot(repmat(1954.5+maxlocEESC-1+period,1,2),[1e-9,5e-9],'LineStyle','--','LineWidth',2,'color','k')

set(gca,'xtick',1955:5:2100,'xticklabel',1960:5:2100,'fontsize',18);
xlabel('Year','fontsize',20);
ylabel('ppbv','fontsize',20);
title('EESC','fontsize',22);
% lh = legend('EESC','EESC without CH3CCl3','EESC without CH3CCl3 and HCFCs');
% set(lh,'box','off','fontsize',18,'location','southeast');
xlim([1968 2018])

annotation('textbox',[.2 .35 .1 .1],'String',[num2str(diffallpast./(period/10)),' ppbv/decade'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',14,...
        'EdgeColor','none','fontweight','bold','color',cbrew(1,:))
annotation('textbox',[.375 .35 .1 .1],'String',[num2str(diffallfuture./(period/10)),' ppbv/decade'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',14,...
    'EdgeColor','none','fontweight','bold','color',cbrew(1,:))

annotation('textbox',[.2 .325 .1 .1],'String',[num2str(diffallpast2./(period/10)),' ppbv/decade'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',14,...
        'EdgeColor','none','fontweight','bold','color',cbrew(2,:))
annotation('textbox',[.375 .325 .1 .1],'String',[num2str(diffallfuture2./(period/10)),' ppbv/decade'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',14,...
    'EdgeColor','none','fontweight','bold','color',cbrew(2,:))

annotation('textbox',[.2 .3 .1 .1],'String',[num2str(diffallpast3./(period/10)),' ppbv/decade'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',14,...
        'EdgeColor','none','fontweight','bold','color',cbrew(3,:))
annotation('textbox',[.375 .3 .1 .1],'String',[num2str(diffallfuture3./(period/10)),' ppbv/decade'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',14,...
    'EdgeColor','none','fontweight','bold','color',cbrew(3,:))

annotation('textbox',[.2 .375 .1 .1],'String','1988-1998','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',14,...
        'EdgeColor','none','fontweight','bold','color','k')
annotation('textbox',[.375 .375 .1 .1],'String','1998-2008','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',14,...
    'EdgeColor','none','fontweight','bold','color','k')

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/EESCwork/',num2str(period),'_methylchloroform.pdf'];

export_fig(filename,'-pdf');
    
