% EESC curve
%Calculate EESC

clear all

%% constants
timeseries = 1765:2500;
alpha = 65; %For polar latitudes
air = 28.97; %(g mol-1)

en = 0;
plot_EESC = 0;
cbrew = cbrewer('qual','Set1',10);
cbrew2 = cbrewer('qual','Set1',17);

%% read in WMO 2011 file
filedir = '/Volumes/ExternalOne/work/data/rcp6.0_table_wmo2011';
fid = fopen(filedir); 
while ~en
    line = fgetl(fid);
    if strcmp(line(1:10),'      Year')
        en = 1;
    end
end
data = fscanf(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',[21 Inf])';
ODS_data = data(1:end,5:21);

%% Extract EESC components
% The equation for each species is as follows:
% (no. atoms)*MMR*(air/species)*(fractional halogen release)*(Br scaling factor if needed)
% Fractional halon release and Br scaling factor obtained from: 
% Newman et al., (2007) https://doi.org/10.5194/acp-7-4537-2007

ODS_data_EESC(:,1) = 3*ODS_data(:,1)*(air./137.37)*.99;         %CFC-11                 (CCl3F) *3
ODS_data_EESC(:,2) = 2*ODS_data(:,2)*(air./120.91)*.86;         %CFC-12                 (CCl2F2) *2
ODS_data_EESC(:,3) = 3*ODS_data(:,3)*(air./187.37)*.90;         %CFC-113                (CCl2FCClF2) *3 
ODS_data_EESC(:,4) = 2*ODS_data(:,4)*(air./170.92)*.40;         %CFC-114                (CClF2CClF2)
ODS_data_EESC(:,5) = 1*ODS_data(:,5)*(air./154.47)*.15;         %CFC-115                (CClF2CF3)
ODS_data_EESC(:,6) = 4*ODS_data(:,6)*(air./153.82)*1;           %Carbon tetrachlorine   (CCl4)
ODS_data_EESC(:,7) = 3*ODS_data(:,7)*(air./133.4)*.99;          %Methyl chloroform      (CH3CCl3)
ODS_data_EESC(:,8) = 1*ODS_data(:,8)*(air./86.47)*.41;          %HCFC-22                (CHClF3)
ODS_data_EESC(:,9) = 2*ODS_data(:,9)*(air./116.95)*.90;         %HCFC-141b              (CH3CCl2F)
ODS_data_EESC(:,10) = 1*ODS_data(:,10)*(air./100.5)*.29;        %HCFC-142b              (CH3CClF2)
ODS_data_EESC(:,11) = 1*ODS_data(:,11)*(air./165.36)*1;         %Halon-1211             (CBrClF2) %for Cl
ODS_data_EESC(:,12) = 1*ODS_data(:,11)*(air./165.36)*1*alpha;   %Halon-1211             (CBrClF2) %for Br
ODS_data_EESC(:,13) = 2*ODS_data(:,12)*(air./209.82)*1*alpha;   %Halon-1202             (CBr2F2)
ODS_data_EESC(:,14) = 1*ODS_data(:,13)*(air./148.91)*.8*alpha;  %Halon-1301             (CBrF3)
ODS_data_EESC(:,15) = 2*ODS_data(:,14)*(air./259.8)*1*alpha;    %Halon-2402             (CBrF2CBrF2)
ODS_data_EESC(:,16) = 1*ODS_data(:,15)*(air./94.94)*.99*alpha;  %Methyl bromide         (CH3Br)
ODS_data_EESC(:,17) = 1*ODS_data(:,16)*(air./50.49)*.91;        %Methyl chloride        (CH3Cl)
%ODS_data_EESC(:,18) = 1*ODS_data(:,17)*.1;                     %HFC134A                (CH2FCF3)
%+((.39e-9*(air./165.36))*.99); %(VSL Bromine); 
EESC = sum(ODS_data_EESC,2)*1e9; 
EESCnomc = sum(ODS_data_EESC(:,[1:6,8:end]),2)*1e9; 
EESC_recent = EESC(timeseries >= 1955 & timeseries <= 2100);
EESCnomc_recent = EESCnomc(timeseries >= 1955 & timeseries <= 2100);

EESC_11_12_113 = sum(ODS_data_EESC(:,1:3),2)*1e9; 
EESC_11_12_113 = EESC_11_12_113(timeseries >= 1955 & timeseries <= 2100);

EESC_Cl= sum(ODS_data_EESC(:,[1,2,3,6,7,8,17]),2)*1e9; 
EESC_Cl= EESC_Cl(timeseries >= 1955 & timeseries <= 2100);

timeperiod = 1954:2099;
value1980 = EESC_11_12_113(timeperiod == 1975);
valuebelow = EESC_11_12_113 < value1980;
belowindex = find(valuebelow == 1);
tp = timeperiod(123+5);

yearval = 1987;
yearvaltoplot = yearval+5;
value1990 = EESC_recent(timeperiod == yearval);
valuebelow = EESC_recent < value1990;
belowindex = find(valuebelow == 1);
tpind = find(belowindex > yearval-1954+5);
tpEESC = timeperiod(belowindex(tpind(1))+5);

value1990nomc = EESCnomc_recent(timeperiod == yearval);
valuebelow = EESCnomc_recent < value1990nomc;
belowindex = find(valuebelow == 1);
tpindnomc = find(belowindex > yearval-1954+5);
tpEESCnomc = timeperiod(belowindex(tpindnomc(1))+5);

%% take cumulative

EESCcs = cumsum(ODS_data_EESC,2)*1e9;
EESCcs = EESCcs(timeseries >= 1955 & timeseries <= 2100,:);

%% plotting
figure;
set(gcf,'color','white','position',[100 100 1000 700]);
%ph(1) = plot(1955:2100,EESC_recent,'LineWidth',3,'color',cbrew(1,:));
hold on
ph(1) = plot(1960:2105,EESC_recent,'LineWidth',3,'color',cbrew(1,:));


ph(2) = plot(1960:2105,EESCnomc_recent,'LineWidth',3,'color',cbrew(2,:));
plot([yearvaltoplot+.5,yearvaltoplot+.5],[0 5],'k--');
plot([tpEESC+.5,tpEESC+.5],[0 value1990],'--','color',cbrew(1,:));
plot([tpEESCnomc+.5,tpEESCnomc+.5],[0 value1990nomc],'--','color',cbrew(2,:));
plot([yearvaltoplot+.5,tpEESC+.5],[value1990,value1990],'--','color',cbrew(1,:));
plot([yearvaltoplot+.5,tpEESCnomc+.5],[value1990nomc,value1990nomc],'--','color',cbrew(2,:));
set(gca,'xtick',1950.5:10:2109.5,'xticklabel',1950:10:2110,'fontsize',18);

xlabel('Year','fontsize',20);
ylabel('ppbv','fontsize',20);
title('EESC','fontsize',22);
lh = legend('South polar EESC','South polar EESC - no CH3CCl3');
set(lh,'box','off','fontsize',18);
box on

%%
figure;
set(gcf,'color','white','position',[100 100 1000 700]);
ph2(1) = plot(1960:2105,EESC_11_12_113,'LineWidth',3,'color',cbrew(3,:));
hold on
plot([1980.5,1980.5],[0 2.2],'k--');
plot([tp+.5,tp+.5],[0 2.2],'k--');
plot([1954.5,2099.5],[value1980,value1980],'k--');
set(gca,'xtick',1950.5:10:2109.5,'xticklabel',1950:10:2110,'fontsize',18);

xlabel('Year','fontsize',20);
ylabel('ppbv','fontsize',20);
title('EESC: CFC-11, CFC-12, CFC-113','fontsize',22);
lh = legend('Polar EESC');
set(lh,'box','off','fontsize',18);

figure;
set(gcf,'color','white','position',[100 100 1000 700]);
ph3(1) = plot(1960:2105,ODS_data_EESC(timeseries >= 1955 & timeseries <= 2100,7)*1e9,'LineWidth',3,'color',cbrew(3,:));

set(gca,'xtick',1950.5:10:2109.5,'xticklabel',1950:10:2110,'fontsize',18);
xlim([1940 2120]);

xlabel('Year','fontsize',20);
ylabel('ppbv','fontsize',20);
title('EESC: Methyl Chloroform','fontsize',22);
lh = legend('Polar Methyl Chloroform');
set(lh,'box','off','fontsize',18);

%%

figure;
set(gcf,'color','white','position',[100 100 1000 700]);
ph3(1) = plot(1960:2105,ODS_data_EESC(timeseries >= 1955 & timeseries <= 2100,2)*1e9,'LineWidth',3,'color',cbrew(3,:));

set(gca,'xtick',1950.5:10:2109.5,'xticklabel',1950:10:2110,'fontsize',18);
xlim([1940 2120]);

xlabel('Year','fontsize',20);
ylabel('ppbv','fontsize',20);
title('EESC: Methyl Chloroform','fontsize',22);
lh = legend('Polar Methyl Chloroform');
set(lh,'box','off','fontsize',18);

%% plot cumulative
% figure;
% set(gcf,'color','white','position',[100 100 1000 700]);
% 
% for i = 1:size(EESCcs,2)
%     ph(i) = plot(1954.5:2099.5,EESCcs(:,i),'LineWidth',3,'color',cbrew(i,:));
%     hold on
% end
