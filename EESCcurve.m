% EESC curve
%Calculate EESC

clear all

%% constants
timeseries = 1765:2500;
alpha = 65; %For polar latitudes
air = 28.97; %(g mol-1)
en = 0;
plot_EESC = 0;
%cbrew = cbrewer('qual','Set1',10);

%% read in WMO 2011 file
filedir = '/Volumes/MyBook/work/data/rcp6.0_table_wmo2011';
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
%+((.39e-9*(air./165.36))*.99); %(VSL Bromine); 
EESC = sum(ODS_data_EESC,2)*1e9; 
EESC_recent = EESC(timeseries >= 1955 & timeseries <= 2100);

%% plotting
figure;
set(gcf,'color','white','position',[100 100 1000 700]);
ph(1) = plot(1954.5:2099.5,EESC_recent,'LineWidth',3,'color',[0.8941 0.1020 0.1098]); %cbrew(1,:));
hold on
ph(2) = plot(1960:2105,EESC_recent,'LineWidth',3,'color',[0.2157 0.4941 0.7216]); %cbrew(2,:));
set(gca,'xtick',1950.5:10:2109.5,'xticklabel',1950:10:2110,'fontsize',18);

xlabel('Year','fontsize',20);
ylabel('ppbv','fontsize',20);
title('EESC','fontsize',22);
lh = legend('Tropical EESC','South polar EESC');
set(lh,'box','off','fontsize',18);

