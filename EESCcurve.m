% EESC curve
%Calculate EESC
clear all

%% constants

alpha = 65; %For polar latitudes
en = 0;
plot_EESC = 0;
cbrew = cbrewer('qual','Set1',10);
cbrew2 = cbrewer('qual','Set1',17);
air = 28.97; %(g mol-1)
SP2013 = 1;
%% read in SPARC 2013 file
if SP2013
    filedir = '/Volumes/ExternalOne/work/data/SPARC2013.txt';
    fid = fopen(filedir); 
    while ~en
        line = fgetl(fid);
        if strcmp(line(1:4),'Year')
            en = 1;
        end
    end
    ODS_data2 = fscanf(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',[17 Inf])';
    ODSyears2 = ODS_data2(:,1);
    ODSyears = 1850:ODSyears2(end);
    ODS_data2 = ODS_data2(:,2:end)*1e-12; % put from ppt to mixing ratio
    ODS_data = [repmat(ODS_data2(1,:),100,1);ODS_data2]; %replicate 1st year for inverse gaussian function

else
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
    ODS_data = data(1:end,5:20);
    ODSyears = data(:,1)-.5;

    ODS_data = ODS_data(86:356,:);
    ODSyears = ODSyears(86:356); 

    conversionfactor = [(air./137.37),(air./120.91),(air./187.37),(air./170.92),(air./154.47),...
        (air./153.82),(air./133.4),(air./86.47),(air./116.95),(air./100.5),...
        (air./165.36),(air./209.82),(air./148.91),(air./259.8),(air./94.94),(air./50.49)];

    ODS_data = ODS_data.*conversionfactor;
end

%% Mean release times
%Mean release times for mean age of air of 5.5 years from Table 2 of Engel et al. 2018 (https://www.atmos-chem-phys.net/18/601/2018/acp-18-601-2018.pdf)

mrtd(1) = 5.5;     %CFC-11 
mrtd(2) = 5.9;     %CFC-12                
mrtd(3) = 5.8;     %CFC-113             
mrtd(4) = 8.3;     %CFC-114               
mrtd(5) = 10.1;     %CFC-115              
mrtd(6) = 5.5;    %Carbon tetrachlorine
mrtd(7) = 5.6;     %Methyl chloroform
mrtd(8) = 7.0;     %HCFC-22 
mrtd(9) = 5.8;     %HCFC-141b
mrtd(10) = 6.5;     %HCFC-142b
mrtd(11) = 5.5;     %Halon-1211
mrtd(12) = 5.5;     %Halon-1211
mrtd(13) = 6.2;     %Halon-1202
mrtd(14) = 5.5;     %Halon-1301
mrtd(15) = 5.5;     %Halon-2402
mrtd(16) = 5.8;     %Methyl bromide

%% Construct inverse Gaussian function
%constructing the age spectrum as an inverse gaussian function following Eq. 9 from Waugh and Hall, (2002) (https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2000RG000101)
% instead of age of air, we use the mean relese times for each molecule
% separately.

% Here, I am using Engel et al.'s delta caculation of lambda = delta^2/mrtd = .7
% this implies delta = sqrt(aoa*.7)

% Newman et al. uses the delta definition of half the aoa. This is also
% used by WMO 2018. In that case delta = 5.5/2 for Newman and mrtd/2 for
% WMO (2018) (although they were not explicit, they may have also used
% 5.5/2)

inversegaussianfunction = zeros(length(ODSyears),size(mrtd,2));
for j = 1:size(mrtd,2)
    %aoa = mrtd(j);
    aoa = 5.5;
    %delta = sqrt(.7*aoa);
    %delta = mrtd(j)/2;
    delta = aoa/2;
    
    count = 1;
    for t = 0:1:length(ODSyears)-1
        inversegaussianfunction(count,j) = 1/(2*delta*sqrt(pi*(t/aoa)^3))*exp(-aoa^2*(t/aoa-1)^2./(4*delta^2*(t/aoa)));
        count = count+1;
    end
end
createfig('medium','on')
hold on
plot(inversegaussianfunction)
set(gca,'xtick',1:1:171,'xticklabel',0:1:171);
title('Age spectrum','fontsize',22)
ylabel('Fraction of mixing ratio','fontsize',20)
xlabel('Transit time from entering the stratosphere (years)','fontsize',20)
xlim([0 20]);
grid on


test = inversegaussianfunction(:,1);

%testweightedmean = test/

%% Extract EESC components
% The equation for each species is as follows:
% (no. atoms)*MMR*(air/species)*(fractional halogen release)*(Br scaling factor if needed)
% Fractional halon release and obtained from: 
% Engel et al. 2018 (https://www.atmos-chem-phys.net/18/601/2018/acp-18-601-2018.pdf)

% % fractional release values: (Engel et al. (2018))
% fr(1) = .99;    %CFC-11                
% fr(2) = .87;    %CFC-12                
% fr(3) = .91;    %CFC-113              
% fr(4) = .41;    %CFC-114               
% fr(5) = .2;     %CFC-115               
% fr(6) = 1;      %Carbon tetrachlorine  
% fr(7) = .99;    %Methyl chloroform      
% fr(8) = .44;    %HCFC-22              
% fr(9) = .90;    %HCFC-141b              
% fr(10) = .65;   %HCFC-142b           
% fr(11) = 1;     %Halon-1211            
% fr(12) = 1;     %Halon-1202           
% fr(13) = .83;   %Halon-1301            
% fr(14) = 1;     %Halon-2402             
% fr(15) = .99;   %Methyl bromide        
% fr(16) = .91;   %Methyl chloride       

% fractional release values: (Newman et al. (2007))
fr(1) = .99;    %CFC-11                 .99 (Newman)
fr(2) = .86;    %CFC-12                 .86 (Newman)
fr(3) = .90;    %CFC-113                .90 (Newman)
fr(4) = .40;    %CFC-114                .40 (Newman)
fr(5) = .15;     %CFC-115                .15 (Newman)
fr(6) = 1;      %Carbon tetrachlorine   1   (Newman)
fr(7) = .99;    %Methyl chloroform      .99 (Newman)
fr(8) = .41;    %HCFC-22                .41 (Newman)
fr(9) = .90;    %HCFC-141b              .90 (Newman)
fr(10) = .65;   %HCFC-142b              .29 (Newman)
fr(11) = 1;     %Halon-1211             1   (Newman)
fr(12) = 1;     %Halon-1202             1   (Newman)
fr(13) = .8;   %Halon-1301             .8  (Newman)
fr(14) = 1;     %Halon-2402             1   (Newman)
fr(15) = .99;   %Methyl bromide         .99 (Newman)
fr(16) = .91;   %Methyl chloride        .91 (Newman)


%% 
clearvars ODS_data_EESC;
ODS_data_EESC = zeros(length(ODSyears),size(mrtd,2));
for i = 1:length(ODSyears)
    if i <= 20
        endint = 1;
        endint2 = 1:i;
    else
        endint = i-19;
        endint2 = 1:20;
    end
    ODS_data_EESC(i,1) = 3*fr(1)*nansum(ODS_data(endint:i,1).*flipud(inversegaussianfunction(endint2,1)));              %CFC-11                 (CCl3F) *3 (air./137.37)
    ODS_data_EESC(i,2) = 2*fr(2)*nansum(ODS_data(endint:i,2).*flipud(inversegaussianfunction(endint2,2)));              %CFC-12                 (CCl2F2) *2 (air./120.91)
    ODS_data_EESC(i,3) = 3*fr(3)*nansum(ODS_data(endint:i,3).*flipud(inversegaussianfunction(endint2,3)));              %CFC-113                (CCl2FCClF2) *3 (air./187.37)
    ODS_data_EESC(i,4) = 2*fr(4)*nansum(ODS_data(endint:i,4).*flipud(inversegaussianfunction(endint2,4)));              %CFC-114                (CClF2CClF2) (air./170.92)
    ODS_data_EESC(i,5) = 1*fr(5)*nansum(ODS_data(endint:i,5).*flipud(inversegaussianfunction(endint2,5)));              %CFC-115                (CClF2CF3) (air./154.47)
    ODS_data_EESC(i,6) = 4*fr(6)*nansum(ODS_data(endint:i,6).*flipud(inversegaussianfunction(endint2,6)));              %Carbon tetrachlorine   (CCl4) (air./153.82)
    ODS_data_EESC(i,7) = 3*fr(7)*nansum(ODS_data(endint:i,7).*flipud(inversegaussianfunction(endint2,7)));              %Methyl chloroform      (CH3CCl3) (air./133.4)
    ODS_data_EESC(i,8) = 1*fr(8)*nansum(ODS_data(endint:i,8).*flipud(inversegaussianfunction(endint2,8)));              %HCFC-22                (CHClF3) (air./86.47)
    ODS_data_EESC(i,9) = 2*fr(9)*nansum(ODS_data(endint:i,9).*flipud(inversegaussianfunction(endint2,9)));              %HCFC-141b              (CH3CCl2F) (air./116.95)
    ODS_data_EESC(i,10) = 1*fr(10)*nansum(ODS_data(endint:i,10).*flipud(inversegaussianfunction(endint2,10)));          %HCFC-142b              (CH3CClF2) (air./100.5)
    ODS_data_EESC(i,11) = 1*fr(11)*nansum(ODS_data(endint:i,11).*flipud(inversegaussianfunction(endint2,11)));          %Halon-1211             (CBrClF2) %for Cl (air./165.36)
    
    ODS_data_EESC(i,12) = 1*fr(11)*alpha*nansum(ODS_data(endint:i,11).*flipud(inversegaussianfunction(endint2,11)));    %Halon-1211             (CBrClF2) %for Br (air./165.36)
    ODS_data_EESC(i,13) = 2*fr(12)*alpha*nansum(ODS_data(endint:i,12).*flipud(inversegaussianfunction(endint2,12)));    %Halon-1202             (CBr2F2) (air./209.82)
    ODS_data_EESC(i,14) = 1*fr(13)*alpha*nansum(ODS_data(endint:i,13).*flipud(inversegaussianfunction(endint2,13)));    %Halon-1301             (CBrF3) (air./148.91)
    ODS_data_EESC(i,15) = 2*fr(14)*alpha*nansum(ODS_data(endint:i,14).*flipud(inversegaussianfunction(endint2,14)));    %Halon-2402             (CBrF2CBrF2) (air./259.8)
    ODS_data_EESC(i,16) = 1*fr(15)*alpha*nansum(ODS_data(endint:i,15).*flipud(inversegaussianfunction(endint2,15)));    %Methyl bromide         (CH3Br) (air./94.94)
    
    ODS_data_EESC(i,17) = 1*fr(16)*nansum(ODS_data(endint:i,16).*flipud(inversegaussianfunction(endint2,16)));          %Methyl chloride        (CH3Cl) (air./50.49)
    
end

% calculate EESC by summing ODS components
EESC = nansum(ODS_data_EESC,2)*1e9; %ppb

fig = figure;
set(fig,'color','white','position',[100 100 1000 700]);
hold on
plot(EESC(111:end-40),'LineWidth',2);

grid on
set(gca,'ytick',1:.1:4.5,'xtick',1:20:200,'xticklabel',1960:20:2080,'fontsize',18);
ylim([1 4.5]);
xlim([0 130]);
ylabel('ppb','fontsize',20);
xlabel('year','fontsize',18);
title('EESC','fontsize',22);
box on
maxval = max(EESC);
val1980 = EESC(131);

