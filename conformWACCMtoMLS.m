% Conform WACCM data to MLS (and perhaps CALIPSO later)
clear all
variable = 'O3';
pres = 150;
lats = [-90 90];
runs = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
daytitle = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};
MLS_compare = 1;
latlon = 0;
%% Reading in WACCM data
if latlon
    commondir = '/Volumes/My Book for Mac/work/data/WACCM/';
else
    commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
end
[NumberDensity, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes, Longitudes,datedata] = ReadWACCMvertical(variable,'monthly',commondir,0);

commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[~, SAD_SULFC, ~, ~, ~, ~, ~, ~, ~] = ReadWACCMvertical('SAD_SULFC','monthly',commondir,0);

[~, Extinction, ~, ~, ~, ~, ~, ~] = ReadWACCMvertical('Extinction','monthly',commondir,0);

Extinction.CCMI = zeros(size(Extinction.MAM));
Extinction.CCMI (Extinction.CCMI == 0) = NaN;
Extinction.Chemonly = zeros(size(Extinction.MAM));
Extinction.Chemonly (Extinction.Chemonly == 0) = NaN;
Extinction.ChemonlyfixedSSTs = zeros(size(Extinction.MAM));
Extinction.ChemonlyfixedSSTs (Extinction.ChemonlyfixedSSTs == 0) = NaN;
Extinction.Chemonlynoleap = zeros(size(Extinction.MAM));
Extinction.Chemonlynoleap (Extinction.Chemonlynoleap == 0) = NaN;

[STS_ND, STS_MolConc, ~, STS_Pressure, ~, ~, STS_Latitudes, STS_datedata] = ...
    ReadWACCMvertical('HNO3_STS','monthly',commondir,0);

[HNO3GAS_ND, HNO3GAS_MolConc, ~, HNO3GAS_Pressure, ~, ~, HNO3GAS_Latitudes, HNO3GAS_datedata] = ...
    ReadWACCMvertical('HNO3_GAS','monthly',commondir,0);

[~, SAD_ICE, ~, ~, ~, ~, ~, ~] = ...
    ReadWACCMvertical('SAD_ICE','monthly',commondir,0);

%% MLS data
addpath('/Users/kanestone/work/projects/MLS/code/');
% a priori
%[MLSpressure,O3ap_monthmean,tempap_monthmean] = importMLSapriori;
[MLS,MLSpresindex,~,MLSpressure] = importMLSmonthly(pres,'globe');
O3ap_monthmean = permute(MLS.O3ap_monthmean,[1,2,4,3]);
tempap_monthmean = permute(MLS.tempap_monthmean,[1,2,4,3]);
%[~,MLSpresindex] = min(abs(MLSpressure - pres));
% AK
filename = '/Users/kanestone/work/projects/MLS/AK/MLS-Aura_L2AK-O3-LAT70N_v04-2x_0000d000.txt';
filename_temp = '/Users/kanestone/work/projects/MLS/AK/MLS-Aura_L2AK-Temperature-LAT70N_v04-2x_0000d000.txt';
startRow = 23;

% For more information, see the TEXTSCAN documentation.
formatSpec = '%13s%13s%13s%13s%13s%s%[^\n\r]';

fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

fileID_temp = fopen(filename_temp,'r');
dataArray_temp = textscan(fileID_temp, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

for i = 1:length(dataArray)-1
    AKdata(:,i) = str2double(dataArray{:,i});
    AKdata_temp(:,i) = str2double(dataArray_temp{:,i});
end

AKdata1 = reshape(AKdata',60,55);
AKdatahold = AKdata1;
AKdatahold (isnan(AKdata1)) = [];
AKdata2 = reshape(AKdatahold,55,55);

AKdata1_temp = reshape(AKdata_temp',60,55);
AKdatahold_temp = AKdata1_temp;
AKdatahold_temp (isnan(AKdata1_temp)) = [];
AKdata2_temp = reshape(AKdatahold_temp,55,55);
%% 
if latlon
    for i = 1:length(runs)       
        for j = 1:size(MolConc.(runs{i}),1)
            for k = 1:size(MolConc.(runs{i}),2)
                for l = 1:size(MolConc.(runs{i}),4)
                    MolConcMLSgrid.(runs{i})(j,k,:,l) = interp1(squeeze(log(Pressure.(runs{i})(j,k,:,l)./100)),...
                        squeeze(MolConc.(runs{i})(j,k,:,l)),log(MLSpressure),'linear','extrap');
                    MolConcMLSgridconvolve = squeeze(O3ap_monthmean(j,k,:,l))+...
                    AKdata2*(squeeze(MolClatMLSperiod.(runs{i})(j,k,:,l))-squeeze(O3ap_monthmean(j,k,:,l)));
                end
            end
        end
    end
end

%% comparing interpolation methods

% piecewise least squares
WACCM_to_MLS_interp = interp1(log(Pressure.MAM(10,:,202)./100),...
    MolConc.MAM(10,:,202),log(MLSpressure),'linear','extrap');
%WACCM_to_MLS_lsq = lsq_lut_piecewise(log(Pressure.MAM(1,:,202)./100)',...
%    double(MolConc.MAM(1,:,202)'),flipud(double(log(MLSpressure))));
WACCM_to_MLS_interp150 = interp1(log(Pressure.MAM(10,:,202)./100),...
    MolConc.MAM(10,:,202),log(150),'linear','extrap');

plot(WACCM_to_MLS_interp,log(MLSpressure))
hold on
plot(MolConc.MAM(10,:,202),log(Pressure.MAM(10,:,202)./100))
scatter(WACCM_to_MLS_interp150,log(150),'d')

%% Conforming data
for i = 1:length(runs)       
    %interp onto MLS grid
    for j = 1:size(MolConc.(runs{i}),1)
        for k = 1:size(MolConc.(runs{i}),3)
            MolConcMLSgrid.(runs{i})(j,:,k) = interp1(log(Pressure.(runs{i})(j,:,k)./100),...
                MolConc.(runs{i})(j,:,k),log(MLSpressure),'linear','extrap');
            TemperatureMLSgrid.(runs{i})(j,:,k) = interp1(log(Pressure.(runs{i})(j,:,k)./100),...
                Temperature.(runs{i}).T(j,:,k),log(MLSpressure),'linear','extrap');
            NDMLSgrid.(runs{i})(j,:,k) = vmr2conc(MolConcMLSgrid.(runs{i})(j,:,k),...
                TemperatureMLSgrid.(runs{i})(j,:,k),MLSpressure'*100,'O3','conc');                                        
            SAD_SULFCMLSgrid.(runs{i})(j,:,k) = interp1(log(Pressure.(runs{i})(j,:,k)./100),...
                SAD_SULFC.(runs{i})(j,:,k),log(MLSpressure),'linear','extrap');
            STS_MolConcMLSgrid.(runs{i})(j,:,k) = interp1(log(Pressure.(runs{i})(j,:,k)./100),...
                STS_MolConc.(runs{i})(j,:,k),log(MLSpressure),'linear','extrap');
            HNO3GAS_MolConcMLSgrid.(runs{i})(j,:,k) = interp1(log(Pressure.(runs{i})(j,:,k)./100),...
                HNO3GAS_MolConc.(runs{i})(j,:,k),log(MLSpressure),'linear','extrap');                          
            SAD_ICEMLSgrid.(runs{i})(j,:,k) = interp1(log(Pressure.(runs{i})(j,:,k)./100),...
                SAD_ICE.(runs{i})(j,:,k),log(MLSpressure),'linear','extrap');                          
        end                        
        for k = 1:size(Extinction.(runs{i}),3)
            ExtinctionMLSgrid.(runs{i})(j,:,k) = interp1(log(Pressure.(runs{i})(j,:,k)./100),...
                Extinction.(runs{i})(j,:,k),log(MLSpressure),'linear','extrap');
        end
    end

    %latforweight = [-90 -85; -85 -80; -80 -75; -75 -70; -70 -65; -65 -60; -60 -55; -55 -50;...
    %    -50 -45; -45 -40; -40 -35; -35 -30; -30 -25];    
    latindex = find(Latitudes >= lats(1) & Latitudes < lats(2));
    MolClat.(runs{i}) = MolConcMLSgrid.(runs{i})(latindex,:,:);
    MolClat.(runs{i}) = permute(MolClat.(runs{i}),[3,2,1]);
    Templat.(runs{i}) = TemperatureMLSgrid.(runs{i})(latindex,:,:);
    Templat.(runs{i}) = permute(Templat.(runs{i}),[3,2,1]);
    NDlat.(runs{i}) = NDMLSgrid.(runs{i})(latindex,:,:);
    NDlat.(runs{i}) = permute(NDlat.(runs{i}),[3,2,1]);
    SAD_SULFClat.(runs{i}) = SAD_SULFCMLSgrid.(runs{i})(latindex,:,:);
    SAD_SULFClat.(runs{i}) = permute(SAD_SULFClat.(runs{i}),[3,2,1]);                
    STSlat.(runs{i}) = STS_MolConcMLSgrid.(runs{i})(latindex,:,:);
    STSlat.(runs{i}) = permute(STSlat.(runs{i}) ,[3,2,1]);
    HNO3GASlat.(runs{i}) = HNO3GAS_MolConcMLSgrid.(runs{i})(latindex,:,:);
    HNO3GASlat.(runs{i}) = permute(HNO3GASlat.(runs{i}) ,[3,2,1]);
    SAD_ICElat.(runs{i}) = SAD_ICEMLSgrid.(runs{i})(latindex,:,:);
    SAD_ICElat.(runs{i}) = permute(SAD_ICElat.(runs{i}) ,[3,2,1]);
    Extinctionlat.(runs{i}) = ExtinctionMLSgrid.(runs{i})(latindex,:,:);
    Extinctionlat.(runs{i}) = permute(Extinctionlat.(runs{i}) ,[3,2,1]);
    
%     for j = 1:size(latforweight,1)
%         latsforweightindex = find(Latitudes >= latforweight(j,1) & Latitudes < latforweight(j,2));
%         for k = 1:size(MolConcMLSgrid.(runs{i}),2)
%             MolClat.(runs{i})(:,k,j) = weightedaverage(squeeze(MolConcMLSgrid.(runs{i})(latsforweightindex,k,:)),...
%                 Latitudes(latsforweightindex));
%             Templat.(runs{i})(:,k,j) = weightedaverage(squeeze(TemperatureMLSgrid.(runs{i})(latsforweightindex,k,:)),...
%                 Latitudes(latsforweightindex));
%             NDlat.(runs{i})(:,k,j) = weightedaverage(squeeze(NDMLSgrid.(runs{i})(latsforweightindex,k,:)),...
%                 Latitudes(latsforweightindex));
%             SAD_SULFC_lat.(runs{i})(:,k,j) = weightedaverage(squeeze(SAD_SULFCMLSgrid.(runs{i})(latsforweightindex,k,:)),...
%                 Latitudes(latsforweightindex));
%         end                                                
%     end                            

    %extract same time period as MLS
    %MLS years = 2004:2016
    MolClat2.(runs{i}) = MolClat.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
        datedata.(runs{i}).date <= 20170101,:,:);
    Templat2.(runs{i}) = Templat.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
        datedata.(runs{i}).date <= 20170101,:,:);
    NDlat2.(runs{i}) = NDlat.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
        datedata.(runs{i}).date <= 20170101,:,:);
    SAD_SULFClat2.(runs{i}) = SAD_SULFClat.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
        datedata.(runs{i}).date <= 20170101,:,:);
    STSlat2.(runs{i}) = STSlat.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
        datedata.(runs{i}).date <= 20170101,:,:);
    HNO3GASlat2.(runs{i}) = HNO3GASlat.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
        datedata.(runs{i}).date <= 20170101,:,:);    
     SAD_ICElat2.(runs{i}) = SAD_ICElat.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
            datedata.(runs{i}).date <= 20170101,:,:);        
    
    monthend = num2str(datedata.(runs{i}).date(end));
    yearend = num2str(datedata.(runs{i}).date(end));
    yearend = str2double(yearend(1:4));
    monthend = str2double(monthend(5:6))-1;
    yearappend = 2016-yearend;
    monthappend = (yearappend*12)+12-monthend;
    MolClat2.(runs{i}) = cat(1,MolClat2.(runs{i}),zeros(monthappend,...
        size(MolClat2.(runs{i}),2),size(MolClat2.(runs{i}),3)));
    Templat2.(runs{i}) = cat(1,Templat2.(runs{i}),zeros(monthappend,...
        size(MolClat2.(runs{i}),2),size(MolClat2.(runs{i}),3)));
    NDlat2.(runs{i}) = cat(1,NDlat2.(runs{i}),zeros(monthappend,...
        size(MolClat2.(runs{i}),2),size(MolClat2.(runs{i}),3)));
    SAD_SULFClat2.(runs{i}) = cat(1,SAD_SULFClat2.(runs{i}),zeros(monthappend,...
        size(SAD_SULFClat2.(runs{i}),2),size(SAD_SULFClat2.(runs{i}),3)));
    STSlat2.(runs{i}) = cat(1,STSlat2.(runs{i}),zeros(monthappend,...
        size(STSlat2.(runs{i}),2),size(STSlat2.(runs{i}),3)));
    HNO3GASlat2.(runs{i}) = cat(1,HNO3GASlat2.(runs{i}),zeros(monthappend,...
        size(HNO3GASlat2.(runs{i}),2),size(HNO3GASlat2.(runs{i}),3)));
    SAD_ICElat2.(runs{i}) = cat(1,SAD_ICElat2.(runs{i}),zeros(monthappend,...
        size(SAD_ICElat2.(runs{i}),2),size(SAD_ICElat2.(runs{i}),3)));
    Extinctionlat.(runs{i}) = cat(1,Extinctionlat.(runs{i}),zeros(3,...
        size(Extinctionlat.(runs{i}),2),size(Extinctionlat.(runs{i}),3)));    
    MolClat2.(runs{i}) (MolClat2.(runs{i}) == 0) = NaN;
    Templat2.(runs{i}) (Templat2.(runs{i}) == 0) = NaN;
    NDlat2.(runs{i}) (NDlat2.(runs{i}) == 0) = NaN;
    Extinctionlat.(runs{i}) (Extinctionlat.(runs{i}) == 0) = NaN;
    %SAD_SULFC_lat2.(runs{i}) (SAD_SULFC_lat2.(runs{i}) == 0) = NaN;
    for j = 1:12
        %for k = 1:17
            MolClatMLSperiod.(runs{i})(:,j,:,:) = MolClat2.(runs{i})(j:12:end,:,:);
            TemplatMLSperiod.(runs{i})(:,j,:,:) = Templat2.(runs{i})(j:12:end,:,:);
            NDlatMLSperiod.(runs{i})(:,j,:,:) = NDlat2.(runs{i})(j:12:end,:,:);
            SAD_SULFClatMLSperiod.(runs{i})(:,j,:,:) = SAD_SULFClat2.(runs{i})(j:12:end,:,:);
            STSlatMLSperiod.(runs{i})(:,j,:,:) = STSlat2.(runs{i})(j:12:end,:,:);
            HNO3GASlatMLSperiod.(runs{i})(:,j,:,:) = HNO3GASlat2.(runs{i})(j:12:end,:,:);
            SAD_ICElatMLSperiod.(runs{i})(:,j,:,:) = SAD_ICElat2.(runs{i})(j:12:end,:,:);
            Extinctionrearrange.(runs{i})(:,j,:,:) = Extinctionlat.(runs{i})(j:12:end,:,:);
        %end
    end
    % convolve with MLS averaging kernels

    for j = 1:size(MolClatMLSperiod.(runs{i}),1)
        for k = 1:size(MolClatMLSperiod.(runs{i}),2)
            for l = 1:size(MolClatMLSperiod.(runs{i}),4)
                WACCM_MolC_convolve.(runs{i})(j,k,:,l) = squeeze(O3ap_monthmean(j,k,:,l))+...
                    AKdata2*(squeeze(MolClatMLSperiod.(runs{i})(j,k,:,l))-squeeze(O3ap_monthmean(j,k,:,l)));
                WACCM_Temp_convolve.(runs{i})(j,k,:,l) = squeeze(tempap_monthmean(j,k,:,l))+...
                    AKdata2_temp*(squeeze(TemplatMLSperiod.(runs{i})(j,k,:,l))-squeeze(tempap_monthmean(j,k,:,l)));                        
                WACCM_ND_convolve.(runs{i})(j,k,:,l) = vmr2conc(squeeze(WACCM_MolC_convolve.(runs{i})(j,k,:,l)),...
                    squeeze(WACCM_Temp_convolve.(runs{i})(j,k,:,l)),MLSpressure*100,'O3','conc');
            end
        end
    end

    %comparing different grid interpolations
%     prestick = [1000 500 200 100 50 20 10 5 2 1];
%     logprestick = log(prestick);
%     figure;
%     plot(squeeze(NDlat.(runs{i})(26*12+10,:,2)),log(MLSpressure))            
%     hold on
%     plot(squeeze(NDlat2.(runs{i})(22,:,2)),log(MLSpressure))
%     plot(squeeze(NDlatMLSperiod.(runs{i})(2,10,:,2)),log(MLSpressure))
%     plot(squeeze(WACCM_ND_convolve.(runs{i})(2,10,:,2)),log(MLSpressure))
%     plot(squeeze(NDMLSgrid.(runs{i})(5,:,26*12+10)),log(MLSpressure))
%     plot(squeeze(NumberDensity.(runs{i})(5,:,26*12+10)),log(squeeze(Pressure.(runs{i})(6,:,27*12+10))./100))
%     set(gca,'Ydir','reverse','ytick',fliplr(logprestick ),'yticklabel',fliplr(prestick))
%     ylim([log(.5) log(1000)])
%     legend('Lat weighted','Lat weighted2','Lat weighted MLS period','Convolved','MLS grid','WACCM grid')
%     figure;
%     plot(squeeze(MolClat.(runs{i})(26*12+10,:,2)),log(MLSpressure))
%     hold on
%     plot(squeeze(MolClatMLSperiod.(runs{i})(2,10,:,2)),log(MLSpressure))
%     plot(squeeze(WACCM_MolC_convolve.(runs{i})(2,10,:,2)),log(MLSpressure))
%     plot(squeeze(MolConcMLSgrid.(runs{i})(5,:,26*12+10)),log(MLSpressure))
%     plot(squeeze(MolConc.(runs{i})(5,:,26*12+10)),log(squeeze(Pressure.(runs{i})(5,:,27*12+10))./100))            
%     set(gca,'Ydir','reverse','ytick',fliplr(logprestick ),'yticklabel',fliplr(prestick))
%     ylim([log(.5) log(1000)])
%     legend('Lat weighted','Lat weighted MLS period','Convolved','MLS grid','WACCM grid')

    %save('/Users/kanestone/work/projects/WACCM/MatlabOutput/WACCM_MLS_compare.mat','WACCM_MolC_convolve');
    %extracting height
    WACCM_MolC_convolve_atpres(i,:,:,:) = WACCM_MolC_convolve.(runs{i})(:,:,MLSpresindex,:);
    
end

save('/Users/kanestone/work/projects/WACCM/MatlabOutput/WACCM_MLS_compare_90S90N.mat',...
    'WACCM_MolC_convolve','WACCM_ND_convolve','SAD_SULFClatMLSperiod','STSlatMLSperiod',...
    'HNO3GASlatMLSperiod','Extinctionrearrange','WACCM_Temp_convolve','TemplatMLSperiod','SAD_ICElatMLSperiod');