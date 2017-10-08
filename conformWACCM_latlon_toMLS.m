% conform lat lon height time to MLS
clear all

variable = 'O3';
pres = 150;
lats = [-90 -26];
runs = {'MAM','VCMAM'};
daytitle = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};
MLS_compare = 1;
latlon = 1;
%% Reading in WACCM data
if latlon
    commondir = '/Volumes/My Book for Mac/work/data/WACCM/';
else
    commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
end
if exist('/Users/kanestone/work/projects/WACCM/output/latlon/MolConc.mat','file') ~= 2;
    [NumberDensity, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes, Longitudes,datedata] = ReadWACCMvertical(variable,'monthly',commondir,1);
    save('/Users/kanestone/work/projects/WACCM/output/latlon/MolConc.mat','MolConc','Pressure','-v7.3');
else
    load('/Users/kanestone/work/projects/WACCM/output/latlon/MolConc.mat');   
    MolConc.VCMAM = MolConc.MAM;
    MolConc.MAM = MolConc.CCMI;   
    MolConc = rmfield(MolConc,'CCMI');
    Pressure.MAM = Pressure.CCMI;
    Pressure.VCMAM = Pressure.MAM;
    Pressure = rmfield(Pressure,'CCMI');
end

%% MLS data
addpath('/Users/kanestone/work/projects/MLS/code/');
% a priori
[MLSpressure,O3ap_monthmean,tempap_monthmean] = importMLSapriori;
[~,MLSpresindex] = min(abs(MLSpressure - pres));

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

%%interpolate to MLS pressure
MolConcMLSgrid.MAM = zeros(size(MolConc.MAM,1),size(MolConc.MAM,2),55,size(MolConc.MAM,4));
MolConcMLSgrid.VCMAM = zeros(size(MolConc.MAM,1),size(MolConc.MAM,2),55,size(MolConc.MAM,4));

%%
if exist('/Users/kanestone/work/projects/WACCM/output/latlon/MolConcMLSgrid.mat','file') ~= 2;
    if latlon
        for i = 1:length(runs)       
            for j = 1:size(MolConc.(runs{i}),1)
                tic;
                j
                for k = 1:size(MolConc.(runs{i}),2)
                    for l = 1:size(MolConc.(runs{i}),4)
                        MolConcMLSgrid.(runs{i})(j,k,:,l) = interp1(squeeze(log(Pressure.(runs{i})(j,k,:,l)./100)),...
                            squeeze(MolConc.(runs{i})(j,k,:,l)),log(MLSpressure),'linear','extrap');
                        %MolConcMLSgridconvolve = squeeze(O3ap_monthmean(j,k,:,l))+...
                        %AKdata2*(squeeze(MolClatMLSperiod.(runs{i})(j,k,:,l))-squeeze(O3ap_monthmean(j,k,:,l)));                    
                    end
                end
                toc;
            end
        end
    end
else
    load('/Users/kanestone/work/projects/WACCM/output/latlon/MolConcMLSgrid.mat');
end
MolConcMLSgrid.MAM = MolConcMLSgrid.MAM(:,1:34,:,5*12+1:end);
MolConcMLSgrid.VCMAM = MolConcMLSgrid.VCMAM(:,1:34,:,5*12+1:end);
szWAC = size(MolConcMLSgrid.MAM);
MolConcMLSgrid.MAM = cat(4,MolConcMLSgrid.MAM,zeros([szWAC(1:3),6]));
MolConcMLSgrid.VCMAM = cat(4,MolConcMLSgrid.VCMAM,zeros([szWAC(1:3),6]));
MolConcMLSgrid.MAM (MolConcMLSgrid.MAM == 0) = NaN;
MolConcMLSgrid.VCMAM (MolConcMLSgrid.VCMAM == 0) = NaN;
% rearrange a priori to [lon lat height time]
count = 1;
for i = 1:13;
    for j = 1:12;
        apriorireshape(count,:,:) = O3ap_monthmean(i,j,:,:);
        count = count+1;
    end
end
apriorireshape = repmat(apriorireshape,[1,1,1,144]);
apriorireshape = permute(apriorireshape,[4,3,2,1]);

%%
for i = 1:length(runs)       
    for j = 1:size(MolConcMLSgrid.(runs{i}),1)
        tic;
        j
        for k = 1:size(MolConcMLSgrid.(runs{i}),2)
            for l = 1:size(MolConcMLSgrid.(runs{i}),4)
                WACCM_MolC_convolve.(runs{i})(j,k,:,l) = squeeze(apriorireshape(j,k,:,l))+...
                    AKdata2*(squeeze(MolConcMLSgrid.(runs{i})(j,k,:,l))-squeeze(apriorireshape(j,k,:,l)));
            end
        end
    end
end
save('/Users/kanestone/work/projects/WACCM/output/latlon/WACCM_MolC_convolve.mat','WACCM_MolC_convolve','MLSpressure');
