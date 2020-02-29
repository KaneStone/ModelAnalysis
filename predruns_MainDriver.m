% predruns maindriver

clear variables
%close all
clc

%% Read in model (user inputs located in function)
inputs = predruns_inputs;

%% 
if strcmp(inputs.var,'TS')
    [surfacedata,tozdata] = predruns_ReadinModel(inputs);
elseif strcmp(inputs.var,'T')
    [~,tozdata] = predruns_ReadinModel(inputs);
    directory = '/Volumes/ExternalOne/work/data/predruns/T/highCl/raw/';
    files = dir([directory,'*.nc']);
    [surfacedata.highCl.years,dataLevel,datacomp,surfacedata.highCl.dataMonthArrange,longitude,latitude] = predruns_ReadIn4d(directory,files,inputs.var,[1995,2024],500,inputs.detrend);
    surfacedata.highCl.data.lon = longitude;
    surfacedata.highCl.data.lat = latitude;
end

% obtaining number and names of fields 
fields = fieldnames(surfacedata);
nofields = length(fields);

%% Calculating the upper and lower percentiles of toz based on user inputs defined below
% user inputs

for i = 1:nofields
    nofiles = length(tozdata.(fields{i}).data);
    [pct.(fields{i}),tozextract.(fields{i}),Eachyear.(fields{i})] = predruns_varPercentiles(...
        tozdata.(fields{i}).toz_composite.montharrange,tozdata.(fields{i}).dataMonthArrange,...
        inputs.tozmonth,inputs.percentile,nofiles,inputs.varmonth);
end

longitude = surfacedata.(fields{1}).data(1).lon;
latitude = surfacedata.(fields{1}).data(1).lat;

%% Calculate correlations and upper and lower temperature differences between TS and toz
if inputs.removeENSO
    filename = ['/Volumes/ExternalOne/work/data/predruns/output/data/',inputs.ClLevel{1},'_',inputs.var,'_ninoremoved_','diffsandcorrs_',...
        monthnames(inputs.varmonthtomean,1,'short'),num2str(inputs.timeperiodvar(1)),'-',...
        num2str(inputs.timeperiodvar(2)),num2str(abs(inputs.lats(1))),'-',...
        num2str(abs(inputs.lats(2))),'_',num2str(inputs.detrend)];
else
    filename = ['/Volumes/ExternalOne/work/data/predruns/output/data/',inputs.ClLevel{1},'_',inputs.var,'_diffsandcorrs_',...
        monthnames(inputs.varmonthtomean,1,'short'),num2str(inputs.timeperiodvar(1)),'-',...
        num2str(inputs.timeperiodvar(2)),num2str(abs(inputs.lats(1))),'-',...
        num2str(abs(inputs.lats(2))),'_',num2str(inputs.detrend)];
end
if ~exist([filename,'.mat'],'file')
    for i = 1:nofields
        [dataVarMonthAve,dataVarMonth,NINO34all,differences,correlations] = ...
            predruns_diffandcorr_percentiles(surfacedata.(fields{i}).dataMonthArrange,...
            tozdata.(fields{i}).dataMonthArrange,latitude,longitude,Eachyear,pct,fields{i},inputs);
    end
    save(filename,'dataVarMonthAve','dataVarMonth','differences','correlations')
else
    load(filename);
end
%% plot correlations and percentile differences
for i = 1:nofields
    predruns_plotcorranddiff(correlations,differences,inputs,longitude,latitude)    
end

%% compare to ERA-Interim (plot 1)
if inputs.compareERA
    [corr2dcorr,diff2dcorr,ERAdata,BodekerTCO] = predruns_compareERA([1980,2016],correlations,differences,inputs,longitude,latitude);
end

%%
if inputs.lastfiveyears
    %% plot time series compare
    predruns_latsfiveyears(surfacedata.(fields{i}).dataMonthArrange,dataVarMonthAve,tozdata.highCl.dataMonthArrange,inputs,latitude,longitude);

    %% root mean squared analysis
    predruns_fiveyearpred(dataVarMonthAve,dataVarMonth,tozdata.highCl.dataMonthArrange,inputs,latitude,longitude,pct,differences,correlations);
    
end

if inputs.lastfiveyears_obs
    %% plot time series compare of observations
    predruns_latsfiveyears_obs(ERAdata.t2m,BodekerTCO,inputs,ERAdata.latitude,ERAdata.longitude);
end

%% Take of ensemble average and perform leave one out metric

%predruns_leaveoneout_empirical(dataVarMonthAve,dataVarMonth,tozdata.highCl.dataMonthArrange,inputs,latitude,longitude,pct);
predruns_leaveoneout_empcomp(dataVarMonthAve,dataVarMonth,tozdata.(inputs.ClLevel{1}).dataMonthArrange,inputs,latitude,longitude,pct,differences.indmonths.individual,differences.composite,surfacedata.(fields{i}).dataMonthArrange);
%predruns_leaveoneout_empcomp_composite(dataVarMonthAve,dataVarMonth,tozdata.highCl.dataMonthArrange,inputs,latitude,longitude,pct,differences.indmonths.individual,differences.composite,surfacedata.(fields{i}).dataMonthArrange);
%predruns_leaveoneout(dataVarMonthAve,dataVarMonth,tozdata.highCl.dataMonthArrange,inputs,latitude,longitude,pct,differences.indmonths.individual,correlations.indmonths.individual);
%predruns_RMSE(dataVarMonthAve,dataVarMonth,tozdata.highCl.dataMonthArrange,inputs,latitude,longitude,pct,differences.indmonths.individual,correlations.indmonths.individual);
%predruns_leaveoneout_test(dataVarMonthAve,dataVarMonth,tozdata.highCl.dataMonthArrange,inputs,latitude,longitude,pct);

%% obs leave one out

Observations = load(['/Volumes/ExternalOne/work/data/predruns/output/data/obs/','obs_perc_diff',monthnames(inputs.varmonthtomean,1,1),'_and_',monthnames(inputs.varmonth,1,1),'.mat']);
[obs.GSS] = predruns_obsleaveoneout(Observations,longitude,latitude,inputs);
%[obs.GSS] = predruns_obsleaveoneout_test(Observations,longitude,latitude,inputs);

