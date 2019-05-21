% predruns SSW coincide with TCO

%% Read in model (user inputs located in function)
inputs.ClLevel = 'lowGHG';
%inputs.timeperiodtoz = [1955,1979];%[1995,2024];
inputs.timeperiodtoz = [1995,2024];
inputs.lats = [63,90];
inputs.detrend_ozone = 0;

%% Read in ozone and Uwind
vardirectory = ['/Volumes/ExternalOne/work/data/predruns/U/',inputs.ClLevel,'/SSW/'];
varfiles = dir([vardirectory,'*.nc']);
tozdirectory = ['/Volumes/ExternalOne/work/data/predruns/','toz','/',inputs.ClLevel,'/'];
tozfiles = dir([tozdirectory,'*.nc']);

%%
for i = 1:length(varfiles)
    % Read in surface temperature or other similar variable        

    [~,U(i),~] = Read_in_netcdf([vardirectory,varfiles(i).name]);
    U(i).U = squeeze(U(i).U);
    % Read in TOZ or other similar variable    
    
end

[tozdata.(inputs.ClLevel).data,tozdata.(inputs.ClLevel).years,tozdata.(inputs.ClLevel).varweighted,...
        tozdata.(inputs.ClLevel).toz_composite,tozdata.(inputs.ClLevel).dataMonthArrange] = ...
        predruns_ReadInlayer_areaaverage(tozdirectory,tozfiles,'toz',inputs.timeperiodtoz,inputs.lats,...
        inputs.detrend_ozone);

%%
Uallmean = squeeze(nanmean(cat(3,U(:).U),1));

%% plot test
mon = 3;
mon2 = 3;
for i = 1:length(varfiles)
%     figure;
%     yyaxis left
%     plot(Uallmean(mon:12:end,i))
%     hold on
%     yyaxis right
%     plot(tozdata.highCl.varweighted(i,mon:12:end))  
    corrs(i) = corr(Uallmean(mon:12:end,i),tozdata.(inputs.ClLevel).varweighted(i,mon2:12:end)');
end

hold on
plot(corrs);

%%
for j = 1:12
    for i = 1:length(varfiles)   
        corrs2(j,i) = corr(Uallmean(j:12:end,i),tozdata.highCl.varweighted(i,j:12:end)');
    end
end
corrsall = nanmean(corrs2,2);
figure;
plot(corrsall);