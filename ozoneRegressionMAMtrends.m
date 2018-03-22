

%% User inputs

timeperiod = [1999,2014];

%% ozone regression - MAM only

SWOOSH.level = [316.22775;261.01572;215.44347;177.82794;146.77992;121.15276;100;82.540421;...
    68.129204;56.234131;46.415890;38.311867;31.622776;26.101572;21.544348;17.782795;14.677993;...
    12.115276;10;8.2540417;6.8129206;5.6234131;4.6415887;3.8311868;3.1622777;2.6101573;2.1544347;...
    1.7782794;1.4677993;1.2115277;1];

%% import highCl
highCllevel = [SWOOSH.level;[.9,.8,.7,.6,.5,.4,.3,.2,.1]'];

%% Read in WACCM

[~,MAM.O3,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/O3/otherMAM/O3_species_f.e11.FWTREFC1SD.f19.f19.ccmi34.79-15.nohetice.mfp.mam.004fg.cam.h0zm.nc');
[~,MAM.T,~] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/O3/otherMAM/T_species_f.e11.FWTREFC1SD.f19.f19.ccmi34.79-15.nohetice.mfp.mam.004fg.cam.h0zm.nc');
MAM.timeperiod.years = 1979:2014;
timeperiodmonths_temp = repmat(1979:2014,12,1);
MAM.timeperiod.months = timeperiodmonths_temp(:);

dateind = MAM.timeperiod.years >= timeperiod(1) & MAM.timeperiod.years <= timeperiod(2);
dateind_months = MAM.timeperiod.months >= timeperiod(1) & MAM.timeperiod.months <= timeperiod(2);

clearvars timeperiodmonths_temp
%% Creating pressure

hyam = permute(repmat(MAM.O3.hyam,[1,size(MAM.O3.O3,1),size(MAM.O3.O3,3)]),[2,1,3]);
hybm = permute(repmat(MAM.O3.hybm,[1,size(MAM.O3.O3,1),size(MAM.O3.O3,3)]),[2,1,3]);
PS = permute(repmat(MAM.O3.PS,[1,1,size(MAM.O3.hyam,1)]),[1,3,2]); 

MAM.Pressure = 100000 .* hyam + hybm .* PS;        

%% Interpolate onto regular pressure

for i = 1:size(MAM.O3.O3,1)
    if i ~= size(MAM.O3.O3,1)
        O3regpres(i,:,:) = intRegPres(squeeze(MAM.O3.O3(i,:,:)),squeeze(MAM.Pressure(i,:,:))./100);
    else
        [O3regpres(i,:,:),repres] = intRegPres(squeeze(MAM.O3.O3(i,:,:)),squeeze(MAM.Pressure(i,:,:))./100);
    end    
end

%% construct yearly average and extract time period
O3regpres_timeextract = O3regpres(:,:,dateind_months);
for i = 1:12
    O3montharrange(:,:,i,:) = O3regpres_timeextract(:,:,i:12:end);
    O3montharrange_percent(:,:,i,:) = (squeeze(O3montharrange(:,:,i,:)) - nanmean(O3montharrange(:,:,i,:),4))./nanmean(O3montharrange(:,:,i,:),4)*100;
end

O3_montharrange_forreg = permute(O3montharrange_percent(:,:,:),[3,1,2]);
O3yearlyaverage = squeeze(nanmean(O3montharrange,3));
O3yearlyaverage_percent = (O3yearlyaverage - nanmean(O3yearlyaverage,3))./nanmean(O3yearlyaverage,3)*100;
O3yearlyaverage_percent_3y = (O3yearlyaverage - nanmean(O3yearlyaverage(:,:,1:3),3))./nanmean(O3yearlyaverage(:,:,1:3),3)*100;

O3yearlyaverage_foreg = permute(O3yearlyaverage_percent,[3,1,2]);
O3yearlyaverage_foreg_3y = permute(O3yearlyaverage_percent_3y,[3,1,2]);

%% constructing percent anomalies

%% take linear trends
regressfunctions = [ones(sum(dateind),1),[1:sum(dateind)]'];
regressfunctions_months = [ones(sum(dateind_months),1),[1:sum(dateind_months)]'];
for i = 1:size(O3yearlyaverage_foreg,2)
    for j = 1:size(O3yearlyaverage_foreg,3)
        
        b.yearlyaverage(i,j,:) = regress(O3yearlyaverage_foreg(:,i,j),regressfunctions);
        b.monthlytimeseries(i,j,:) = regress(O3_montharrange_forreg(:,i,j),regressfunctions_months);
        b.yearlyaverage_3y(i,j,:) = regress(O3yearlyaverage_foreg_3y(:,i,j),regressfunctions);
        
    end
end

%% plot linear regression
btoplot = permute(cat(3,b.yearlyaverage(:,:,2)*10,b.yearlyaverage_3y(:,:,2)*10,b.monthlytimeseries(:,:,2)*120),[3,1,2]);

prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
    5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
logprestick = log(prestick);

titles = {['Yearly average,',' normalized to entire time series, ', num2str(timeperiod(1)),char(8211),num2str(timeperiod(2))],...
    ['Yearly average,',' normalized to first 3 years, ', num2str(timeperiod(1)),char(8211),num2str(timeperiod(2))],...
    ['Monthly time series, normalized to entire time series, ', num2str(timeperiod(1)),char(8211),num2str(timeperiod(2))]}; 

mtit = {'MAM'};

subplotmaps(btoplot,MAM.O3.lat,log(repres),{'div','RdBu'},1,[],16,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
    [-10 10],12,-90:30:90,-90:30:90,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(1) log(300)],1,'-',0,'');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/WilkaShahSolomonPaper/','MAM_trends'];

export_fig(filename,'-pdf');



%% take regression of SD WACCM (chem-only and MAM)
% 
% [bchemonly,predictorschemonly,O3Anomalychemonly] = ...
%     ozoneRegressionTrends(SDWaccmData(2).O3,SDWaccmData(2).U,SDWaccmData(2).solar,SDWaccmData(2).SPE,...
%     SDWaccmData(2).NINO34,SDWaccmData(2).HFS,SDWaccmData(2).HFN,SDWaccmData(2).NO2,SDtimeperiod,1,0,1);%SDWaccmData(2).SPEintpres 

%% 