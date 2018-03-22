% CanESM2 percentiles
% Read in and analyse all zonal wave data
clear all
directory = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/';
files = dir([directory,'*.mat']);
% user inputs

CanESM2years = 1950:2100;
%% reading in data
for i = 1:length(files)
    nameofdata{i} = files(i).name(1:end-4);
    if i < length(files)
        temp = load([directory,files(i).name]);
        data.(nameofdata{i}) = temp.(nameofdata{i});
    else
        WACCM_predruns = load([directory,files(i).name]);
    end
    nameofdatafortitles{i} = nameofdata{i}; 
    nameofdatafortitles{i} (nameofdatafortitles{i} == '_') = ' ';
end

%files,longitudes,pressures,years
load('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/lonvert/CanESM2vert3060S.mat');

%% variation in minimum and maximum longitudes
minall = cat(3,data.Can2ESM.minlongitude_m(51:100).m);
minallspring = squeeze(nanmean(minall(9:11,4,:)));
minallstd = std(minallspring);



%% find percentiles
names = {'highcl','highGHGlowcl','lowcl'};
lats = [-50,-55,-60,-65,-70,-75,-80,-85];
latind = 4;
perc = 20;
timeperiod = [1995,2024;2050,2079;1955,1979];    
for i = 1:size(timeperiod,1)
    timelength = timeperiod(i,2) - timeperiod(i,1)+1;
    CanESM2timeextract.(names{i}) = CanESM2verens.spring(1:50,:,:,CanESM2years >= timeperiod(i,1) & CanESM2years <= timeperiod(i,2));
    CanESM24060_composite.(names{i}) = reshape(permute(CanESM2timeextract.(names{i}),[2,3,4,1]),...
        [size(CanESM2timeextract.(names{i}),2),size(CanESM2timeextract.(names{i}),3),size(CanESM2timeextract.(names{i}),4)*size(CanESM2timeextract.(names{i}),1)]); 
    amplitude = cat(3,data.Can2ESM.amplitude_m_spr(51:100).m);
    amplitude_timeextract.(names{i}) = amplitude(:,CanESM2years >= timeperiod(i,1) & CanESM2years <= timeperiod(i,2),:);
    amplitude_composite.(names{i}) = reshape(amplitude_timeextract.(names{i}),[8,50*timelength]);
    percentiles.(names{i}).low = prctile(amplitude_composite.(names{i})(latind,:),perc);
    percentiles.(names{i}).high  = prctile(amplitude_composite.(names{i})(latind,:),100-perc);
    prctile_ind.(names{i}).low = find(amplitude_composite.(names{i})(latind,:) <= percentiles.(names{i}).low); 
    prctile_ind.(names{i}).high = find(amplitude_composite.(names{i})(latind,:) >= percentiles.(names{i}).high);
    
    differences(i,:,:) = nanmean(CanESM24060_composite.(names{i})(:,:,prctile_ind.(names{i}).high),3) - ...
        nanmean(CanESM24060_composite.(names{i})(:,:,prctile_ind.(names{i}).low),3);    
    for j = 1:size(CanESM24060_composite.(names{i}),1)
        for k = 1:size(CanESM24060_composite.(names{i}),2)
            [h(i,j,k)] = ttest2(CanESM24060_composite.(names{i})(j,k,prctile_ind.(names{i}).high),CanESM24060_composite.(names{i})(j,k,prctile_ind.(names{i}).low));            
        end
    end
end

names = {'highcl','highGHGlowcl','lowcl'};
lats = [-50,-55,-60,-65,-70,-75,-80,-85];
latind = 4;
perc = 20;
timeperiod = [1995,2024;2050,2079;1955,1979];    
for i = 1:size(timeperiod,1)
    timelength = timeperiod(i,2) - timeperiod(i,1)+1;
    CanESM2timeextract.(names{i}) = CanESM2verens.spring(51:100,:,:,CanESM2years >= timeperiod(i,1) & CanESM2years <= timeperiod(i,2));
    CanESM24060_composite.(names{i}) = reshape(permute(CanESM2timeextract.(names{i}),[2,3,4,1]),...
        [size(CanESM2timeextract.(names{i}),2),size(CanESM2timeextract.(names{i}),3),size(CanESM2timeextract.(names{i}),4)*size(CanESM2timeextract.(names{i}),1)]); 
    amplitude = cat(3,data.Can2ESM.amplitude_m_spr(102:2:200).m);
    amplitude_timeextract.(names{i}) = amplitude(:,CanESM2years >= timeperiod(i,1) & CanESM2years <= timeperiod(i,2),:);
    amplitude_composite.(names{i}) = reshape(amplitude_timeextract.(names{i}),[8,50*timelength]);
    percentiles.(names{i}).low = prctile(amplitude_composite.(names{i})(latind,:),perc);
    percentiles.(names{i}).high  = prctile(amplitude_composite.(names{i})(latind,:),100-perc);
    prctile_ind.(names{i}).low = find(amplitude_composite.(names{i})(latind,:) <= percentiles.(names{i}).low); 
    prctile_ind.(names{i}).high = find(amplitude_composite.(names{i})(latind,:) >= percentiles.(names{i}).high);

    differencesstrat(i,:,:) = nanmean(CanESM24060_composite.(names{i})(:,:,prctile_ind.(names{i}).high),3) - ...
        nanmean(CanESM24060_composite.(names{i})(:,:,prctile_ind.(names{i}).low),3);
    for j = 1:size(CanESM24060_composite.(names{i}),1)
        for k = 1:size(CanESM24060_composite.(names{i}),2)
            [hstrat(i,j,k)] = ttest2(CanESM24060_composite.(names{i})(j,k,prctile_ind.(names{i}).high),CanESM24060_composite.(names{i})(j,k,prctile_ind.(names{i}).low));            
        end
    end
    
end
h (h == 1) = .05;
hstrat (hstrat == 1) = .05;
h = cat(1,h(1,:,:),h(3,:,:),h(2,:,:));
hstrat = cat(1,hstrat(1,:,:),hstrat(3,:,:),hstrat(2,:,:));
differences2 = cat(1,differences(1,:,:),differences(3,:,:),differences(2,:,:));
differences2strat = cat(1,differencesstrat(1,:,:),differencesstrat(3,:,:),differencesstrat(2,:,:));
h(:,:,10:end) = 0;
hstrat(:,:,10:end) = 0;
%% plotting

%rCanESM2 = reshape(rCanESM2,[1,size(rCanESM2)]);
cbrew = cbrewer('div','RdBu',16);         

prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];

presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
    5,[],[],2,1};


logprestick = log(prestick);

contourtitle = {['30-60{\circ}S T diff of upper and lower ',num2str(perc),'th percentiles in ',num2str(abs(lats(latind))),'{\circ}S',' 50 hPa amplitude']};       

contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','2050-2079 (high GHGs - low chlorine)'};       

subplotmaps(differences2,data.Can2ESM.longitude,log(CanESM2verens.lev./100),{'div','RdBu'},1,h,18,contourtitle2,'Longitude','Pressure (hPa)','Kelvin','on',...
    [-2,2],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/',...
    'CanESM2_Tempdiff_',num2str(abs(lats(latind))),'S_',num2str(perc),'_perc_spring'];
export_fig(filename,'-pdf');

% strat
contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','2050-2079 (low chlorine)'};       

subplotmaps(differences2strat,data.Can2ESM.longitude,log(CanESM2verens.lev./100),{'div','RdBu'},1,hstrat,18,contourtitle2,'Longitude','Pressure (hPa)','Kelvin','on',...
    [-2,2],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/',...
    'CanESM2strat_Tempdiff_',num2str(abs(lats(latind))),'S_',num2str(perc),'_perc_spring'];
export_fig(filename,'-pdf');



