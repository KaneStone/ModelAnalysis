%predruns_ampdifontrop

% plot toz and temperature longitude lines
clear all
%% Read in temperature and toz
tozvar = 'toz';
lats = [-90 -60];
detrend_ozone = 0;
sig = .05;

%% Read in highcl Temperature at 50hPa
var = 'T';
ClLevel = 'highCl';
timeperiodhigh = [1995,2024];%[1955,1975]

vardirectory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/50hPa/'];
varfilespast = dir([vardirectory,'*.nc']);
[data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl]...
    = predruns_ReadInlayer(vardirectory,varfilespast,var,timeperiodhigh,lats,1);

%% Read in lowGHG Temperature at 50hPa
var = 'T';
ClLevel = 'lowGHG';
timeperiodhigh = [1995,2024];%[1955,1975]

vardirectory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/50hPa/'];
varfilesGHG = dir([vardirectory,'*.nc']);
[data.lowGHG,years.lowGHG,composite.lowGHG,dataMonthArrange.lowGHG]...
    = predruns_ReadInlayer(vardirectory,varfilesGHG,var,timeperiodhigh,lats,1);

%% Read in lowcl Temperature at 50hPa
var = 'T';
ClLevel = 'lowCl';
timeperiodhigh = [1955,1979];%[1955,1975]

vardirectory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/50hPa/'];
varfilespast = dir([vardirectory,'*.nc']);
[data.lowcl,years.lowcl,composite.lowcl,dataMonthArrange.lowcl]...
    = predruns_ReadInlayer(vardirectory,varfilespast,var,timeperiodhigh,lats,1);

fieldnames = fields(data);

latitudes = data.lowcl.lat;
longitudes = data.lowcl.lon;

%% Read in zonal average 4060S data

wa4060S.highcl = load('/Volumes/MyBook/work/data/predruns/output/zonalanalysis/highCl_standard_4060Swa.mat');
wa4060S.lowcl = load('/Volumes/MyBook/work/data/predruns/output/zonalanalysis/lowCl_standard_4060Swa.mat');
wa4060S.lowGHG = load('/Volumes/MyBook/work/data/predruns/output/zonalanalysis/lowGHG_standard_4060Swa.mat');
regpres = wa4060S.lowcl.regpres;

%% create composite for correlations
month = 10;
for i = 1:length(fieldnames)
    % month arrange first
    for j = 1:length(wa4060S.(fieldnames{i}).standardwa)
        for k = 1:12
            montharrange.(fieldnames{i})(j).T(:,:,k,:) = wa4060S.(fieldnames{i}).standardwa(j).T(:,:,k:12:end);
        end
    end
    wa4060Scomposite.(fieldnames{i}).months = cat(4,montharrange.(fieldnames{i})(:).T);    
    wa4060Scomposite.(fieldnames{i}).spring = squeeze(nanmean(wa4060Scomposite.(fieldnames{i}).months(:,:,9:11,:),3));    
end

%% extracting data to plot
perc = 20;
lats = [-50,-55,-60,-65,-70,-75,-80,-85];
latindforperc = 4;
for i = 1:length(fieldnames)    
    for m = 1:12
        % extracting latitude and month from toz and temperature data
        for l = 1:length(lats)
            [~,latind] = min(abs(lats(l)-latitudes));
            for j = 1:length(data.(fieldnames{i}))                                
                var_extract.(fieldnames{i})(j).d(l,:,m,:) = squeeze(data.(fieldnames{i})(j).(var)(:,latind,1,m:12:end));
            end
        end
    end
    for j = 1:length(data.(fieldnames{i}))
        var_extract.(fieldnames{i})(j).spring = squeeze(nanmean(var_extract.(fieldnames{i})(j).d(:,:,9:11,:),3));
    end
    var_extract_composite.(fieldnames{i}).spring = cat(3,var_extract.(fieldnames{i})(:).spring);
    amplitude_composite.(fieldnames{i}).spring = squeeze(max(var_extract_composite.(fieldnames{i}).spring(latindforperc,:,:),[],2) - ...
        min(var_extract_composite.(fieldnames{i}).spring(latindforperc,:,:),[],2));
    percentiles.(fieldnames{i}).spring.low = prctile(amplitude_composite.(fieldnames{i}).spring,perc);
    percentiles.(fieldnames{i}).spring.high = prctile(amplitude_composite.(fieldnames{i}).spring,100-perc);
    prctile_ind.(fieldnames{i}).spring.low = find(amplitude_composite.(fieldnames{i}).spring <= percentiles.(fieldnames{i}).spring.low); 
    prctile_ind.(fieldnames{i}).spring.high = find(amplitude_composite.(fieldnames{i}).spring >= percentiles.(fieldnames{i}).spring.high); 
    
    differences.spring(i,:,:) = nanmean(wa4060Scomposite.(fieldnames{i}).spring(:,:,prctile_ind.(fieldnames{i}).spring.high),3) - ...
        nanmean(wa4060Scomposite.(fieldnames{i}).spring(:,:,prctile_ind.(fieldnames{i}).spring.low),3);
    
    
end

%% take differences between highGHG and lowGHG (highcl) (need to normalize)
%plot([nanmean(amplitude_composite.lowGHG.spring),nanmean(amplitude_composite.highcl.spring),nanmean(amplitude_composite.lowcl.spring)])
clresponse = reshape(nanmean(wa4060Scomposite.lowGHG.spring,3) - nanmean(wa4060Scomposite.lowcl.spring,3),[1,144,52]);

differences2(1,:,:) = nanmean(wa4060Scomposite.lowGHG.spring(:,:,prctile_ind.lowGHG.spring.high),3) - ...
        nanmean(wa4060Scomposite.lowcl.spring(:,:,prctile_ind.lowcl.spring.low),3);

%% find differences
cbrew = cbrewer('div','RdBu',16);         

prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];

presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
    5,[],[],2,1};


logprestick = log(prestick);

contourtitle = {['40-60{\circ}S T diff of upper and lower ',num2str(perc),'th percentiles in ',num2str(abs(lats(latindforperc))),'{\circ}S',' 50 hPa amplitude']};       
contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','1995-2024 (high chlorine - low GHG)'};       

differences.spring (differences.spring <= -2) = -2;

subplotmaps(differences.spring,longitudes,log(regpres),{'div','RdBu'},1,[],18,contourtitle2,'Longitude','Pressure (hPa)','Kelvin','on',...
    [-2,2],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/','Tempdiff_',num2str(abs(lats(latindforperc))),'S_',num2str(perc),'_perc_spring'];
export_fig(filename,'-pdf');

%%

subplotmaps(differences2,longitudes,log(regpres),{'div','RdBu'},1,[],18,contourtitle2,'Longitude','','Correlation','on',...
    [-2,2],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/clresponse_spring'];
export_fig(filename,'-pdf');

subplotmaps(clresponse,longitudes,log(regpres),{'div','RdBu'},1,[],18,contourtitle2,'Longitude','','Correlation','on',...
    [-2,2],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');