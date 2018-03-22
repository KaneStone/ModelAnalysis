% LENS percentiles
% Read in and analyse all zonal wave data
clear all
directory = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/';
files = dir([directory,'*.mat']);
% user inputs

LENSyears = 1920:2100;
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
load('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/lonvert/LENSvert3060S.mat');

%% variation in minimum and maximum longitudes
minall = cat(3,data.LENS.minlongitude_m(:).m);
minallCan = cat(3,data.Can2ESM.minlongitude_m(51:100).m);
for i = 1:181
    minallspringstd(i) = std(squeeze(nanmean(minall(9:11,4,i:181:end))));    
end
for i = 1:150
    minallspringstdCan(i) = std(squeeze(nanmean(minallCan(9:11,4,i:181:end))));
end

predmaxlong = std(squeeze(nanmean(WACCM_predruns.var_maxlong.lowGHG(4,9:11,:),2)));
predminlong = std(squeeze(nanmean(WACCM_predruns.var_minlong.lowGHG(4,9:11,:),2)));

% figure;
% plot(1920:2100,minallspringstd);
% hold on
% plot(1951:2100,minallspringstdCan);

%% find percentiles
names = {'highcl','highGHGlowcl','lowcl'};
lats = [-50,-55,-60,-65,-70,-75,-80,-85];
latind = 4;
perc = 20;
timeperiod = [1995,2024;2050,2079;1955,1979];    
for i = 1:size(timeperiod,1)
    timelength = timeperiod(i,2) - timeperiod(i,1)+1;
    LENStimeextract.(names{i}) = LENSverens.spring(:,:,:,LENSyears >= timeperiod(i,1) & LENSyears <= timeperiod(i,2));
    LENS4060_composite.(names{i}) = reshape(permute(LENStimeextract.(names{i}),[2,3,4,1]),...
        [size(LENStimeextract.(names{i}),2),size(LENStimeextract.(names{i}),3),size(LENStimeextract.(names{i}),4)*size(LENStimeextract.(names{i}),1)]); 
    amplitude = cat(3,data.LENS.amplitude_m_spr(:).m);
    amplitude_timeextract.(names{i}) = amplitude(:,LENSyears >= timeperiod(i,1) & LENSyears <= timeperiod(i,2),:);
    amplitude_composite.(names{i}) = reshape(amplitude_timeextract.(names{i}),[8,30*timelength]);
    percentiles.(names{i}).low = prctile(amplitude_composite.(names{i})(latind,:),perc);
    percentiles.(names{i}).high  = prctile(amplitude_composite.(names{i})(latind,:),100-perc);
    prctile_ind.(names{i}).low = find(amplitude_composite.(names{i})(latind,:) <= percentiles.(names{i}).low); 
    prctile_ind.(names{i}).high = find(amplitude_composite.(names{i})(latind,:) >= percentiles.(names{i}).high);

    differences(i,:,:) = nanmean(LENS4060_composite.(names{i})(:,:,prctile_ind.(names{i}).high),3) - ...
        nanmean(LENS4060_composite.(names{i})(:,:,prctile_ind.(names{i}).low),3);      
    for j = 1:size(LENS4060_composite.(names{i}),1)
        for k = 1:size(LENS4060_composite.(names{i}),2)
            [h(i,j,k)] = ttest2(LENS4060_composite.(names{i})(j,k,prctile_ind.(names{i}).high),LENS4060_composite.(names{i})(j,k,prctile_ind.(names{i}).low));            
        end
    end
end
h (h == 1) = .05;
% %% two-sampled t-test
% for i = 1:length(longitude)-1
%     for j = 1:length(latitude)        
%         [highcl.h(j,i),highcl.p(j,i),highcl.ci(j,i,:),~] = ttest2(squeeze(varextract.highcl.lowind(:,j,i)),squeeze(varextract.highcl.highind(:,j,i)),'Alpha',.05,'Tail','both');
%         [lowcl.h(j,i),lowcl.p(j,i),lowcl.ci(j,i,:),~] = ttest2(squeeze(varextract.lowcl.lowind(:,j,i)),squeeze(varextract.lowcl.highind(:,j,i)),'Alpha',.05,'Tail','both');
%     end
% end
differences2 = cat(1,differences(1,:,:),differences(3,:,:),differences(2,:,:));
h = cat(1,h(1,:,:),h(3,:,:),h(2,:,:));

%% plotting
h(:,:,1:14) = 0;
%rLENS = reshape(rLENS,[1,size(rLENS)]);
cbrew = cbrewer('div','RdBu',16);         

prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];

presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
    5,[],[],2,1};


logprestick = log(prestick);

contourtitle = {['30-60{\circ}S T diff of upper and lower ',num2str(perc),'th percentiles in ',num2str(abs(lats(latind))),'{\circ}S',' 50 hPa amplitude']};       
contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','2050-2079 (high GHGs - low chlorine)'};       

differences2 (differences2 <= -2) = -2;

subplotmaps(differences2,data.LENS.longitude,log(LENSverens.lev),{'div','RdBu'},1,h,18,contourtitle2,'Longitude','Pressure (hPa)','Kelvin','on',...
    [-2,2],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/',...
    'LENS_Tempdiff_',num2str(abs(lats(latind))),'S_',num2str(perc),'_perc_spring'];
export_fig(filename,'-pdf');


