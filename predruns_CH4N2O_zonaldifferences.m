% Read in manipulate CH4 and N20
clear variables


%% Read in CH4 and N2O emissions
% constants
timeseries = 1765:2500;
alpha = 65; %For polar latitudes
air = 28.97; %(g mol-1)

en = 0;
plot_EESC = 0;
cbrew = cbrewer('qual','Set1',10);
cbrew2 = cbrewer('qual','Set1',17);

% read in WMO 2011 file
filedir = '/Volumes/ExternalOne/work/data/rcp6.0_table_wmo2011';
fid = fopen(filedir); 
while ~en
    line = fgetl(fid);
    if strcmp(line(1:10),'      Year')
        en = 1;
    end
end

dataemis = fscanf(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',[21 Inf])';
dataemis(:,1) = dataemis(:,1) - .5;
emissions = dataemis(1:end,3:4);
emissionslow = emissions(dataemis >= 1955 & dataemis <= 1979,:);
emissionshigh = emissions(dataemis >= 1995 & dataemis <= 2024,:);

%% Reading in data
ens = {'highCl','lowCl'};
var = {'CH4','N2O'};

if exist('/Volumes/ExternalOne/work/data/predruns/output/CH4N2O/regpres_detrend.mat')    
    load('/Volumes/ExternalOne/work/data/predruns/output/CH4N2O/regpres_detrend.mat');
else
    for j = 1:length(ens)
        for k = 1:length(var)

            directory = ['/Volumes/ExternalOne/work/data/predruns/',var{k},'/',ens{j},'/zonalmean/'];
            files = dir([directory,'*.nc']);

            for i = 1:length(files)

                [~,data(i).(ens{j}).(var{k}),~] = Read_in_netcdf([directory,files(i).name]);

                % calculating pressure

                Pressure(i).(ens{j}).(var{k}) = 100000*repmat(data(i).(ens{j}).(var{k}).hyam,[1,size(data(i).(ens{j}).(var{k}).PS)])...
                    + repmat(data(i).(ens{j}).(var{k}).hybm,[1,size(data(i).(ens{j}).(var{k}).PS)])...
                    .* permute(repmat(data(i).(ens{j}).(var{k}).PS,[1,1,66]),[3,1,2]);

                % putting onto regular pressure
                for l = 1:size(data(i).(ens{j}).(var{k}).(var{k}),1)
                    dataregpres.(ens{j}).(var{k})(i,l,:,:) = intRegPres(squeeze(data(i).(ens{j}).(var{k}).(var{k})(l,:,:)),squeeze(Pressure(i).(ens{j}).(var{k})(:,l,:)));
                    for m = 1:12
                        dataregpres_monarr.(ens{j}).(var{k})(i,l,:,:,m) = dataregpres.(ens{j}).(var{k})(i,l,:,m:12:end);
                        % take monthly mean
                        dataregpres_yearmean.(ens{j}).(var{k}) = nanmean(dataregpres_monarr.(ens{j}).(var{k}),5);
                        % regress against CH4 or N2O
                        %dataregpres_monarr_detrend.(ens{j}).(var{k})(i,l,:,:,m) = detrend(squeeze(dataregpres_monarr.(ens{j}).(var{k})(i,l,:,:,m))')+nanmean(squeeze(dataregpres_monarr.(ens{j}).(var{k})(i,l,:,:,m))');
                    end
                end            
            end
        end
    end
    WAClat = data(1).highCl.CH4.lat;
    save('/Volumes/ExternalOne/work/data/predruns/output/CH4N2O/regpres_detrend.mat','dataregpres_yearmean','WAClat');
end

%% normalising to lowcl global average
dataregpres_yearmean.highCl.CH4 = permute(dataregpres_yearmean.highCl.CH4,[1,3,2,4]);
dataregpres_yearmean.lowCl.CH4 = permute(dataregpres_yearmean.lowCl.CH4,[1,3,2,4]);
dataregpres_yearmean.highCl.N2O = permute(dataregpres_yearmean.highCl.N2O,[1,3,2,4]);
dataregpres_yearmean.lowCl.N2O = permute(dataregpres_yearmean.lowCl.N2O,[1,3,2,4]);

%%
CH4highclnormalized =  nanmean(permute(dataregpres_yearmean.highCl.CH4 - nanmean(dataregpres_yearmean.highCl.CH4(:,:,:),3) + nanmean(dataregpres_yearmean.lowCl.CH4(:,:,:),3),[3,2,1,4]),4);
CH4lowclnormalized = nanmean(permute(dataregpres_yearmean.lowCl.CH4,[3,2,1,4]),4);
CH4difference = (nanmean(CH4highclnormalized,3) - nanmean(CH4lowclnormalized,3))./nanmean(CH4lowclnormalized,3)*100;

N2Ohighclnormalized =  nanmean(permute(dataregpres_yearmean.highCl.N2O - nanmean(dataregpres_yearmean.highCl.N2O(:,:,:),3) + nanmean(dataregpres_yearmean.lowCl.N2O(:,:,:),3),[3,2,1,4]),4);
N2Olowclnormalized = nanmean(permute(dataregpres_yearmean.lowCl.N2O,[3,2,1,4]),4);
N2Odifference = (nanmean(N2Ohighclnormalized,3) - nanmean(N2Olowclnormalized,3))./nanmean(N2Olowclnormalized,3)*100;


%% put together surface concentration
% CH4highclsurface = squeeze(nanmean(dataregpres_yearmean.highCl.CH4(:,5,48,:),2));
% CH4lowclsurface = squeeze(nanmean(dataregpres_yearmean.lowCl.CH4(:,5,48,:),2));
% 
% N2Ohighclsurface = squeeze(nanmean(dataregpres_yearmean.highCl.N2O(:,5,48,:),2));
% N2Olowclsurface = squeeze(nanmean(dataregpres_yearmean.lowCl.N2O(:,5,48,:),2));
lattoind = 48;

CH4highclsurface = squeeze(dataregpres_yearmean.highCl.CH4(:,1,lattoind,:));
CH4lowclsurface = squeeze(dataregpres_yearmean.lowCl.CH4(:,1,lattoind,:));

N2Ohighclsurface = squeeze(dataregpres_yearmean.highCl.N2O(:,1,lattoind,:));
N2Olowclsurface = squeeze(dataregpres_yearmean.lowCl.N2O(:,1,lattoind,:));


combinedCH4 = [CH4lowclsurface,CH4highclsurface]*1e6;
combinedN2O = [CH4lowclsurface,CH4highclsurface]*1e6;

%% if regressing
regressing = 1;
if regressing
    %% put together ensembles and regressing
    CH4_dataregpres_yearmean = cat(4,dataregpres_yearmean.lowCl.CH4,dataregpres_yearmean.highCl.CH4);
    N2O_dataregpres_yearmean = cat(4,dataregpres_yearmean.lowCl.N2O,dataregpres_yearmean.highCl.N2O);

    emissionsCH4 = [emissionslow(:,1);emissionshigh(:,1)];
    emissionsN2O = [emissionslow(:,2);emissionshigh(:,2)];

    %%
    for i  = 1:size(CH4_dataregpres_yearmean,1)
        for j = 1:size(CH4_dataregpres_yearmean,2)
            for k = 1:size(CH4_dataregpres_yearmean,3)
%                 b.CH4(i,j,k,:) = regress(squeeze(CH4_dataregpres_yearmean(i,j,k,:)),[ones(size(CH4_dataregpres_yearmean,4),1),emissionsCH4]);
%                 b.N2O(i,j,k,:) = regress(squeeze(N2O_dataregpres_yearmean(i,j,k,:)),[ones(size(N2O_dataregpres_yearmean,4),1),emissionsN2O]);
                
                b.CH4(i,j,k,:) = regress(squeeze(CH4_dataregpres_yearmean(i,j,k,:)),[ones(size(CH4_dataregpres_yearmean,4),1),combinedCH4(i,:)']);
                b.N2O(i,j,k,:) = regress(squeeze(N2O_dataregpres_yearmean(i,j,k,:)),[ones(size(N2O_dataregpres_yearmean,4),1),combinedN2O(i,:)']);
                
                % removing trend from each ensemble e
                CH4detrend(i,j,k,:) = squeeze(CH4_dataregpres_yearmean(i,j,k,:)) - b.CH4(i,j,k,2).*combinedCH4(i,:)';% - b.CH4(i,j,k,1);% 
                N2Odetrend(i,j,k,:) = squeeze(N2O_dataregpres_yearmean(i,j,k,:)) - b.N2O(i,j,k,2).*combinedN2O(i,:)';% - b.N2O(i,j,k,1);%
            end
        end
    end

    CH4_dataregpres_yearmean_globalmean = permute(CH4detrend,[1,4,2,3]);
    N2O_dataregpres_yearmean_globalmean = permute(N2Odetrend,[1,4,2,3]);
    CH4_dataregpres_yearmean_globalmean2 = permute(CH4_dataregpres_yearmean_globalmean - nanmean(CH4_dataregpres_yearmean_globalmean(:,:,:),3),[1,3,4,2]);
    N2O_dataregpres_yearmean_globalmean2 = permute(N2O_dataregpres_yearmean_globalmean - nanmean(N2O_dataregpres_yearmean_globalmean(:,:,:),3),[1,3,4,2]);
    CH4_dataregpres_yearmean_globalmean = CH4detrend;
    N2O_dataregpres_yearmean_globalmean = N2Odetrend;
end
%% taking ensemble average
CH4ens = squeeze(nanmean(CH4_dataregpres_yearmean_globalmean,1));
N2Oens = squeeze(nanmean(N2O_dataregpres_yearmean_globalmean,1));

%% extracting relative years
clearvars CH4enstoplot N2Oenstoplot
ind = 0;
timeperiod = [1955:1979,1995:2024];
if ind
    lowyear = 1960;
    titextlow = num2str(lowyear);
    highyear = 2015;
    titexthigh = num2str(highyear);
    index1975 = find(timeperiod == lowyear);
index2015 = find(timeperiod == highyear);
else    
    titextlow2 = [];
    titexthigh2 = [];
    titextlow = 'LowCl';
    titexthigh = 'HighCl';
end

percent = 1;


if ind
    CH4enstoplot(1,:,:) = CH4ens(:,:,index2015);
    CH4enstoplot(2,:,:) = CH4ens(:,:,index1975);
    CH4enstoplot(3,:,:) = (CH4enstoplot(1,:,:) - CH4enstoplot(2,:,:))./CH4enstoplot(2,:,:)*100;

    N2Oenstoplot(1,:,:) = N2Oens(:,:,index2015);
    N2Oenstoplot(2,:,:) = N2Oens(:,:,index1975);
    N2Oenstoplot(3,:,:) = (N2Oenstoplot(1,:,:) - N2Oenstoplot(2,:,:))./N2Oenstoplot(2,:,:)*100;
else
    CH4enstoplot(1,:,:) = nanmean(CH4ens(:,:,26:end),3);
    CH4enstoplot(2,:,:) = nanmean(CH4ens(:,:,1:25),3);
    CH4enstoplot(3,:,:) = (CH4enstoplot(1,:,:) - CH4enstoplot(2,:,:))./abs(CH4enstoplot(2,:,:))*100;
    
%     CH4enstoplot(1,:,:) = nanmean(CH4highclnormalized,3);
%     CH4enstoplot(1,:,:) = nanmean(CH4lowclnormalized,3);
%     CH4enstoplot(3,:,:) = (nanmean(CH4highclnormalized,3) - nanmean(CH4lowclnormalized,3))./abs(nanmean(CH4lowclnormalized,3)).*100;

    N2Oenstoplot(1,:,:) = nanmean(N2Oens(:,:,26:end),3);
    N2Oenstoplot(2,:,:) = nanmean(N2Oens(:,:,1:25),3);
    N2Oenstoplot(3,:,:) = (N2Oenstoplot(1,:,:) - N2Oenstoplot(2,:,:))./abs(N2Oenstoplot(2,:,:))*100;


% N2Oenstoplot(1,:,:) = nanmean(N2Ohighclnormalized,3);
% N2Oenstoplot(1,:,:) = nanmean(N2Olowclnormalized,3);
% N2Oenstoplot(3,:,:) = (nanmean(N2Ohighclnormalized,3) - nanmean(N2Olowclnormalized,3))./abs(nanmean(N2Olowclnormalized,3)).*100;end
end
if ~percent 
    N2Oenstoplot = N2Oenstoplot.*1e6;
    CH4enstoplot = CH4enstoplot.*1e6;
else
    
    CH4enstoplot(1:2,:,:) = CH4enstoplot(1:2,:,:).*1e6;
    N2Oenstoplot(1:2,:,:) = N2Oenstoplot(1:2,:,:).*1e6;
end

CH4enstoplot = permute(CH4enstoplot,[1,3,2]);
N2Oenstoplot = permute(N2Oenstoplot,[1,3,2]);

regpres = [1000;975;950;925;900;875;850;825;800;775;750;700;650;600;550;500;450;400;350;...
        300;250;225;200;175;150;125;100;70;50;30;20;10;7;5;3;2;1;0.7;0.5;0.3;0.2;0.1;.07;.05;...
        .03;.02;.01;.007;.005;.003;.002;.001];

%% plotting all three

prestick = regpres;
    presticklabel = {1000,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],300,[],[],[],[],[],[],100,[],[],30,[],10,[],[],...
        3,[],1,[],[],.3,[],.1,[],[],.03,[],.01,[],[],.003,[],.001};
    logprestick = log(prestick);

titles = {['HighCl ',titexthigh2],['LowCl ',titextlow2]};
mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions and global yearly means'};
%mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions'};

fig = subplotmaps(CH4enstoplot(1:2,:,:),WAClat,log(regpres),{'div','RdBu'},1,[],14,titles,'Latitude','Pressure (hPa)','ppm','on',...
    [-.1 .1],22,-90:30:90,-90:30:90,...
    flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'-',0,'');


filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/CH4_',titexthigh,'_',titextlow,'.eps',];
%export_fig(filename,'-eps');
print(filename,'-depsc');
mtit = {'Ensemble mean N2O concentrations normalized to N2O emissions and global yearly means'};
%mtit = {'Ensemble mean N2O concentrations normalized to N2O emissions'};

fig = subplotmaps(N2Oenstoplot(1:2,:,:),WAClat,log(regpres),{'div','RdBu'},1,[],14,titles,'Latitude','Pressure (hPa)','ppm','on',...
    [-.1 .1],22,-90:30:90,-90:30:90,...
    flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'-',0,'');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/N2O_',titexthigh,'_',titextlow,'.eps',];
%export_fig(filename,'-eps');
print(filename,'-depsc');


%% plot differences
titles = {'CH4','N2O'};
mtit = {'Ensemble mean HighCl - LowCl differences'};
fig = subplotmaps(cat(1,CH4enstoplot(3,:,:),N2Oenstoplot(3,:,:)),WAClat,log(regpres),{'div','RdBu'},1,[],14,titles,'Latitude','Pressure (hPa)','percent','on',...
    [-20 20],22,-90:30:90,-90:30:90,...
    flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'none',0,'');


filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/Diffs_',titexthigh,'_',titextlow,'.eps',];
%export_fig(filename,'-eps');
print(filename,'-depsc');





