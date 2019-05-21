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


%% calculate and plot raw ensemble averages
calcraw = 1
if calcraw
    highclCH4 = permute(dataregpres_yearmean.highCl.CH4,[2,3,1,4]);
    highclCH4 = nanmean(highclCH4(:,:,:),3)*1e6;

    lowclCH4 = permute(dataregpres_yearmean.lowCl.CH4,[2,3,1,4]);
    lowclCH4 = nanmean(lowclCH4(:,:,:),3)*1e6;

    difference = (highclCH4 - lowclCH4)./lowclCH4*100;

    regpres = [1000;975;950;925;900;875;850;825;800;775;750;700;650;600;550;500;450;400;350;...
            300;250;225;200;175;150;125;100;70;50;30;20;10;7;5;3;2;1;0.7;0.5;0.3;0.2;0.1;.07;.05;...
            .03;.02;.01;.007;.005;.003;.002;.001];

    prestick = regpres;
        presticklabel = {1000,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],300,[],[],[],[],[],[],100,[],[],30,[],10,[],[],...
            3,[],1,[],[],.3,[],.1,[],[],.03,[],.01,[],[],.003,[],.001};
        logprestick = log(prestick);

    titles = {['HighCl '],['LowCl ']};
    mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions and global yearly means'};
    %mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions'};

    fig = subplotmaps(permute(cat(3,highclCH4,lowclCH4),[3,1,2]),WAClat,log(regpres),{'seq','YlOrBr'},0,[],14,titles,'Latitude','Pressure (hPa)','ppm','on',...
        [0 .6],21,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/CH4_','rawdata.eps',];

    print(filename,'-depsc');

    % plot differences
    titles = {'CH4','N2O'};
    mtit = {'Raw data, Ensemble mean HighCl - LowCl differences'};
    fig = subplotmaps(permute(repmat(difference,[1,1,1]),[3,1,2]),WAClat,log(regpres),{'div','RdBu'},1,[],14,titles,'Latitude','Pressure (hPa)','percent','on',...
        [-30 30],22,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'none',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/Diffs_raw','.eps',];
    %export_fig(filename,'-eps');
    print(filename,'-depsc');
end

%% normalized highcl lat height to lowcl global mean.
norm = 1;
if norm
    highclCH4 = permute(dataregpres_yearmean.highCl.CH4,[2,3,1,4]);
    lowclCH4 = permute(dataregpres_yearmean.lowCl.CH4,[2,3,1,4]);
    
    highclCH4norm = highclCH4 - nanmean(highclCH4(:)) + nanmean(lowclCH4(:));
    highclCH4norm = nanmean(highclCH4norm(:,:,:),3)*1e6;
    lowclCH4norm = nanmean(lowclCH4(:,:,:),3)*1e6;
    
    difference = (highclCH4norm - lowclCH4norm)./lowclCH4norm*100;

    regpres = [1000;975;950;925;900;875;850;825;800;775;750;700;650;600;550;500;450;400;350;...
            300;250;225;200;175;150;125;100;70;50;30;20;10;7;5;3;2;1;0.7;0.5;0.3;0.2;0.1;.07;.05;...
            .03;.02;.01;.007;.005;.003;.002;.001];

    prestick = regpres;
        presticklabel = {1000,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],300,[],[],[],[],[],[],100,[],[],30,[],10,[],[],...
            3,[],1,[],[],.3,[],.1,[],[],.03,[],.01,[],[],.003,[],.001};
        logprestick = log(prestick);

    titles = {['HighCl '],['LowCl ']};
    mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions and global yearly means'};
    %mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions'};

    fig = subplotmaps(permute(cat(3,highclCH4norm,lowclCH4norm),[3,1,2]),WAClat,log(regpres),{'seq','YlOrBr'},0,[],14,titles,'Latitude','Pressure (hPa)','ppm','on',...
        [0 .6],21,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/CH4_','normdata.eps',];

    print(filename,'-depsc');

    % plot differences
    titles = {'CH4','N2O'};
    mtit = {'Normalized to gliobal mean, Ensemble mean HighCl - LowCl differences'};
    fig = subplotmaps(permute(repmat(difference,[1,1,1]),[3,1,2]),WAClat,log(regpres),{'div','RdBu'},1,[],14,titles,'Latitude','Pressure (hPa)','percent','on',...
        [-100 100],22,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'none',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/Diffs_norm','.eps',];
    %export_fig(filename,'-eps');
    print(filename,'-depsc');
end

%% normalized highcl lat height to lowcl meridional mean.
norm2 = 1;
if norm2
    highclCH4 = permute(dataregpres_yearmean.highCl.CH4,[3,2,1,4]);
    lowclCH4 = permute(dataregpres_yearmean.lowCl.CH4,[3,2,1,4]);
    
    highclCH4norm = highclCH4 - nanmean(highclCH4(:,:),2) + nanmean(lowclCH4(:,:),2);
    highclCH4norm = permute(nanmean(highclCH4norm(:,:,:),3)*1e6,[2,1]);
    lowclCH4norm = permute(nanmean(lowclCH4(:,:,:),3)*1e6,[2,1]);
    
    difference = (highclCH4norm - lowclCH4norm)./lowclCH4norm*100;

    regpres = [1000;975;950;925;900;875;850;825;800;775;750;700;650;600;550;500;450;400;350;...
            300;250;225;200;175;150;125;100;70;50;30;20;10;7;5;3;2;1;0.7;0.5;0.3;0.2;0.1;.07;.05;...
            .03;.02;.01;.007;.005;.003;.002;.001];

    prestick = regpres;
        presticklabel = {1000,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],300,[],[],[],[],[],[],100,[],[],30,[],10,[],[],...
            3,[],1,[],[],.3,[],.1,[],[],.03,[],.01,[],[],.003,[],.001};
        logprestick = log(prestick);

    titles = {['HighCl '],['LowCl ']};
    mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions and global yearly means'};
    %mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions'};

    fig = subplotmaps(permute(cat(3,highclCH4norm,lowclCH4norm),[3,1,2]),WAClat,log(regpres),{'seq','YlOrBr'},0,[],14,titles,'Latitude','Pressure (hPa)','ppm','on',...
        [0 .6],21,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/CH4_','normmeddata.eps',];

    print(filename,'-depsc');

    % plot differences
    titles = {'CH4','N2O'};
    mtit = {'Normalized to individual pressure global mean, Ensemble mean HighCl - LowCl differences'};
    fig = subplotmaps(permute(repmat(difference,[1,1,1]),[3,1,2]),WAClat,log(regpres),{'div','RdBu'},1,[],14,titles,'Latitude','Pressure (hPa)','percent','on',...
        [-30 30],22,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'none',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/Diffs_normmed','.eps',];
    %export_fig(filename,'-eps');
    print(filename,'-depsc');
end

%% remove CH4 trend through regression
reg = 1;
if reg
    CH4_dataregpres_yearmean = squeeze(nanmean(cat(4,dataregpres_yearmean.lowCl.CH4,dataregpres_yearmean.highCl.CH4)));
    CH4_dataregpres_yearmean2 = squeeze(cat(4,dataregpres_yearmean.lowCl.CH4,dataregpres_yearmean.highCl.CH4));
    CH4predictor = squeeze(nanmean(CH4_dataregpres_yearmean(:,1,:)))*1e6;
    CH4predictor(:,2) = ones(length(CH4predictor),1);

    for i = 1:size(CH4_dataregpres_yearmean2,1)
        for j = 1:size(CH4_dataregpres_yearmean2,2)
            for k = 1:size(CH4_dataregpres_yearmean2,3)
                b = regress(squeeze(CH4_dataregpres_yearmean2(i,j,k,:)),CH4predictor);
                CH4regressed(i,j,k,:) = squeeze(CH4_dataregpres_yearmean2(i,j,k,:)) - b(1).*CH4predictor(:,1);%- b(2).*CH4predictor(:,2);
            end
        end
    end

    CH4mean = squeeze(nanmean(CH4regressed));
    CH4highcl = nanmean(CH4mean(:,:,26:end),3);
    CH4lowcl = nanmean(CH4mean(:,:,1:25),3);
    difference = (CH4highcl - CH4lowcl)./abs(CH4lowcl)*100;

    regpres = [1000;975;950;925;900;875;850;825;800;775;750;700;650;600;550;500;450;400;350;...
                300;250;225;200;175;150;125;100;70;50;30;20;10;7;5;3;2;1;0.7;0.5;0.3;0.2;0.1;.07;.05;...
                .03;.02;.01;.007;.005;.003;.002;.001];

    prestick = regpres;
        presticklabel = {1000,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],300,[],[],[],[],[],[],100,[],[],30,[],10,[],[],...
            3,[],1,[],[],.3,[],.1,[],[],.03,[],.01,[],[],.003,[],.001};
        logprestick = log(prestick);

    titles = {['HighCl '],['LowCl ']};
    mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions and global yearly means'};
    %mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions'};

    fig = subplotmaps(permute(cat(3,CH4highcl,CH4lowcl),[3,1,2]),WAClat,log(regpres),{'seq','YlOrBr'},0,[],14,titles,'Latitude','Pressure (hPa)','ppm','on',...
        [0 .6],21,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/CH4_','regresseddata.eps',];

    print(filename,'-depsc');

    % plot differences
    titles = {'CH4','N2O'};
    mtit = {'Normalized through linear trend line regression, Ensemble mean HighCl - LowCl differences'};
    fig = subplotmaps(permute(repmat(difference,[1,1,1]),[3,1,2]),WAClat,log(regpres),{'div','RdBu'},1,[],14,titles,'Latitude','Pressure (hPa)','percent','on',...
        [-30 30],22,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'none',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/Diffs_regressed','.eps',];
    %export_fig(filename,'-eps');
    print(filename,'-depsc');
end

%%
det = 1;
if det


    CH4_dataregpres_yearmean = squeeze(nanmean(cat(4,dataregpres_yearmean.lowCl.CH4,dataregpres_yearmean.highCl.CH4)));
    CH4_dataregpres_yearmean2 = squeeze(cat(4,dataregpres_yearmean.lowCl.CH4,dataregpres_yearmean.highCl.CH4));
    CH4predictor = [1:25,41:70]';
    CH4predictor(:,2) = ones(length(CH4predictor),1);

    for i = 1:size(CH4_dataregpres_yearmean2,1)
        for j = 1:size(CH4_dataregpres_yearmean2,2)
            for k = 1:size(CH4_dataregpres_yearmean2,3)
                b(i,j,k,:) = regress(squeeze(CH4_dataregpres_yearmean2(i,j,k,26:end)),CH4predictor(26:end,:));
                b2(i,j,k,:) = regress(squeeze(CH4_dataregpres_yearmean2(i,j,k,1:25)),CH4predictor(1:25,:));
                %CH4regressed(i,j,k,:) = squeeze(CH4_dataregpres_yearmean2(i,j,k,:)) - b(1).*CH4predictor(:,1);%- b(2).*CH4predictor(:,2);
                CH4regressed_highcl(i,j,k,:) = squeeze(CH4_dataregpres_yearmean2(i,j,k,26:end)) - b(i,j,k,1).*CH4predictor(26:end,1);%- b(2).*CH4predictor(:,2);
                CH4regressed_lowcl(i,j,k,:) = squeeze(CH4_dataregpres_yearmean2(i,j,k,1:25)) - b2(i,j,k,1).*CH4predictor(1:25,1);%- b(2).*CH4predictor(:,2);
            end
        end
    end

%     CH4mean = squeeze(nanmean(CH4regressed));
%     CH4highcl = nanmean(CH4mean(:,:,26:end),3);
%     CH4lowcl = nanmean(CH4mean(:,:,1:25),3);
%     difference = (CH4highcl - CH4lowcl)./abs(CH4lowcl)*100;
    
    CH4meanhighcl = squeeze(nanmean(CH4regressed_highcl));
    CH4meanlowcl = squeeze(nanmean(CH4regressed_lowcl));
    CH4highcl = nanmean(CH4mean,3);
    CH4lowcl = nanmean(CH4mean,3);
    difference = (CH4highcl - CH4lowcl)./abs(CH4lowcl)*100;
    

    regpres = [1000;975;950;925;900;875;850;825;800;775;750;700;650;600;550;500;450;400;350;...
                300;250;225;200;175;150;125;100;70;50;30;20;10;7;5;3;2;1;0.7;0.5;0.3;0.2;0.1;.07;.05;...
                .03;.02;.01;.007;.005;.003;.002;.001];

    prestick = regpres;
        presticklabel = {1000,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],300,[],[],[],[],[],[],100,[],[],30,[],10,[],[],...
            3,[],1,[],[],.3,[],.1,[],[],.03,[],.01,[],[],.003,[],.001};
        logprestick = log(prestick);

    titles = {['HighCl '],['LowCl ']};
    mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions and global yearly means'};
    %mtit = {'Ensemble mean CH4 concentrations normalized to CH4 emissions'};

    fig = subplotmaps(permute(cat(3,CH4highcl,CH4lowcl),[3,1,2]),WAClat,log(regpres),{'seq','YlOrBr'},0,[],14,titles,'Latitude','Pressure (hPa)','ppm','on',...
        [0 .6],21,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'-',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/CH4_','detrenddata.eps',];

    print(filename,'-depsc');

    % plot differences
    titles = {'CH4','N2O'};
    mtit = {'Ensemble mean HighCl - LowCl differences'};
    fig = subplotmaps(permute(repmat(difference,[1,1,1]),[3,1,2]),WAClat,log(regpres),{'div','RdBu'},1,[],14,titles,'Latitude','Pressure (hPa)','percent','on',...
        [-5 5],22,-90:30:90,-90:30:90,...
        flipud(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.01) log(1000)],1,'none',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/CH4N2O/Diffs_detrend','.eps',];
    %export_fig(filename,'-eps');
    print(filename,'-depsc');
end