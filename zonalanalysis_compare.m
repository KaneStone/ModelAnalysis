% Read in and analyse all zonal wave data
clear all
directory = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/';
files = dir([directory,'*.mat']);
plottingens = 0;
plottingens2 = 0;
% user inputs

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

%% 50hPa temperature amplitude
coorwithozone = 0;
if coorwithozone
    %% read in LENS prescribed ozone

    [~,LENSozone,~] = Read_in_netcdf('/Volumes/MyBook/work/data/LENS/O3/9030S_50hPa_b.e11.BRCP85C5CNBDRD.f09_g16.034.cam.h0.ozone.192001-210012.nc');

    %% manipulate
    lat = [-50:-5:-85];
    LENSozone.years = 1920:2100;
    for j = 1:length(lat)
            [~,latind] = min(abs(lat(j) - LENSozone.lat));
        for i = 1:12            
            LENSozonemonthly(:,j,:,i) = squeeze(LENSozone.ozone(:,latind,i:12:end));
        end
    end


    %%
    latind = 1;
    LENSozonespring = nanmean(LENSozonemonthly(:,:,:,9:11),4);
    LENSozonespringens = repmat(LENSozonespring,1,1,30);

    LENSozonespringens_amplitude = squeeze(max(LENSozonespringens,[],1) - min(LENSozonespringens,[],1));
    LENSozonespring_amplitude = squeeze(max(LENSozonespring,[],1) - min(LENSozonespring,[],1));
    %% smoothing temperature
    for i = 1:30    
        for j = 1:12
            for k = 1:8
                LENStempamp_smooth(i,j,k,:) = smooth(squeeze(data.LENS.amplitude_m(i).m(j,k,:)),10,'moving');
            end
        end
    end

    %% temperature
    for i = 1:30
        temp = squeeze(nanmean(LENStempamp_smooth(i,9:11,:,:),2));    
        temp2(i,:,:) = temp;
        rm(i) = corr(detrend(squeeze(temp(latind,40:75))'),detrend(LENSozonespring_amplitude(latind,40:75)'));
        rm3(i) = corr(detrend(squeeze(temp(latind,6:35))'),detrend(LENSozonespring_amplitude(latind,6:35)'));
        rm2(i) = corr(detrend(squeeze(temp(latind,81:170))'),detrend(LENSozonespring_amplitude(latind,81:170)'));
        %plot(squeeze(temp(4,:)))
        hold on
    end

    pastcombine = squeeze(permute(temp2(:,latind,40:75),[2,1,3]))';
    pastcombine_sinvec = pastcombine(:);

    pastcombine2 = squeeze(permute(temp2(:,latind,81:175),[2,1,3]))';
    pastcombine_sinvec2 = pastcombine(:);

    for i = 1:30
        pastcombine_detrend(:,i) = detrend(pastcombine(:,i));
        pastcombine_detrend2(:,i) = detrend(pastcombine2(:,i));
    end

    [r4,p4] = corr(pastcombine_detrend(:),repmat(LENSozonespring_amplitude(latind,40:75),[1,30])');
    [r5,p5] = corr(pastcombine_detrend2(:),repmat(LENSozonespring_amplitude(latind,81:175),[1,30])');


    plot(pastcombine_detrend(:,1),LENSozonespring_amplitude(latind,40:75));

    figure;
    plot(detrend(LENSozonespring_amplitude(latind,40:75)')*1e7)
    hold on
    plot(detrend(squeeze(nanmean(temp2(:,latind,40:75),1))))
    
end
%% data quality 
data.Can2ESM.ensyears = [];
data.Can2ESM.ensyears(1).y = data.Can2ESM.years(1).y(1:852);  
data.Can2ESM.ensyears(2).y = data.Can2ESM.years(2).y;  
data.Can2ESM.ensyears(3).y = data.Can2ESM.years(1).y(1:852);  
data.Can2ESM.ensyears(4).y = data.Can2ESM.years(2).y;  

%% taking linear trends
timeperiods.(nameofdata{1}) = [1960,2001,2100];
timeperiods.(nameofdata{2}) = [1960,2001,2100];
timeperiods.(nameofdata{3}) = [1979,2001,2016];
timeperiods.(nameofdata{4}) = [1960,2001,2100];
timeperiods.(nameofdata{5}) = [1960,2001,2100];

for i = 1:length(nameofdata)-1
    if strcmp(nameofdata{i},'ERA')        
        ln = 1;
        trendvariablestoplot = {'maxlongitude_spr','minlongitude_spr','amplitude_spr'};
    else
        trendvariablestoplot = {'maxlongitude_spr2','minlongitude_spr2','amplitude_spr2'};
        fieldstrends.(nameofdata{i}) = fieldnames(data.(nameofdata{i}).(trendvariablestoplot{1}));
        % removing WACCM RCP85 for now
        if strcmp(nameofdata{i},'WACCM_CCMI')        
            fieldstrends.WACCM_CCMI(2) = [];
            data.(nameofdata{i}).ensyears = [];
            data.(nameofdata{i}).ensyears(1).y = data.(nameofdata{i}).years(1).y;
            data.(nameofdata{i}).ensyears(2).y = data.(nameofdata{i}).years(7).y;
            data.(nameofdata{i}).ensyears(3).y = data.(nameofdata{i}).years(10).y;
            data.(nameofdata{i}).ensyears(4).y = data.(nameofdata{i}).years(13).y;
        end        
        ln = length(fieldstrends.(nameofdata{i}));            
    end    
    for l = 1:ln             
        if strcmp(nameofdata{i},'ERA')        
            yearstouse(i,l).y = data.(nameofdata{i}).years.y(1:12:end);
            dateindstart_past = find(yearstouse(i,l).y == timeperiods.(nameofdata{i})(1));
            dateindend_past = find(yearstouse(i,l).y == timeperiods.(nameofdata{i})(2)-1);
            dateindstart_future = find(yearstouse(i,l).y == timeperiods.(nameofdata{i})(2));
            dateindend_future = find(yearstouse(i,l).y == timeperiods.(nameofdata{i})(3));
        else            
            yearstouse(i,l).y = data.(nameofdata{i}).ensyears(l).y(1:12:end);
            dateindstart_past = find(yearstouse(i,l).y == timeperiods.(nameofdata{i})(1));
            dateindend_past = find(yearstouse(i,l).y == timeperiods.(nameofdata{i})(2)-1);
            dateindstart_future = find(yearstouse(i,l).y == timeperiods.(nameofdata{i})(2));
            dateindend_future = find(yearstouse(i,l).y,1,'last');
        end        
        for j = 1:length(trendvariablestoplot)
            count = 1;
            for k = 2:2:length(data.(nameofdata{i}).latitude)                
                if strcmp(nameofdata{i},'ERA')        
                    tempdata = data.(nameofdata{i}).(trendvariablestoplot{j})(k,dateindstart_past:dateindend_past);
                    [bpast(i,l,j,count,:),bpastint(i,l,j,count,:,:)] = regress(tempdata',[ones(1,length(tempdata));1:length(tempdata)]');
                    tempdata2 = data.(nameofdata{i}).(trendvariablestoplot{j})(k,dateindstart_future:dateindend_future);
                    [bfuture(i,l,j,count,:),bfutureint(i,l,j,count,:,:)] = regress(tempdata2',[ones(1,length(tempdata2));1:length(tempdata2)]');
                else
                    tempdata = data.(nameofdata{i}).(trendvariablestoplot{j}).(fieldstrends.(nameofdata{i}){l})(k,dateindstart_past:dateindend_past);
                    [bpast(i,l,j,count,:),bfutureint(i,l,j,count,:,:)] = regress(tempdata',[ones(1,length(tempdata));1:length(tempdata)]');
                    tempdata2 = data.(nameofdata{i}).(trendvariablestoplot{j}).(fieldstrends.(nameofdata{i}){l})(k,dateindstart_future:dateindend_future);
                    [bfuture(i,l,j,count,:),bfutureint(i,l,j,count,:,:)] = regress(tempdata2',[ones(1,length(tempdata2));1:length(tempdata2)]');
                end
                count = count+1;
            end
        end
    end
end

%% put together REFC2 type simulations
% historical
plot_trends = 0;
if plot_trends
orderpast = {'ACCESS-CCM REFC2','ACCESS-CCM REFC1','CanESM2 Historical','ERA-Interim','LENS RCP8.5','WACCM-CCMI REF-C2','WACCM-CCM REF-C1'};
bREFC2 = permute(cat(4,squeeze(bpast(1,2,:,:,:)),squeeze(bpast(1,1,:,:,:)),squeeze(bpast(2,2,:,:,:)),...
squeeze(bpast(3,1,:,:,:)),squeeze(bpast(4,1,:,:,:)),squeeze(bpast(5,1,:,:,:)),squeeze(bpast(5,4,:,:,:))),[4,1,2,3]);

bREFC2future = permute(cat(4,squeeze(bfuture(1,2,:,:,:)),squeeze(bfuture(1,1,:,:,:)),squeeze(bfuture(2,2,:,:,:)),...
squeeze(bfuture(3,1,:,:,:)),squeeze(bfuture(4,1,:,:,:)),squeeze(bfuture(5,1,:,:,:)),squeeze(bfuture(5,4,:,:,:))),[4,1,2,3]);

orderstratozonepast = {'ACCESS-CCM lowGHG','CanESM2 StratOzone','WACCM-CCMI lowGHG'};
bstratozone = permute(cat(4,squeeze(bpast(1,3,:,:,:)),squeeze(bpast(2,4,:,:,:)),squeeze(bpast(5,2,:,:,:))),[4,1,2,3]);
bstratozonefuture = permute(cat(4,squeeze(bfuture(1,3,:,:,:)),squeeze(bfuture(2,4,:,:,:)),squeeze(bfuture(5,2,:,:,:))),[4,1,2,3]);

orderGHGonlypast = {'ACCESS-CCM lowGHG','WACCM-CCMI lowGHG'};
bGHGonly = permute(cat(4,squeeze(bpast(1,4,:,:,:)),squeeze(bpast(5,3,:,:,:))),[4,1,2,3]);
bGHGonlyfuture = permute(cat(4,squeeze(bfuture(1,4,:,:,:)),squeeze(bfuture(5,3,:,:,:))),[4,1,2,3]);

%% plot smoothed amplitudes
latind = 4;
smoothlength = 10;
cbrew = cbrewer('qual','Set1',10);
createfig('large','on');

fieldtoplot = 'amplitude_spr2';
ERAfieldtoplot = 'amplitude_spr';
fi = 1;

if strcmp(fieldtoplot,'amplitude_spr2')
    xlab = 'Temperature';
else
    xlab = 'Longitude ({\circ}E';
end

subplot(3,1,1);

smootheddatatemp = smooth(data.ACCESS.(fieldtoplot).REFC2(latind,1:length(yearstouse(1,2).y)),smoothlength,'moving');
ph(1) = plot(yearstouse(1,2).y(smoothlength/2+1:end-smoothlength/2),(smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2))),'LineWidth',2);
hold on
    
smootheddatatemp = smooth(data.Can2ESM.(fieldtoplot).Historical(latind,1:length(yearstouse(2,2).y)),smoothlength,'moving');
%ph(1) = plot(yearstouse(2,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
ph(2) = plot(yearstouse(2,2).y(smoothlength/2+1:end-smoothlength/2),(smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2))),'LineWidth',2);

smootheddatatemp = smooth(data.WACCM_CCMI.(fieldtoplot).REFC2(latind,1:length(yearstouse(5,1).y)),smoothlength,'moving');
%ph(4) = plot(yearstouse(5,1).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
ph(3) = plot(yearstouse(5,1).y(smoothlength/2+1:end-smoothlength/2),(smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2))),'LineWidth',2);

smootheddatatemp = smooth(data.ERA.(ERAfieldtoplot)(latind,1:length(yearstouse(3,1).y)-1),smoothlength,'moving');
%ph(2) = plot(yearstouse(3,1).y(smoothlength/2+1:end-smoothlength/2-1),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
ph(4) = plot(yearstouse(3,1).y(smoothlength/2+1:end-smoothlength/2-1),(smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2))),'LineWidth',2);

smootheddatatemp = smooth(data.LENS.(fieldtoplot).All(latind,1:length(yearstouse(4,1).y)),smoothlength,'moving');
%ph(3) = plot(yearstouse(4,1).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
ph(5) = plot(yearstouse(4,1).y(smoothlength/2+1:end-smoothlength/2),(smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2))),'LineWidth',2);
fsize = 18;
xlim([1960 2101])
title('All forcings','fontsize',fsize+4);
set(lh,'fontsize',18','location','SouthWest');

ylabel(xlab,'fontsize',fsize+2)
xlabel('Year','fontsize',fsize+2)
lh = legend(ph,'ACCESS-CCM','CanESM2','WACCM-CCMI','ERA-Interim','LENS');
set(lh,'box','off','location','SouthWest','fontsize',fsize-2);
set(gca,'fontsize',fsize);

subplot(3,1,2);

smootheddatatemp = smooth(data.ACCESS.(fieldtoplot).lowGHG(latind,1:length(yearstouse(1,2).y)),smoothlength,'moving');
ph(1) = plot(yearstouse(1,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2)),'LineWidth',2);
hold on

smootheddatatemp = smooth(data.Can2ESM.(fieldtoplot).StratOzone(latind,1:length(yearstouse(2,2).y)),smoothlength,'moving');
ph(2) = plot(yearstouse(2,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2)),'LineWidth',2);

smootheddatatemp = smooth(data.WACCM_CCMI.(fieldtoplot).lowGHG(latind,1:length(yearstouse(5,3).y)),smoothlength,'moving');
ph(3) = plot(yearstouse(5,3).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2)),'LineWidth',2);
title('Stratospheric ozone forcing only','fontsize',fsize+4);
set(gca,'fontsize',fsize);
ylabel(xlab,'fontsize',fsize+2)
xlim([1960 2101])
subplot(3,1,3);

smootheddatatemp = smooth(data.ACCESS.(fieldtoplot).lowODS(latind,1:length(yearstouse(1,2).y)),smoothlength,'moving');
ph(1) = plot(yearstouse(1,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2)),'LineWidth',2);
hold on

smootheddatatemp =smooth(data.Can2ESM.(fieldtoplot).Historical(latind,1:length(yearstouse(2,2).y)) - data.Can2ESM.(fieldtoplot).StratOzone(latind,1:length(yearstouse(2,2).y)),smoothlength,'moving');
ph(2) = plot(yearstouse(2,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2)),'LineWidth',2);

smootheddatatemp = smooth(data.WACCM_CCMI.(fieldtoplot).lowODS(latind,1:length(yearstouse(5,3).y)),smoothlength,'moving');
ph(2) = plot(yearstouse(5,3).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2) - nanmean(smootheddatatemp(smoothlength/2+1:end-smoothlength/2)),'LineWidth',2);
title('GHG forcing only','fontsize',fsize+4);
set(gca,'fontsize',fsize);
ylabel(xlab,'fontsize',fsize+2)
xlabel('Year','fontsize',fsize+2)
xlim([1960 2101])

annotation('textbox',[0 .99 1 0],'String','Springtime zonal wav - 1 amplitude (normalized by mean)','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,...
        'EdgeColor','none','fontweight','bold')

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/Amplitudechanges'];
    export_fig(filename,'-pdf');    
    
%%
    fsize = 18;
    fieldtoplot = 'amplitude_spr2';
    ERAfieldtoplot = 'amplitude_spr';
    fi = 1;

    if strcmp(fieldtoplot,'amplitude_spr2')
        xlab = 'Temperature';
    else
        xlab = 'Longitude ({\circ}E';
    end

    % plot (REF-C2)
    latind = 4;
    latind2 = 3;
    smoothlength = 10;
    cbrew = cbrewer('qual','Set1',10);
    createfig('medium','on');
    smootheddatatemp = smooth(data.ACCESS.(fieldtoplot).REFC2(latind,1:length(yearstouse(1,2).y)),smoothlength,'moving');
    ph(1) = plot(yearstouse(1,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    hold on
    plot([1960:2000],bREFC2(1,fi,latind2,1)+bREFC2(1,fi,latind2,2)*[1:41],'k','LineWidth',3);
    plot([2001:2100],bREFC2future(1,fi,latind2,1)+bREFC2future(1,fi,latind2,2)*[1:100],'k','LineWidth',3);

    smootheddatatemp = smooth(data.ACCESS.(fieldtoplot).REFC1(latind,1:length(yearstouse(1,1).y)),smoothlength,'moving');
    ph(2) = plot(yearstouse(1,1).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    plot([1960:2000],bREFC2(2,fi,latind2,1)+bREFC2(2,fi,latind2,2)*[1:41],'k','LineWidth',3);
    %plot([2001:2010],bREFC2future(2,fi,latind2,1)+bREFC2future(2,fi,latind2,2)*[1:10],'k','LineWidth',3);

    smootheddatatemp = smooth(data.Can2ESM.(fieldtoplot).Historical(latind,1:length(yearstouse(2,2).y)),smoothlength,'moving');
    ph(3) = plot(yearstouse(2,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    plot([1960:2000],bREFC2(3,fi,latind2,1)+bREFC2(3,fi,latind2,2)*[1:41],'k','LineWidth',3);
    plot([2001:2100],bREFC2future(3,fi,latind2,1)+bREFC2future(3,fi,latind2,2)*[1:100],'k','LineWidth',3);

    smootheddatatemp = smooth(data.ERA.(ERAfieldtoplot)(latind,1:length(yearstouse(3,1).y)-1),smoothlength,'moving');
    ph(4) = plot(yearstouse(3,1).y(smoothlength/2+1:end-smoothlength/2-1),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    plot([1979:2000],bREFC2(4,fi,latind2,1)+bREFC2(4,fi,latind2,2)*[1:22],'k','LineWidth',3);


    smootheddatatemp = smooth(data.LENS.(fieldtoplot).All(latind,1:length(yearstouse(4,1).y)),smoothlength,'moving');
    ph(5) = plot(yearstouse(4,1).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    plot([1960:2000],bREFC2(5,fi,latind2,1)+bREFC2(5,fi,latind2,2)*[1:41],'k','LineWidth',3);
    plot([2001:2100],bREFC2future(5,fi,latind2,1)+bREFC2future(5,fi,latind2,2)*[1:100],'k','LineWidth',3);

    smootheddatatemp = smooth(data.WACCM_CCMI.(fieldtoplot).REFC2(latind,1:length(yearstouse(5,1).y)),smoothlength,'moving');
    ph(6) = plot(yearstouse(5,1).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    plot([1960:2000],bREFC2(6,fi,latind2,1)+bREFC2(6,fi,latind2,2)*[1:41],'k','LineWidth',3);
    plot([2001:2100],bREFC2future(6,fi,latind2,1)+bREFC2future(6,fi,latind2,2)*[1:100],'k','LineWidth',3);

    smootheddatatemp = smooth(data.WACCM_CCMI.(fieldtoplot).REFC1(latind,1:length(yearstouse(5,4).y)),smoothlength,'moving');
    ph(7) = plot(yearstouse(5,4).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    plot([1960:2000],bREFC2(7,fi,latind2,1)+bREFC2(7,fi,latind2,2)*[1:41],'k','LineWidth',3);

    title([num2str(data.(nameofdata{i}).latitude(latind)),'{\circ}N, ','All-forcing simulations'],'fontsize',fsize+4);
    ylabel(xlab,'fontsize',fsize+2)
    ylabel('Year','fontsize',fsize+2)
    lh = legend(ph,orderpast);
    set(lh,'box','off');
    set(gca,'fontsize',fsize);
    if contains(fieldtoplot,'amplitude')
        forfilename = 'Amplitude';
    elseif contains(fieldtoplot,'max')
        forfilename = 'Maximum Longitude';
    elseif contains(fieldtoplot,'min')
        forfilename = 'Minimum Longitude';
    end
    title([num2str(data.(nameofdata{i}).latitude(latind)),'{\circ}N, ',forfilename,' All-forcing simulations'],'fontsize',fsize+4);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/Allforcings_spring_',num2str(abs(data.ERA.latitude(latind))),'S_',forfilename];
    export_fig(filename,'-pdf');
    clearvars ph



    % strat only
    createfig('medium','on');
    smootheddatatemp = smooth(data.ACCESS.(fieldtoplot).lowGHG(latind,1:length(yearstouse(1,2).y)),smoothlength,'moving');
    ph(1) = plot(yearstouse(1,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    hold on
    plot([1960:2000],bstratozone(1,fi,latind2,1)+bstratozone(1,fi,latind2,2)*[1:41],'k','LineWidth',3);
    plot([2001:2100],bstratozonefuture(1,fi,latind2,1)+bstratozonefuture(1,fi,latind2,2)*[1:100],'k','LineWidth',3);

    smootheddatatemp = smooth(data.Can2ESM.(fieldtoplot).StratOzone(latind,1:length(yearstouse(2,2).y)),smoothlength,'moving');
    ph(2) = plot(yearstouse(2,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    plot([1960:2000],bstratozone(2,fi,latind2,1)+bstratozone(2,fi,latind2,2)*[1:41],'k','LineWidth',3);
    plot([2001:2100],bstratozonefuture(2,fi,latind2,1)+bstratozonefuture(2,fi,latind2,2)*[1:100],'k','LineWidth',3);

    smootheddatatemp = smooth(data.WACCM_CCMI.(fieldtoplot).lowGHG(latind,1:length(yearstouse(5,3).y)),smoothlength,'moving');
    ph(3) = plot(yearstouse(5,3).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    plot([1960:2000],bstratozone(3,fi,latind2,1)+bstratozone(3,fi,latind2,2)*[1:41],'k','LineWidth',3);
    plot([2001:2100],bstratozonefuture(3,fi,latind2,1)+bstratozonefuture(3,fi,latind2,2)*[1:100],'k','LineWidth',3);


    ylabel(xlab,'fontsize',fsize+2)
    ylabel('Year','fontsize',fsize+2)
    lh = legend(ph,orderstratozonepast);
    set(lh,'box','off');
    set(gca,'fontsize',fsize);
    if contains(fieldtoplot,'amplitude')
        forfilename = 'Amplitude';
    elseif contains(fieldtoplot,'max')
        forfilename = 'Maximum Longitude';
    elseif contains(fieldtoplot,'min')
        forfilename = 'Minimum Longitude';
    end
    title([num2str(data.(nameofdata{i}).latitude(latind)),'{\circ}N, ',forfilename,' Stratospheric ozone simulations'],'fontsize',fsize+4);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/StratOzone_spring_',num2str(abs(data.ERA.latitude(latind))),'S_',forfilename];
    export_fig(filename,'-pdf');
    clearvars ph



    % GHG only
    createfig('medium','on');
    smootheddatatemp = smooth(data.ACCESS.(fieldtoplot).lowODS(latind,1:length(yearstouse(1,2).y)),smoothlength,'moving');
    ph(1) = plot(yearstouse(1,2).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    hold on
    plot([1960:2000],bGHGonly(1,fi,latind2,1)+bGHGonly(1,fi,latind2,2)*[1:41],'k','LineWidth',3);
    plot([2001:2100],bGHGonlyfuture(1,fi,latind2,1)+bGHGonlyfuture(1,fi,latind2,2)*[1:100],'k','LineWidth',3);

    smootheddatatemp = smooth(data.WACCM_CCMI.(fieldtoplot).lowODS(latind,1:length(yearstouse(5,3).y)),smoothlength,'moving');
    ph(2) = plot(yearstouse(5,3).y(smoothlength/2+1:end-smoothlength/2),smootheddatatemp(smoothlength/2+1:end-smoothlength/2),'LineWidth',2);
    plot([1960:2000],bGHGonly(2,fi,latind2,1)+bGHGonly(2,fi,latind2,2)*[1:41],'k','LineWidth',3);
    plot([2001:2100],bGHGonlyfuture(2,fi,latind2,1)+bGHGonlyfuture(2,fi,latind2,2)*[1:100],'k','LineWidth',3);


    ylabel(xlab,'fontsize',fsize+2)
    ylabel('Year','fontsize',fsize+2)
    lh = legend(ph,orderGHGonlypast);
    set(lh,'box','off');
    set(gca,'fontsize',fsize);
    if contains(fieldtoplot,'amplitude')
        forfilename = 'Amplitude';
    elseif contains(fieldtoplot,'max')
        forfilename = 'Maximum Longitude';
    elseif contains(fieldtoplot,'min')
        forfilename = 'Minimum Longitude';
    end
    title([num2str(data.(nameofdata{i}).latitude(latind)),'{\circ}N, ',forfilename,' GHG only simulations'],'fontsize',fsize+4);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/GHGonly_spring_',num2str(abs(data.ERA.latitude(latind))),'S_',forfilename];
    export_fig(filename,'-pdf');
    clearvars ph
end
    %% plotting
    if plottingens
        variabletoplot = 'amplitude_spr2'; %'minlongitude_spr'; %amplitude_spr
        if contains(variabletoplot,'amplitude')
            xtit = 'Temperature (K)';
        else
            xtit = 'Longitude ({\circ}E)'; 
        end
        fsize = 18;
        cbrew = cbrewer('qual','Set1',10);
        smoothlength = 20;
        for i = length(nameofdata)-1
            if strcmp(nameofdata{i},'ERA')
                variabletoplot2 = variabletoplot(1:end-1);
            else
                fields.(nameofdata{i}) = fieldnames(data.(nameofdata{i}).(variabletoplot));
            end

            for j = 2:2:length(data.(nameofdata{i}).latitude)
                createfig('medium','on');        
                if strcmp(nameofdata{i},'ERA')
                    plottingdate = data.(nameofdata{i}).years.y(1):data.(nameofdata{i}).years.y(end-2);
                    smootheddatatemp = smooth(data.(nameofdata{i}).(variabletoplot2)(j,1:length(plottingdate)),smoothlength,'moving');
                    ph(k) = plot(plottingdate((smoothlength/2+1):length(smootheddatatemp(1:end-smoothlength/2))),...
                        smootheddatatemp(smoothlength/2+1:end-smoothlength/2));
                    set(ph(k),'LineWidth',2,'color',cbrew(k,:));   
                    title([num2str(data.(nameofdata{i}).latitude(j)),'N, ',nameofdatafortitles{i}]);
                else
                    for k = 1:length(fields.(nameofdata{i}))
                        plottingdate = data.(nameofdata{i}).ensyears(k).y(1):data.(nameofdata{i}).ensyears(k).y(end);
                        smootheddatatemp = smooth(data.(nameofdata{i}).(variabletoplot).(fields.(nameofdata{i}){k})(j,1:length(plottingdate)),smoothlength,'moving');
                        ph(k) = plot(plottingdate((smoothlength/2+1):length(smootheddatatemp(1:end-smoothlength/2))),...
                            smootheddatatemp(smoothlength/2+1:end-smoothlength/2));
                        set(ph(k),'LineWidth',2,'color',cbrew(k,:));
                        title([num2str(data.(nameofdata{i}).latitude(j)),'N, ',nameofdatafortitles{i}]);

                        hold on
                    end
                    lh = legend([ph(:)],fields.(nameofdata{i}));
                    set(lh,'fontsize',fsize,'box','off');            
                end                
                set(gca,'fontsize',fsize);
                ylabel(xtit,'fontsize',fsize+2);
                xlabel('Year','fontsize',fsize+2);

                clearvars ph
                filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/modelSummary/',...
                    nameofdata{i},'_',num2str(data.(nameofdata{i}).latitude(j)),'N_',variabletoplot,'.pdf'];
                export_fig(filename,'-pdf');
            end
        end
    end

%%
if plottingens2           
    variabletoplot = 'minlongitude_m_spr'; %'minlongitude_spr'; %amplitude_spr
    if contains(variabletoplot,'amplitude')
        xtit = 'Temperature (K)';
    else
        xtit = 'Longitude ({\circ}E)'; 
    end
    fsize = 18;
    cbrew = cbrewer('qual','Set1',10);
    smoothlength = 20;
    for i = length(nameofdata)-1
        % creating ensembles from individuals
        if strcmp(nameofdata{i},'WACCM_CCMI')
            
        elseif strcmp(nameofdata{i},'ACCESS_CCMI')
            
        elseif strcmp(nameofdata{i},'Can2ESM')
            %data.Can2ESM.(variabletoplot) = orderfields(data.Can2ESM.(variabletoplot));
            
        end
                        
        if strcmp(nameofdata{i},'ERA')

        else
            fields.(nameofdata{i}) = fieldnames(data.(nameofdata{i}).(variabletoplot));            
        end

        for j = 1:2:length(data.(nameofdata{i}).latitude)
            createfig('medium','on');        
            if strcmp(nameofdata{i},'ERA')
                plottingdate = data.(nameofdata{i}).years.y(1):data.(nameofdata{i}).years.y(end-2);
                smootheddatatemp = smooth(data.(nameofdata{i}).(variabletoplot)(j,1:length(plottingdate)),smoothlength,'moving');
                ph(k) = plot(plottingdate((smoothlength/2+1):length(smootheddatatemp(1:end-smoothlength/2))),...
                    smootheddatatemp(smoothlength/2+1:end-smoothlength/2));
                set(ph(k),'LineWidth',2,'color',cbrew(k,:));   
                title([num2str(data.(nameofdata{i}).latitude(j)),'N, ',nameofdatafortitles{i}]);
            else
                for k = 1:length(fields.(nameofdata{i}))
                    plottingdate = data.(nameofdata{i}).ensyears(k).y(1):data.(nameofdata{i}).ensyears(k).y(end);
                    smootheddatatemp = smooth(data.(nameofdata{i}).(variabletoplot).(fields.(nameofdata{i}){k})(j,1:length(plottingdate)),smoothlength,'moving');
                    ph(k) = plot(plottingdate((smoothlength/2+1):length(smootheddatatemp(1:end-smoothlength/2))),...
                        smootheddatatemp(smoothlength/2+1:end-smoothlength/2));
                    set(ph(k),'LineWidth',2,'color',cbrew(k,:));
                    title([num2str(data.(nameofdata{i}).latitude(j)),'N, ',nameofdatafortitles{i}]);

                    hold on
                end
                lh = legend([ph(:)],fields.(nameofdata{i}));
                set(lh,'fontsize',fsize,'box','off');            
            end                
            set(gca,'fontsize',fsize);
            ylabel(xtit,'fontsize',fsize+2);
            xlabel('Year','fontsize',fsize+2);

            clearvars ph
            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/modelSummary/',...
                nameofdata{i},'_',num2str(data.(nameofdata{i}).latitude(j)),'N_',variabletoplot,'.pdf'];
            export_fig(filename,'-pdf');
        end
    end
end
    close all
    plin = 0;
    %% plotting individual amplitudes
if plin
    cannames = fieldnames(data.WACCM_CCMI.amplitude_m_spr);
    figure;
    initial = 10;
    toplot = 'minlongitude_m_spr';
    for i = initial:initial+2%length(cannames)
        plot(squeeze(data.WACCM_CCMI.(toplot)(i).m(5,:)),'b')
        hold on
        
    end
    testmean = nanmean([squeeze(data.WACCM_CCMI.(toplot)(initial).m(5,:));...
            squeeze(data.WACCM_CCMI.(toplot)(initial+1).m(5,:));...
            squeeze(data.WACCM_CCMI.(toplot)(initial+2).m(5,:))]);
    testmeansmooth = smooth(testmean,10,'moving');
        
    plot(testmean,'r','LineWidth',2);
    figure;
    plot(testmeansmooth,'k','LineWidth',2);
end    

    %plot(data.WACCM_CCMI.amplitude_spr2.lowODS(5,:),'k');

%% correlations between amplitude and 4060profiles
readvertdata = 0;
LENScorr = 0;
if LENScorr
    if readvertdata
        LENSverticaldir = '/Volumes/MyBook/work/data/LENS/3060S/';
        LENSverfiles = dir([LENSverticaldir,'*.nc']);

        for i = 1:length(LENSverfiles)
            [~,LENSverdata(i),~] = Read_in_netcdf([LENSverticaldir,LENSverfiles(i).name]);
            for j = 1:12
                LENSverens.months(i,:,:,j,:) = LENSverdata(i).T(:,:,j:12:end);  
            end
            LENSverens.spring(i,:,:,:) = squeeze(nanmean(LENSverens.months(i,:,:,9:11,:),4));  
        end
        LENScombine = cat(4,LENSverdata(:).T);
        LENSverens.mean = nanmean(LENScombine,4);
        for i = 1:12
            LENSverens.mean_months(:,:,i,:) = LENSverens.mean(:,:,i:12:end);  
        end
        LENSverens.mean_spring = squeeze(nanmean(LENSverens.mean_months(:,:,[9,10,11],:),3));
        LENSverens.lev = LENSverdata(1).lev;
        save('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/lonvert/LENSvert3060S.mat','LENSverens','-v7.3');
    else
        load('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/lonvert/LENSvert3060S.mat');
    end

    %% talking correlations
    clearvars rLENS rLENS2
    dateind = 76:106;
    dateind2 = 36:60;
    dateind3 = 131:161;
    oct = 0;
    for i = 1:size(LENSverens.mean_spring,1)
        tic;
        for j = 1:size(LENSverens.mean_spring,2)
            if oct
                temp = cat(4,data.LENS.amplitude_m(:).m);
                temp = detrend(squeeze(temp(10,4,dateind,:)));
                temp = temp(:);
                %LENSverens.october = squeeze(LENSverens.months(:,:,:,10,:));
                temp2 = detrend(squeeze(LENSverens.months(:,i,j,10,dateind))');
                temp2 = temp2(:);
                octext = 'oct';
            else
                temp = cat(3,data.LENS.amplitude_m_spr(:).m);
                templowcl = detrend(squeeze(temp(4,dateind2,:)));
                tempozre = detrend(squeeze(temp(4,dateind3,:)));
                temp = detrend(squeeze(temp(4,dateind,:)));                
                temp = temp(:);
                templowcl = templowcl(:);
                tempozre = tempozre(:);
                temp2 = detrend(squeeze(LENSverens.spring(:,i,j,dateind))');
                temp2 = temp2(:);
                temp2lowcl = detrend(squeeze(LENSverens.spring(:,i,j,dateind2))');
                temp2lowcl = temp2lowcl(:);
                temp2ozre = detrend(squeeze(LENSverens.spring(:,i,j,dateind3))');
                temp2ozre = temp2ozre(:);
                
                octext = 'spr';                
            end                                    
            %rLENS(i,j) = corr(detrend(data.LENS.amplitude_spr2.All(4,41:81))',detrend(squeeze(LENSverens.mean_spring(i,j,41:81))));
            [rLENS(1,i,j),pLENS(1,i,j)] = corr(temp,temp2);
            [rLENS(2,i,j),pLENS(2,i,j)] = corr(templowcl,temp2lowcl);
            [rLENS(3,i,j),pLENS(3,i,j)] = corr(tempozre,temp2ozre);
        end
        toc;
    end

    for i = 1:size(LENSverens.spring,2)
        for j = 1:size(LENSverens.spring,3)
            count = 1;
            for k = [1:1]
                rLENS2(count,i,j) = corr(detrend(data.LENS.amplitude_m_spr(k).m(4,dateind))',detrend(squeeze(LENSverens.spring(k,i,j,dateind))));
                count = count+1;
            end
        end
    end

    pLENS(:,:,1:14) = 1;
    %rLENS = reshape(rLENS,[1,size(rLENS)]);
    cbrew = cbrewer('div','RdBu',16);         

        prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];

        presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
            5,[],[],2,1};


        logprestick = log(prestick);

        contourtitle = {'Correlation between 65{\circ}S Amplitude and 30-60{\circ}S temperature'};       
        contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','2050-2079 (high GHGs - low chlorine)'};       

        %%
        subplotmaps(rLENS,data.LENS.longitude,log(LENSverens.lev),{'div','RdBu'},1,pLENS,18,contourtitle2,'Longitude','','Correlation','on',...
            [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');
        
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/LENSamp4060Scorr_',octext];
        export_fig(filename,'-pdf');
%         for i = 1:size(rLENS2,1)
%             subplotmaps(rLENS2(i,:,:),data.LENS.longitude,log(LENSverens.lev),{'div','RdBu'},1,[],12,contourtitle,'Longitude ({\circ})E','Pressure (hPa)','Correlation','on',...
%                 [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),'',1,[0 360],[log(1) log(1000)],1,'none',0,'none');
%             filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/LENSamp4060Scorr_ens_',sprintf('%02d',i)];
%             export_fig(filename,'-png');
%         end        
        
end

%%
readvertdata = 1;
CanESM2corr = 1;
if CanESM2corr
    if readvertdata
        CanESM2verticaldir = '/Volumes/MyBook/work/data/CanESM2/3060S/';
        CanESM2verfiles = dir([CanESM2verticaldir,'*.nc']);

        for i = 1:length(CanESM2verfiles)
            [~,CanESM2verdata(i),~] = Read_in_netcdf([CanESM2verticaldir,CanESM2verfiles(i).name]);
            for j = 1:12
                CanESM2verens.months(i,:,:,j,:) = CanESM2verdata(i).ta(:,:,j:12:end);  
            end
            CanESM2verens.spring(i,:,:,:) = squeeze(nanmean(CanESM2verens.months(i,:,:,9:11,:),4));  
        end
        CanESM2combine.hist = cat(4,CanESM2verdata(1:50).ta);
        CanESM2combine.strat = cat(4,CanESM2verdata(51:100).ta);
        CanESM2verens.meanhist = nanmean(CanESM2combine.hist,4);
        CanESM2verens.meanstrat = nanmean(CanESM2combine.strat,4);
        for i = 1:12
            CanESM2verens.meanhist_months(:,:,i,:) = CanESM2verens.meanhist(:,:,i:12:end);  
            CanESM2verens.meanstrat_months(:,:,i,:) = CanESM2verens.meanstrat(:,:,i:12:end);  
        end
        CanESM2verens.meanhist_spring = squeeze(nanmean(CanESM2verens.meanhist_months(:,:,[9,10,11],:),3));
        CanESM2verens.meanstrat_spring = squeeze(nanmean(CanESM2verens.meanstrat_months(:,:,[9,10,11],:),3));
        CanESM2verens.lev = CanESM2verdata(1).plev;
        save('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/lonvert/CanESM2vert3060S.mat','CanESM2verens','-v7.3');
    else
        load('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/lonvert/CanESM2vert3060S.mat');
    end

    %% talking correlations
    clearvars rCanESM2 rCanESM22
    dateind = 46:76;
    dateind2 = 6:30;
    dateind3 = 101:130;
    for i = 1:size(CanESM2verens.meanhist_spring,1)
        tic;
        for j = 1:size(CanESM2verens.meanhist_spring,2)
            
            if oct
                temp = cat(4,data.Can2ESM.amplitude_m(51:100).m);
                temp = detrend(squeeze(temp(10,4,dateind,:)));                
                temp = temp(:);
                %CanESM2verens.october = squeeze(CanESM2verens.months(:,:,:,10,:));
                temp2 = detrend(squeeze(CanESM2verens.months(1:50,i,j,10,dateind))');
                temp2 = temp2(:);
                octext = 'oct';
            else
                temp = cat(3,data.Can2ESM.amplitude_m_spr(51:100).m);                
                tempstrat = cat(3,data.Can2ESM.amplitude_m_spr(102:2:end).m);
                
                templowcl = detrend(squeeze(temp(4,dateind2,:)));
                temphighGHG = detrend(squeeze(temp(4,dateind3,:)));
                temp = detrend(squeeze(temp(4,dateind,:)));
                
                tempstratlowcl = detrend(squeeze(tempstrat(4,dateind2,:)));
                tempstrathighGHG = detrend(squeeze(tempstrat(4,dateind3,:)));
                tempstrat = detrend(squeeze(tempstrat(4,dateind,:)));
                
                temp = temp(:);
                templowcl = templowcl(:);
                temphighGHG = temphighGHG(:);
                tempstrat = tempstrat(:);
                tempstratlowcl = tempstratlowcl(:);
                tempstrathighGHG = tempstrathighGHG(:);
                
                temp2 = detrend(squeeze(CanESM2verens.spring(1:50,i,j,dateind))');
                temp2lowcl = detrend(squeeze(CanESM2verens.spring(1:50,i,j,dateind2))');
                temp2highGHG = detrend(squeeze(CanESM2verens.spring(1:50,i,j,dateind3))');
                temp2strat = detrend(squeeze(CanESM2verens.spring(51:end,i,j,dateind))');
                temp2stratlowcl = detrend(squeeze(CanESM2verens.spring(51:end,i,j,dateind2))');
                temp2strathighGHG = detrend(squeeze(CanESM2verens.spring(51:end,i,j,dateind3))');
                temp2 = temp2(:);
                temp2lowcl = temp2lowcl(:);
                temp2highGHG = temp2highGHG(:);
                temp2strat = temp2strat(:); 
                temp2stratlowcl = temp2stratlowcl(:); 
                temp2strathighGHG = temp2strathighGHG(:); 
                octext = 'spr';
                
            end   
                        
            %rCanESM2(i,j) = corr(detrend(data.Can2ESM.amplitude_spr2.Historical(4,11:51))',detrend(squeeze(CanESM2verens.meanhist_spring(i,j,11:51))));
            [rCanESM2(1,i,j),pCanESM2(1,i,j)] = corr(temp,temp2);
            [rCanESM2(2,i,j),pCanESM2(2,i,j)] = corr(templowcl,temp2lowcl);
            [rCanESM2(3,i,j),pCanESM2(3,i,j)] = corr(temphighGHG,temp2highGHG);
            
            [rCanESMstrat(1,i,j),pCanESMstrat(1,i,j)] = corr(tempstrat,temp2strat);
            [rCanESMstrat(2,i,j),pCanESMstrat(2,i,j)] = corr(tempstratlowcl,temp2stratlowcl);
            [rCanESMstrat(3,i,j),pCanESMstrat(3,i,j)] = corr(tempstrathighGHG,temp2strathighGHG);
        end
        toc;
    end

    for i = 1:size(CanESM2verens.spring,2)
        for j = 1:size(CanESM2verens.spring,3)
            count = 1;
            for k = [1:1]
                rCanESM22(count,i,j) = corr(detrend(data.Can2ESM.amplitude_m_spr(50+k).m(4,dateind))',detrend(squeeze(CanESM2verens.spring(k,i,j,dateind))));
                count = count+1;
            end
        end
    end


    %rCanESM2 = reshape(rCanESM2,[1,size(rCanESM2)]);
    %%
    pCanESM2(:,:,10:end) = 1;
    pCanESMstrat(:,:,10:end) = 1;
    cbrew = cbrewer('div','RdBu',16);         

    prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];

    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
        5,[],[],2,1};


    logprestick = log(prestick);

    contourtitle = {'Correlation between 65{\circ}S 50 hPa amplitude and 30-60{\circ}S temperature'};       
    contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','2050-2079 (high GHG - low chlorine)'};       

    subplotmaps(rCanESM2,data.Can2ESM.longitude,log(CanESM2verens.lev/100),{'div','RdBu'},1,pCanESM2,18,contourtitle2,'Longitude','','Correlation','on',...
        [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/CanESM2amp3060Scorr_',octext];
    export_fig(filename,'-pdf');
    
    contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','2050-2079 (low chlorine)'};       
    
    subplotmaps(rCanESMstrat,data.Can2ESM.longitude,log(CanESM2verens.lev/100),{'div','RdBu'},1,pCanESMstrat,18,contourtitle2,'Longitude','','Correlation','on',...
        [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/Strat_CanESM2amp3060Scorr_',octext];
    export_fig(filename,'-pdf');
    
    
%     for i = 1:size(rLENS22,1)
%         subplotmaps(rCanESM22(i,:,:),data.Can2ESM.longitude,log(CanESM2verens.lev/100),{'div','RdBu'},1,[],12,'','Longitude ({\circ})E','Pressure (hPa)','Correlation','on',...
%             [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');
% 
%         filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/CanSM2amp4060Scorr_ens_',sprintf('%02d',i)];
%         export_fig(filename,'-png');
%     end
end
