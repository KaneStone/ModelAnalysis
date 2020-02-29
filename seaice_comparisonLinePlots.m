function [] = seaice_comparisonLinePlots(seaicemodel,observations,tozdata,tozobs,Seasplot,...
    modlat,modlon,obslat,obslon,inputs,pct)

seaiceobs = observations.montharrange(:,observations.timeperiod >= observations.toz.timeperiod(1) & observations.timeperiod <= observations.toz.timeperiod(2),:,:);

% plots sea ice extent for ensemble and observations
colors = cbrewer('qual','Set1',10);
colors2 = cbrewer('qual','Set2',7);
%% Read in area grid cells
areacell = gridcellarea(modlat,modlon);
areacellobs = gridcellarea(obslat,obslon);

%% calculating sea ice extent for different regions

for i = 1:size(Seasplot.lon,1)
    latind = modlat >= Seasplot.lat(i,1) & modlat <= Seasplot.lat(i,2);
    latindobs = obslat >= Seasplot.lat(i,1) & obslat <= Seasplot.lat(i,2);
    if i == 3
        lonind = modlon >= Seasplot.lon(i,1) | modlon < Seasplot.lon(i,2)-360;
        lonindobs = obslon >= Seasplot.lon(i,1) | obslon < Seasplot.lon(i,2)-360;
    else
        lonind = modlon >= Seasplot.lon(i,1) & modlon < Seasplot.lon(i,2);
        lonindobs = obslon >= Seasplot.lon(i,1) & obslon < Seasplot.lon(i,2);
    end
    
    % model
    for j = 1:9
        temp = permute(areacell(latind,lonind).*permute(squeeze(seaicemodel(j,:,:,latind,lonind)),[3,4,1,2]),[3,4,1,2]);
        seaiceextent.model.ind(:,:,j,i) = (detrend(sum(temp(:,:,:),3)')+nanmean(sum(temp(:,:,:),3)',1))';
        
    end
    
    %obervations
    temp1 = permute(areacellobs(latindobs,lonindobs).*permute(squeeze(seaiceobs(:,:,lonindobs,latindobs)),[4,3,1,2]),[3,4,1,2]);
    seaiceextent.obs(:,:,i) = sum(temp1(:,:,:),3);
    
    
end
seaiceextent.model.indanomaly = seaiceextent.model.ind - nanmedian(seaiceextent.model.ind,2);
seaiceextent.obsanomaly = seaiceextent.obs - nanmedian(seaiceextent.obs,2);

%%
tozobs_anomaly = tozobs - nanmean(tozobs);
tozdata_anomaly = tozdata.highCl.dataMonthArrange - nanmedian(tozdata.highCl.dataMonthArrange,3);





%% testing gaussian
clearvars testresult
testing = 0;
if testing
    
    test = squeeze(seaiceextent.model.indanomaly(10,:,:,1));
    test2 = test(:);
    histogram(test2)
    
    for j = 1:1000
        for i = 1:9
            x1 = randsample(30,6);
            x2 = randsample(30,6);

            testextracthigh(i) = nanmean(test(x1,i));
            testextractlow(i) = nanmean(test(x2,i));
            testdifference(i) = testextracthigh(i) - testextractlow(i);
        end
        testresult(j) = nanmean(testdifference);
    end
end


%% plotting

mon = 10;
sea = 1;
for i = sea
    createfig('medium','on');
    for j = 1:10
        subplot(5,2,j)
        if j < 10
            plot(seaiceextent.model.indanomaly(mon,:,j,i));
            ylim([-4e5 4e5]);
            yyaxis right
            plot(squeeze(tozdata_anomaly(j,3,:)));
            ylim([-75 75]);
            hold on
            plot([0 31],[0 0],'--k');
            r(i,j) = corr(squeeze(seaiceextent.model.indanomaly(mon,:,j,i))',squeeze(squeeze(tozdata_anomaly(j,3,:))));
        else
            plot(seaiceextent.obsanomaly(mon,:,i));
            ylim([-4e5 4e5]);
            yyaxis right
            plot(squeeze(tozobs_anomaly(1,:)));
            ylim([-75 75]);
            hold on
            plot([0 38],[0 0],'--k');
            r(i,j) = corr(squeeze(seaiceextent.obsanomaly(mon,:,i))',squeeze(tozobs_anomaly(1,:))');
        end
    end
end
    
% extracting only extremes and plotting
for i = 1:9
    ice_extremes(:,:,i,:) = seaiceextent.model.indanomaly(:,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],i,:);
    ice_extremeshigh(:,:,i,:) = nanmean(seaiceextent.model.indanomaly(:,pct.highCl.ind.highind(i,:),i,:),2);
    ice_extremeslow(:,:,i,:) = nanmean(seaiceextent.model.indanomaly(:,pct.highCl.ind.lowind(i,:),i,:),2);
    toz_extremes(i,:,:) = tozdata_anomaly(i,:,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)]);
end
iceall = seaiceextent.model.indanomaly;
save('/Volumes/ExternalOne/work/data/predruns/output/ice/predtesting.mat','ice_extremes','toz_extremes','iceall');

icehighmean = squeeze(nanmean(ice_extremeshigh,3));
icelowmean = squeeze(nanmean(ice_extremeslow,3));
icediff = icehighmean - icelowmean;

obs_extremeslow = seaiceextent.obsanomaly(:,[observations.toz.highind,observations.toz.lowind],:);
obs_tozextremes = tozobs_anomaly(:,[observations.toz.highind,observations.toz.lowind]);


for i = sea
    createfig('medium','on');
    for j = 1:10
        subplot(5,2,j)
        if j < 10
            plot(ice_extremes(mon,:,j,i),'-o');
            ylim([-4e5 4e5]);
            yyaxis right
            plot(squeeze(toz_extremes(j,3,:)),'-o');
            ylim([-75 75]);
            hold on
            plot([0 13],[0 0],'--k');
            rext(i,j) = corr(squeeze(ice_extremes(mon,:,j,i))',squeeze(squeeze(toz_extremes(j,3,:))));
        else
            plot(obs_extremeslow(mon,:,i),'-o');
            ylim([-4e5 4e5]);
            yyaxis right
            plot(squeeze(obs_tozextremes(1,:)),'-o');
            ylim([-75 75]);
            hold on
            plot([0 16],[0 0],'--k');
            rext(i,j) = corr(squeeze(obs_extremeslow(mon,:,i))',squeeze(obs_tozextremes(1,:))');
        end
    end
end

end
