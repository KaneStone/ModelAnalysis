function [] = predruns_leaveoneout(surfacedata,surfacedataind,tozdata,inputs,latitude,longitude,pct,differences,correlations)

% over 1995-2014. Leave one year out of the regression and see if I can
% predict the remaining years anomaly

% Do this for the ensemble mean and produce a map metric of predictability
% for the ensemble average.

% Using the high signal to noise plots, extract regions near these places
% for each 

% keep leaving years out

% Take regression of ensemble

% DON'T predefine anything!!!

%% Read in ERA-Interim
% Observations = load(['/Volumes/ExternalOne/work/data/predruns/output/data/obs/','obs_perc_diff',monthnames(inputs.varmonthtomean,1,1),'_and_',monthnames(inputs.varmonth,1,1),'.mat']);
% [obs.GSS] = predruns_obsleaveoneout(Observations,longitude,latitude,inputs);
%% rearrange data into composite

surfacedata_composite = permute(surfacedata,[3,4,2,1]);
surfacedata_composite = surfacedata_composite(:,:,:); 

for i = 1:9
    precondtoz(:,i) = squeeze(tozdata(i,inputs.tozmonth,:));%-nanmean(squeeze(tozdata(i,inputs.tozmonth,:)));
    precondtoz_pct(:,i) = squeeze(tozdata(i,inputs.tozmonth,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)])) -...
        nanmean(squeeze(tozdata(i,inputs.tozmonth,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)])));            
end

tozdata_composite = precondtoz(:);
tozdata_composite_pct = precondtoz_pct(:);

surfacedata = double(surfacedata);
surfacedataind = double(surfacedataind);

%% each individual separately 

%% testing
%close all
testing = 0;
if testing
    %sd1 = squeeze(surfacedata(1,:,79,31));
    %sd1 = double(squeeze(surfacedata(1,:,79,35)));
    %sd1 = squeeze(surfacedata(1,:,77,43));
    sd1 = squeeze(surfacedata(1,:,80,79));
    %sd1 = squeeze(surfacedata(1,:,36,79));
    %sd1 = squeeze(surfacedata(1,:,84,44));
    b = regress(sd1',[precondtoz(:,1)]);
    modelpred = b(1)*precondtoz(:,1);%+b(2);

    % leaveout
    for i = 1:length(sd1)
        sd12 = 1:length(sd1);
        sd12 (sd12 == i) = [];
        ozoneind = precondtoz(sd12,1);%-nanmean(precondtoz(sd12,1));        

        %sumpos = sum(
        ozoneind_left = precondtoz(i,1);%-nanmean(precondtoz(sd12,1));
        %sd123(i,:) = sd1(sd12)-nanmean(sd1(sd12),'double');
        b = regress(sd1(sd12)',[ones(size(sd1(sd12)))',[1:length(sd1(sd12))]']);
        sd123(i,:) = sd1(sd12) - b(1) - b(2)*[1:length(sd1(sd12))];
        bloo(i,:) = regress(sd123(i,:)',[ones(size(ozoneind)),ozoneind]);
        modelpredloo(i) = bloo(i,2)*ozoneind_left+bloo(i,1);
        r(i) = corr(ozoneind,sd123(i,:)');
        %modelpredloo(i) = (bloo(i,1)*ozoneind_left);
        acdata(i) = sd1(i)-nanmean(sd1(sd12));
        acdata3(i) = sd1(i) - b(1) - b(2)*i;
    end
    acdata2 = detrend(acdata);
    modelpredloo2 = detrend(modelpredloo);
    figure
    hold on
    plot(acdata)
    plot(acdata2)
    plot(acdata3)
    plot(modelpredloo2,'k--')
    plot(modelpredloo,'k')
    figure    
    plot(r)


% plot(squeeze(modelprediction_ind(1,:,80,30)),'k')
% hold on
% plot(squeeze(actualdata_ind(1,:,80,30))-nanmeGSScorrecttestan(squeeze(actualdata_ind(1,:,80,30))),'k--')

%%


    modelpredictionign_ind = sign(modelpredloo);
    modelpredictionign_ind2 = sign(modelpred');
    actualdatasign_ind =  sign(acdata3);

    modelpredictionresults_ind = modelpredictionign_ind;
    modelpredictionresults_ind (modelpredictionresults_ind ~= actualdatasign_ind) = -2;
    modelpredictionresults_ind (modelpredictionresults_ind == actualdatasign_ind) = 1;
    modelpredictionresults_ind (modelpredictionresults_ind == -2) = 0;

    modelpredictionresults_ind2 = modelpredictionign_ind2;
    modelpredictionresults_ind2 (modelpredictionresults_ind2 ~= actualdatasign_ind) = -2;
    modelpredictionresults_ind2 (modelpredictionresults_ind2 == actualdatasign_ind) = 1;
    modelpredictionresults_ind2 (modelpredictionresults_ind2 == -2) = 0;
    
    abs(sum(sign(bloo(:,2)))) ~= length(modelpredloo);    
    

    GSStest = (sum(modelpredictionresults_ind,2) - size(modelpredictionresults_ind,2)./2)./...
        (size(modelpredictionresults_ind,2)-size(modelpredictionresults_ind,2)./2)*100;  
%     if abs(sum(sign(bloo(:,2)))) ~= length(modelpredloo)   
%         GSStest = 0;
%     end
   

    GSStest2 = (sum(modelpredictionresults_ind2,2) - size(modelpredictionresults_ind2,2)./2)./...
        (size(modelpredictionresults_ind2,2)-size(modelpredictionresults_ind2,2)./2)*100;  

end
%%

for k = 1:size(surfacedata,1)
    count = 1;    
    allin = 1:size(precondtoz,1);
    
    %Ozone informaton
    ozone_anom(k,:) = precondtoz(:,k) - nanmean(precondtoz(:,k));    
    surface_anom(k,:,:,:) = surfacedata(k,:,:,:) - nanmean(surfacedata(k,:,:,:),2);
    
    for j = 1:length(latitude)   
        detrendedsurface = detrend(squeeze(surface_anom(k,:,j,:)));
        ball = [ozone_anom(k,:)',ones(size(ozone_anom(k,:)'))]\detrend(squeeze(surface_anom(k,:,j,:)));
        ballpct = [ozone_anom(k,[pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)])',...
            ones(size(ozone_anom(k,[pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)])'))]\detrendedsurface([pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)],:);
        modelpred_all(k,:,j,:) = ozone_anom(k,:).*ball(1,:)'+ball(2,:)';
        modelpred_allpct(k,:,j,:) = ozone_anom(k,[pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)]).*ballpct(1,:)'+ballpct(2,:)';
        actualdata_all(k,:,j,:) = detrend(squeeze(surface_anom(k,:,j,:)));
        actualdata_allpct(k,:,j,:) = detrendedsurface([pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)],:);
        for l = 1:size(surfacedataind,3)
            detrendedsurfacemonth = detrend(squeeze(surfacedataind(k,:,l,j,:)));
            ballmonth = [ozone_anom(k,:)',ones(size(ozone_anom(k,:)'))]\detrendedsurfacemonth;
            ballpctmonth = [ozone_anom(k,[pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)])',...
            ones(size(ozone_anom(k,[pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)])'))]\detrendedsurfacemonth([pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)],:);
            modelpred_allmonth(k,:,l,j,:) = ozone_anom(k,:).*ballmonth(1,:)'+ballmonth(2,:)';
            modelpred_allpctmonth(k,:,l,j,:) = ozone_anom(k,[pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)]).*ballpctmonth(1,:)'+ballpctmonth(2,:)';
            actualdata_allmonth(k,:,l,j,:) = detrendedsurfacemonth;
            actualdata_allpctmonth(k,:,l,j,:) = detrendedsurfacemonth([pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)],:);
        end
    end
    
    for i = 1:length(allin)
        leaveout = allin;        
        leaveout (leaveout == count) = [];         
        
        %predictors_anom = [leaveout',ones(size(leaveout))'];
        predictors_anom = [leaveout',ones(size(leaveout))'];
        %ozone_ind_b = precondtoz(leaveout,k)\predictors_anom;        
        ozone_ind_b = predictors_anom\precondtoz(leaveout,k);        
        ozone_ind = precondtoz(leaveout,k) - ozone_ind_b(1).*leaveout';% - ozone_ind_b(2);
        ozoneleft = precondtoz(count,k) - ozone_ind_b(1).*count;% - ozone_ind_b(2);
        
%         ozone_ind = precondtoz(leaveout,k); 
%         ozoneleft = precondtoz(count,k); 
        
        predictors_ind = [ozone_ind,ones(size(ozone_ind))];        
        
        for j = 1:length(latitude)   
            leaveout_mean_ind = nanmean(squeeze(surfacedata(k,leaveout,j,:)),'double');            
            %anomalydata_ind = squeeze(surfacedata(k,leaveout,j,:)) - leaveout_mean_ind;            
            b_anom = predictors_anom\squeeze(surfacedata(k,leaveout,j,:));
            if j == 87
                abc= 1;
            end
            %anomalydata_ind = squeeze(surfacedata(k,leaveout,j,:))-b_anom(1,:).*[1:length(ozone_ind)]' - b_anom(2,:);
            anomalydata_ind = squeeze(surfacedata(k,leaveout,j,:))-b_anom(1,:).*leaveout' - b_anom(2,:);
            med(k,i,j,:) = median(anomalydata_ind);
            %anomalydata_ind = anomalydata_ind - nanmedian(anomalydata_ind);
            r(k,i,j,:) = corr(anomalydata_ind,ozone_ind);
            b_ind(k,i,j,:,:) = predictors_ind\anomalydata_ind;              
            modelprediction_ind(k,i,j,:) = ozoneleft.*squeeze(b_ind(k,i,j,1,:))+squeeze(b_ind(k,i,j,2,:));
            actualdata_ind(k,i,j,:) = squeeze(surfacedata(k,count,j,:)) - b_anom(2,:)' - b_anom(1,:)'*count;% - nanmedian(anomalydata_ind)';                                         
            
            % individual months
            for l = 1:size(surfacedataind,3)
                leaveout_mean_ind_months = nanmean(squeeze(surfacedataind(k,leaveout,l,j,:)),'double');  
                b_anom_months = predictors_anom\squeeze(surfacedataind(k,leaveout,l,j,:));
                
                %anomalydata_ind_months = squeeze(surfacedataind(k,leaveout,l,j,:))-leaveout_mean_ind_months;
                %anomalydata_ind_months = squeeze(surfacedataind(k,leaveout,l,j,:))-b_anom_months(1,:).*[1:length(ozone_ind)]' - b_anom_months(2,:);
                anomalydata_ind_months = squeeze(surfacedataind(k,leaveout,l,j,:))-b_anom_months(1,:).*leaveout' - b_anom_months(2,:);
                medmonths(k,i,l,j,:) = nanmedian(anomalydata_ind_months);
                anomalydata_ind_months = anomalydata_ind_months;% - nanmedian(anomalydata_ind_months);
                r_months(k,i,l,j,:) = corr(anomalydata_ind_months,ozone_ind);
                b_ind_months(k,i,l,j,:,:) = predictors_ind\anomalydata_ind_months;                  
                actualdata_ind_months(k,i,l,j,:) = squeeze(surfacedataind(k,count,l,j,:)) - b_anom_months(2,:)' - b_anom_months(1,:)'*count;% - nanmedian(anomalydata_ind_months)';
                modelprediction_ind_months(k,i,l,j,:) = ozoneleft.*squeeze(b_ind_months(k,i,l,j,1,:))+squeeze(b_ind_months(k,i,l,j,2,:));         
            end
            
        end                            
        count = count+1;
    end
    
    
    %percentile
    allin_pct = 1:size(precondtoz_pct,1);
    surfacedata_pct(k,:,:,:) = surfacedata(k,[pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)],:,:);
    surfacedataind_pct(k,:,:,:,:) = surfacedataind(k,[pct.highCl.ind.highind(k,:),pct.highCl.ind.lowind(k,:)],:,:,:);
    for i = 1:length(allin_pct)   
        leaveout_pct = allin_pct;
        leaveout_pct (leaveout_pct == i) = []; 
        ozone_ind_pct = precondtoz_pct(leaveout_pct,k);
        predictors_ind_pct = [ozone_ind_pct ones(size(ozone_ind_pct))];
        predictors_anom_pct = [[1:length(ozone_ind_pct)]',ones(size(leaveout_pct))'];
        for j = 1:length(latitude)   
            leaveout_mean_ind_pct = nanmean(squeeze(surfacedata_pct(k,leaveout_pct,j,:)),'double');
            b_anom_pct = predictors_anom_pct\squeeze(surfacedata_pct(k,leaveout_pct,j,:));
            %anomalydata_ind_pct = squeeze(surfacedata_pct(k,leaveout_pct,j,:))-b_anom_pct(1,:).*[1:length(ozone_ind_pct)]' - b_anom_pct(2,:);
            anomalydata_ind_pct = squeeze(surfacedata_pct(k,leaveout_pct,j,:)) - leaveout_mean_ind_pct;            
            anomalydata_ind_pct = anomalydata_ind_pct;% - nanmedian(anomalydata_ind_pct);
            r_pct(k,i,j,:) = corr(anomalydata_ind_pct,ozone_ind_pct);
            b_ind_pct(k,i,j,:,:) = predictors_ind_pct\anomalydata_ind_pct;  
            %actualdata_ind_pct(k,i,j,:) = squeeze(surfacedata_pct(k,i,j,:)) - b_anom_pct(2,:)' - b_anom_pct(1,:)'*i;   
            actualdata_ind_pct(k,i,j,:) = squeeze(surfacedata_pct(k,i,j,:)) - leaveout_mean_ind_pct';
            modelprediction_ind_pct(k,i,j,:) = precondtoz_pct(i,k).*squeeze(b_ind_pct(k,i,j,1,:)) + squeeze(b_ind_pct(k,i,j,2,:));          % squeeze(b(i,:,2,:)) + 
            for l = 1:size(surfacedataind_pct,3)
                leaveout_mean_ind_pct_months = nanmean(squeeze(surfacedataind_pct(k,leaveout_pct,l,j,:)),'double');
                %
                b_anom_pct_months = predictors_anom_pct\squeeze(surfacedataind_pct(k,leaveout_pct,l,j,:));
                
                %anomalydata_ind_pct_months = squeeze(surfacedataind(k,leaveout_pct,l,j,:))-b_anom_pct_months(1,:).*[1:length(ozone_ind_pct)]' - b_anom_pct_months(2,:);
                anomalydata_ind_pct_months = squeeze(surfacedataind_pct(k,leaveout_pct,l,j,:)) - leaveout_mean_ind_pct_months;
                medmonths_pct(k,i,l,j,:) = nanmedian(anomalydata_ind_pct_months);
                r_pct_months(k,i,l,j,:) = corr(anomalydata_ind_pct_months,ozone_ind_pct);
                b_ind_pct_months(k,i,l,j,:,:) = predictors_ind_pct\anomalydata_ind_pct_months;  
                %actualdata_ind_pct_months(k,i,l,j,:) = squeeze(surfacedataind_pct(k,i,l,j,:)) - b_anom_pct_months(2,:)' - b_anom_pct_months(1,:)'*i;   
                actualdata_ind_pct_months(k,i,l,j,:) = squeeze(surfacedataind_pct(k,i,l,j,:)) - leaveout_mean_ind_pct_months';
                modelprediction_ind_pct_months(k,i,l,j,:) = precondtoz_pct(i,k).*squeeze(b_ind_pct_months(k,i,l,j,1,:)) + squeeze(b_ind_pct_months(k,i,l,j,2,:));          % squeeze(b(i,:,2,:)) + 
            end
            
        end
        
        
    end
end

%%

%notcons = squeeze(sum(sign(squeeze(b_ind(:,:,:,1,:))),2));
notcons = squeeze(sum(r,2));
notcons (abs(notcons) <= size(r,2)*.1) = 0;
notcons (abs(notcons) > size(r,2)*.1) = 1;
% notcons (abs(notcons) < size(b_ind,2)) = 0;
% notcons (abs(notcons) == size(b_ind,2)) = 1;

notcons_months = squeeze(sum(r_months,2));
notcons_months (abs(notcons_months) <= size(r_months,2)*.1) = 0;
notcons_months (abs(notcons_months) > size(r_months,2)*.1) = 1;

notcons_pct = squeeze(sum(r_pct,2));
notcons_pct (abs(notcons_pct) <= size(r_pct,2)*.2) = 0;
notcons_pct (abs(notcons_pct) > size(r_pct,2)*.2) = 1;

notcons_pct_months = squeeze(sum(r_pct_months,2));
notcons_pct_months (abs(notcons_pct_months) <= size(r_pct_months,2)*.2) = 0;
notcons_pct_months (abs(notcons_pct_months) > size(r_pct_months,2)*.2) = 1;

ozone_anom = repmat(ozone_anom,[1,1,size(surface_anom,3),size(surface_anom,4)]);

%% calculating cheating method
modelpred_all = permute(modelpred_all,[1,4,3,2]);
[GSS.cheat,Heidke.cheat] =  predruns_calcHeidke(modelpred_all,actualdata_all,0);

modelpred_allpct = permute(modelpred_allpct,[1,4,3,2]);
[GSS.cheatpct,Heidke.cheatpct] =  predruns_calcHeidke(modelpred_allpct,actualdata_allpct,0);

modelpred_allmonth = permute(modelpred_allmonth,[1,5,3,4,2]);
[GSS.cheatmonth,Heidke.cheatmonth] =  predruns_calcHeidke(modelpred_allmonth,actualdata_allmonth,0);

modelpred_allpctmonth = permute(modelpred_allpctmonth,[1,5,3,4,2]);
[GSS.cheatpctmonth,Heidke.cheatpctmonth] =  predruns_calcHeidke(modelpred_allpctmonth,actualdata_allpctmonth,0);

%% Calculate skill score all
%[GSS.all,Heidke.all] =  predruns_calcHeidke(modelprediction_ind,actualdata_ind,0);

[GSS.all,Heidke.all] =  predruns_calcHeidke(modelprediction_ind,actualdata_ind-med,0);
%[GSS.all,Heidke.all] =  predruns_calcHeidke(modelprediction_ind,actualdata_ind-med,0);

%% testing bootsrapping
mp = permute(modelprediction_ind,[3,4,1,2]);
%mp = permute(modelprediction_ind-med,[3,4,1,2]);
mp = mp(:,:,:);
ad = permute(actualdata_ind,[3,4,1,2]);
%ad = permute(actualdata_ind-med,[3,4,1,2]);
ad = ad(:,:,:);

mpm = permute(modelprediction_ind_months-medmonths,[3,4,5,1,2]);
mpm = mpm(:,:,:,:);
adm = permute(actualdata_ind_months-medmonths,[3,4,5,1,2]);
adm = adm(:,:,:,:);

%%

for i = 1:size(mp,1)        
    for j = 1:size(mp,2)  
        HSS.composite(i,j) = predruns_calcHeidke_forbootstrap(squeeze(mp(i,j,:)),squeeze(ad(i,j,:)));
    end
end

%%
if ~exist('/Volumes/ExternalOne/work/data/predruns/output/HSS/Allperc.mat','file')
    bootstat = zeros(size(mp,1),size(mp,2),500);
    percentiles.eighty = zeros(size(mp,1),size(mp,2));
    percentiles.ninety = zeros(size(mp,1),size(mp,2));
    percentiles.ninetyfive = zeros(size(mp,1),size(mp,2));
    bootstatm = zeros(5,size(mp,1),size(mp,2),500);
    percentilesm.eighty = zeros(5,size(mp,1),size(mp,2));
    percentilesm.ninety = zeros(5,size(mp,1),size(mp,2));
    percentilesm.ninetyfive = zeros(5,size(mp,1),size(mp,2));
    for i = 1:size(mp,1)        
        tic;
        for j = 1:size(mp,2)       
            mpe = squeeze(mp(i,j,:));
            ade = squeeze(ad(i,j,:));
            bootstat(i,j,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap(mpe,ade),mpe);
            percentiles.eighty(i,j) = prctile(squeeze(bootstat(i,j,:)),80);
            percentiles.ninty(i,j) = prctile(squeeze(bootstat(i,j,:)),90);
            percentiles.ninetyfive(i,j) = prctile(squeeze(bootstat(i,j,:)),95);
            for k = 1:5
                mpme = squeeze(mpm(k,i,j,:));
                adme = squeeze(adm(k,i,j,:));
                bootstatm(k,i,j,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap(mpme,adme),mpme);
                percentilesm.eighty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),80);
                percentilesm.ninty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),90);
                percentilesm.ninetyfive(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),95);
            end
        end
        toc;
    end
    save('/Volumes/ExternalOne/work/data/predruns/output/HSS/Allperc.mat','bootstat','bootstatm','percentiles','percentilesm');
else
    allsig = load('/Volumes/ExternalOne/work/data/predruns/output/HSS/Allperc.mat');
end

%%
p = zeros(size(GSS.all.mean));
p (GSS.all.mean < reshape(allsig.percentiles.ninetyfive,[1,size(allsig.percentiles.ninetyfive)])) = -1;
%p = repmat(p,[2,1,1]);
%%

[GSS.all2,Heidke.all2] =  predruns_calcHeidke(modelprediction_ind,actualdata_ind,0);
% GSS.all.ind(:,1,:,:) = squeeze(GSS.all.ind).*notcons;
% GSS.all.mean(1,:,:) = squeeze(nanmean(GSS.all.ind));

%% extracting upper and lower 20th percentiles
for i = 1:9
    modelprediction_ind_pct2(i,:,:,:) = cat(2,modelprediction_ind(i,pct.highCl.ind.highind(i,:),:,:),modelprediction_ind(i,pct.highCl.ind.lowind(i,:),:,:));
    actualdata_ind_pct2(i,:,:,:) = cat(2,actualdata_ind(i,pct.highCl.ind.highind(i,:),:,:),actualdata_ind(i,pct.highCl.ind.lowind(i,:),:,:));
    med2(i,:,:,:) = cat(2,med(i,pct.highCl.ind.highind(i,:),:,:),med(i,pct.highCl.ind.lowind(i,:),:,:));
    
    modelprediction_ind_pct_months2(i,:,:,:,:) = cat(2,modelprediction_ind_months(i,pct.highCl.ind.highind(i,:),:,:,:),...
        modelprediction_ind_months(i,pct.highCl.ind.lowind(i,:),:,:,:));
    actualdata_ind_pct_months2(i,:,:,:,:) = cat(2,actualdata_ind_months(i,pct.highCl.ind.highind(i,:),:,:,:),...
        actualdata_ind_months(i,pct.highCl.ind.lowind(i,:),:,:,:));
    medmonths_pct2(i,:,:,:,:) = cat(2,medmonths(i,pct.highCl.ind.highind(i,:),:,:,:),...
        medmonths(i,pct.highCl.ind.lowind(i,:),:,:,:));    
end


%% Calculate skill score percentage

%[GSS.pct,Heidke.pct] =  predruns_calcHeidke(modelprediction_ind_pct2-med2,actualdata_ind_pct2-med2,0);
[GSS.pct,Heidke.pct] =  predruns_calcHeidke(modelprediction_ind_pct2,actualdata_ind_pct2-med2,0);
GSS.pct.ind(:,1,:,:) = squeeze(GSS.pct.ind);
GSS.pct.mean(1,:,:) = squeeze(nanmean(GSS.pct.ind));


%% percentile percentile
mp = permute(modelprediction_ind_pct2,[3,4,1,2]);
%mp = permute(modelprediction_ind_pct2-med2,[3,4,1,2]);
mp = mp(:,:,:);
ad = permute(actualdata_ind_pct2-med2,[3,4,1,2]);
ad = ad(:,:,:);

mpm = permute(modelprediction_ind_pct_months2,[3,4,5,1,2]);
%mpm = permute(modelprediction_ind_pct_months2-medmonths_pct2,[3,4,5,1,2]);
mpm = mpm(:,:,:,:);
adm = permute(actualdata_ind_pct_months2-medmonths_pct2,[3,4,5,1,2]);
adm = adm(:,:,:,:);

%%
tic;
if ~exist('/Volumes/ExternalOne/work/data/predruns/output/HSS/Pctperc.mat','file')
    bootstat = zeros(size(mp,1),size(mp,2),500);
    bootstatm = zeros(5,size(mp,1),size(mp,2),500);
        
    percentiles.eighty = zeros(size(mp,1),size(mp,2));
    percentiles.ninety = zeros(size(mp,1),size(mp,2));
    percentiles.ninetyfive = zeros(size(mp,1),size(mp,2));    
    percentilesm.eighty = zeros(5,size(mp,1),size(mp,2));
    percentilesm.ninety = zeros(5,size(mp,1),size(mp,2));
    percentilesm.ninetyfive = zeros(5,size(mp,1),size(mp,2));
    
    for i = 1:size(mp,1)     
        tic;
        for j = 1:size(mp,2)       
            mpe = squeeze(mp(i,j,:));
            ade = squeeze(ad(i,j,:));
            bootstat(i,j,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap(mpe,ade),mpe);            
            
            percentiles.eighty(i,j) = prctile(squeeze(bootstat(i,j,:)),80);
            percentiles.ninty(i,j) = prctile(squeeze(bootstat(i,j,:)),90);
            percentiles.ninetyfive(i,j) = prctile(squeeze(bootstat(i,j,:)),95);
            
            for k = 1:5
                mpme = squeeze(mpm(k,i,j,:));
                adme = squeeze(adm(k,i,j,:));
                bootstatm(k,i,j,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap(mpme,adme),mpme);                
                
                percentilesm.eighty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),80);
                percentilesm.ninty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),90);
                percentilesm.ninetyfive(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),95);
                
            end
        end
        toc;
    end
    save('/Volumes/ExternalOne/work/data/predruns/output/HSS/Pctperc.mat','bootstat','bootstatm','percentiles','percentilesm');
else
    pctsig = load('/Volumes/ExternalOne/work/data/predruns/output/HSS/Pctperc.mat');
end

%%

p2 = zeros(size(GSS.pct.mean));
p2 (GSS.pct.mean < reshape(pctsig.percentiles.ninetyfive,[1,size(pctsig.percentiles.ninetyfive)])) = -1;


%% normalizing to average correlation
% rmean = nanmean(r,2);
% GSS.all.ind (abs(rmean) < .2) = 0;
% GSS.all.mean(1,:,:) = squeeze(nanmean(GSS.all.ind,1));
% rpct_mean = nanmean(r_pct,2);
% GSS.pct.ind (abs(rmean) < .2) = 0;
% GSS.pct.mean(1,:,:) = squeeze(nanmean(GSS.pct.ind,1));
%% Calculate skill score all months

%[GSS.monthsall,Heidke.monthsall] =  predruns_calcHeidke(modelprediction_ind_months-medmonths,actualdata_ind_months-medmonths,0);
[GSS.monthsall,Heidke.monthsall] =  predruns_calcHeidke(modelprediction_ind_months,actualdata_ind_months-medmonths,0);

%% Calculate skill score percentage months

%[GSS.monthspct,Heidke.monthspct] =  predruns_calcHeidke(modelprediction_ind_pct_months2-medmonths_pct2,actualdata_ind_pct_months2-medmonths_pct2,0);
[GSS.monthspct,Heidke.monthspct] =  predruns_calcHeidke(modelprediction_ind_pct_months2,actualdata_ind_pct_months2-medmonths_pct2,0);
   
%% plot     
titles = {'Ensemble Heidke skill score (leave one out)','Ensemble 20th percentile Heidke skill score (leave one out)','Observed percentile general skill score (leave one out)'};
% subplotmaps(permute(cat(1,reshape(squeeze(GSS.cheatmonth.mean(1,2,:,:)),[1,96,144]),reshape(squeeze(GSS.cheatpctmonth.mean(1,2,:,:)),[1,96,144])),[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,titles,'Longitude','Latitude','','on',...
%         [-80,80],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
toplot1 = cat(1,GSS.all.mean,GSS.pct.mean);
toplotp = permute(cat(1,p,p2),[1,3,2]);
% subplotmaps(permute(toplot1,[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,titles,'Longitude','Latitude','','on',...
%          [-40,40],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'-',1,'Miller Cylindrical');
%toplotp(1,:,:) = -1;
subplotmaps(permute(toplot1,[1,3,2]),longitude,latitude,{'seq','YlGnBu'},0,toplotp,16,titles,'Longitude','Latitude','HSS','on',...
    [0,50],11,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
hold on
set(gcf,'Renderer','Painters');
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/Ensmean_GSSothercolor_nodetrend',monthnames(inputs.varmonth,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'.eps'];
print(filename,'-depsc');        
           
%     
%% plot skill lines    

%Find areas of large difference and large GSS
% difference above 3 degrees and GSS above 4
% areas_lons = [290,320;30,120;240,290;30,120];%lons (Greenland,Russia,America,Asia)
% areas_lats = [60,80;50,75;30,60;15,45];%lats

diffforplot = differences;
diffforplot = circshift(diffforplot,[0,144/2,0,0]);

cbrewdiff = cbrewer('div','RdBu',21);
cbrewdiff = flipud(cbrewdiff);
cbrewdiff = [cbrewdiff(1:10,:);cbrewdiff(12:end,:)];

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([4,1,10,8],:);
cbrewqual3 = cbrewqual2./2;

xticklab = {'March','April','May','June','July'};

lstyle = {':','-','-.','--','-'};    
fsize = 16;
Mks = {'d','s','o','v'};



%%
same = 1;
% corrtoplot = 0;
difftoplot = 1;
corrtoplot = 0;
 areas_lons = [80,180;30,80;240,290;60,120];%lons (Greenland,East Russia, West Russia,America,Asia)
areas_lats = [60,80;60,80;35,60;25,45];%lats

for i = 1:2
    [GSSmonthtoplot_areamean(i).g,GSSmonthtoplot_areamean_pct(i).g,GSSmonthtoplot_areamean_pctsame(i).g,...
        forerrorextract(i).g,forerrorextractpct(i).g,lattoplot(i).g,lontoplot(i).g,lattoplot_pct(i).g,lontoplot_pct(i).g] = ...
        predruns_calculateHSStoplot(same,difftoplot,corrtoplot,latitude,longitude,differences,...
        modelprediction_ind_months,actualdata_ind_months,medmonths,modelprediction_ind_pct_months2,...
        actualdata_ind_pct_months2,medmonths_pct2,GSS,areas_lons,areas_lats,r_months,allsig,pctsig);
    same = 0;
end
%%

if corrtoplot
    corrtoplotext = 'corr';
elseif difftoplot
    corrtoplotext = 'diff';
else
    corrtoplotext = 'GSS';
end

if same
    sameext = 'same';
else
    sameext = '';
end


titles = {'Ensemble composite','Ensemble composite 20^t^h percentiles','Ensemble member','Ensemble member 20^t^h percentiles'}; 
labels = {'a','b','c','d'};
for mon = 1:5
    count = 1;
    fig = figure;
    set(fig,'position',[100 1 1000 984],'color','white','Visible','on');
    for j = 1:2
        areas_lons2 = areas_lons;
        areas_lats2 = areas_lats;
               
        sp(count) = subplot(3,2,count);
        sppos = get(sp(count),'position');
        if count == 3
            set(sp(count),'position',[sppos(1)-.025,sppos(2)+.020,sppos(3:4)]);            
        else
            set(sp(count),'position',[sppos(1)-.025,sppos(2:4)]);
        end
        sppos = get(sp(count),'position');
        %toplot = squeeze(nanmean(GSSmonthtoplot,2));
        toplot = squeeze(nanmean(GSSmonthtoplot_areamean(j).g(mon,:,:,:),3));
        totimes = [-3,-1,1,3];
        for i = 1:size(areas_lats,1)   

            eh1(i) = errorbar([1:5]+(totimes(i)./20),toplot(:,i),forerrorextract(j).g(mon,:,i),...
                'linewidth',2,'LineStyle','none','color',cbrewqual2(i,:));
            hold on
            phl(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'MarkerSize',8,'MarkerEdgeColor',cbrewqual3(i,:),'MarkerFaceColor',cbrewqual2(i,:),'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));

            phl2(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},...
                'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'MarkerSize',10);
        end
        uistack(eh1,'bottom')
        uistack(phl2,'top')
        
        annotation('textbox',sppos,'String',labels{count},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+4,... % 
        'EdgeColor','none','fontweight','bold');    
        
        plot([0,6],[0,0],'linewidth',2,'LineStyle','--','color','k');
        set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
        if count == 3
            xlabel('Month','fontsize',fsize+2);
        end
        xlim([.5,5.5]);
        ylim([-20 80]);
        set(gca,'ytick',-20:10:100,'yticklabel',-20:10:100);
        ylabel('HSS','fontsize',fsize+2);
        title(titles{count},'fontsize',fsize+4)
        %March TCO - surface temperature ensemble mean HSS
        if count == 1
            lh = legend(phl,'Eastern Russia','Western Russia','Northern America','Asia');                    
            set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
        end
        
        count = count+1;
        
        sp(count) = subplot(3,2,count);
        sppos = get(sp(count),'position');
        if count == 4            
            set(sp(count),'position',[sppos(1)-.075,sppos(2)+.020,sppos(3:4)]);
        else
            set(sp(count),'position',[sppos(1)-.075,sppos(2:4)]);
        end
        
        sppos = get(sp(count),'position');
        
        if same
            toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pctsame(j).g(mon,:,:,:),3));
        else
            toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pct(j).g(mon,:,:,:),3));
        end
        for i = 1:size(areas_lats,1)   
            ehp1(i) = errorbar([1:5]+(totimes(i)./20),toplot(:,i),squeeze(forerrorextractpct(j).g(mon,:,i)),...
                'linewidth',2,'LineStyle','none','color',cbrewqual2(i,:));
            hold on
            php(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));    
            php2(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},...
                'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'MarkerSize',10);
        end
        uistack(ehp1,'bottom')
        uistack(php2,'top')
        
        annotation('textbox',sppos,'String',labels{count},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+4,... % 
        'EdgeColor','none','fontweight','bold');    
        
        plot([0,6],[0,0],'linewidth',2,'LineStyle','--','color','k');
        set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
        if count == 4
            xlabel('Month','fontsize',fsize+2);
        end
        %ylabel('HSS','fontsize',fsize+2);
        title(titles{count},'fontsize',fsize+4)
        xlim([.5,5.5]);
        ylim([-20 80]);
        set(gca,'ytick',-20:10:100,'yticklabel',-20:10:100);
        % plotting latitudes
        for i = 1:numel(lontoplot(j).g)
            if lontoplot(j).g(i) > 180
                lontoplot(j).g(i) = lontoplot(j).g(i) - 360;
            end
        end

        for i = 1:numel(lontoplot_pct(j).g)
            if lontoplot_pct(j).g(i) > 180
                lontoplot_pct(j).g(i) = lontoplot_pct(j).g(i) - 360;
            end
        end

        for i = 1:numel(areas_lons2)
            if areas_lons2(i) > 180
                areas_lons2(i) = areas_lons2(i) - 360;
            end
        end
        count = count+1;
    end

sp2 = subplot(3,1,3);

% fig = figure;
% set(fig,'position',[100 100 1000 400],'color','white');
m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
[~,h] = m_contourf(longitude-180,latitude,squeeze(nanmean(diffforplot(:,:,:,mon),1))',-5:.5:5,'LineStyle','none');
hold on
m_coast('color','k','LineWidth',1);
m_grid('ytick',[0:15:90],'xtick',[-180:60:360],'XaxisLocation','bottom','fontsize',fsize);  
caxis([-5 5]);
ch = colorbar;
set(ch,'YTick',[-5:1:5],'fontsize',fsize)%,
set(get(ch,'ylabel'),'string','Temperature','fontsize',fsize+2)
colormap(cbrewdiff);

ylab = ylabel('Latitude','fontsize',fsize+2);
ylabpos = get(ylab,'position');
set(ylab,'position',[ylabpos(1)-.25,ylabpos(2),ylabpos(3)]);
xlabel('Longitude','fontsize',fsize+2);
set(sp2,'position',[.14 .075 .65 .32]);
sp2pos = get(sp2,'position');

annotation('textbox',[sp2pos(1),sp2pos(2)-.04,sp2pos(3:4)],'String','e','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+4,... % 
        'EdgeColor','none','fontweight','bold');    

%cbarrow
if corrtoplot
    title(['Location of maximum ensemble March TCO and ',monthnames(mon+2,0,0), ' TS differences'],'fontsize',fsize+2);
elseif difftoplot
    title(['Locations of ensemble March TCO and ',monthnames(mon+2,0,0), ' surface temperature differences'],'fontsize',fsize+2);
    
else
    title(['Location of maximum ensemble member March TCO and ',monthnames(mon+2,0,0), ' HSS'],'fontsize',fsize+2);
end
%a = 1
box on
hold on
for i = 1:size(areas_lats,1)
    %m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
    m_plot([areas_lons2(i,1),areas_lons2(i,2),areas_lons2(i,2),areas_lons2(i,1),areas_lons2(i,1)],[areas_lats(i,1),areas_lats(i,1),areas_lats(i,2),areas_lats(i,2),areas_lats(i,1)],...
        'LineStyle','-','LineWidth',4,'color',cbrewqual2(i,:),'LineStyle',lstyle{i})
    mp(i) = m_plot(lontoplot(2).g(mon,:,i)',lattoplot(2).g(mon,:,i)','LineStyle','none','Marker',Mks{i},'MarkerSize',12,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',1);
    mp2(i) = m_plot(lontoplot(1).g(mon,:,i)',lattoplot(1).g(mon,:,i)','LineStyle','none','Marker','p','MarkerSize',20,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor','k','LineWidth',2);
end
%axes.SortMethod='ChildOrder';

annotation('textbox',[.01 .995 .925 0],'String','March TCO - surface temperature ensemble HSSs','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');    

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/Ind_',corrtoplotext,'_','mean_Locations_withpct_test','_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),monthnames(mon+2,0,0),'both'];
print(filename,'-depsc');
end
%onstruct table of individual GSS differenes and correlations for april for
%those locations.
%table = 


%% BACKUP
% for mon = 1:5
% areas_lons2 = areas_lons;
% areas_lats2 = areas_lats;
% 
% fig = figure;
% set(fig,'position',[100 100 1000 800],'color','white','Visible','on');
% 
% sp(1) = subplot(2,2,1);
% sppos = get(sp(1),'position');
% set(sp(1),'position',[sppos(1)-.025,sppos(2:4)]);
% get(sp(1),'position')
% %toplot = squeeze(nanmean(GSSmonthtoplot,2));
% toplot = squeeze(nanmean(GSSmonthtoplot_areamean,2));
% totimes = [-3,-1,1,3];
% for i = 1:size(areas_lats,1)   
%     
%     eh1(i) = errorbar([1:5]+(totimes(i)./20),toplot(:,i),forerrorextract(:,i),...
%         'linewidth',2,'LineStyle','none','color',cbrewqual2(i,:));
%     hold on
%     phl(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'MarkerSize',8,'MarkerEdgeColor',cbrewqual3(i,:),'MarkerFaceColor',cbrewqual2(i,:),'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
%     
%     phl2(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},...
%         'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'MarkerSize',10);
% end
% uistack(eh1,'bottom')
% uistack(phl2,'top')
% plot([0,6],[0,0],'linewidth',2,'LineStyle','--','color','k');
% set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
% xlabel('Month','fontsize',fsize+2);
% xlim([.5,5.5]);
% ylim([-20 80]);
% set(gca,'ytick',-20:10:100,'yticklabel',-20:10:100);
% ylabel('HSS','fontsize',fsize+2);
% title('Ensemble','fontsize',fsize+4)
% %March TCO - surface temperature ensemble mean HSS
% lh = legend(phl,'Eastern Russia','Western Russia','Northern America','Asia');
% %lh2 = legend(phl2,'','','','');
% set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
% 
% sp(2) = subplot(2,2,2);
% sppos = get(sp(2),'position');
% set(sp(2),'position',[sppos(1)-.075,sppos(2:4)]);
% if same
%     toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pctsame,2));
% else
%     toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pct,2));
% end
% for i = 1:size(areas_lats,1)   
%     ehp1(i) = errorbar([1:5]+(totimes(i)./20),toplot(:,i),forerrorextractpct(:,i),...
%         'linewidth',2,'LineStyle','none','color',cbrewqual2(i,:));
%     hold on
%     php(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));    
%     php2(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},...
%         'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'MarkerSize',10);
% end
% uistack(ehp1,'bottom')
% uistack(php2,'top')
% plot([0,6],[0,0],'linewidth',2,'LineStyle','--','color','k');
% set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
% xlabel('Month','fontsize',fsize+2);
% %ylabel('HSS','fontsize',fsize+2);
% title('Ensemble 20^t^h percentile','fontsize',fsize+4)
% xlim([.5,5.5]);
% ylim([-20 80]);
% set(gca,'ytick',-20:10:100,'yticklabel',-20:10:100);
% % plotting latitudes
% for i = 1:numel(lontoplot)
%     if lontoplot(i) > 180
%         lontoplot(i) = lontoplot(i) - 360;
%     end
% end
% 
% for i = 1:numel(lontoplot_pct)
%     if lontoplot_pct(i) > 180
%         lontoplot_pct(i) = lontoplot_pct(i) - 360;
%     end
% end
% 
% for i = 1:numel(areas_lons2)
%     if areas_lons2(i) > 180
%         areas_lons2(i) = areas_lons2(i) - 360;
%     end
% end
% 
% sp2 = subplot(2,1,2);
% 
% % fig = figure;
% % set(fig,'position',[100 100 1000 400],'color','white');
% m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
% [~,h] = m_contourf(longitude-180,latitude,squeeze(nanmean(diffforplot(:,:,:,mon),1))',-5:.5:5,'LineStyle','none');
% hold on
% m_coast('color','k','LineWidth',1);
% m_grid('ytick',[0:15:90],'xtick',[-180:60:360],'XaxisLocation','bottom','fontsize',fsize);  
% caxis([-5 5]);
% ch = colorbar;
% set(ch,'YTick',[-5:1:5],'fontsize',fsize)%,
% set(get(ch,'ylabel'),'string','Temperature','fontsize',fsize+2)
% colormap(cbrewdiff);
% 
% ylab = ylabel('Latitude','fontsize',fsize+2);
% ylabpos = get(ylab,'position');
% set(ylab,'position',[ylabpos(1)-.25,ylabpos(2),ylabpos(3)]);
% xlabel('Longitude','fontsize',fsize+2);
% set(sp2,'position',[.14 .18 .65 .32]);
% %cbarrow
% if corrtoplot
%     title(['Location of maximum ensemble March TCO and ',monthnames(mon+2,0,0), ' TS differences'],'fontsize',fsize+2);
% elseif difftoplot
%     title(['Location of maximum ensemble March TCO and ',monthnames(mon+2,0,0), ' TS differences'],'fontsize',fsize+2);
%     
% else
%     title(['Location of maximum ensemble member March TCO and ',monthnames(mon+2,0,0), ' HSS'],'fontsize',fsize+2);
% end
% 
% box on
% hold on
% for i = 1:size(areas_lats,1)
%     %m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
%     m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker',Mks{i},'MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
%     m_plot([areas_lons2(i,1),areas_lons2(i,2),areas_lons2(i,2),areas_lons2(i,1),areas_lons2(i,1)],[areas_lats(i,1),areas_lats(i,1),areas_lats(i,2),areas_lats(i,2),areas_lats(i,1)],...
%         'LineStyle','-','LineWidth',4,'color',cbrewqual2(i,:),'LineStyle',lstyle{i})
% end
% axes.SortMethod='ChildOrder';
% 
% annotation('textbox',[.01 1 .925 0],'String','March TCO - surface temperature ensemble HSS','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
%         'EdgeColor','none','fontweight','bold');    
% 
% filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/Ind_',corrtoplotext,'_','mean_Locations_withpct_test','_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),monthnames(mon+2,0,0),sameext];
% print(filename,'-depsc');
% end


end
