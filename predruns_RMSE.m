function [] = predruns_RMSE(surfacedata,surfacedataind,tozdata,inputs,latitude,longitude,pct,differences,correlations)

% over 1995-2014. Leave one year RMSE cross validation out of the regression and see if I can
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
        ozone_ind = precondtoz(leaveout,k);        
        predictors_ind = [ozone_ind,ones(size(ozone_ind))];        
        predictors_anom = [[1:length(ozone_ind)]',ones(size(leaveout))'];
        ozoneleft = precondtoz(count,k);
        for j = 1:length(latitude)   
            leaveout_mean_ind = nanmean(squeeze(surfacedata(k,leaveout,j,:)),'double');            
            %anomalydata_ind = squeeze(surfacedata(k,leaveout,j,:)) - leaveout_mean_ind;            
            b_anom = predictors_anom\squeeze(surfacedata(k,leaveout,j,:));
            if j == 87
                abc= 1;
            end
            %anomalydata_ind = squeeze(surfacedata(k,leaveout,j,:))-b_anom(1,:).*[1:length(ozone_ind)]' - b_anom(2,:);
            anomalydata_ind = squeeze(surfacedata(k,leaveout,j,:))-b_anom(1,:).*[1:length(ozone_ind)]';
            med(k,i,j,:) = median(anomalydata_ind);
            %anomalydata_ind = anomalydata_ind - nanmedian(anomalydata_ind);
            r(k,i,j,:) = corr(anomalydata_ind,ozone_ind);
            b_ind(k,i,j,:,:) = predictors_ind\anomalydata_ind;              
            modelprediction_ind(k,i,j,:) = ozoneleft.*squeeze(b_ind(k,i,j,1,:))+squeeze(b_ind(k,i,j,2,:));
            %actualdata_ind(k,i,j,:) = squeeze(surfacedata(k,count,j,:)) - b_anom(2,:)' - b_anom(1,:)'*count;% - nanmedian(anomalydata_ind)';                                         
            actualdata_ind(k,i,j,:) = squeeze(surfacedata(k,count,j,:)) - b_anom(1,:)'*count;% - nanmedian(anomalydata_ind)';                                         
            
            % individual months
            for l = 1:size(surfacedataind,3)
                leaveout_mean_ind_months = nanmean(squeeze(surfacedataind(k,leaveout,l,j,:)),'double');  
                b_anom_months = predictors_anom\squeeze(surfacedataind(k,leaveout,l,j,:));
                
                %anomalydata_ind_months = squeeze(surfacedataind(k,leaveout,l,j,:))-leaveout_mean_ind_months;
                anomalydata_ind_months = squeeze(surfacedataind(k,leaveout,l,j,:))-b_anom_months(1,:).*[1:length(ozone_ind)]';% - b_anom_months(2,:);
                medmonths(k,i,l,j,:) = nanmedian(anomalydata_ind_months);
                %anomalydata_ind_months = anomalydata_ind_months;% - nanmedian(anomalydata_ind_months);
                r_months(k,i,l,j,:) = corr(anomalydata_ind_months,ozone_ind);
                b_ind_months(k,i,l,j,:,:) = predictors_ind\anomalydata_ind_months;                  
                actualdata_ind_months(k,i,l,j,:) = squeeze(surfacedataind(k,count,l,j,:)) - b_anom_months(1,:)'*count;% - nanmedian(anomalydata_ind_months)';
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
            anomalydata_ind_pct = squeeze(surfacedata_pct(k,leaveout_pct,j,:));
            r_pct(k,i,j,:) = corr(anomalydata_ind_pct,ozone_ind_pct);
            b_ind_pct(k,i,j,:,:) = predictors_ind_pct\anomalydata_ind_pct;  
            %actualdata_ind_pct(k,i,j,:) = squeeze(surfacedata_pct(k,i,j,:)) - b_anom_pct(2,:)' - b_anom_pct(1,:)'*i;   
            actualdata_ind_pct(k,i,j,:) = squeeze(surfacedata_pct(k,i,j,:));% - leaveout_mean_ind_pct';
            modelprediction_ind_pct(k,i,j,:) = precondtoz_pct(i,k).*squeeze(b_ind_pct(k,i,j,1,:)) + squeeze(b_ind_pct(k,i,j,2,:));          % squeeze(b(i,:,2,:)) + 
            for l = 1:size(surfacedataind_pct,3)
                leaveout_mean_ind_pct_months = nanmean(squeeze(surfacedataind_pct(k,leaveout_pct,l,j,:)),'double');
                %anomalydata_ind_pct_months = squeeze(surfacedataind_pct(k,leaveout_pct,l,j,:)) - leaveout_mean_ind_pct_months;
                b_anom_pct_months = predictors_anom_pct\squeeze(surfacedataind_pct(k,leaveout_pct,l,j,:));
                
                anomalydata_ind_pct_months = squeeze(surfacedataind(k,leaveout_pct,l,j,:));
                medmonths_pct(k,i,l,j,:) = nanmedian(anomalydata_ind_pct_months);
                r_pct_months(k,i,l,j,:) = corr(anomalydata_ind_pct_months,ozone_ind_pct);
                b_ind_pct_months(k,i,l,j,:,:) = predictors_ind_pct\anomalydata_ind_pct_months;  
                actualdata_ind_pct_months(k,i,l,j,:) = squeeze(surfacedataind_pct(k,i,l,j,:));% - b_anom_pct_months(2,:)' - b_anom_pct_months(1,:)'*i;   
                modelprediction_ind_pct_months(k,i,l,j,:) = precondtoz_pct(i,k).*squeeze(b_ind_pct_months(k,i,l,j,1,:)) + squeeze(b_ind_pct_months(k,i,l,j,2,:));          % squeeze(b(i,:,2,:)) + 
            end
            
        end
        
        
    end
end

%% calculate RMSE

qnt25 = prctile(actualdata_ind,25,2);
qnt75 = prctile(actualdata_ind,75,2);
maxval = max(actualdata_ind,[],2);
minval = min(actualdata_ind,[],2);
datamean = nanmean(actualdata_ind,2);

maxvalmonths = max(actualdata_ind_months,[],2);
minvalmonths = min(actualdata_ind_months,[],2);

RMSE = sqrt(nanmean((actualdata_ind - modelprediction_ind).^2,2))./(qnt75-qnt25);
RMSEmaxmin = sqrt(nanmean((actualdata_ind - modelprediction_ind).^2,2))./(maxval-minval);
RMSEmaxmin_months = sqrt(nanmean((actualdata_ind_months - modelprediction_ind_months).^2,2))./(maxvalmonths-minvalmonths);
RMSE_raw = sqrt(nanmean((actualdata_ind - modelprediction_ind).^2,2));
RMSE_mean = sqrt(nanmean((actualdata_ind - modelprediction_ind).^2,2))./abs(datamean);

RMSEmaxmin_mean = reshape(nanmean(RMSEmaxmin),[1,size(RMSEmaxmin,3),size(RMSEmaxmin,4)]);
RMSEmaxmin_months_mean = reshape(nanmean(RMSEmaxmin_months),[1,size(RMSEmaxmin_months,3),size(RMSEmaxmin_months,4),size(RMSEmaxmin_months,5)]);

%% calculate RMSE pct
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

RMSEmaxmin_pct = sqrt(nanmean((actualdata_ind_pct2 - modelprediction_ind_pct2).^2,2))./(maxval-minval);
RMSEmaxmin_months_pct = sqrt(nanmean((actualdata_ind_pct_months2 - modelprediction_ind_pct_months2).^2,2))./(maxvalmonths-minvalmonths);

RMSEmaxmin_mean_pct = reshape(nanmean(RMSEmaxmin_pct),[1,size(RMSEmaxmin_pct,3),size(RMSEmaxmin_pct,4)]);
RMSEmaxmin_months_mean_pct = reshape(nanmean(RMSEmaxmin_months_pct),[1,size(RMSEmaxmin_months_pct,3),size(RMSEmaxmin_months_pct,4),size(RMSEmaxmin_months_pct,5)]);

%%
% qnt25_pct = prctile(actualdata_ind_pct,25,2);
% qnt75_pct = prctile(actualdata_ind_pct,75,2);
% maxval_pct = max(actualdata_ind_pct,[],2);
% minval_pct = min(actualdata_ind_pct,[],2);
% datamean_pct = nanmean(actualdata_ind_pct,2);
% 
% maxvalmonths_pct = max(actualdata_ind_pct_months,[],2);
% minvalmonths_pct = min(actualdata_ind_pct_months,[],2);
% 
% RMSE_pct = sqrt(nanmean((actualdata_ind_pct - modelprediction_ind_pct).^2,2))./(qnt75_pct-qnt25_pct);
% RMSEmaxmin_pct = sqrt(nanmean((actualdata_ind_pct - modelprediction_ind_pct).^2,2))./(maxval_pct-minval_pct);
% RMSEmaxmin_months_pct = sqrt(nanmean((actualdata_ind_pct_months - modelprediction_ind_pct_months).^2,2))./(maxvalmonths_pct-minvalmonths_pct);
% 
% RMSEmaxmin_mean_pct = reshape(nanmean(RMSEmaxmin_pct),[1,size(RMSEmaxmin_pct,3),size(RMSEmaxmin_pct,4)]);
% RMSEmaxmin_months_mean_pct = reshape(nanmean(RMSEmaxmin_months_pct),[1,size(RMSEmaxmin_months_pct,3),size(RMSEmaxmin_months_pct,4),size(RMSEmaxmin_months_pct,5)]);

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
[GSS.all,Heidke.all] =  predruns_calcHeidke(modelprediction_ind-med,actualdata_ind-med,0);
[GSS.all2,Heidke.all2] =  predruns_calcHeidke(modelprediction_ind,actualdata_ind,0);
% GSS.all.ind(:,1,:,:) = squeeze(GSS.all.ind).*notcons;
% GSS.all.mean(1,:,:) = squeeze(nanmean(GSS.all.ind));

% normalizing to average correlation
rmean = nanmean(r,2);
GSS.all.ind (abs(rmean) < .2) = 0;
GSS.all.mean(1,:,:) = squeeze(nanmean(GSS.all.ind,1));
%% Calculate skill score percentage

[GSS.pct,Heidke.pct] =  predruns_calcHeidke(modelprediction_ind_pct,actualdata_ind_pct,0);
GSS.pct.ind(:,1,:,:) = squeeze(GSS.pct.ind).*notcons_pct;
GSS.pct.mean(1,:,:) = squeeze(nanmean(GSS.pct.ind));

%% Calculate skill score all months

[GSS.monthsall,Heidke.monthsall] =  predruns_calcHeidke(modelprediction_ind_months-medmonths,actualdata_ind_months-medmonths,0);
[GSS.monthsall2,Heidke.monthsall2] =  predruns_calcHeidke(modelprediction_ind_months,actualdata_ind_months,0);

rmonthmean = nanmean(r_months,2);
GSS.monthsall.ind (abs(rmonthmean) < .2) = 0;
GSS.monthsall.mean(1,:,:,:) = squeeze(nanmean(GSS.monthsall.ind,1));
%GSS.monthsall.ind(:,1,:,:,:) = squeeze(GSS.monthsall.ind).*notcons_months;
%GSS.monthsall.mean(1,:,:,:) = squeeze(nanmean(GSS.monthsall.ind));
%% Calculate skill score percentage months

[GSS.monthspct,Heidke.monthspct] =  predruns_calcHeidke(modelprediction_ind_pct_months-medmonths_pct,actualdata_ind_pct_months-medmonths_pct,0);

   
%% plot     
titles = {'Heidke skill score (leave one out) ensemble mean','Percentile Heidke skill score (leave one out) ensemble mean','Observed percentile general skill score (leave one out)'};
% subplotmaps(permute(cat(1,reshape(squeeze(GSS.cheatmonth.mean(1,2,:,:)),[1,96,144]),reshape(squeeze(GSS.cheatpctmonth.mean(1,2,:,:)),[1,96,144])),[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,titles,'Longitude','Latitude','','on',...
%         [-80,80],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
toplot1 = GSS.all.mean;
subplotmaps(permute(toplot1,[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,titles,'Longitude','Latitude','','on',...
         [-80,80],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
    
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/Ensmean_GSS_',monthnames(inputs.varmonth,1,1)];
%export_fig(filename,'-pdf');        

% toplot2 = GSS.cheat.mean;
% subplotmaps(permute(toplot2,[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,titles,'Longitude','Latitude','','on',...
%          [-80,80],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

     %%
titles = {'Cross validation NRMSE (leave one out) ensemble mean','20th Percentile cross validation NRMSE (leave one out) ensemble mean','Observed percentile general skill score (leave one out)'};
toplot3 = reshape(squeeze(nanmean(RMSEmaxmin)),[1,size(squeeze(nanmean(RMSEmaxmin)))]);
toplot31 = reshape(squeeze(nanmean(RMSEmaxmin_pct)),[1,size(squeeze(nanmean(RMSEmaxmin_pct)))]);
toplot5 = cat(1,toplot3,toplot31);
subplotmaps(permute(toplot5,[1,3,2]),longitude,latitude,{'seq','YlGnBu'},0,[],16,titles,'Longitude','Latitude','','on',...
         [.2,.3],11,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');     
     
% toplot4 = reshape(squeeze(nanmean(RMSE)),[1,size(squeeze(nanmean(RMSEmaxmin)))]);
% subplotmaps(permute(toplot4,[1,3,2]),longitude,latitude,{'seq','YlOrBr'},0,[],16,titles,'Longitude','Latitude','','on',...
%          [.4,1.4],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');          
     
%export_fig(filename,'-pdf');        
%print(filename,'-depsc');

% subplotmaps(permute(squeeze(cat(1,GSS.all.ind(1,1,:,:),GSS.all.ind(3,1,:,:),GSS.all.ind(5,1,:,:))),[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,{'No. 1','No. 3','No. 5'},'Longitude','Latitude','','on',...
%         [-80,80],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
%     
% filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/Ind_GSS_',monthnames(inputs.varmonth,1,1)];
%     export_fig(filename,'-png');            
%     
% subplotmaps(permute(squeeze(cat(1,GSS.pct.ind(1,1,:,:),GSS.pct.ind(3,1,:,:),GSS.pct.ind(5,1,:,:))),[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,{'No. 1','No. 3','No. 5'},'Longitude','Latitude','','on',...
%     [-80,80],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
%     
% filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/Indpct_GSS_',monthnames(inputs.varmonth,1,1)];
%     export_fig(filename,'-png');            
%     
%     
% subplotmaps(permute(cat(1,Heidke.all.mean,Heidke.pct.mean),[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,{'Heidke skill score (leave three out) mean','Percentile Heidke skill score (leave one out) mean'},'Longitude','Latitude','','on',...
%         [-1,1],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
%     
% filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/Ensmean_HSS_',monthnames(inputs.varmonth,1,1)];
%     export_fig(filename,'-png');            
%     
%% plot skill lines    
    
%% plot line plots
fsize = 16;
% obslats = [70,67,50,30];%Greenland,Russia,USA,India
% obslons = [310,70,260,78];
complats = [67,65,65,48,30];
complons = [305,140,60,238,62];

% complats = [68,67,52,36];
% complons = [310,70,255,75];

% for i = 1:length(obslats)
%     [~,latind(i)] = min(abs(obslats(i)-ERAdata.latitude));
%     [~,lonind(i)] = min(abs(obslons(i)-ERAdata.longitude));
% end
for i = 1:length(complats)
    [~,modlatind(i)] = min(abs(complats(i)-latitude));
    [~,modlonind(i)] = min(abs(complons(i)-longitude));
end

% lnames = {['Greenland, (',num2str(abs(obslons(1))),'{\circ}','W, ',num2str(abs(obslats(1))),'{\circ}','N)'],...
%     ['Russia, (',num2str(abs(obslons(2))),'{\circ}','W, ',num2str(abs(obslats(2))),'{\circ}','N)'],...
%     ['North America, (',num2str(abs(obslons(3))),'{\circ}','W, ',num2str(abs(obslats(3))),'{\circ}','N)'],...
%     ['Himalayas, (',num2str(abs(obslons(4))),'{\circ}','W, ',num2str(abs(obslats(4))),'{\circ}','N)']};
lnamesmod = {['Greenland, (',num2str(abs(complons(1))),'{\circ}','W, ',num2str(abs(complats(1))),'{\circ}','N)'],...
    ['Eastern Russia, (',num2str(abs(complons(2))),'{\circ}','W, ',num2str(abs(complats(2))),'{\circ}','N)'],...
    ['Western Russia, (',num2str(abs(complons(3))),'{\circ}','W, ',num2str(abs(complats(3))),'{\circ}','N)'],...
    ['North America, (',num2str(abs(complons(4))),'{\circ}','W, ',num2str(abs(complats(4))),'{\circ}','N)'],...
    ['Asia, (',num2str(abs(complons(5))),'{\circ}','W, ',num2str(abs(complats(5))),'{\circ}','N)']};

xticklab = {'March','April','May','June','July'};

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([3,4,10,8,1],:);

fig = figure;
set(fig,'position',[100 100 1000 700],'color','white');
lstyle = {':','-','-.','--','-'};        
% sp(1) = subplot(1,2,1);
% sp_pos(1,:) = get(sp(1),'position');
% for i = 1:length(obslats)    
%     pho(i) = plot(differences_ind(:,lonind(i),latind(i)),'o','linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
%     hold on
% end
% plot([0 5],[0,0],'--k','lineWidth',2);
% xlim([.5,5.5]);
% ylim([-7 7]);
% set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
% xlabel('Month','fontsize',fsize+2);
% ylabel('Temperature difference (K)','fontsize',fsize+2);
% title('Observed temperature difference','fontsize',fsize+4)
% lh = legend(pho,lnames);
% set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
% sp(2) = subplot(1,2,2);
monthstd = std(GSS.cheatmonth.ind,0,1);
totimes = [-1.5,-.5,.5,1.5];
for i = 1:length(complats)          
    %plot(squeeze(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:))','linewidth',1,'LineStyle',lstyle{i},'color',cbrewqual2(i,:)./1.5);
%     ph(i) = plot([1:5]+(totimes(i)./20),squeeze(GSS.monthsall.mean(1,:,modlatind(i),modlonind(i))),'o',...
%         'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
    
    ph(i) = plot([1:5],squeeze(RMSEmaxmin_months_mean(1,:,modlatind(i),modlonind(i))),'o',...
        'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
    
    hold on
%     errorbar([1:5]+(totimes(i)./20),squeeze(GSS.monthsall.mean(1,:,modlatind(i),modlonind(i),:)),...
%         squeeze(std(GSS.monthsall.ind(:,:,:,modlatind(i),modlonind(i)),0,1)),...
%         'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
end
% sp_pos(2,:) = get(sp(2),'position');
% set(sp(2),'position',[sp_pos(2,1)-.05,sp_pos(2,2),sp_pos(2,3),sp_pos(1,4)]);
% set(sp(1),'position',[sp_pos(1,1),sp_pos(1,2),sp_pos(1,3),sp_pos(1,4)]);
%plot(meandiff,'linewidth',3);
hold on
plot([0 6],[0,0],'--k','lineWidth',2);
xlim([.5,5.5]);
ylim([0 .4]);
set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
xlabel('Month','fontsize',fsize+2);
ylabel('HSS','fontsize',fsize+2);
title('March TCO - surface temperature ensemble mean skill lines','fontsize',fsize+4)
lh = legend(ph,lnamesmod);
set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
axes.SortMethod='ChildOrder';
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/North_skillines_RMSE'];
export_fig(filename,'-pdf');

%% 
for mon = 1:5
same = 1;
corrtoplot = 1;

if same
    sameext = 'same';
else
    sameext = '';
end
if corrtoplot
    corrtoplotext = 'corr';
else
    corrtoplotext = 'GSS';
end
areas_lons = [300,320;80,180;30,80;240,290;60,120];%lons (Greenland,East Russia, West Russia,America,Asia)
areas_lats = [60,80;60,80;60,80;35,60;25,45];%lats

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([3,4,1,10,8],:);

r_monthsmean = squeeze(nanmean(r_months,1));
for i = 1:size(areas_lons,1)
    lats = latitude > areas_lats(i,1) & latitude < areas_lats(i,2);
    lons = longitude > areas_lons(i,1) & longitude < areas_lons(i,2);
    latextract = latitude(lats);
    lonextract = longitude(lons);
    [latmesh,lonmesh] = meshgrid(latextract,lonextract);
    
    latmesh = latmesh';
    lonmesh = lonmesh';
    for k = 1:size(differences,1)    
        
%         mult = squeeze(GSS.monthsall.ind(k,1,mon,lats,lons));
%         mult_pct = squeeze(GSS.monthspct.ind(k,1,mon,lats,lons));
        
%         mult = squeeze(GSS.all.mean(1,lats,lons));        
%         mult_pct = squeeze(GSS.pct.mean(1,lats,lons));        
        
        if corrtoplot
            if same
                mult = squeeze(nanmean(r_monthsmean(:,mon,lats,lons),1));        
                mult_pct = squeeze(nanmean(r_monthsmean(:,mon,lats,lons),1));        
            else
                mult = squeeze(nanmean(r_months(k,:,mon,lats,lons),2));        
                mult_pct = squeeze(nanmean(r_months(k,:,mon,lats,lons),2));        
            end
        else
            if ~same
                mult = squeeze(RMSEmaxmin_months(k,1,mon,lats,lons));
                mult_pct = squeeze(RMSEmaxmin_months_pct(k,1,mon,lats,lons));
            else
                mult = squeeze(nanmean(RMSEmaxmin_months(:,1,mon,lats,lons),1));
                mult_pct = squeeze(nanmean(RMSEmaxmin_months(:,1,mon,lats,lons),1));
%                 mult = squeeze(GSS.all.mean(1,lats,lons));        
%                 mult_pct = squeeze(GSS.pct.mean(1,lats,lons));        
            end
        end
        %mult = squeeze(nanmean(r_monthsmean(:,mon,lats,lons),1));        
        %mult = squeeze(nanmean(r(k,:,lats,lons),2));        
        %mult = squeeze(GSS.cheatpct.ind(k,1,lats,lons));
        [maxval(i,k),maxind] = min(abs(mult(:)));
        [maxval_pct(i,k),maxind_pct] = min(abs(mult_pct(:)));
        %[~,maxind2] = max(mult2(:));
        lattoplot(k,i) = latmesh(maxind);
        lontoplot(k,i) = lonmesh(maxind);
        
        lattoplot_pct(k,i) = latmesh(maxind_pct);
        lontoplot_pct(k,i) = lonmesh(maxind_pct);
        
        %[maxrow(k,i),maxcol(k,i)] = find(mult == max(mult(:)));
        %GSSmonthtoplot2 = squeeze(GSS.monthsall.ind(k,1,:,lats,lons));
        GSSmonthtoplot2 = squeeze(RMSEmaxmin_months(k,1,:,lats,lons));
        
        GSSmonthtoplot3 = squeeze(RMSEmaxmin_months_pct(k,1,:,lats,lons));
        
        if maxind <= size(mult,1)
            if maxind ~= 1
                GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                    GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                                
                
                GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                    GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                                
                
            else
                GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),...
                    GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                            
                
                GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),...
                    GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                            
            end
        elseif maxind == size(mult,1)+1
            GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
                GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        
            
            GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
                GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        
            
        elseif maxind == numel(mult) - size(mult,1)
            GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)),...
                GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        
            
            GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)),...
                GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                    
            
        elseif maxind > numel(mult) - size(mult,1)
            if maxind ~= numel(mult)
                GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                    GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                
                
                GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                    GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                
                
            else                
                GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind-1),...
                    GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                
                
                GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind-1),...
                    GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                
            end
        else
            GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
                GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        
            
            GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
                GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        
            
        end
        
        %pct
        if maxind_pct <= size(mult_pct,1)
            if maxind_pct ~= 1               
                GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                    GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
            else                
            
                GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),...
                    GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
            end
        elseif maxind_pct == size(mult_pct,1)+1            
            
            GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
                GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
            
        elseif maxind_pct == numel(mult_pct) - size(mult_pct,1)
            
            GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)),...
                GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
            
        elseif maxind_pct > numel(mult_pct) - size(mult_pct,1)
            if maxind_pct ~= numel(mult_pct)                
                
                GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                    GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
                
            else                
                
                GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct-1),...
                    GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
            end
        else
           
            GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
                GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
            
        end
        
        
    end
end


areas_lons2 = areas_lons;
areas_lats2 = areas_lats;

fig = figure;
set(fig,'position',[100 100 1000 1200],'color','white','Visible','off');

subplot(3,2,1)
%toplot = squeeze(nanmean(GSSmonthtoplot,2));
toplot = squeeze(nanmean(GSSmonthtoplot_areamean,2));
for i = 1:5
    plot(toplot(:,i),'o','linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
    hold on
end
set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
xlabel('Month','fontsize',fsize+2);
xlim([.5,5.5]);
ylim([.1 .3]);
ylabel('NRMSE','fontsize',fsize+2);
title('Ensemble mean','fontsize',fsize+4)
%March TCO - surface temperature ensemble mean HSS
lh = legend('Greenland','Eastern Russia','Western Russia','North America','Asia');
set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')

subplot(3,2,2)
if same
    toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pctsame,2));
else
    toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pct,2));
end
for i = 1:5    
    plot(toplot(:,i),'o','linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
    hold on
end
set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
xlabel('Month','fontsize',fsize+2);
%ylabel('HSS','fontsize',fsize+2);
title('Ensemble mean 20th percentile','fontsize',fsize+4)
xlim([.5,5.5]);
ylim([.1 .3]);
% plotting latitudes
for i = 1:numel(lontoplot)
    if lontoplot(i) > 180
        lontoplot(i) = lontoplot(i) - 360;
    end
end

for i = 1:numel(lontoplot_pct)
    if lontoplot_pct(i) > 180
        lontoplot_pct(i) = lontoplot_pct(i) - 360;
    end
end

for i = 1:numel(areas_lons2)
    if areas_lons2(i) > 180
        areas_lons2(i) = areas_lons2(i) - 360;
    end
end

subplot(3,1,2);
cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([3,4,1,10,8],:);
cbrewqual3 = cbrewqual2./2;
% fig = figure;
% set(fig,'position',[100 100 1000 400],'color','white');
m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
m_coast('color','k','LineWidth',1);
m_grid('ytick',[0:15:90],'xtick',[-180:60:360],'XaxisLocation','bottom','fontsize',fsize);            
ylab = ylabel('Latitude','fontsize',fsize+2);
ylabpos = get(ylab,'position');
set(ylab,'position',[ylabpos(1)-.2,ylabpos(2),ylabpos(3)]);
xlabel('Longitude','fontsize',fsize+2);
if corrtoplot
    title(['Location of maximum ensemble average March TCO and ',monthnames(mon+2,0,0), ' TS correlation'],'fontsize',fsize+2);
else
    title(['Location of minimum ensemble member March TCO and ',monthnames(mon+2,0,0), ' NRMSE'],'fontsize',fsize+2);
end

box on
hold on
for i = 1:5
    %m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
    m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
    m_plot([areas_lons2(i,1),areas_lons2(i,2),areas_lons2(i,2),areas_lons2(i,1),areas_lons2(i,1)],[areas_lats(i,1),areas_lats(i,1),areas_lats(i,2),areas_lats(i,2),areas_lats(i,1)],...
        'LineStyle','-','LineWidth',4,'color',cbrewqual2(i,:),'LineStyle',lstyle{i})
end
axes.SortMethod='ChildOrder';

% pct
if ~same && ~corrtoplot
    subplot(3,1,3);
    m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
    m_coast('color','k','LineWidth',1);
    m_grid('ytick',[0:15:90],'xtick',[-180:60:360],'XaxisLocation','bottom','fontsize',fsize);            
    ylab = ylabel('Latitude','fontsize',fsize+2);
    ylabpos = get(ylab,'position');
    set(ylab,'position',[ylabpos(1)-.2,ylabpos(2),ylabpos(3)]);
    xlabel('Longitude','fontsize',fsize+2);
    if corrtoplot
        title(['Location of maximum ensemble member March TCO and ',monthnames(mon+2,0,0), ' TS correlation'],'fontsize',fsize+2);
    else
        title(['Location of minimum ensemble member March TCO and ',monthnames(mon+2,0,0), ' 20th percentile NRMSE'],'fontsize',fsize+2);
    end
    box on
    hold on
    for i = 1:5
        %m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
        m_plot(lontoplot_pct(:,i)',lattoplot_pct(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
        m_plot([areas_lons2(i,1),areas_lons2(i,2),areas_lons2(i,2),areas_lons2(i,1),areas_lons2(i,1)],[areas_lats(i,1),areas_lats(i,1),areas_lats(i,2),areas_lats(i,2),areas_lats(i,1)],...
            'LineStyle','-','LineWidth',4,'color',cbrewqual2(i,:),'LineStyle',lstyle{i})
    end
    axes.SortMethod='ChildOrder';
end

annotation('textbox',[.01 .98 1 0],'String','March TCO - surface temperature ensemble mean NRMSE','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');    

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/NRMSE_Ind_',corrtoplotext,'_','mean_Locations_withpct_test','_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),monthnames(mon+2,0,0),sameext];
print(filename,'-depsc');
end











end
