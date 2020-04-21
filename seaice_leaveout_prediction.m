function [] = seaice_leaveout_prediction(surfacedata,tozdata,inputs,latitude,longitude,differences,observations,Seasplot)




%% STEPS
% 1. Leave three out of both the predictor and predictand time series
% 2. Construct prediction times series of ozone by first removing the linear trend 
% 3. Calculate sea ice extent or snow depth 

% over 1995-2014. Leave one year out of the regression and see if I can
% predict the remaining years anomaly

% find average temperature anomaly of upper 50th percentile of ozone
% find average temperature anomaly of lower 50th percentile of ozone
% two cases that can be used as a binary prediction (shift in pdf)

if inputs.lats(1) < 0
    ext = 1;
else
    ext = 0;
end

%pcttouse = [10,90];
%pcttouse = [20,80];
pcttouse = [23,77];

%% calculate sea ice extent or snow depth time series.
obs.data = permute(observations.montharrange,[1,2,4,3]);
if strcmp(inputs.var,'ICEFRAC')
    icefrac = .15;
    areacell = gridcellarea(latitude,longitude);
    areacellobs = gridcellarea(observations.latitude,observations.longitude);
    noseas = size(Seasplot.lon,1);
        
    modeldata = surfacedata.highCl.dataMonthArrange;
    modeldata (modeldata < inputs.fraction) = 0;
    modeldata (modeldata >= inputs.fraction) = 1;
    obs.data (obs.data < inputs.fraction) = 0;
    obs.data (obs.data >= inputs.fraction) = 1;
else
    noseas = size(Seasplot.mod.lon,1);
    modeldata = surfacedata.highCl.dataMonthArrange;
end
        
%% rearrange data into composite

for i = 1:size(tozdata.highCl.dataMonthArrange,1)
    precondtoz(:,i) = squeeze(tozdata.highCl.dataMonthArrange(i,inputs.tozmonth,:));%-nanmean(squeeze(tozdata(i,inputs.tozmonth,:)));    
end

tozdata_composite = precondtoz(:);

%surfacedata = double(surfacedata.highCl.dataMonthArrange);


%% observations
obs.allin = 1:size(observations.toz.zm,2);
obs.toz = observations.toz.zm;
obs.upperpctind.u = find(obs.toz >= prctile(obs.toz,80)); 
obs.lowerpctind.l = find(obs.toz <= prctile(obs.toz,20)); 
count = 1;
laglength = 3;
for i = 1:length(obs.allin)
    
    obs.leaveout = obs.allin;
    if i == 1
        obs.leaveout (obs.leaveout == count | obs.leaveout == count+1) = [];             
    elseif i == length(obs.allin)
        obs.leaveout (obs.leaveout == count | obs.leaveout == count-1) = []; 
    else
        obs.leaveout (obs.leaveout == count-1 | obs.leaveout == count | obs.leaveout == count+1) = [];         
    end


    obs.leaveoutenso = [obs.leaveout,obs.leaveout(end)+ext];
    obs.predictors_anom = [obs.leaveout',ones(size(obs.leaveout))'];        

    obs.ozone_ind_b = obs.predictors_anom\obs.toz(obs.leaveout)';        
    obs.ozone_anomaly = obs.toz(obs.leaveout) - obs.ozone_ind_b(1).*obs.leaveout';% - ozone_ind_b(2); 

    obs.ozone_pct(i) = prctile(obs.toz(obs.leaveout),50);
    obs.ozone_pct1(i) = prctile(obs.toz(obs.leaveout),pcttouse(1));
    obs.ozone_pct2(i) = prctile(obs.toz(obs.leaveout),pcttouse(2));

    obs.ozonelowerind(i).m = find(obs.toz(obs.leaveout) < obs.ozone_pct1(i));
    obs.ozoneupperind(i).m = find(obs.toz(obs.leaveout) > obs.ozone_pct2(i));                      

    obs.ozone_left(i) = obs.toz(count) - obs.ozone_pct(i);

    obs.predsign(i) = sign(obs.ozone_left(i));            

    for j = 1:noseas
        if strcmp(inputs.var,'ICEFRAC')
            obslatind = observations.latitude >= Seasplot.lat(j,1) & observations.latitude <= Seasplot.lat(j,2);
            if j == 3
                obslonind = observations.longitude >= Seasplot.lon(j,1) | observations.longitude < Seasplot.lon(j,2)-360;             
            else
                obslonind = observations.longitude >= Seasplot.lon(j,1) & observations.longitude < Seasplot.lon(j,2);             
            end            
            obs.dataextent(j).a = permute(areacellobs(obslatind,obslonind).*permute(obs.data(:,:,obslatind,obslonind),[3,4,1,2]),[3,4,1,2]);
            obs.dataextentfinal(j,:,:) = nansum(obs.dataextent(j).a(:,:,:),3);
        else
            obslatind = observations.latitude >= Seasplot.lat(j,1) & observations.latitude <= Seasplot.lat(j,2);
            obslonind = observations.longitude >= Seasplot.lon(j,1) & observations.longitude < Seasplot.lon(j,2);    
            obs.dataextent(j).a = obs.data(:,:,obslatind,obslonind);            
            obs.dataextentfinal(j,:,:) = nansum(obs.dataextent(j).a(:,:,:),3);            
            
            %obsdifftemp = permute(repmat(permute(nanmean(differences.observations.difference([2,3,4],obslonind,obslatind),1),[3,2,1]),[1,1,size(obs.dataextent(j).a,1),size(obs.dataextent(j).a,2)]),[3,4,1,2]);            
            obsdifftemp = permute(repmat(permute(nanmean(differences.observations.difference(1,obslonind,obslatind),1),[3,2,1]),[1,1,size(obs.dataextent(j).a,1),size(obs.dataextent(j).a,2)]),[3,4,1,2]);            
            if i < 3
                obs.dataextent(j).a (obsdifftemp > 0) = NaN;            
            else
                obs.dataextent(j).a (obsdifftemp < 0) = NaN;            
            end
            %obs.dataextentfinal(j,:,:) = nanmean(obs.dataextent(j).a(:,:,:),3);            

%             moddiff = differences.difference.ens;
%             moddiffseason = repmat(nanmean(differences.difference.ens(:,:,[2,3,4]),3),[1,1,length(inputs.varmonth)]);                                        
            
        end   
        
         for m = 1:size(obs.data,1)                
             if inputs.removeENSO
                for lag = 1:laglength                                                            
                    ENSOyearind = reshape(observations.ENSO,12,length(observations.ENSO)./12);
                    ENSOleft = ENSOyearind(:,count);
                    ENSOyearind = ENSOyearind(:,obs.leaveout);
                    ENSOyearind = ENSOyearind(:);
                    
                    if m == 1 || m == 2
                        m2 = 3;
                    else
                        m2 = m;
                    end
                    ensopredictors = [squeeze(ENSOyearind(m2-lag+1:12:end,:)),ones(size(obs.leaveout))'];
                    benso(lag,:,:) = ensopredictors\squeeze(obs.dataextentfinal(j,m,obs.leaveout));

%                     ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonth(m)-lag+1:12:end,:))',ones(size(alldata,3)-ext,1)];
%                     benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedataind(k,:,m,j,:));                    
                end
                %find max lag
                [~,bensomax_ind] = max(abs(benso),[],1);                
                %[~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   

                bensomax = benso(bensomax_ind(1),1);
                    %bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
                obs_leaveoutenso = squeeze(obs.dataextentfinal(j,m,obs.leaveout)) - bensomax.*squeeze(ENSOyearind(m2-bensomax_ind(1)+1:12:end,:));   
                obs_leftenso = squeeze(obs.dataextentfinal(j,m,count)) - bensomax.*squeeze(ENSOleft(m2-bensomax_ind(1)+1));   

                obs.b_months = obs.predictors_anom\squeeze(obs_leaveoutenso);
                %seaiceextent_leaveout_detrend(j,k,m,1:length(leaveout)) = squeeze(seaiceextent(j,k,m,leaveout)) - leaveout'.*b_months(1,:);

                obs_leaveout_detrend = squeeze(obs.dataextentfinal(j,m,obs.leaveout)) - obs.leaveout'.*obs.b_months(1,:);                    
                obs_left_detrend = obs_leftenso - count.*obs.b_months(1,:);        
          
                clearvars obs_leaveoutenso

            
            else
                 % now detrend sea ice extent data from the left out                      

                obs.b_months = obs.predictors_anom\squeeze(obs.dataextentfinal(j,m,obs.leaveout));
                %seaiceextent_leaveout_detrend(j,k,m,1:length(leaveout)) = squeeze(seaiceextent(j,k,m,leaveout)) - leaveout'.*b_months(1,:);
                obs_leaveout_detrend = squeeze(obs.dataextentfinal(j,m,obs.leaveout)) - obs.leaveout'.*obs.b_months(1,:);                    
                obs_left_detrend = squeeze(obs.dataextentfinal(j,m,count)) - count.*obs.b_months(1,:);        

            end
            % calculating anomaly on data point left out
            obs.left_months(i,j,m) = obs_left_detrend - nanmedian(obs_leaveout_detrend);

            % calculating anomaly of the training dataset
            obs.anomaly_months = obs_leaveout_detrend - nanmedian(obs_leaveout_detrend);

            % extracting mean of ozone upper 20th percentile in training dataset
            obs.anomaly_upper_months = squeeze(nanmedian(obs.anomaly_months(obs.ozoneupperind(i).m)));

            % extracting mean of ozone lower 20th percentile in training dataset
            obs.anomaly_lower_months = squeeze(nanmedian(obs.anomaly_months(obs.ozonelowerind(i).m)));                        

            % calculating the sign of the difference in extremes
            obs.signchange_months(i,j,m) = sign(obs.anomaly_upper_months-obs.anomaly_lower_months);

            % difference greater than 10 percent

            % calculating the sign of the difference in the data point left out
            obs.datasign_months(i,j,m) = sign(obs.left_months(i,j,m));

            % deciding sign of prediction
            if obs.signchange_months(i,j,m) == 1 %&& abs(rcond) > .1
                obs.pred_months(i,j,m) = 1;
            elseif obs.signchange_months(i,j,m) == -1 %&& abs(rcond) > .1
                obs.pred_months(i,j,m) = -1;        
            elseif obs.signchange_months(i,j,m) == 0 
                obs.pred_months(i,j,m) = 0;        
            else
                obs.pred_months(i,j,m) = 0;        
            end

            if obs.pred_months(i,j,m) == 1
                obs.predsign_months2(i,j,m) = obs.predsign(i).*1;
                if obs.predsign(i) == obs.datasign_months(i,j,m)
                    obs.correct_months(i,j,m) = 1;
                else
                    obs.correct_months(i,j,m) = 0;
                end
            elseif obs.pred_months(i,j,m) == -1
                obs.predsign_months2(i,j,m) = obs.predsign(i).*-1;
                if obs.predsign(i) == obs.datasign_months(i,j,m)
                    obs.correct_months(i,j,m) = 0;
                else
                    obs.correct_months(i,j,m) = 1;
                end
            elseif obs.pred_months(i,j,m) == 0
                obs.predsign_months2(i,j,m) = obs.predsign(i).*0;
                if obs.predsign(i) == obs.datasign_months(i,j,m)
                    obs.correct_months(i,j,m) = NaN;
                else
                    obs.correct_months(i,j,m) = NaN;
                end
            end   
        end
                                
    end
    count = count+1;
    clearvars obs_leaveout_detrend
end

%% conditions
obs.condition1 = repmat(abs(sum(obs.signchange_months,1)),[37,1,1]);
obs.condition1 (obs.condition1 < 1) = NaN; %23
obs.condition1 (obs.condition1 >= 1) = 1;  %23
obs.correct_months2 = obs.correct_months.*obs.condition1;
obs.predsign_months3 = obs.predsign_months2.*obs.condition1;
obs.datasign_months3 = obs.datasign_months.*obs.condition1;

obs.correct_months_pct = obs.correct_months2([obs.upperpctind.u,obs.lowerpctind.l],:,:);

obs.predsignpct_months = obs.predsign_months3([obs.upperpctind.u,obs.lowerpctind.l],:,:);
obs.datasignpct_months = obs.datasign_months3([obs.upperpctind.u,obs.lowerpctind.l],:,:);

obs.GSS.monthsall = (nansum(obs.correct_months2,1) - sum(~isnan(obs.correct_months2),1)./2)./...
    (sum(~isnan(obs.correct_months2),1)-sum(~isnan(obs.correct_months2),1)./2)*100;  
obs.GSS.monthspct = (nansum(obs.correct_months_pct,1) - sum(~isnan(obs.correct_months_pct),1)./2)./...
    (sum(~isnan(obs.correct_months_pct),1)-sum(~isnan(obs.correct_months_pct),1)./2)*100;  

%% calculate significance
obs_mpm = permute(obs.predsignpct_months,[3,2,1]);
obs_adm = permute(obs.datasignpct_months,[3,2,1]);
 for j = 1:size(obs_mpm,2)       

    for k = 1:12
        obs_mpme = squeeze(obs_mpm(k,j,:));
        obs_adme = squeeze(obs_adm(k,j,:));
        obs_bootstatm(k,j,:) = bootstrp(1000, @(obs_mpme) predruns_calcHeidke_forbootstrap_compemp(obs_mpme,obs_adme,1),obs_mpme);
        obs_percentilesm.eighty(k,j) = prctile(squeeze(obs_bootstatm(k,j,:)),80);
        obs_percentilesm.ninety(k,j) = prctile(squeeze(obs_bootstatm(k,j,:)),90);
        obs_percentilesm.ninetyfive(k,j) = prctile(squeeze(obs_bootstatm(k,j,:)),95);
    end
end
 obs_percentilesm.ninetyfive (obs_percentilesm.ninetyfive == -100) = NaN;
 obs_percentilesm.ninety (obs_percentilesm.ninety == -100) = NaN;
 obs_percentilesm.eighty (obs_percentilesm.eighty == -100) = NaN;

%% each individual separately 

%%

for k = 1:size(surfacedata.highCl.dataMonthArrange,1)
    tic;
    allin = 1:size(precondtoz,1);
    
    tozautocorr = precondtoz(:,k);
    upperpctind(k).u = find(tozautocorr >= prctile(tozautocorr,80));
    lowerpctind(k).l = find(tozautocorr <= prctile(tozautocorr,20));
    count = 1;
    for i = 1:length(allin)
        leaveout = allin;                                    
        if i == 1
            leaveout (leaveout == count | leaveout == count+1) = [];             
        elseif i == 30
            leaveout (leaveout == count | leaveout == count-1) = []; 
        else
            leaveout (leaveout == count-1 | leaveout == count | leaveout == count+1) = [];         
        end

        leaveoutenso = [leaveout,leaveout(end)+ext];

        predictors_anom = [leaveout',ones(size(leaveout))'];        

        ozone_ind_b = predictors_anom\precondtoz(leaveout,k);        
        ozone_anomaly = precondtoz(leaveout,k) - ozone_ind_b(1).*leaveout';% - ozone_ind_b(2); 
       
        % calculating 50th percentile 
        ozone_pct(k,i) = prctile(precondtoz(leaveout,k),50);
        ozone_pct1(k,i) = prctile(precondtoz(leaveout,k),pcttouse(1));
        ozone_pct2(k,i) = prctile(precondtoz(leaveout,k),pcttouse(2));

        ozonelowerind(k,i).m = find(precondtoz(leaveout,k) < ozone_pct1(k,i));
        ozoneupperind(k,i).m = find(precondtoz(leaveout,k) > ozone_pct2(k,i));                      

        ozone_left(k,i) = precondtoz(count,k) - ozone_pct(k,i);
        
        predsign(k,i) = sign(ozone_left(k,i));        
              
       

        for j = 1:noseas
            
            if strcmp(inputs.var,'ICEFRAC')
                latind = latitude >= Seasplot.lat(j,1) & latitude <= Seasplot.lat(j,2);
                if j == 3
                    lonind = longitude >= Seasplot.lon(j,1) | longitude < Seasplot.lon(j,2)-360;
                 %   lonindobs = observations.longitude >= Seasplot.lon(j,1) | observations.longitude < Seasplot.lon(j,2)-360;
                else
                    lonind = longitude >= Seasplot.lon(j,1) & longitude < Seasplot.lon(j,2);
                  %  lonindobs = observations.longitude >= Seasplot.lon(j,1) & observations.longitude < Seasplot.lon(j,2);
                end            
                dataextent(j).a = permute(areacell(latind,lonind).*permute(modeldata(:,:,:,latind,lonind),[4,5,1,2,3]),[3,4,5,1,2]);
                dataextentfinal(j,:,:,:) = nansum(dataextent(j).a(:,:,:,:),4);
            else
                latind = latitude >= Seasplot.lat(j,1) & latitude <= Seasplot.lat(j,2);
                lonind = longitude >= Seasplot.lon(j,1) & longitude < Seasplot.lon(j,2);
                dataextent(j).a = modeldata(:,:,:,latind,lonind);
                dataextentfinal(j,:,:,:) = nanmean(dataextent(j).a(:,:,:,:),4);
                
                %difftemp = permute(repmat(nansum(differences.difference.ens(latind,lonind,[2,3,4]),3),[1,1,size(dataextent(j).a,1),size(dataextent(j).a,2),size(dataextent(j).a,3)]),[3,4,5,1,2]);
                difftemp = permute(repmat(nansum(differences.difference.ens(latind,lonind,1),3),[1,1,size(dataextent(j).a,1),size(dataextent(j).a,2),size(dataextent(j).a,3)]),[3,4,5,1,2]);
                
                if i < 3
                    dataextent(j).a (difftemp > 0) = NaN;            
                else
                    dataextent(j).a (difftemp < -.025) = NaN;            
                end
                %dataextentfinal(j,:,:,:) = nanmean(dataextent(j).a(:,:,:,:),4);         

            end          

            for m = 1:size(surfacedata.highCl.dataMonthArrange,2)                                
                %remove enso
                % remove ENSO

                
                if inputs.removeENSO
                    for lag = 1:laglength                                                            
                        ENSOyearind = reshape(surfacedata.highCl.ENSO(k,:),12,size(surfacedata.highCl.ENSO,2)./12);
                        ENSOleft = ENSOyearind(:,count);
                        ENSOyearind = ENSOyearind(:,leaveout);
                        ENSOyearind = ENSOyearind(:);

                        if m == 1 || m == 2
                            m2 = 3;
                        else
                            m2 = m;
                        end
                        ensopredictors = [squeeze(ENSOyearind(m2-lag+1:12:end,:)),ones(size(leaveout))'];
                        benso(lag,:,:) = ensopredictors\squeeze(dataextentfinal(j,k,m,leaveout));

%                     ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonth(m)-lag+1:12:end,:))',ones(size(alldata,3)-ext,1)];
%                     benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedataind(k,:,m,j,:));                    
                    end
                    %find max lag
                    [~,bensomax_ind] = max(abs(benso),[],1);                
                    %[~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   

                    bensomax = benso(bensomax_ind(1),1);
                        %bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
                    leaveoutenso = squeeze(dataextentfinal(j,k,m,leaveout)) - bensomax.*squeeze(ENSOyearind(m2-bensomax_ind(1)+1:12:end,:));   
                    leftenso = squeeze(dataextentfinal(j,k,m,count)) - bensomax.*squeeze(ENSOleft(m2-bensomax_ind(1)+1));   

                    b_months = predictors_anom\squeeze(leaveoutenso);
                    %seaiceextent_leaveout_detrend(j,k,m,1:length(leaveout)) = squeeze(seaiceextent(j,k,m,leaveout)) - leaveout'.*b_months(1,:);

                    leaveout_detrend = squeeze(dataextentfinal(j,k,m,leaveout)) - leaveout'.*b_months(1,:);                    
                    left_detrend = leftenso - count.*b_months(1,:);        

                    clearvars obs_leaveoutenso   
                else
                     % now detrend sea ice extent data from the left out                      

                    b_months = predictors_anom\squeeze(dataextentfinal(j,k,m,leaveout));
                    %dataextentfinal_leaveout_detrend(j,k,m,1:length(leaveout)) = squeeze(dataextentfinal(j,k,m,leaveout)) - leaveout'.*b_months(1,:);
                    leaveout_detrend = squeeze(dataextentfinal(j,k,m,leaveout)) - leaveout'.*b_months(1,:);
                    left_detrend = squeeze(dataextentfinal(j,k,m,count)) - count.*b_months(1,:);                       
                end

%                 if m == 10 && j == 1
%                     abc = 1;
%                 end
                
                % calculating anomaly on data point left out
                left_months(k,i,j,m) = left_detrend - nanmedian(leaveout_detrend);
                
                % calculating anomaly of the training dataset
                anomaly_months = leaveout_detrend - nanmedian(leaveout_detrend);
                
                % extracting mean of ozone upper 20th percentile in training dataset
                anomaly_upper_months = squeeze(nanmean(anomaly_months(ozoneupperind(k,i).m)));
                
                % extracting mean of ozone lower 20th percentile in training dataset
                anomaly_lower_months = squeeze(nanmean(anomaly_months(ozonelowerind(k,i).m)));                        
                
%                 % calculating r condition
%                 rcond = corr(anomaly_months,ozone_anomaly);
                
                % calculating the sign of the difference in extremes
                signchange_months(k,i,j,m) = sign(anomaly_upper_months-anomaly_lower_months);

                % difference greater than 10 percent
                
                % calculating the sign of the difference in the data point left out
                datasign_months(k,i,j,m) = sign(left_months(k,i,j,m));
                
                % deciding sign of prediction
                if signchange_months(k,i,j,m) == 1 %&& abs(rcond) > .1
                    pred_months(k,i,j,m) = 1;
                elseif signchange_months(k,i,j,m) == -1 %&& abs(rcond) > .1
                    pred_months(k,i,j,m) = -1;        
                elseif signchange_months(k,i,j,m) == 0 
                    pred_months(k,i,j,m) = 0;        
                else
                    pred_months(k,i,j,m) = 0;        
                end

                if pred_months(k,i,j,m) == 1
                    predsign_months2(k,i,j,m) = predsign(k,i).*1;
                    if predsign(k,i) == datasign_months(k,i,j,m)
                        correct_months(k,i,j,m) = 1;
                    else
                        correct_months(k,i,j,m) = 0;
                    end
                elseif pred_months(k,i,j,m) == -1
                    predsign_months2(k,i,j,m) = predsign(k,i).*-1;
                    if predsign(k,i) == datasign_months(k,i,j,m)
                        correct_months(k,i,j,m) = 0;
                    else
                        correct_months(k,i,j,m) = 1;
                    end
                elseif pred_months(k,i,j,m) == 0
                    predsign_months2(k,i,j,m) = predsign(k,i).*0;
                    if predsign(k,i) == datasign_months(k,i,j,m)
                        correct_months(k,i,j,m) = NaN;
                    else
                        correct_months(k,i,j,m) = NaN;
                    end
                end                

            end

        end    
        count = count+1;
        clearvars leaveout_detrend
    end   
    toc;
end
%% condition where if the sign change is not 90% the same, the prediction is thrown out.
condition1 = repmat(abs(sum(signchange_months,2)),[1,30,1,1]);
condition1 (condition1 < 1) = NaN; %18
condition1 (condition1 >= 1) = 1;  %18
correct_months2 = correct_months.*condition1;
predsign_months3 = predsign_months2.*condition1;
datasign_months3 = datasign_months.*condition1;

% correct_months2 = correct_months;
% predsign_months3 = predsign_months2;
% datasign_months3 = datasign_months;

%%

% extracting pct
for i = 1:size(tozdata.highCl.dataMonthArrange,1)

    correct_months_pct(i,:,:,:) = correct_months2(i,[upperpctind(i).u;lowerpctind(i).l],:,:);

    predsignpct_months(i,:,:,:) = predsign_months3(i,[upperpctind(i).u;lowerpctind(i).l],:,:);
    datasignpct_months(i,:,:,:) = datasign_months3(i,[upperpctind(i).u;lowerpctind(i).l],:,:);

end

%%
month = 10;
sea = 1;
mem = 9;
noleavout = load('/Volumes/ExternalOne/work/data/predruns/output/ice/predtesting.mat');
noleavout.test = squeeze(noleavout.ice_extremes(month,:,:,sea))';
noleavout.testall = squeeze(noleavout.iceall(month,:,:,sea))';
test = squeeze(left_months(:,:,sea,month));
testpred1 = squeeze(signchange_months(:,:,sea,month));
testpred = squeeze(pred_months(:,:,sea,month));
testpred2 = squeeze(predsign_months2(:,:,sea,month));
for i = 1:9
    testpct(i,:) = test(i,[upperpctind(i).u;lowerpctind(i).l]);
    testpredpct(i,:) = testpred(i,[upperpctind(i).u;lowerpctind(i).l]);
    testpredpct1(i,:) = testpred1(i,[upperpctind(i).u;lowerpctind(i).l]);
    testpredpct2(i,:) = testpred2(i,[upperpctind(i).u;lowerpctind(i).l]);
    ozonetest(i,:) = ozone_left(i,[upperpctind(i).u;lowerpctind(i).l]);
    predsignpct(i,:) = predsign(i,[upperpctind(i).u;lowerpctind(i).l]);
end

% figure; plot(testpct(mem,:),'-o'); hold on; plot(noleavout.test(mem,:),'-ro'); yyaxis right; plot(ozonetest(mem,:),'-o');
% figure; plot(test(mem,:),'-o'); yyaxis right; plot(ozone_left(mem,:),'-o'); hold on; plot([0 31],[0 0],'k--');
% figure; plot(noleavout.testall(mem,:),'-o'); yyaxis right; plot(precondtoz(:,mem) - nanmedian(precondtoz(:,mem)),'-ro'); hold on; plot([0 31],[0 0],'k--');
%yyaxis right; hold on; plot(testpredpct2(9,:),'-ro');


%%

GSS.monthsall.ind = (nansum(correct_months2,2) - sum(~isnan(correct_months2),2)./2)./...
    (sum(~isnan(correct_months2),2)-sum(~isnan(correct_months2),2)./2)*100;  

% GSS.monthspcttest.ind = (nansum(correct_months_pct,2) - size(correct_months_pct,2)./2)./...
%     (size(correct_months_pct,2)-size(correct_months_pct,2)./2)*100;  
GSS.monthspct.ind = (nansum(correct_months_pct,2) - sum(~isnan(correct_months_pct),2)./2)./...
    (sum(~isnan(correct_months_pct),2)-sum(~isnan(correct_months_pct),2)./2)*100;  
GSS.monthsall.mean(1,:,:) = squeeze(nanmean(GSS.monthsall.ind));
GSS.monthspct.mean(1,:,:) = squeeze(nanmean(GSS.monthspct.ind));


%% temp plot
% seanames = {'BKS','LEC and OB','GS'};
% createfig('medium','on')
% plot(3:12,squeeze(GSS.monthspct.mean(1,:,3:12))','LineWidth',2)
% ylabel('HSS','fontsize',20);
% xlabel('Month','fontsize',20);
% set(gca,'fontsize',18);
% title('HSSs during ozone extremes','fontsize',22);
% lh = legend(seanames);
% set(lh,'fontsize',20,'location','southeast','box','off');
%  filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
% filename = [filedir,'Seaice_extendHSS_LINEPLOT_from_',...
%     monthnames(inputs.tozmonth,1,'long'),'_',inputs.obstouse,'_Arcticozoneextremes_over_',...
%     num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
% export_fig(filename,'-pdf');
%% bootsrapping

mpm = permute(predsignpct_months,[4,3,1,2]);
mpm = mpm(:,:,:);
adm = permute(datasignpct_months,[4,3,1,2]);
adm = adm(:,:,:);




bootstatm = zeros(12,size(mpm,2),1000);
percentilesm.eighty = zeros(12,size(mpm,2));
percentilesm.ninety = zeros(12,size(mpm,2));
percentilesm.ninetyfive = zeros(12,size(mpm,2));
%for i = 1:size(mpm,3)        
    tic;
    for j = 1:size(mpm,2)       

        for k = 1:12       
            mpme = squeeze(mpm(k,j,:));
            adme = squeeze(adm(k,j,:));
            mpme (isnan(mpme)) = [];
            adme (isnan(adme)) = [];
            bootstatm(k,j,:) = bootstrp(1000, @(mpme) predruns_calcHeidke_forbootstrap_compemp(mpme,adme,0),mpme);
            percentilesm.eighty(k,j) = prctile(squeeze(bootstatm(k,j,:)),80);
            percentilesm.ninety(k,j) = prctile(squeeze(bootstatm(k,j,:)),90);
            percentilesm.ninetyfive(k,j) = prctile(squeeze(bootstatm(k,j,:)),95);
        end
    end
    toc;
    %end
 
%p2 (GSS.monthspct.mean < reshape(allsigpct.percentilesm.ninetyfive,[1,size(allsigpct.percentilesm.ninetyfive)])) = -1;


%% plot     


%% plot skill lines    

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([4,1,8,10],:);
cbrewqual3 = cbrewqual2./2;

if inputs.lats(1) < 0
    xticklab = {'November','December','January','February','March'};
else    
    xticklab = {'Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
end

lstyle = {':','-','--'};    
fsize = 16;
Mks = {'d','s','o'};
xshift = [-.15,0,.15];
%% 
if strcmp(inputs.var,'ICEFRAC')
    xin = 3:12;   
else
    xin = 3:7;
end

figure;
set(gcf,'position',[100 100 1400 350],'color','white');

plot([2,13],[0,0],'--k','LineWidth',2);
for j = 1:2
    sp(j) = subplot(1,2,j);
    sppos(j,:) = get(sp(j),'position');
    for i = 1:size(GSS.monthspct.mean,2)
        hold on
        if j == 1
            ph(i) = plot([xin]+xshift(i),squeeze(obs.GSS.monthspct(1,i,xin)),'LineStyle',lstyle{i},'color',cbrewqual2(i,:),...
                'LineWidth',3,'Marker',Mks{i},'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),...
                 'MarkerSize',12);
            hold on
            errorbar([xin]+xshift(i),squeeze(obs.GSS.monthspct(1,i,xin)),obs_percentilesm.ninetyfive(xin,i),'LineStyle',lstyle{i},'color',cbrewqual2(i,:),...
                'LineWidth',3,'Marker',Mks{i},'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),...
                'MarkerSize',15);
            title('Observations, ozone extremes','fontsize',fsize+4)
            if strcmp(inputs.var,'ICEFRAC')
                ylim([-100 120]);   
                ylab = [-100:20:120];
            else
                ylim([-100 130]);   
                ylab = -100:50:100;
            end
        else
            errorbar([xin]+xshift(i),squeeze(GSS.monthspct.mean(1,i,xin)),percentilesm.ninetyfive(xin,i),'LineStyle',lstyle{i},'color',cbrewqual2(i,:),...
                'LineWidth',3,'Marker',Mks{i},'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),...
                'MarkerSize',15);
            %xlabel('Month','fontsize',fsize+2);
            title('Model composite, ozone extremes','fontsize',fsize+4)
            set(sp(j),'position',[sppos(j,1)-.05,sppos(j,2),sppos(2,3),sppos(1,4)-.01]);
           if strcmp(inputs.var,'ICEFRAC')
                ylim([-50 70]);    
                ylab = [-50:10:70];
            else
                ylim([-50 50]);    
                ylab = [-50:10:50];
            end
        end
        xlim([xin(1)-1,xin(end)+1])
        box on
        %ylabel('HSS','fontsize',fsize+2);
        set(gca,'fontsize',fsize+2,'xtick',xin,'xticklabel',xticklab,'ytick',ylab);
    end
    if j == 1
        if strcmp(inputs.var,'ICEFRAC')
            lh = legend(ph,'BK','OBCLEB','GN');    
        else
            lh = legend(ph,'Northern Russia','Eastern Russia','South East Alaska');    
        end
        
        set(lh,'box','off','fontsize',fsize+2,'location','SouthWest');
    end
end

% filename
 filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
filename = [filedir,'Prediction_',num2str(pcttouse(1)),'-',num2str(pcttouse(2)),inputs.var,'_LINEPLOT_from_',...
    monthnames(inputs.tozmonth,1,'long'),'_',inputs.obstouse,'_Arcticozoneextremes_over_',...
    num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
export_fig(filename,'-pdf');


end