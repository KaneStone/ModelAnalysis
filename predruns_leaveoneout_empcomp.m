function [] = predruns_leaveoneout_empcomp(surfacedata,surfacedataind,tozdata,inputs,latitude,longitude,pct,differences,diffcomp,alldata)

% over 1995-2014. Leave one year out of the regression and see if I can
% predict the remaining years anomaly

% Do this for the ensemble mean and produce a map metric of predictability
% for the ensemble average.

% Using the high signal to noise plots, extract regions near these places
% for each 

% keep leaving years out

% find average temperature anomaly of upper 50th percentile of ozone
% find average temperature anomaly of lower 50th percentile of ozone
% two cases that can be used as a binary prediction (shift in pdf)

if inputs.lats(1) < 0
    ext = 1;
else
    ext = 0;
end

%% rearrange data into composite

%% rearrange data into composite

surfacedata_composite = permute(surfacedata,[3,4,2,1]);
surfacedata_composite = surfacedata_composite(:,:,:); 

for i = 1:size(tozdata,1)
    precondtoz(:,i) = squeeze(tozdata(i,inputs.tozmonth,:));%-nanmean(squeeze(tozdata(i,inputs.tozmonth,:)));    
end

tozdata_composite = precondtoz(:);

surfacedata = double(surfacedata);
surfacedataind = double(surfacedataind);

%% each individual separately 

%% testing
%close all
testing = 0;

if testing
    %sd1 = squeeze(surfacedata(1,:,79,31));
    
    sd1 = squeeze(surfacedata(1,:,88,106));
    
    %sd1 = squeeze(surfacedata(1,:,87,22));
    %sd1 = squeeze(surfacedata(1,:,80,40));
    %sd1 = squeeze(surfacedata(1,:,87,30));
    %sd1 = double(squeeze(surfacedata(1,:,79,35)));
    %sd1 = squeeze(surfacedata(1,:,77,44));
    %sd1 = squeeze(surfacedata(1,:,90,60));
    %sd1 = squeeze(surfacedata(1,:,85,45));           
    %sd1 = squeeze(surfacedata(1,:,14,105));           
    %sd1 = squeeze(surfacedata(1,:,80,79));
    % leaveout
    
    for i = 1:length(sd1)
        
        sd12 = 1:length(sd1);
        sd12 (sd12 == i) = [];
        
        predictors_anom = [sd12',ones(size(sd12))'];
        %ozone_ind_b = precondtoz(leaveout,k)\predictors_anom;        
%         ozone_ind_b = predictors_anom\precondtoz(sd12,1);        
        ozone_pct_test(i) = prctile(precondtoz(sd12,1),50);
        
        ozonelowerind_test = find(precondtoz(sd12,1) < ozone_pct_test(i));
        ozoneupperind_test = find(precondtoz(sd12,1) > ozone_pct_test(i));                      
                     
        ozone_left_test(i) = precondtoz(i,1) - ozone_pct_test(i);
        
        %temp_anomaly = detrend(sd1(sd12));
        b = predictors_anom\sd1(sd12)';
        temp_anomaly_test(i,:) = sd1(sd12) - b(1).*sd12;% - b(2);
        temp_left_test(i) = sd1(i) - b(1).*i;% - b(2);
        
        temp_left_test(i) = temp_left_test(i) - nanmedian(temp_anomaly_test(i,:));
        temp_anomaly_test(i,:) = temp_anomaly_test(i,:) - nanmedian(temp_anomaly_test(i,:));
        
        temp_anomaly_upper_test(i) = nanmean(temp_anomaly_test(i,ozoneupperind_test),2);
        temp_anomaly_lower_test(i) = nanmean(temp_anomaly_test(i,ozonelowerind_test),2);
        
        signlower_test(i) = sign(temp_anomaly_lower_test(i));
        signupper_test(i) = sign(temp_anomaly_upper_test(i));        
        
        predsign_test(i) = sign(ozone_left_test(i));
        
        datasign_test(i) = sign(temp_left_test(i));
        
        if signlower_test(i) < 0 
            pred_test(i) = 1;
        elseif signlower_test(i) > 0 
            pred_test(i) = -1;        
        end
        
        if pred_test(i) == 1
            if predsign_test(i) == datasign_test(i)
               correct_test(i) = 1;
            else
                correct_test(i) = 0;
            end
        elseif pred_test(i) == -1
            if predsign_test(i) == datasign_test(i)
               correct_test(i) = 0;
            else
                correct_test(i) = 1;
            end
        end
        
        %r(i) = corr(temp_anomaly_test(i,:)',ozone_anomaly_test(i,:)');        
    end

        figure;
        plot(temp_left_test)
        hold on
        %plot(ozone_left_test)


    GSScorrecttest = (sum(correct_test,2) - size(correct_test,2)./2)./...
        (size(correct_test,2)-size(correct_test,2)./2)*100;  
    
    GSScorrecttest_pct = (sum(correct_test([pct.highCl.ind.lowind(1,:),pct.highCl.ind.highind(1,:)]),2) - size(correct_test([pct.highCl.ind.lowind(1,:),pct.highCl.ind.highind(1,:)]),2)./2)./...
        (size(correct_test([pct.highCl.ind.lowind(1,:),pct.highCl.ind.highind(1,:)]),2)-size(correct_test([pct.highCl.ind.lowind(1,:),pct.highCl.ind.highind(1,:)]),2)./2)*100;  

    

end
%%
% pred = zeros(9,28,96,144);
% pred_months = zeros(9,28,96,5,144);
% pred2 = zeros(9,28,96,144);
% correct = zeros(9,28,96,144);
% correct2 = zeros(9,28,96,144);
% correct_months = zeros(9,28,96,5,144);


pred = zeros(size(surfacedata,1),size(surfacedata,2),96,144); %THIS NEEDS TO BE CHANGED TO SIZE OF DATASET
pred_months = zeros(size(surfacedata,1),size(surfacedata,2),96,5,144);
pred2 = zeros(size(surfacedata,1),size(surfacedata,2),96,144);
correct = zeros(size(surfacedata,1),size(surfacedata,2),96,144);
correct2 = zeros(size(surfacedata,1),size(surfacedata,2),96,144);
correct_months = zeros(size(surfacedata,1),size(surfacedata,2),96,5,144);

for k = 1:size(surfacedata,1)
    tic;
    allin = 1:size(precondtoz,1);
    
    %tozautocorr = precondtoz(2:end-1,k);
    tozautocorr = precondtoz(:,k);
    upperpctind(k).u = find(tozautocorr >= prctile(tozautocorr,80));
    lowerpctind(k).l = find(tozautocorr <= prctile(tozautocorr,20));
    count = 1;
    for i = 1:length(allin)
        leaveout = allin;        
        %leaveout (leaveout == i) = [];                         
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
        
        %ozone_anomaly = repmat(ozone_anomaly,[1,length(longitude)]);                            
        %ozone_left(k,i) = precondtoz(count,1) - ozone_ind_b(1).*count;% - ozone_ind_b(2);    
        
        % calculating 50th percentile 
        ozone_pct(k,i) = prctile(precondtoz(leaveout,k),50);
        ozone_pct1(k,i) = prctile(precondtoz(leaveout,k),20);
        ozone_pct2(k,i) = prctile(precondtoz(leaveout,k),80);
        
%         ozone_pct(k,i) = prctile(ozone_anomaly,50);
%         ozone_pct1(k,i) = prctile(ozone_anomaly,20);
%         ozone_pct2(k,i) = prctile(ozone_anomaly,80);
        
        %ozone_pct(k,i) = prctile(ozone_anomaly,50);
        
        ozonelowerind(k,i).m = find(precondtoz(leaveout,k) <= ozone_pct1(k,i));
        ozoneupperind(k,i).m = find(precondtoz(leaveout,k) >= ozone_pct2(k,i));                      
        
%         ozonelowerind(k,i).m = find(ozone_anomaly <= ozone_pct(k,i));
%         ozoneupperind(k,i).m = find(ozone_anomaly >= ozone_pct(k,i));                      
        
        ozone_left(k,i) = precondtoz(count,k) - ozone_pct(k,i);
        %ozone_left(k,i) = ozone_left(k,i) - ozone_pct(k,i);
        %ozone_left(k,i) = precondtoz(i,1) - ozone_ind_b(1).*i - ozone_ind_b(2) - ozone_pct(k,i);    
        predsign(k,i) = sign(ozone_left(k,i));        
%         ozone_left = ozone_left - nanmedian(ozone_anomaly);
%         ozone_anomaly = ozone_anomaly - nanmedian(ozone_anomaly);
        
        %ozone_anomalysign = sign(ozone_anomaly);
                
        ENSOyearind = predruns_removeENSOforpred(alldata(:,:,leaveoutenso,:,:),latitude,longitude,0);
        ENSOyearind = ENSOyearind(:,:);
        ENSOall = predruns_removeENSOforpred(alldata,latitude,longitude,0);
        ENSOall = ENSOall(:,:); 
        laglength = 3;        
        
        for j = 1:length(latitude)
            
            % remove ENSO
            if inputs.removeENSO
                for lag = 1:laglength
                    %ensopredictors = [squeeze(ENSOyearind(k,inputs.varmonthtomean-lag+1,:)),ones(size(leaveout))'];
                    ensopredictors = [squeeze(ENSOyearind(k,inputs.varmonthtomean-lag+1:12:end-ext*6))',ones(1,length(leaveout))'];
                    benso(lag,:,:) = ensopredictors\squeeze(surfacedata(k,leaveout,j,:));

                    %ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonthtomean-lag+1,:)),ones(size(ENSOall,3),1)];
                    ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonthtomean-lag+1:12:end-ext*6))',ones(size(alldata,3)-ext,1)];
                    benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedata(k,:,j,:));

                end
                %find max lag
                [~,bensomax_ind] = max(abs(benso),[],1);                
                [~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   
                
                for li = 1:size(longitude)
                    bensomax(li) = benso(bensomax_ind(1,1,li),1,li);
                    bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
                end

                %removing enso
                if i == 2
                    abc = 1;
                end
                for li = 1:size(longitude)
                    temp_anomaly2(:,li) = squeeze(surfacedata(k,leaveout,j,li)) - (squeeze(bensomax(li)).*squeeze(ENSOyearind(k,inputs.varmonthtomean-bensomax_ind(1,1,li)+1)))';
                    temp_left2(li,1) = squeeze(surfacedata(k,count,j,li)) - (squeeze(bensomax_alldata(li)).*squeeze(ENSOall(k,inputs.varmonthtomean-bensomax_ind_alldata(1,1,li)+1))');                   
                end                
                temp_left = temp_left2;
                temp_anomaly = temp_anomaly2;
                clearvars temp_left2 temp_anomaly2
                % empirical model % remove trend                                                                        
    %             b = predictors_anom\squeeze(surfacedata(k,leaveout,j,:));
    %             temp_anomaly = squeeze(surfacedata(k,leaveout,j,:)) - leaveout'.*b(1,:);% - b(2,:); %looks good
    %             temp_left = squeeze(surfacedata(k,count,j,:))' - count.*b(1,:);% - b(2,:); %looks good            

                b = predictors_anom\temp_anomaly;
                temp_anomaly = temp_anomaly - leaveout'.*b(1,:);% - b(2,:); %looks good
                temp_left = temp_left' - count.*b(1,:);% - b(2,:); %looks good      
            else
                b = predictors_anom\squeeze(surfacedata(k,leaveout,j,:));
                temp_anomaly = squeeze(surfacedata(k,leaveout,j,:)) - leaveout'.*b(1,:);% - b(2,:); %looks good
                temp_left = squeeze(surfacedata(k,count,j,:))' - count.*b(1,:);% - b(2,:); %looks good      
            end
%             if j == 69 %lonind = 30
%                 abc = 1;
%             end
            temp_left = squeeze(temp_left) - squeeze(nanmedian(temp_anomaly));
            temp_anomaly = squeeze(temp_anomaly) - squeeze(nanmedian(temp_anomaly));

            temp_anomaly_upper = squeeze(nanmean(temp_anomaly(ozoneupperind(k,i).m,:)));
            temp_anomaly_lower = squeeze(nanmean(temp_anomaly(ozonelowerind(k,i).m,:)));            
                        
%             if j == 84 && i > 1 && i < 30
%                 templow(k,i,:) = temp_anomaly(ozoneupperind(k,i).m,29);
%                 temphigh(k,i,:) = temp_anomaly(ozonelowerind(k,i).m,29);
%             end
            
            signchange = sign(temp_anomaly_upper-temp_anomaly_lower);
            
            datasign(k,i,j,:) = sign(temp_left);
                  
%             if j == 88 && i == 30
%                 abc = 1;
%             end
            
            for l = 1:length(squeeze(datasign(k,i,j,:)))
%                 if signlower(k,i,j,l) < 0 
%                     pred(k,i,j,l) = 1;
%                 elseif signlower(k,i,j,l) > 0 
%                     pred(k,i,j,l) = -1;        
%                 end
                
                if signchange(l) == 1
                    pred(k,i,j,l) = 1;
                else
                    pred(k,i,j,l) = -1;
                end
                
%                 if pred(k,i,j,l) == 1
%                     if predsign(k,i) == datasign(k,i,j,l)
%                         correct(k,i,j,l) = 1;
%                     else
%                         correct(k,i,j,l) = 0;
%                     end
%                 elseif pred(k,i,j,l) == -1
%                     if predsign(k,i) == datasign(k,i,j,l)
%                         correct(k,i,j,l) = 0;
%                     else
%                         correct(k,i,j,l) = 1;
%                     end
%                 end
                
                if pred(k,i,j,l) == 1
                    presign2(k,i,j,l) = predsign(k,i);
                    if predsign(k,i) == datasign(k,i,j,l)
                        correct(k,i,j,l) = 1;
                    else
                        correct(k,i,j,l) = 0;
                    end
                elseif pred(k,i,j,l) == -1
                    presign2(k,i,j,l) = predsign(k,i).*-1;
                    if predsign(k,i) == datasign(k,i,j,l)
                        correct(k,i,j,l) = 0;
                    else
                        correct(k,i,j,l) = 1;
                    end
                end
                
            end  
            
            for m = 1:size(surfacedataind,3)                                
                %remove enso
                % remove ENSO
                for lag = 1:laglength
                    ensopredictors = [squeeze(ENSOyearind(k,inputs.varmonth(m)-lag+1:12:end-ext*6,:))',ones(size(leaveout))'];
                    benso(lag,:,:) = ensopredictors\squeeze(surfacedataind(k,leaveout,m,j,:));

                    ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonth(m)-lag+1:12:end-ext*6,:))',ones(size(alldata,3)-ext,1)];
                    benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedataind(k,:,m,j,:));

                end
                
                %find max lag
                [~,bensomax_ind] = max(abs(benso),[],1);                
                [~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   
                
                for li = 1:size(longitude)
                    bensomax(li) = benso(bensomax_ind(1,1,li),1,li);
                    bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
                end

                %removing enso
                for li = 1:size(longitude)
                    temp_anomaly_months2(:,li) = squeeze(surfacedataind(k,leaveout,m,j,li)) - (squeeze(bensomax(li)).*squeeze(ENSOyearind(k,inputs.varmonth(m)-bensomax_ind(1,1,li)+1)))';
                    temp_left_months2(li,1) = squeeze(surfacedataind(k,count,m,j,li)) - (squeeze(bensomax_alldata(li)).*squeeze(ENSOall(k,inputs.varmonth(m)-bensomax_ind_alldata(1,1,li)+1))');
                end
                temp_left_months = temp_left_months2;
                temp_anomaly_months = temp_anomaly_months2;
                clearvars temp_left_months2 temp_anomaly_months2
                
                
%                 %find max lag
%                 [bensomax,bensomax_ind] = max(abs(benso),[],1);            
%                 [~,bensomax_ind] = max(abs(benso),[],1);            
%                 bensomax = squeeze(benso(bensomax_ind,1,:));
%                 [~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);            
%                 bensomax_alldata = squeeze(benso_alldata(bensomax_ind_alldata,1,:));
%                 
%                 
%                 %removing enso
%                 temp_anomaly_months = squeeze(temp_anomaly_months) - (squeeze(bensomax(1,1,:)).*squeeze(ENSOyearind(k,inputs.varmonth(m)-bensomax_ind(1,1,:)+1,:)))';
%                 temp_left_months = squeeze(temp_left_months) - (squeeze(bensomax_alldata(1,1,:)).*squeeze(ENSOall(k,inputs.varmonth(m)-bensomax_ind_alldata(1,1,:)+1,count))')';
                

                b_months = predictors_anom\temp_anomaly_months;
                temp_anomaly_months = temp_anomaly_months - leaveout'.*b_months(1,:);
                temp_left_months = temp_left_months' - count.*b_months(1,:);% - b(2,:); %looks good        
%                 temp_anomaly_months = squeeze(surfacedataind(k,leaveout,m,j,:)) - leaveout'.*b_months(1,:);% - b(2,:); %looks good
%                 temp_left_months = squeeze(surfacedataind(k,count,m,j,:))' - count.*b_months(1,:);% - b(2,:); %looks good        

                temp_left_months = squeeze(temp_left_months) - squeeze(nanmedian(temp_anomaly_months,1));
                temp_anomaly_months = squeeze(temp_anomaly_months)' - squeeze(nanmedian(temp_anomaly_months,1))';
                temp_anomaly_upper_months = squeeze(nanmean(temp_anomaly_months(:,ozoneupperind(k,i).m),2));
                temp_anomaly_lower_months = squeeze(nanmean(temp_anomaly_months(:,ozonelowerind(k,i).m),2));                        
                
                signchange_months = sign(temp_anomaly_upper_months-temp_anomaly_lower_months);
                
                datasign_months(k,i,j,m,:) = sign(temp_left_months);
                
                for l = 1:length(squeeze(datasign_months(k,i,j,m,:)))
                    if signchange_months(l) == 1 
                        pred_months(k,i,j,m,l) = 1;
                    elseif signchange_months(l) == -1 
                        pred_months(k,i,j,m,l) = -1;        
                    end

                    if pred_months(k,i,j,m,l) == 1
                        predsign_months2(k,i,j,m,l) = predsign(k,i).*1;
                        if predsign(k,i) == datasign_months(k,i,j,m,l)
                            correct_months(k,i,j,m,l) = 1;
                        else
                            correct_months(k,i,j,m,l) = 0;
                        end
                    elseif pred_months(k,i,j,m,l) == -1
                        predsign_months2(k,i,j,m,l) = predsign(k,i).*-1;
                        if predsign(k,i) == datasign_months(k,i,j,m,l)
                            correct_months(k,i,j,m,l) = 0;
                        else
                            correct_months(k,i,j,m,l) = 1;
                        end
                    end
                end  
                
            end
            
        end    
        count = count+1;
    end   
    toc;
end
%% emp

% extracting pct
for i = 1:size(tozdata,1)
%     correctpct(i,:,:,:) = correct(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:);   
%     correct_months_pct(i,:,:,:,:) = correct_months(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:,:);
%         
%     predsignpct(i,:,:,:) = presign2(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:);
%     datasignpct(i,:,:,:) = datasign(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:);
%     predsignpct_months(i,:,:,:,:) = predsign_months2(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:,:);
%     datasignpct_months(i,:,:,:,:) = datasign_months(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:,:);

    correctpct(i,:,:,:) = correct(i,[upperpctind(i).u,lowerpctind(i).l],:,:);   
    correct_months_pct(i,:,:,:,:) = correct_months(i,[upperpctind(i).u,lowerpctind(i).l],:,:,:);
        
    predsignpct(i,:,:,:) = presign2(i,[upperpctind(i).u,lowerpctind(i).l],:,:);
    datasignpct(i,:,:,:) = datasign(i,[upperpctind(i).u,lowerpctind(i).l],:,:);
    predsignpct_months(i,:,:,:,:) = predsign_months2(i,[upperpctind(i).u,lowerpctind(i).l],:,:,:);
    datasignpct_months(i,:,:,:,:) = datasign_months(i,[upperpctind(i).u,lowerpctind(i).l],:,:,:);

end

%%

GSS.all.ind = (sum(correct,2) - size(correct,2)./2)./...
    (size(correct,2)-size(correct,2)./2)*100;  
GSS.pct.ind = (sum(correctpct,2) - size(correctpct,2)./2)./...
    (size(correctpct,2)-size(correctpct,2)./2)*100;  
GSS.all.mean(1,:,:) = squeeze(nanmean(GSS.all.ind));
GSS.pct.mean(1,:,:) = squeeze(nanmean(GSS.pct.ind));

GSS.monthsall.ind = (sum(correct_months,2) - size(correct_months,2)./2)./...
    (size(correct_months,2)-size(correct_months,2)./2)*100;  
GSS.monthspct.ind = (sum(correct_months_pct,2) - size(correct_months_pct,2)./2)./...
    (size(correct_months_pct,2)-size(correct_months_pct,2)./2)*100;  
GSS.monthsall.mean(1,:,:,:) = squeeze(nanmean(GSS.monthsall.ind));
GSS.monthspct.mean(1,:,:,:) = squeeze(nanmean(GSS.monthspct.ind));
  
%% testing bootsrapping
mp = permute(presign2,[3,4,1,2]);
%mp = permute(modelprediction_ind-med,[3,4,1,2]);
mp = mp(:,:,:);
ad = permute(datasign,[3,4,1,2]);
%ad = permute(actualdata_ind-med,[3,4,1,2]);
ad = ad(:,:,:);

mpm = permute(predsign_months2,[4,3,5,1,2]);
mpm = mpm(:,:,:,:);
adm = permute(datasign_months,[4,3,5,1,2]);
adm = adm(:,:,:,:);


if ~exist(['/Volumes/ExternalOne/work/data/predruns/output/HSS/',inputs.ClLevel{1},'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'_empcomp_Allperc.mat'],'file')
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
            bootstat(i,j,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap_compemp(mpe,ade),mpe);
            percentiles.eighty(i,j) = prctile(squeeze(bootstat(i,j,:)),80);
            percentiles.ninty(i,j) = prctile(squeeze(bootstat(i,j,:)),90);
            percentiles.ninetyfive(i,j) = prctile(squeeze(bootstat(i,j,:)),95);
            for k = 1:5
                mpme = squeeze(mpm(k,i,j,:));
                adme = squeeze(adm(k,i,j,:));
                bootstatm(k,i,j,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap_compemp(mpme,adme),mpme);
                percentilesm.eighty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),80);
                percentilesm.ninty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),90);
                percentilesm.ninetyfive(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),95);
            end
        end
        toc;
    end
    save(['/Volumes/ExternalOne/work/data/predruns/output/HSS/',inputs.ClLevel{1},'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'_empcomp_Allperc.mat'],'bootstat','bootstatm','percentiles','percentilesm');
%     allsig.percentiles  = percentiles;
%     allsig.percentiles  = percentilesm;
%     allsig.percentiles  = bootstat;
%     allsig.percentiles  = bootstatm;

end
allsig = load(['/Volumes/ExternalOne/work/data/predruns/output/HSS/',inputs.ClLevel{1},'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'_','empcomp_Allperc.mat']);
p = zeros(size(GSS.all.mean));
p (GSS.all.mean < reshape(allsig.percentiles.ninetyfive,[1,size(allsig.percentiles.ninetyfive)])) = -1;

%%

mp = permute(predsignpct,[3,4,1,2]);
%mp = permute(modelprediction_ind-med,[3,4,1,2]);
mp = mp(:,:,:);
ad = permute(datasignpct,[3,4,1,2]);
%ad = permute(actualdata_ind-med,[3,4,1,2]);
ad = ad(:,:,:);

mpm = permute(predsignpct_months,[4,3,5,1,2]);
mpm = mpm(:,:,:,:);
adm = permute(datasignpct_months,[4,3,5,1,2]);
adm = adm(:,:,:,:);


if ~exist(['/Volumes/ExternalOne/work/data/predruns/output/HSS/',inputs.ClLevel{1},'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'_empcomp_Pctperc.mat'],'file')
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
            bootstat(i,j,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap_compemp(mpe,ade),mpe);
            percentiles.eighty(i,j) = prctile(squeeze(bootstat(i,j,:)),80);
            percentiles.ninty(i,j) = prctile(squeeze(bootstat(i,j,:)),90);
            percentiles.ninetyfive(i,j) = prctile(squeeze(bootstat(i,j,:)),95);
            for k = 1:5
                mpme = squeeze(mpm(k,i,j,:));
                adme = squeeze(adm(k,i,j,:));
                bootstatm(k,i,j,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap_compemp(mpme,adme),mpme);
                percentilesm.eighty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),80);
                percentilesm.ninty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),90);
                percentilesm.ninetyfive(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),95);
            end
        end
        toc;
    end
    save(['/Volumes/ExternalOne/work/data/predruns/output/HSS/',inputs.ClLevel{1},'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'_empcomp_Pctperc.mat'],'bootstat','bootstatm','percentiles','percentilesm');
end
allsigpct = load(['/Volumes/ExternalOne/work/data/predruns/output/HSS/',inputs.ClLevel{1},'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'_empcomp_Pctperc.mat']);
p2 = zeros(size(GSS.all.mean));
p2 (GSS.pct.mean < reshape(allsigpct.percentiles.ninetyfive,[1,size(allsigpct.percentiles.ninetyfive)])) = -1;


%% plot     

if inputs.lats(1) < 0
    ylims = [-90 0];
    diffclims = [-3 3];
else
    ylims = [0 90];
    diffclims = [-5 5];
end

titles = {'Ensemble composite HSSs ','Ensemble composite HSSs during ozone extremes','Observed percentile general skill score (leave one out)'};
toplotp = permute(cat(1,p,p2),[1,3,2]);
[fig,fh] = subplotmaps(permute(cat(1,GSS.all.mean,GSS.pct.mean),[1,3,2]),longitude,latitude,{'seq','YlGnBu'},0,toplotp,22,titles,'Longitude','Latitude','HSS','on',...
        [0,50],11,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[ylims(1):15:ylims(2)],[ylims(1):15:ylims(2)],{''},1,[0 360],ylims,0,'none',1,'Miller Cylindrical');         
 
% subplotmaps(permute(cat(1,GSS.all2.mean,GSS.pct2.mean),[1,3,2]),longitude,latitude,{'seq','YlGnBu'},0,[],16,titles,'Longitude','Latitude','','on',...
%         [0,50],11,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');             
    
child = get(fig,'children');
axes1 = child(4);
axes2 = child(2);

axes(axes1);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','c','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',24,... % 
    'EdgeColor','none','fontweight','bold');    

axes(axes2);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','d','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',24,... % 
    'EdgeColor','none','fontweight','bold');    


set(gcf,'Renderer','Painters');
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/CompEmp_Ensmean_GSSothercolor_nodetrend_upto80_',...
    inputs.ClLevel{1},'_',monthnames(inputs.varmonth,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',inputs.ClLevel{1},'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'.eps'];
print(filename,'-depsc');            
    
%% plot skill lines    
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

if inputs.lats(1) < 0
    xticklab = {'November','December','January','February','March'};
else    
    xticklab = {'March','April','May','June','July'};
end

lstyle = {':','-','-.','--','-'};    
fsize = 16;
Mks = {'d','s','o','v'};

%%
same = 1;
% corrtoplot = 0;
difftoplot = 1;
corrtoplot = 0;
%  areas_lons = [80,180;30,80;240,290;60,120];%lons (Greenland,East Russia, West Russia,America,Asia)
% areas_lats = [60,80;60,80;35,60;25,45];%lats

if inputs.lats(1) < 0
    areas_lons = [115,155;10,40;285,320;280,300]; %lons (East Russia, West Russia,America,Asia)
    areas_lats = [-40,-15;-35,-15;-40,-15;-55,-40]; %lats
else
    areas_lons = [85,180;30,85;240,290;30,120]; %lons (East Russia, West Russia,America,Asia)
    areas_lats = [55,80;55,80;40,65;20,45]; %lats
end

for j = 1:2
    for mon = 1:5

    if same
        sameext = 'same';
    else
        sameext = '';
    end
    if corrtoplot
        corrtoplotext = 'corr';
    elseif difftoplot
        corrtoplotext = 'diff';
    else
        corrtoplotext = 'GSS';
    end
    
        for i = 1:size(areas_lons,1)
            lats = latitude > areas_lats(i,1) & latitude < areas_lats(i,2);
            lons = longitude > areas_lons(i,1) & longitude < areas_lons(i,2);
            latextract = latitude(lats);
            lonextract = longitude(lons);
            [latmesh,lonmesh] = meshgrid(latextract,lonextract);

            latmesh = latmesh';
            lonmesh = lonmesh';
            for k = 1:size(differences,1)           

                if corrtoplot               
                elseif difftoplot            
                    diff = permute(differences,[1,4,3,2]);           
                    diff2 = permute(diffcomp,[1,3,2]);
                    if same                
                        mult = squeeze(nanmean(diff(:,mon,lats,lons),1));        
                        %mult = squeeze(nanmean(diff2(:,lats,lons),1));        
                        mult_pct = squeeze(nanmean(diff(:,mon,lats,lons),1));        
                    else
                        mult = squeeze(diff(k,mon,lats,lons));        
                        mult_pct = squeeze(diff(k,mon,lats,lons));        
                    end
                else                
                end    
                if i == 1 && inputs.lats(1) < 0
                    [maxval(i,k),maxind] = max(mult(:));
                    [maxval_pct(i,k),maxind_pct] = max(abs(mult_pct(:)));        
                else
                    [maxval(i,k),maxind] = max(abs(mult(:)));
                    [maxval_pct(i,k),maxind_pct] = max(abs(mult_pct(:)));        
                end                
                lattoplot(j).g(mon,k,i) = latmesh(maxind);
                lontoplot(j).g(mon,k,i) = lonmesh(maxind);

                lattoplot_pct(j).g(mon,k,i) = latmesh(maxind_pct);
                lontoplot_pct(j).g(mon,k,i) = lonmesh(maxind_pct);

                if k == size(differences,1) && same
                    forerror = squeeze(allsig.percentilesm.ninetyfive(:,lats,lons));
                    forerrorextract(j).g(mon,:,i) = forerror(:,maxind);

                    forerrorpct = squeeze(allsigpct.percentilesm.ninetyfive(:,lats,lons));
                    forerrorextractpct(j).g(mon,:,i) = forerrorpct(:,maxind);  

                    modelmonthsforerror_extract = [];
                    datamonthsforerror_extract = [];
                    %medmonthsforerror_extract = [];
                    modelmonthsforerror_pct_extract = [];
                    datamonthsforerror_pct_extract = [];
                    %medmonthsforerror_pct_extract = [];
                elseif ~same

                    %modelmonthsforerror = squeeze(modelprediction_ind_months(k,:,:,lats,lons));
                    modelmonthsforerror = permute(squeeze(predsign_months2(k,:,lats,:,lons)),[1,3,2,4]);
                    modelmonthsforerror_extract(i,:,k,:) = modelmonthsforerror(:,:,maxind)';
                    %datamonthsforerror = squeeze(actualdata_ind_months(k,:,:,lats,lons));
                    datamonthsforerror = permute(squeeze(datasign_months(k,:,lats,:,lons)),[1,3,2,4]);
                    datamonthsforerror_extract(i,:,k,:) = datamonthsforerror(:,:,maxind)';
                    %medmonthsforerror = squeeze(medmonths(k,:,:,lats,lons));
                    %medmonthsforerror_extract(i,:,k,:) = medmonthsforerror(:,:,maxind)';

                    modelmonthsforerror_pct = permute(squeeze(predsignpct_months(k,:,lats,:,lons)),[1,3,2,4]);
                    modelmonthsforerror_pct_extract(i,:,k,:) = modelmonthsforerror_pct(:,:,maxind)';
                    datamonthsforerror_pct = permute(squeeze(datasignpct_months(k,:,lats,:,lons)),[1,3,2,4]);
                    datamonthsforerror_pct_extract(i,:,k,:) = datamonthsforerror_pct(:,:,maxind)';
                    %medmonthsforerror_pct = squeeze(medmonths_pct2(k,:,:,lats,lons));
                    %medmonthsforerror_pct_extract(i,:,k,:) = medmonthsforerror_pct(:,:,maxind)';

                end

                %[maxrow(k,i),maxcol(k,i)] = find(mult == max(mult(:)));
                GSSmonthtoplot2 = permute(squeeze(GSS.monthsall.ind(k,1,lats,:,lons)),[2,1,3]);
                GSSmonthtoplot(:,k,i) = GSSmonthtoplot2(:,maxind);
                GSSmonthtoplot3 = permute(squeeze(GSS.monthspct.ind(k,1,lats,:,lons)),[2,1,3]);
                GSSmonthtoplot4(:,k,i) = GSSmonthtoplot3(:,maxind);
                if maxind <= size(mult,1)
                    if maxind ~= 1
                        GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                            GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                                

                        GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                            GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                                

                    else
                        GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),...
                            GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                            

                        GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),...
                            GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                            
                    end
                elseif maxind == size(mult,1)+1
                    GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                        GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
                        GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        

                    GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                        GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
                        GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        

                elseif maxind == numel(mult) - size(mult,1)
                    GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                        GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)),...
                        GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        

                    GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                        GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)),...
                        GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                    

                elseif maxind > numel(mult) - size(mult,1)
                    if maxind ~= numel(mult)
                        GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                            GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                

                        GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                            GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                

                    else                
                        GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind-1),...
                            GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                

                        GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind-1),...
                            GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                
                    end
                else
                    GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                        GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
                        GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        

                    GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                        GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
                        GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        

                end

                %pct
                if maxind_pct <= size(mult_pct,1)
                    if maxind_pct ~= 1               
                        GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                            GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
                    else                

                        GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),...
                            GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
                    end
                elseif maxind_pct == size(mult_pct,1)+1            

                    GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                        GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
                        GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);

                elseif maxind_pct == numel(mult_pct) - size(mult_pct,1)

                    GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                        GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)),...
                        GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);

                elseif maxind_pct > numel(mult_pct) - size(mult_pct,1)
                    if maxind_pct ~= numel(mult_pct)                

                        GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                            GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);

                    else                

                        GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct-1),...
                            GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
                    end
                else

                    GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                        GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
                        GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);            
                end                
            end
        end
    
    

        %% constructing bootstrap error
        if ~same
            mp = modelmonthsforerror_extract(:,:,:);% - medmonthsforerror_extract(:,:,:);
            ap = datamonthsforerror_extract(:,:,:); %- medmonthsforerror_extract(:,:,:);
            mpp = modelmonthsforerror_pct_extract(:,:,:);% - medmonthsforerror_pct_extract(:,:,:);
            app = datamonthsforerror_pct_extract(:,:,:);% - medmonthsforerror_pct_extract(:,:,:);
            for x = 1:size(mp,1)       

                for i = 1:size(mp,2)

                    mpe = squeeze(mp(x,i,:));
                    ade = squeeze(ap(x,i,:));

                    mpme = squeeze(mpp(x,i,:));
                    adme = squeeze(app(x,i,:));


                    bootstatspec(x,i,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap_compemp(mpe,ade),mpe);            
                    forerrorextract(2).g(mon,i,x) = prctile(squeeze(bootstatspec(x,i,:)),95);

                    bootstatmspec(x,i,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap_compemp(mpme,adme),mpme);                
                    forerrorextractpct(2).g(mon,i,x) = prctile(squeeze(bootstatmspec(x,i,:)),95);

                end
            end
        end

    end
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

titles = {'Ensemble composite','Ensemble composite, ozone extremes','Ensemble members','Ensemble members, ozone extremes'}; 
labels = {'c','d','a','b'};
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
            phltemp(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'MarkerSize',8,'MarkerEdgeColor',cbrewqual3(i,:),'MarkerFaceColor',cbrewqual2(i,:),'linewidth',2,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
            phl(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'MarkerSize',8,'MarkerEdgeColor',cbrewqual3(i,:),'MarkerFaceColor',cbrewqual2(i,:),'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));

            phl2(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},...
                'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'MarkerSize',10);
        end
        uistack(eh1,'bottom')
        uistack(phl2,'top')
        plot([0,6],[0,0],'linewidth',2,'LineStyle','--','color','k');
        set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
        if count == 3
            xlabel('Month','fontsize',fsize+2);
        end
        xlim([.5,5.5]);
        ylim([-30 80]);
        set(gca,'ytick',-20:10:100,'yticklabel',-20:10:100);
        ylabel('HSS','fontsize',fsize+2);
        title(titles{count},'fontsize',fsize+4)
        %March TCO - surface temperature ensemble mean HSS
        if count == 1
            lh = legend(phltemp,'Eastern Russia','Western Russia','Northern America','Southern Asia');                    
            set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
        end
        annotation('textbox',sppos,'String',labels{count},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+2,... % 
            'EdgeColor','none','fontweight','bold');    
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
        plot([0,6],[0,0],'linewidth',2,'LineStyle','--','color','k');
        set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
        if count == 4
            xlabel('Month','fontsize',fsize+2);
        end
         annotation('textbox',sppos,'String',labels{count},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+2,... % 
            'EdgeColor','none','fontweight','bold');    
        %ylabel('HSS','fontsize',fsize+2);
        title(titles{count},'fontsize',fsize+4)
        xlim([.5,5.5]);
        ylim([-30 80]);
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
m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',ylims);
[~,h] = m_contourf(longitude-180,latitude,squeeze(nanmean(diffforplot(:,:,:,mon),1))',-5:.5:5,'LineStyle','none');
hold on
m_coast('color','k','LineWidth',1);
m_grid('ytick',[ylims(1):15:ylims(2)],'xtick',[-180:60:360],'XaxisLocation','bottom','fontsize',fsize);  
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
sppos = get(sp2,'position');
annotation('textbox',[sppos(1),sppos(2)-.04,sppos(3:4)],'String','c','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+2,... % 
            'EdgeColor','none','fontweight','bold');    
cbarrow
if corrtoplot
    title(['Location of maximum ensemble ', monthnames(inputs.tozmonth,0,0),' TCO and ',xticklab{mon}, ' TS differences'],'fontsize',fsize+2);
elseif difftoplot
    title(['Locations of ensemble ' monthnames(inputs.tozmonth,0,0) ' TCO and ',xticklab{mon}, ' surface temperature differences'],'fontsize',fsize+2);
    
else
    title(['Location of maximum ensemble member ', monthnames(inputs.tozmonth,0,0),' TCO and ',xticklab{mon}, ' HSS'],'fontsize',fsize+2);
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

annotation('textbox',[.01 .995 .925 0],'String',[monthnames(inputs.tozmonth,0,0),' TCO - surface temperature ensemble composite HSSs'],'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');    

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/CompEmp_Ind_',...
    inputs.ClLevel{1},'_',corrtoplotext,'_','mean_Locations_withpct_test','_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),monthnames(mon+2,0,0),'_',num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'both'];
print(filename,'-depsc');
end
end
