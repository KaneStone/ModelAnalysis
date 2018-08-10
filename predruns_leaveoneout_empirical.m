function [] = predruns_leaveoneout_empirical(surfacedata,surfacedataind,tozdata,inputs,latitude,longitude,pct,correlations)

% over 1995-2014. Leave one year out of the regression and see if I can
% predict the remaining years anomaly

% Do this for the ensemble mean and produce a map metric of predictability
% for the ensemble average.

% Using the high signal to noise plots, extract regions near these places
% for each 

% keep leaving years out

% Take regression of ensemble

% DON'T predefine anything!!!

%% rearrange data into composite

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
testing = 1;

if testing
    %sd1 = squeeze(surfacedata(1,:,79,31));
    %sd1 = squeeze(surfacedata(1,:,87,22));
    %sd1 = squeeze(surfacedata(1,:,80,40));
    %sd1 = squeeze(surfacedata(1,:,87,30));
    %sd1 = double(squeeze(surfacedata(1,:,79,35)));
    %sd1 = squeeze(surfacedata(1,:,77,44));
    %sd1 = squeeze(surfacedata(1,:,90,60));
    %sd1 = squeeze(surfacedata(1,:,85,45));           
    sd1 = squeeze(surfacedata(1,:,14,105));           
    %sd1 = squeeze(surfacedata(1,:,80,79));
    % leaveout
    
    for i = 1:length(sd1)
        
        sd12 = 1:length(sd1);
        sd12 (sd12 == i) = [];
        
        predictors_anom = [sd12',ones(size(sd12))'];
        %ozone_ind_b = precondtoz(leaveout,k)\predictors_anom;        
        ozone_ind_b = predictors_anom\precondtoz(sd12,1);        
        ozone_anomaly_test(i,:) = precondtoz(sd12,1) - ozone_ind_b(1).*sd12';% - ozone_ind_b(2);        
        ozone_left_test(i) = precondtoz(i,1) - ozone_ind_b(1).*i;% - ozone_ind_b(2);    
        
%         ozone_left_test(i) = ozone_left_test(i) - nanmedian(ozone_anomaly_test(i,:));
%         ozone_anomaly_test(i,:) = ozone_anomaly_test(i,:) - nanmedian(ozone_anomaly_test(i,:));
        
        ozone_left_test(i) = ozone_left_test(i) - nanmean(ozone_anomaly_test(i,:));
        ozone_anomaly_test(i,:) = ozone_anomaly_test(i,:) - nanmean(ozone_anomaly_test(i,:));
        
        
%         ozone_anomaly_test = precondtoz(sd12,1);%-nanmean(precondtoz(sd12,1),'double');    
%         b = regress(sd1(sd12)',[ones(size(sd1(sd12)))',[1:length(sd1(sd12))]']);
                        
        %temp_anomaly = detrend(sd1(sd12));
        b = predictors_anom\sd1(sd12)';
        temp_anomaly_test(i,:) = sd1(sd12) - b(1).*sd12;% - b(2);
        temp_left_test(i) = sd1(i) - b(1).*i;% - b(2);
        
        temp_left_test(i) = temp_left_test(i) - nanmedian(temp_anomaly_test(i,:));
        temp_anomaly_test(i,:) = temp_anomaly_test(i,:) - nanmedian(temp_anomaly_test(i,:));
        
        r(i) = corr(temp_anomaly_test(i,:)',ozone_anomaly_test(i,:)');
        %temp_anomaly_test = sd1(sd12) - nanmean(sd1(sd12));
        ozone_anomalysign = sign(ozone_anomaly_test(i,:));
        temp_anomalysign = sign(temp_anomaly_test(i,:));        

        bothpositive = sum((ozone_anomalysign > 0 & temp_anomalysign > 0));
        bothnegative = sum((ozone_anomalysign < 0 & temp_anomalysign < 0));
        bothdifferent = sum((ozone_anomalysign ~= temp_anomalysign));
        bothsame = bothpositive+bothnegative;
        
        if bothsame > bothdifferent
            prediction_test(i) = 1;
        elseif bothdifferent > bothsame
            prediction_test(i) = -1; 
        else
            prediction_test(i) = 0;
        end
                
        sign_temp_left(i) = sign(temp_left_test(i));
        sign_ozone_left(i) = sign(ozone_left_test(i));
        
        if prediction_test(i) == 1 % ozone anomaly = temperature anomaly
            if sign_ozone_left(i) == sign_temp_left(i)
                correct_test(i) = 1;
            else
                correct_test(i) = 0;
            end
        elseif prediction_test(i) == -1
            if sign_ozone_left(i) ~= sign_temp_left(i)
                correct_test(i) = 1;
            else
                correct_test(i) = 0;
            end
        elseif prediction_test(i) == 0
                correct_test(i) = .5;
        end

        i = i+1;
    end
    
%     acdata2 = detrend(acdata);
%     modelpredloo2 = detrend(modelpredloo);
%     figure
%     hold on
%     plot(acdata)
%     plot(acdata2)
%     plot(modelpredloo2,'k--')
%     plot(modelpredloo,'k')
%     figure    
%     plot(r)

figure;
plot(temp_left_test)
hold on
plot(ozone_left_test)
% plot(squeeze(modelprediction_ind(1,:,80,30)),'k')
% hold on
% plot(squeeze(sd1(1,:,80,30)),'k--')


    GSScorrecttest = (sum(correct_test,2) - size(correct_test,2)./2)./...
        (size(correct_test,2)-size(correct_test,2)./2)*100;  

    

end
%%

for k = 1:size(surfacedata,1)
    
    allin = 1:size(precondtoz,1);
    
    for i = 1:length(allin)
        leaveout = allin;        
        leaveout (leaveout == i) = [];                         
        
        predictors_anom = [leaveout',ones(size(leaveout))'];        
        ozone_ind_b = predictors_anom\precondtoz(leaveout,k);        
        ozone_anomaly = precondtoz(leaveout,k) - ozone_ind_b(1).*leaveout' - ozone_ind_b(2); 
        
        ozone_anomaly = repmat(ozone_anomaly,[1,length(longitude)]);                            
        ozone_left(k,i) = precondtoz(i,1) - ozone_ind_b(1).*i - ozone_ind_b(2);    
        
%         ozone_left = ozone_left - nanmedian(ozone_anomaly);
%         ozone_anomaly = ozone_anomaly - nanmedian(ozone_anomaly);
        
        ozone_anomalysign = sign(ozone_anomaly);
                
        for j = 1:length(latitude)                           
            % empirical model
                                                                        
            b = predictors_anom\squeeze(surfacedata(k,leaveout,j,:));
            temp_anomaly = squeeze(surfacedata(k,leaveout,j,:)) - leaveout'.*b(1,:) - b(2,:); %looks good
            temp_left = squeeze(surfacedata(k,i,j,:))' - i.*b(1,:) - b(2,:); %looks good            
            
%             temp_left = temp_left - nanmedian(temp_anomaly);
%             temp_anomaly = temp_anomaly - nanmedian(temp_anomaly);
            
            
            med(k,i,j,:) = median(temp_anomaly);
            %ozonemed(k,i,j) = median(ozone_anomaly);
            %temp_anomaly = detrend(squeeze(surfacedata(k,leaveout,j,:)));
            
            temp_anomalysign = sign(temp_anomaly);
%             for m = 1:size(temp_anomaly,2)
%                 r(m) = corr(temp_anomaly(:,m),ozone_anomaly(:,m));                       
%             end
            
            bothpositive = sum(ozone_anomalysign > 0 & temp_anomalysign > 0);
            bothnegative = sum(ozone_anomalysign < 0 & temp_anomalysign < 0);
            bothdifferent = sum(ozone_anomalysign ~= temp_anomalysign);
            bothsame = bothpositive+bothnegative;
            
            prediction = zeros(size(bothsame));
            %prediction2 = zeros(size(bothsame));
            prediction (bothsame > bothdifferent) = 1;
            prediction (bothdifferent > bothsame) = -1;
            
%             prediction2 (r < 0) = -1;
%             prediction2 (r > 0) = 1;
            
                                
            for l = 1:length(prediction)                                
                if prediction(l) == 1 
                    if sign(ozone_left(k,i)) == sign(temp_left(l))
                        correct(k,i,j,l) = 1;
                    else
                        correct(k,i,j,l) = 0;
                    end
                elseif prediction(l) == -1
                    if sign(ozone_left(k,i)) ~= sign(temp_left(l))
                        correct(k,i,j,l) = 1;
                    else
                        correct(k,i,j,l) = 0;
                    end
                elseif prediction(l) == 0
                    correct(k,i,j,l) = .5;
                end
            end       
            
%             for l = 1:length(prediction2)                                
%                 if prediction2(l) == 1 
%                     if sign(ozone_left) == sign(temp_left(l))
%                         correct2(k,i,j,l) = 1;
%                     else
%                         correct2(k,i,j,l) = 0;
%                     end
%                 elseif prediction2(l) == -1
%                     if sign(ozone_left) ~= sign(temp_left(l))
%                         correct2(k,i,j,l) = 1;
%                     else
%                         correct2(k,i,j,l) = 0;
%                     end
%                 end
%             end       
            
            
        end        
    end
    
       
%     allin = 1:size(precondtoz_pct,1);    
%      for i = 1:length(allin)
%         leaveout = allin;        
%         leaveout (leaveout == i) = [];         
%         ozoneleft = precondtoz_pct(i,k);
%         for j = 1:length(latitude)                           
%             %% empirical model
%             ozone_anomaly = precondtoz_pct(leaveout,k)-nanmean(precondtoz_pct(leaveout,k),'double');    
%             b = [ones(size(ozone_anomaly)),[1:length(ozone_anomaly)]']\squeeze(surfacedata(k,leaveout,j,:));%regress(sd1(sd12)',[ones(size(sd1(sd12)))',[1:length(sd1(sd12))]']);            
%             temp_anomaly = squeeze(surfacedata_pct(k,leaveout,j,:)) - [1:length(ozone_anomaly)]'.*b(2,:) - b(1,:); %looks good
%             %temp_anomaly = detrend(squeeze(surfacedata(k,leaveout,j,:)));
%             ozone_anomaly = repmat(ozone_anomaly,[1,size(temp_anomaly,2)]);            
%             ozone_anomalysign = sign(ozone_anomaly);
%             temp_anomalysign = sign(temp_anomaly);
%             if j == 77
%                 abc = 1;
%             end
%             bothpositive = sum(ozone_anomalysign > 0 & temp_anomalysign > 0);
%             bothnegative = sum(ozone_anomalysign < 0 & temp_anomalysign < 0);
%             bothdifferent = sum(ozone_anomalysign ~= temp_anomalysign);
%             bothsame = bothpositive+bothnegative;
%             prediction = zeros(size(bothsame));
%             prediction (bothsame-3 > bothdifferent) = 1;
%             prediction (bothdifferent-3 > bothsame) = -1;
% %             if bothsame >= bothdifferent
% %                 prediction = 1;
% %             else
% %                 prediction = -1;
% %             end
%             %actualdataall(k,i,j,:) = squeeze(surfacedata_pct(k,i,j,:))' - b(2,:)*i - b(1,:);
%             actualdata = sign(squeeze(surfacedata_pct(k,i,j,:))' - b(2,:)*i - b(1,:));
%             %predictiondataall(k,i,j,:) = repmat(ozoneleft - nanmean(precondtoz_pct(leaveout,k),'double'),[1,length(actualdata)]);
%             predictiondata = repmat(sign(ozoneleft - nanmean(precondtoz_pct(leaveout,k))),[1,length(actualdata)]);
%             %modeldata = sign(precondtoz(i,1) - nanmean(precondtoz(sd12,1)));
%             for l = 1:length(actualdata)
%                 if j == 77 && l == 44
%                     abc = 1;
%                 end
%                 
%                 if prediction(l) == 1 % ozone anomaly = temperature anomaly
%                     if predictiondata(l) == actualdata(l)
%                         correct_pct(k,i,j,l) = 1;
%                     else
%                         correct_pct(k,i,j,l) = 0;
%                     end
%                 elseif prediction(l) == -1
%                     if predictiondata(l) ~= actualdata(l)
%                         correct_pct(k,i,j,l) = 1;
%                     else
%                         correct_pct(k,i,j,l) = 0;
%                     end
%                 elseif prediction(l) == 0
%                     correct_pct(k,i,j,l) = .5;
%                 end
%             end                           
%         end        
%      end        
end
%% emp

% extracting pct
for i = 1:9
    correctpct(i,:,:,:) = correct(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:);
    %correctpct2(i,:,:,:) = correct2(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:);
end

%%

GSS.emp = (sum(correct,2) - size(correct,2)./2)./...
    (size(correct,2)-size(correct,2)./2)*100;  
GSS.emp_pct = (sum(correctpct,2) - size(correctpct,2)./2)./...
    (size(correctpct,2)-size(correctpct,2)./2)*100;  
GSS.empmean(1,:,:) = squeeze(nanmean(GSS.emp));
GSS.empmean_pct(1,:,:) = squeeze(nanmean(GSS.emp_pct));

% GSS.emp2 = (sum(correct2,2) - size(correct2,2)./2)./...
%     (size(correct2,2)-size(correct2,2)./2)*100;  
% GSS.emp_pct2 = (sum(correctpct2,2) - size(correctpct2,2)./2)./...
%     (size(correctpct2,2)-size(correctpct2,2)./2)*100;  
% GSS.empmean2(1,:,:) = squeeze(nanmean(GSS.emp2));
% GSS.empmean_pct2(1,:,:) = squeeze(nanmean(GSS.emp_pct2));

  
%% plot     
titles = {'general skill score (leave three out) mean','Percentile general skill score (leave one out) mean','Observed percentile general skill score (leave one out)'};
subplotmaps(permute(cat(1,GSS.empmean,GSS.empmean_pct),[1,3,2]),longitude,latitude,{'seq','YlGnBu'},0,[],16,titles,'Longitude','Latitude','','on',...
        [0,50],11,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
    
% subplotmaps(permute(cat(1,GSS.empmean2,GSS.empmean_pct2),[1,3,2]),longitude,latitude,{'seq','YlGnBu'},0,[],16,titles,'Longitude','Latitude','','on',...
%     [0,50],11,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
    
    %subplotmaps(permute(squeeze(cat(1,GSS.emp_pct(2,1,:,:),GSS.emp_pct(4,1,:,:),GSS.emp_pct(5,1,:,:))),[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,titles,'Longitude','Latitude','','on',...
    %    [-60,60],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
    
%%    
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/Ensmean_GSS_',monthnames(inputs.varmonth,1,1)];
    export_fig(filename,'-png');        
    
subplotmaps(permute(cat(1,Heidke.all.mean,Heidke.pct.mean),[1,3,2]),longitude,latitude,{'div','RdBu'},1,[],16,{'Heidke skill score (leave three out) mean','Percentile Heidke skill score (leave one out) mean'},'Longitude','Latitude','','on',...
        [-1,1],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');
    
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/Ensmean_HSS_',monthnames(inputs.varmonth,1,1)];
    export_fig(filename,'-png');            
    
%% plot skill lines    
    
end