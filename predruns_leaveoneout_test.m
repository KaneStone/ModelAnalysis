function [] = predruns_leaveoneout_test(surfacedata,surfacedataind,tozdata,inputs,latitude,longitude,pct)


% over 1995-2014. Leave one year out of the regression and see if I can
% predict the remaining years anomaly

% Do this for the ensemble mean and produce a map metric of predictability
% for the ensemble average.

% Using the high signal to noise plots, extract regions near these places
% for each 

% keep leaving years out

% Take regression of ensemble

% DON'T predefine anything!!!

surfacedata_composite = permute(surfacedata,[3,4,2,1]);
surfacedata_composite = surfacedata_composite(:,:,:); 

for i = 1:9
    precondtoz(:,i) = squeeze(tozdata(i,inputs.tozmonth,:));%-nanmean(squeeze(tozdata(i,inputs.tozmonth,:)));
    precondtoz2(:,i) = squeeze(tozdata(i,inputs.tozmonth,:));%-nanmean(squeeze(tozdata(i,inputs.tozmonth,:)));
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
    %sd1 = squeeze(surfacedata(1,:,78,90));
    %sd1 = squeeze(surfacedata(1,:,78,91));
    %sd1 = squeeze(surfacedata(1,:,87,1));
    %sd1 = squeeze(surfacedata(1,:,87,22));
    sd1 = squeeze(surfacedata(1,:,87,30));
    %sd1 = double(squeeze(surfacedata(1,:,79,35)));
    %sd1 = squeeze(surfacedata(1,:,77,43));
    %sd1 = squeeze(surfacedata(1,:,74,39));
    %sd1 = squeeze(surfacedata(1,:,80,40));
    b = regress(sd1',[precondtoz(:,1),ones(size(precondtoz(:,1)))]);
    modelpred = b(1)*precondtoz(:,1)+b(2);

    % leaveout
    for i = 1:length(sd1)
        sd12 = 1:length(sd1);
        sd12 (sd12 == i) = [];
        ozonedetrended = regress(precondtoz(sd12,1),[ones(size(sd12));sd12]');
        %ozoneind = precondtoz(sd12,1)-nanmean(precondtoz(sd12,1));        
        ozoneind = precondtoz(sd12,1) - ozonedetrended(2).*sd12' - ozonedetrended(1);
        
        %ozoneind_left = precondtoz(i,1)-nanmean(precondtoz(sd12,1));
        ozoneind_left = precondtoz(i,1) - ozonedetrended(2).*i' - ozonedetrended(1);
        
        b = regress(sd1(sd12)',[ones(size(sd1(sd12)))',sd12']);        
        temp_loo_detrended(i,:) = sd1(sd12) - b(2)*sd12-b(1);        
        med(i) = nanmedian(temp_loo_detrended(i,:));
        temp_loom_detrended(i,:) = temp_loo_detrended(i,:) - med(i);
        nomedtest(i) = sum(sign(temp_loo_detrended(i,:)));
        medtest(i) = sum(sign(temp_loom_detrended(i,:)));
        
        bloo_med(i,:) = regress(temp_loom_detrended(i,:)',[ones(size(ozoneind)),ozoneind]);
        bloo_nomed(i,:) = regress(temp_loo_detrended(i,:)',[ones(size(ozoneind)),ozoneind]);
        modelpredloo_nomed(i) = bloo_nomed(i,2)*ozoneind_left+bloo_nomed(i,1);%+nanmedian(sd123(i,:));
        modelpredloo_med(i) = bloo_med(i,2)*ozoneind_left+bloo_med(i,1);
        r_nomed(i) = corr(ozoneind,temp_loo_detrended(i,:)');
        r_med(i) = corr(ozoneind,temp_loom_detrended(i,:)');
        %modelpredloo(i) = (bloo(i,1)*ozoneind_left);
        acdata(i) = sd1(i)-nanmean(sd1(sd12));
        acdata_nomed(i) = sd1(i) - b(2)*i - b(1);
        acdata_med(i) = sd1(i) - b(2)*i-med(i) - b(1);
    end
    
    createfig('large','on')
    
    subplot(2,2,1)
    plot(r_nomed,'LineWidth',2)
    hold on
    plot(r_med,'LineWidth',2)    
    title('Correlations (58N, 220E)','fontsize',22); %
    set(gca,'fontsize',20);
    ylabel('r','fontsize',20);
    subplot(2,2,2)
    plot(bloo_nomed(:,2),'-o','LineWidth',2)
    hold on
    plot(bloo_med(:,2),'LineWidth',2)
    
    
    title('Ozone beta functions (58N, 220E)','fontsize',22);
    set(gca,'fontsize',20);
    ylabel('r','fontsize',20);
    subplot(2,2,3)
    plot((modelpred-nanmean(modelpred))*10,'k','LineWidth',2);
    hold on
    plot(modelpredloo_nomed*2,'k--','LineWidth',2)
    plot(acdata_nomed-med,'-o','LineWidth',2)
    lh = legend('Ozone time series (HSS = 0)','Leave-one-out model predictions','Surface temperature data');
    set(lh,'fontsize',14,'box','off','location','NorthWest');
    ylabel('normalized DU or K','fontsize',20);
    xlabel('year','fontsize',20);
    subplot(2,2,4)
    plot(modelpredloo_med*2,'k--','LineWidth',2)
    hold on
    plot(modelpredloo_nomed*2,'k','LineWidth',2)
    
    plot(acdata_med,'-o','LineWidth',2)
    title('Scaled predictions (58N, 220E)','fontsize',22);
    ylabel('normalized DU or K','fontsize',20);
    xlabel('year','fontsize',20);
    
%    plot(acdata4,'-o','LineWidth',2)
    
    set(gca,'fontsize',20);
    
    %export_fig('/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Testing58N220E.pdf','-pdf');
end

% without each climatology normalized to the median
modelpredictionign_ind = sign(modelpredloo_nomed);
actualdatasign_ind =  sign(acdata_nomed);

actualdatasign_ind3 = acdata_nomed; 
actualdatasign_ind3 (actualdatasign_ind3 > med) = 1;
actualdatasign_ind3 (actualdatasign_ind3 < med) = -1;

test2 = modelpredloo_nomed;
test2 (test2 > med) = 1;
test2 (test2 < med) = -1;

modelpredictionresults_ind = modelpredictionign_ind;
modelpredictionresults_ind (modelpredictionresults_ind ~= actualdatasign_ind) = -2;
modelpredictionresults_ind (modelpredictionresults_ind == actualdatasign_ind) = 1;
modelpredictionresults_ind (modelpredictionresults_ind == -2) = 0;

GSStest = (sum(modelpredictionresults_ind,2) - size(modelpredictionresults_ind,2)./2)./...
    (size(modelpredictionresults_ind,2)-size(modelpredictionresults_ind,2)./2)*100;  

% with each climatology normalized to the median
modelpredictionign_ind2 = sign(modelpredloo_nomed-med);
actualdatasign_ind2 =  sign(acdata_med);

modelpredictionresults_ind2 = modelpredictionign_ind2;
modelpredictionresults_ind2 (modelpredictionresults_ind2 ~= actualdatasign_ind2) = -2;
modelpredictionresults_ind2 (modelpredictionresults_ind2 == actualdatasign_ind2) = 1;
modelpredictionresults_ind2 (modelpredictionresults_ind2 == -2) = 0;

GSStest2 = (sum(modelpredictionresults_ind2,2) - size(modelpredictionresults_ind2,2)./2)./...
    (size(modelpredictionresults_ind2,2)-size(modelpredictionresults_ind2,2)./2)*100;  

% with only TS normalized to the median
modelpredictionign_ind3 = sign(modelpredloo_nomed);
actualdatasign_ind3 =  sign(acdata_med);

modelpredictionresults_ind3 = modelpredictionign_ind3;
modelpredictionresults_ind3 (modelpredictionresults_ind3 ~= actualdatasign_ind3) = -2;
modelpredictionresults_ind3 (modelpredictionresults_ind3 == actualdatasign_ind3) = 1;
modelpredictionresults_ind3 (modelpredictionresults_ind3 == -2) = 0;

GSStest3 = (sum(modelpredictionresults_ind3,2) - size(modelpredictionresults_ind3,2)./2)./...
    (size(modelpredictionresults_ind3,2)-size(modelpredictionresults_ind3,2)./2)*100;  

% all data
modelpredictionign_ind4 = sign(modelpred'-nanmean(modelpred));
acdata2 = sign(detrend(sd1));

modelpredictionresults_ind4 = modelpredictionign_ind4;
modelpredictionresults_ind4 (modelpredictionresults_ind4 ~= acdata2) = -2;
modelpredictionresults_ind4 (modelpredictionresults_ind4 == acdata2) = 1;
modelpredictionresults_ind4 (modelpredictionresults_ind4 == -2) = 0;

GSStest4 = (sum(modelpredictionresults_ind4,2) - size(modelpredictionresults_ind4,2)./2)./...
    (size(modelpredictionresults_ind4,2)-size(modelpredictionresults_ind4,2)./2)*100;  

% actualdatasign_indpercent = sum(actualdatasign_ind)./length(actualdatasign_ind);
% actualdatasign_indpercent2 = sum(actualdatasign_ind2)./length(actualdatasign_ind2);

HA_ind = sum(modelpredictionign_ind >= 0 & actualdatasign_ind >= 0,2); % if predicted correctly (hit)
HB_ind = sum(modelpredictionign_ind > actualdatasign_ind,2); % if predictly incorrectly (miss)
HC_ind = sum(modelpredictionign_ind < actualdatasign_ind,2); % if not predicted but occurred
HD_ind = sum(modelpredictionign_ind <= 0 & actualdatasign_ind <= 0,2); % if did not predict and did not occur

HSStest = 2.*(HA_ind.*HD_ind - HB_ind.*HC_ind)./(((HA_ind+HC_ind).*(HC_ind+HD_ind))+((HA_ind+HB_ind).*(HB_ind+HD_ind)))*100; %2(ad-bc)/(a+c)(c+d)+(a+b)(b+d)

HA_ind3 = sum(modelpredictionign_ind3 >= 0 & acdata_med >= 0,2);
HB_ind3 = sum(modelpredictionign_ind3 > acdata_med,2);
HC_ind3 = sum(modelpredictionign_ind3 < acdata_med,2);
HD_ind3 = sum(modelpredictionign_ind3 <= 0 & acdata_med <= 0,2);

HSStest3 = 2.*(HA_ind3.*HD_ind3 - HB_ind3.*HC_ind3)./(((HA_ind3+HC_ind3).*(HC_ind3+HD_ind3))+((HA_ind3+HB_ind3).*(HB_ind3+HD_ind3)))*100; %2(ad-bc)/(a+c)(c+d)+(a+b)(b+d)

HA_ind2 = sum(modelpredictionign_ind2 >= 0 & actualdatasign_ind2 >= 0,2);
HB_ind2 = sum(modelpredictionign_ind2 > actualdatasign_ind2,2);
HC_ind2 = sum(modelpredictionign_ind2 < actualdatasign_ind2,2);
HD_ind2 = sum(modelpredictionign_ind2 <= 0 & actualdatasign_ind2 <= 0,2);

HSStest2 = 2.*(HA_ind2.*HD_ind2 - HB_ind2.*HC_ind2)./(((HA_ind2+HC_ind2).*(HC_ind2+HD_ind2))+((HA_ind2+HB_ind2).*(HB_ind2+HD_ind2)))*100; %2(ad-bc)/(a+c)(c+d)+(a+b)(b+d)
    
    
end
