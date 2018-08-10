function [GSS,HSS] = predruns_calcHeidke(modelprediction,actualdata,obs)
%calculates Heidke

%data needs to be as [no. runs,time,...]
if obs
    modelpredictionign_ind = sign(modelprediction);
    actualdatasign_ind =  sign(actualdata);

    modelpredictionresults_ind = modelpredictionign_ind;
    modelpredictionresults_ind (modelpredictionresults_ind ~= actualdatasign_ind) = -2;
    modelpredictionresults_ind (modelpredictionresults_ind == actualdatasign_ind) = 1;
    modelpredictionresults_ind (modelpredictionresults_ind == -2) = 0;

    GSS.ind = (sum(modelpredictionresults_ind,1) - size(modelpredictionresults_ind,1)./2)./...
        (size(modelpredictionresults_ind,1)-size(modelpredictionresults_ind,1)./2)*100;    

    HA_ind = sum(modelpredictionign_ind >= 0 & actualdatasign_ind >= 0,1);
    HB_ind = sum(modelpredictionign_ind > actualdatasign_ind,1);
    HC_ind = sum(modelpredictionign_ind < actualdatasign_ind,1);
    HD_ind = sum(modelpredictionign_ind <= 0 & actualdatasign_ind <= 0,1);

    HSS.ind = 2.*(HA_ind.*HD_ind - HB_ind.*HC_ind)./(((HA_ind+HC_ind).*(HC_ind+HD_ind))+((HA_ind+HB_ind).*(HB_ind+HD_ind))); %2(ad-bc)/(a+c)(c+d)+(a+b)(b+d)    
else
    modelpredictionign_ind = sign(modelprediction);
    actualdatasign_ind =  sign(actualdata);

    modelpredictionresults_ind = modelpredictionign_ind;
    modelpredictionresults_ind (modelpredictionresults_ind ~= actualdatasign_ind) = -2;
    modelpredictionresults_ind (modelpredictionresults_ind == actualdatasign_ind) = 1;
    modelpredictionresults_ind (modelpredictionresults_ind == -2) = 0;

    GSS.ind = (sum(modelpredictionresults_ind,2) - size(modelpredictionresults_ind,2)./2)./...
        (size(modelpredictionresults_ind,2)-size(modelpredictionresults_ind,2)./2)*100;
    GSS.mean = nanmean(squeeze(GSS.ind),1);

    HA_ind = sum(modelpredictionign_ind >= 0 & actualdatasign_ind >= 0,2);
    HB_ind = sum(modelpredictionign_ind > actualdatasign_ind,2);
    HC_ind = sum(modelpredictionign_ind < actualdatasign_ind,2);
    HD_ind = sum(modelpredictionign_ind <= 0 & actualdatasign_ind <= 0,2);

    HSS.ind = 2.*(HA_ind.*HD_ind - HB_ind.*HC_ind)./(((HA_ind+HC_ind).*(HC_ind+HD_ind))+((HA_ind+HB_ind).*(HB_ind+HD_ind))); %2(ad-bc)/(a+c)(c+d)+(a+b)(b+d)
    HSS.mean = nanmean(squeeze(HSS.ind),1);
end


end