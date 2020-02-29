function [GSS] = predruns_calcHeidke_forbootstrap_compemp(modelpredictionign_ind,actualdatasign_ind,condition1)
%calculates Heidke

%data needs to be as [no. runs,time,...]
% modelpredictionign_ind = sign(modelprediction);
% actualdatasign_ind =  sign(actualdata);

modelpredictionresults_ind = modelpredictionign_ind;
modelpredictionresults_ind (modelpredictionresults_ind ~= actualdatasign_ind) = -2;
modelpredictionresults_ind (modelpredictionresults_ind == actualdatasign_ind) = 1;
modelpredictionresults_ind (modelpredictionresults_ind == -2) = 0;
modelpredictionresults_ind (isnan(modelpredictionign_ind)) = NaN;

if condition1
    GSS = (nansum(modelpredictionresults_ind) - nansum(~isnan(modelpredictionresults_ind))./2)./...
        (nansum(~isnan(modelpredictionresults_ind))-nansum(~isnan(modelpredictionresults_ind))./2)*100;   
else
    GSS = (sum(modelpredictionresults_ind) - length(modelpredictionresults_ind)./2)./...
        (length(modelpredictionresults_ind)-length(modelpredictionresults_ind)./2)*100;   
end


end
