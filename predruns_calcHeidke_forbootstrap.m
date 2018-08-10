function [GSS] = predruns_calcHeidke_forbootstrap(modelprediction,actualdata)
%calculates Heidke

%data needs to be as [no. runs,time,...]
modelpredictionign_ind = sign(modelprediction);
actualdatasign_ind =  sign(actualdata);

modelpredictionresults_ind = modelpredictionign_ind;
modelpredictionresults_ind (modelpredictionresults_ind ~= actualdatasign_ind) = -2;
modelpredictionresults_ind (modelpredictionresults_ind == actualdatasign_ind) = 1;
modelpredictionresults_ind (modelpredictionresults_ind == -2) = 0;

GSS = (sum(modelpredictionresults_ind) - length(modelpredictionresults_ind)./2)./...
    (length(modelpredictionresults_ind)-length(modelpredictionresults_ind)./2)*100;   


end