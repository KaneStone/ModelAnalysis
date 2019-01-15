function [SOI] = predruns_calculateSOI(PSL,PSLall,latitudes,longitudes,indmonth)

Tahitilat = 0; %-17.5
Tahitilon = 214.5; 
Darwinlat = 0; %-12.5
Darwinlon = 131; %131

[~,tahlatindex] = min(abs(latitudes-Tahitilat));
[~,tahlonindex] = min(abs(longitudes-Tahitilon));
[~,darlatindex] = min(abs(latitudes-Darwinlat));
[~,darlonindex] = min(abs(longitudes-Darwinlon));

fieldnames = fields(PSL);
for k = 1:length(fieldnames)
    varcombine.(fieldnames{k}) = cat(1,PSL.(fieldnames{k}).highind,PSL.(fieldnames{k}).lowind);    
    
    standarddeviation.(fieldnames{k}).tahiti = std(varcombine.(fieldnames{k})(:,tahlatindex,tahlonindex),1);
    standarddeviation.(fieldnames{k}).darwin = std(varcombine.(fieldnames{k})(:,darlatindex,darlonindex),1);
    
    average.(fieldnames{k}).tahiti = nanmean(varcombine.(fieldnames{k})(:,tahlatindex,tahlonindex),1);
    average.(fieldnames{k}).darwin = nanmean(varcombine.(fieldnames{k})(:,darlatindex,darlonindex),1);
   
    standardized.(fieldnames{k}).tahiti = (varcombine.(fieldnames{k})(:,tahlatindex,tahlonindex) - ...
        average.(fieldnames{k}).tahiti)./standarddeviation.(fieldnames{k}).tahiti;
    standardized.(fieldnames{k}).darwin = (varcombine.(fieldnames{k})(:,darlatindex,darlonindex) - ...
        average.(fieldnames{k}).darwin)./standarddeviation.(fieldnames{k}).darwin;
    
end

stdforanomaly_tahiti =  std(squeeze(PSLall.highcl.data(tahlonindex,tahlatindex,indmonth:12:end)),1);
averageforanomaly_tahiti =  nanmean(squeeze(PSLall.highcl.data(tahlonindex,tahlatindex,indmonth:12:end)));
stdforanomaly_darwin =  std(squeeze(PSLall.highcl.data(darlonindex,darlatindex,indmonth:12:end)),1);
averageforanomaly_darwin =  nanmean(squeeze(PSLall.highcl.data(darlonindex,darlatindex,indmonth:12:end)));

standforanom_tahiti = (squeeze(PSLall.highcl.data(tahlonindex,tahlatindex,indmonth:12:end)) - averageforanomaly_tahiti)./...
    stdforanomaly_tahiti;
standforanom_darwin = (squeeze(PSLall.highcl.data(darlonindex,darlatindex,indmonth:12:end)) - averageforanomaly_darwin)./...
    stdforanomaly_darwin;

MSD2 = sqrt(sum((standforanom_tahiti - standforanom_darwin).^2)/length(standforanom_darwin)); 

SOI.lowcl = (standardized.lowcl.tahiti - standardized.lowcl.darwin)./MSD2;
SOI.highcl = (standardized.highcl.tahiti - standardized.highcl.darwin)./MSD2;    
    
end
