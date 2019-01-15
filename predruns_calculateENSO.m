function [NINO34,IOD,NH] = predruns_calculateENSO(TS,TSall,latitudes,longitudes,indmonth)

latlimits = [-5 5];
lonlimits = [190 240];
lonlimitsIOD = [55 90];
lonlimitsNH = [355 5];
latlimitsNH = [75 85];
latindex = find(latitudes >= latlimits(1) & latitudes <= latlimits(2));
lonindex = find(longitudes >= lonlimits(1) & longitudes <= lonlimits(2));
lonindexIOD = find(longitudes >= lonlimitsIOD(1) & longitudes <= lonlimitsIOD(2));
if lonlimitsNH(1) > lonlimitsNH(2)
    lonindexNH = find(longitudes >= lonlimitsNH(1) | longitudes <= lonlimitsNH(2));
else
    lonindexNH = find(longitudes >= lonlimitsNH(1) & longitudes <= lonlimitsNH(2));
end
latindexNH = find(latitudes >= latlimitsNH(1) & latitudes <= latlimitsNH(2));

fieldnames = fields(TS);
for k = 1:length(fieldnames)
    varcombine.(fieldnames{k}) = cat(1,TS.(fieldnames{k}).highind,TS.(fieldnames{k}).lowind);    
    TSallextract.(fieldnames{k}) = nanmean(TSall.(fieldnames{k}).data(lonindex,latindex,:),1);
    TSallextract.(fieldnames{k}) = squeeze(nanmean(TSallextract.(fieldnames{k}),2));
    varcombineextract.(fieldnames{k}) = nanmean(varcombine.(fieldnames{k})(:,latindex,lonindex),2);
    varcombineextract.(fieldnames{k}) = squeeze(nanmean(varcombineextract.(fieldnames{k}),3));
    NINO34.(fieldnames{k}) = (varcombineextract.(fieldnames{k}) - nanmean(TSallextract.(fieldnames{k})(indmonth:12:end)))./std(TSallextract.(fieldnames{k})(indmonth:12:end),1);        
    
    TSallextractIOD.(fieldnames{k}) = nanmean(TSall.(fieldnames{k}).data(lonindexIOD,latindex,:),1);
    TSallextractIOD.(fieldnames{k}) = squeeze(nanmean(TSallextractIOD.(fieldnames{k}),2));
    varcombineextractIOD.(fieldnames{k}) = nanmean(varcombine.(fieldnames{k})(:,latindex,lonindexIOD),2);
    varcombineextractIOD.(fieldnames{k}) = squeeze(nanmean(varcombineextractIOD.(fieldnames{k}),3));
    IOD.(fieldnames{k}) = (varcombineextractIOD.(fieldnames{k}) - ...
        nanmean(TSallextractIOD.(fieldnames{k})(indmonth:12:end)))./std(TSallextractIOD.(fieldnames{k})(indmonth:12:end),1);        
    
    TSallextractNH.(fieldnames{k}) = nanmean(TSall.(fieldnames{k}).data(lonindexNH,latindexNH,:),1);
    TSallextractNH.(fieldnames{k}) = squeeze(nanmean(TSallextractNH.(fieldnames{k}),2));
    varcombineextractNH.(fieldnames{k}) = nanmean(varcombine.(fieldnames{k})(:,latindexNH,lonindexNH),2);
    varcombineextractNH.(fieldnames{k}) = squeeze(nanmean(varcombineextractNH.(fieldnames{k}),3));
    NH.(fieldnames{k}) = (varcombineextractNH.(fieldnames{k}) - ...
        nanmean(TSallextractNH.(fieldnames{k})(indmonth:12:end)))./std(TSallextractNH.(fieldnames{k})(indmonth:12:end),1);        
            
end
end
