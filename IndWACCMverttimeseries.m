function [NDfinal,NDdateslatmean_GMH] = IndWACCMverttimeseries(variable,lats,run,remove2002,yearmin,yearmax,altsmean)


szlat = size(lats);
convert_to_ND = 1;

%% Reading in WACCM data
[NumberDensity, MolConc, Pressure, Altitude, GMH, Latitudes, datedata] = ReadWACCMvertical(variable);

% extracting latitude ranges and dates 

%extracting dates

%finding latitudes
if numel(lats) > 1
    for p = 1:szlat(1)
        latindex(p).p = find(Latitudes >= lats(p,1) & Latitudes <= lats(p,2));
    end
else
    [~,latindex] = min(abs(Latitudes - lats));
end

%extracting dates
NDdates = NumberDensity.(run)(:,:,...
    datedata.(run).date >= yearmin & datedata.(run).date < yearmax);
MolConcdates = MolConc.(run)(:,:,...
    datedata.(run).date >= yearmin & datedata.(run).date < yearmax);
datetemp = datedata.(run).date(datedata.(run).date >= yearmin & ...
    datedata.(run).date < yearmax);
Pressuredates = Pressure.(run)(:,:,...
    datedata.(run).date >= yearmin & datedata.(run).date < yearmax);
Altitudedates = Altitude.(run)(:,:,...
    datedata.(run).date >= yearmin & datedata.(run).date < yearmax);
GMHdates = GMH.(run)(:,:,...
    datedata.(run).date >= yearmin & datedata.(run).date < yearmax);

switch remove2002
    case 'all_lats'
        NDdates.(run)(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        MolConcdates.(run)(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        Pressuredates.(run)(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        Altitudedates.(run)(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        GMHdates.(run)(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
    case 'Southern'        
        NDdates.(run)(Latitudes < 0,...
            datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        MolConcdates.(run)(Latitudes < 0,...
            datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        Altitudedates.(run)(Latitudes < 0,...
            datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        GMHdates.(run)(Latitudes < 0,...
            datetemp >= 20021001 & datetemp < 20031001) = NaN; 
end
NDdates = cat(3,NDdates,NDdates(:,:,end-1:end));
MolConcdates = cat(3,MolConcdates,MolConcdates(:,:,end-1:end));
Pressuredates = cat(3,Pressuredates,Pressuredates(:,:,end-1:end))./100;
Altitudedates = cat(3,Altitudedates,Altitudedates(:,:,end-1:end));
GMHdates = cat(3,GMHdates,GMHdates(:,:,end-1:end));

[NDdatesint, ~] = presinterpWACCM(Pressuredates,NDdates);        
[MolConcdatesint, presint] = presinterpWACCM(Pressuredates,MolConcdates);        

%altitude interpolation
%NDdatesAltint = interp1(
for i = 1:size(NDdates,1)
    for j = 1:size(NDdates,3)        
        NDdatesAltint(i,j,:) = interp1(squeeze(Altitudedates(i,:,j)),squeeze(NDdates(i,:,j)),500:1000:64500);        
        NDdatesGMHint(i,j,:) = interp1(squeeze(GMHdates(i,:,j)),squeeze(NDdates(i,:,j)),500:1000:64500);        
    end
end

if latindex(1).p > 1
    NDdateslatmean = squeeze(nanmean(NDdatesAltint(latindex(1).p,:,:)));
    NDdateslatmean_GMH = squeeze(nanmean(NDdatesGMHint(latindex(1).p,:,:)));
end
    
%NDfinal = nanmean(NDdateslatmean(:,altsmean),2);
NDfinal = nanmean(NDdateslatmean_GMH(:,altsmean),2);

end