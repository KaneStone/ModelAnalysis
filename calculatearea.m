function [areadifferencecorrect] = calculatearea(low,high,longitude,latitude,toplot,tactogetherind,tactogetherindlow)

% Barents Kara sea; Laptet, East Siberian, and Chukchi seas, Greenland
% Seas; Sea of Okhotsk and Berling Sea
Seasplot.lon = [25,80;90,235;310,378;100,200];
Seasplot.lat = [65,80;66,80;55,85;40,65];
fraction = .15;
% I want to calculate the area of each curve for each latitude in the
% vector

earthellipsoid = referenceSphere('earth','km');
%area = areaquad(-90,-180,90,180,earthellipsoid)
count = 1;
for k = 1:size(low,1)
    for i = 1:size(Seasplot.lon,1)
        % find latitude vectors for each sea
        %extract lowmean and highmean data
        latind = latitude >= Seasplot.lat(i,1) & latitude <= Seasplot.lat(i,2);
        latbounds(1) = min(latitude(latind));
        latbounds(2) = max(latitude(latind));
        if i == 3
            lowmeantoplot = cat(3,low,low(:,:,1:11));
            highmeantoplot = cat(3,high,high(:,:,1:11));
            lon = [longitude;longitude(1:11)+360];
            lonind = lon >= Seasplot.lon(i,1) & lon < Seasplot.lon(i,2);
            
            lonbounds(1) = min(lon(lonind));
            lonbounds(2) = max(lon(lonind));
            
        else
            lon = longitude;
            lowmeantoplot = low;
            highmeantoplot = high;
            lonind = longitude >= Seasplot.lon(i,1) & longitude < Seasplot.lon(i,2);
            
            lonbounds(1) = min(longitude(lonind));
            lonbounds(2) = max(longitude(lonind));
            
        end

        lowmeanextract = squeeze(lowmeantoplot(k,latind,lonind));
        highmeanextract = squeeze(highmeantoplot(k,latind,lonind));
        if toplot
            fig = createfig('medium','on');
        else
            fig = createfig('medium','off');
        end
        ylim([Seasplot.lat(i,1),Seasplot.lat(i,2)]);
        [cindlow,chl] = contour(lon(lonind),latitude(latind),lowmeanextract,[.15,.15],'r');
        hold on
        [cindhigh,chh] = contour(lon(lonind),latitude(latind),highmeanextract,[.15,.15],'b');        
        cindhigh (cindhigh == .15) = NaN;
        cindhigh(2,isnan(cindhigh(1,:))) = NaN;
        cindhighnanind = [find(isnan(cindhigh(1,:))),length(cindhigh(1,:))+1];
        
        cindlow (cindlow == .15) = NaN;
        cindlow(2,isnan(cindlow(1,:))) = NaN;
        cindlownanind = [find(isnan(cindlow(1,:))),length(cindlow(1,:))+1];
                                
        % using a function to find area based on different cases.
        
        [areadifferencecorrect(k,i)] = findarea(cindhigh,cindlow,cindhighnanind,cindlownanind,lonbounds,latbounds,count,toplot,tactogetherind,tactogetherindlow);

        count = count+1;
        if ~toplot
            close all
        end
    end
end
end