function [yearmeananomaly, yearmean] = TCOanomaly(data,latboundaries,latitudes,anomyear,refyears,months)

% data needs to be the size of [latitude, time;

%yearly latitude weighted averages
W = cosd(latitudes);
W (latitudes <latboundaries(1) | latitudes > latboundaries(2)) = [];
data(latitudes < latboundaries(1) | latitudes > latboundaries(2),:) = []; 

count = 1;
for j = months;    
    anomalymean(:,count) = nanmean(data(:,12*(anomyear(1)-refyears(1))+j:12:12*(anomyear(2)-refyears(1))),2);
    count = count+1;
end

count = min(months);
for i = 1:refyears(2)-refyears(1)
    temp = (data(:,count:count+length(months)-1)-anomalymean)./anomalymean*100;    
    temp2 = (data(:,count:count+length(months)-1));
    W1 = W;
    W1 = repmat(W1,[1,length(months)]);
    W1 (isnan(temp)) = NaN;
    yearmeananomaly(i) = nansum(temp(:).*W1(:))./nansum(W1(:));    
    yearmean(i) = nansum(temp2(:).*W1(:))./nansum(W1(:));    
    count = count+12;
end

end
