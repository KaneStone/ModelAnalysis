function [yearmeananomaly, yearmean] = Verticalanomaly(data,pressure,latboundaries,...
    latitudes,anomyear,refyears,months,yearly_average)

% Vertical Anomaly of zonally averaged data

% data needs to be the size of [latitude, height, time];

%yearly latitude weighted averages
W = cosd(latitudes);
W (latitudes <latboundaries(1) | latitudes > latboundaries(2)) = [];
data(latitudes < latboundaries(1) | latitudes > latboundaries(2),:,:) = []; 
%pressure(latitudes < latboundaries(1) | latitudes > latboundaries(2),:,:) = []; 

yearmeananomaly = [];
yearmean = [];

if yearly_average
    count = 1;
    for j = months;    
        anomalymean(:,:,count) = nanmean(data(:,:,12*(anomyear(1)-refyears(1))+j:12:12*(anomyear(2)-refyears(1))),3);
        count = count+1;
    end
end

count = min(months);
if yearly_average
    for i = 1:refyears(2)-refyears(1)
        temp = (data(:,:,count:count+length(months)-1)-anomalymean)./anomalymean*100;    
        temp2 = (data(:,:,count:count+length(months)-1));
        W1 = W;
        W1 = repmat(W1,[1,length(months)]);    
        for j = 1:length(pressure)
            temp3 = squeeze(temp(:,j,:));
            W2 = W1;
            W2 (isnan(temp3)) = NaN;
            temp4 = squeeze(temp2(:,j,:));
            yearmeananomaly(i,j) = nansum(temp3(:).*W2(:))./nansum(W2(:));    
            yearmean(i,j) = nansum(temp4(:).*W2(:))./nansum(W2(:));            
        end
        count = count+12;
    end
else
    for i = 1:size(data,3)
        %temp = (data(:,:,count:count+length(months)-1)-anomalymean)./anomalymean*100;    
        %temp2 = (data(:,:,count:count+length(months)-1));
        W1 = W;
        %W1 = repmat(W1,[1,length(months)]);    
        for j = 1:length(pressure)
            temp3 = squeeze(data(:,j,i));
            W2 = W1;
            W2 (isnan(temp3)) = NaN;            
            %yearmeananomaly(i,j) = nansum(temp3(:).*W2(:))./nansum(W2(:));    
            yearmean(i,j) = nansum(temp3(:).*W2(:))./nansum(W2(:));     
            
        end
        if isempty(yearmean(i,:))
            yearmean(i,:) = NaN;
        end
        count = count+12;
    end
end
end