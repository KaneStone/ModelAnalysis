function [weightedmean] = weightedaverage(data,latitudes)

% takes the latitude weighted averaged of a dataset over the first
% dimension

% data = [lats,time]

%yearly latitude weighted averages
W = cosd(latitudes);
weightedmean = zeros(1,size(data,2));
for i = 1:size(data,2)
    temp = data(:,i);
    weightedmean(i) = nansum(temp(:).*W(:))./nansum(W(:));       
end

end