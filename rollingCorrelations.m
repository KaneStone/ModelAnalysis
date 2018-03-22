function [correlations,trend] = rollingCorrelations(data,data2,rollingwindow)
% rolling correlation function

% data = [x,y,time] - data1 and data2 have to be the same time length
% data2 = [time]
% rolling window = length of rolling window

%% calculating the number of correlations
corrlength = size(data,3) - rollingwindow;

%% taking rolling correlations
for i = 1:size(data,1)    
    count = 1;
    for j = 1:corrlength
        if j == 17
            a = 1;
        end
        [correlations.r(i,:,j),correlations.p(i,:,j)] = corr(squeeze(data(i,:,count:rollingwindow+count-1))',...
            squeeze(data2(count:rollingwindow+count-1)));
        count = count+1;
    end
    % calculate trends
    for j = 1:size(data,2)        
        [trend.r(i,j,:),trend.int(i,j,:,:)] = regress(squeeze(correlations.r(i,j,:)),...
            [ones(size(correlations.r,3),1),[1:size(correlations.r,3)]']);        
    end
        
end

end
