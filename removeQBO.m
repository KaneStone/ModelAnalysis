function [dataQBOremoved,dataQBOremoved_nolag,b] = removeQBO(data,Uwind)

%% needs to be lagged
%lag by up to 5 months
%dataQBOremoved = zeros(size(data));
data = data(2:end,:,:);
for i = 1:size(data,2)
    for j = 1:size(data,3)
        count = 0;
        for k = 1:2            
            if count >= i 
                b(i,j,k,:) = regress(data(:,i,j),[ones(length(data(:,i,j)),1),squeeze(Uwind(1:end-1,i-count+12,:))]);
                b2(i,j,k,:) = regress(data(:,i,j),[ones(length(data(:,i,j)),1),squeeze(Uwind(1:end-1,i-count+12,1))]);
            else
                b(i,j,k,:) = regress(data(:,i,j),[ones(length(data(:,i,j)),1),squeeze(Uwind(2:end,i-count,:))]);
                b2(i,j,k,:) = regress(data(:,i,j),[ones(length(data(:,i,j)),1),squeeze(Uwind(2:end,i-count,1))]);
            end
            count = count+1;
        end
        % find maximum
        [maxb(i,j,1),maxbind(i,j,1)] = max(abs(squeeze(b2(i,j,:,2))));
        %[maxb(i,j,2),maxbind(i,j,2)] = max(abs(squeeze(b2(i,j,:,3))));
        
        if maxbind(i,j,1) > i 
            Uwindtouse(:,i,1) = Uwind(1:end-1,i-(maxbind(i,j,1)-1)+12,1); 
        else
            Uwindtouse(:,i,1) = Uwind(2:end,i-(maxbind(i,j,1)-1),1);
        end        
        dataQBOremoved(:,i,j) = data(:,i,j) - maxb(i,j,1)*Uwindtouse(:,i,1);% - maxb(i,j,2)*Uwindtouse(:,i,2);       
        dataQBOremoved_nolag(:,i,j) = data(:,i,j) - b(i,j,1,2)*Uwind(2:end,i,1) - b(i,j,1,3)*Uwind(2:end,i,2);
    end
end
end