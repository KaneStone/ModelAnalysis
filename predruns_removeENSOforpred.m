function [NINOmonth] = predruns_removeENSOforpred(dataMonthArrange,lats,lons,obs)

%% calculate ENSO and remove it from the surface temperature data

latlimits = [-5 5];
lonlimits = [190 240];

latindex = lats >= latlimits(1) & lats <= latlimits(2);
lonindex = lons >= lonlimits(1) & lons <= lonlimits(2);
if obs
    for j = 1:12
        NINO_mn(j,:,:) = squeeze(nanmean(dataMonthArrange(j,:,latindex,lonindex),4));
        NINO_mn2(j,:) = squeeze(nanmean(NINO_mn(j,:,:),3));
        NINOmonth(j,:) = (squeeze(NINO_mn2(j,:)) - nanmean(squeeze(NINO_mn2(j,:)),2))./std(squeeze(NINO_mn2(j,:)),1,2);
    end   
else
    for j = 1:12
        NINO_mn(:,j,:,:) = squeeze(nanmean(dataMonthArrange(:,j,:,latindex,lonindex),5));
        NINO_mn2(:,j,:) = squeeze(nanmean(NINO_mn(:,j,:,:),4));
        NINOmonth(:,j,:) = (squeeze(NINO_mn2(:,j,:)) - nanmean(squeeze(NINO_mn2(:,j,:)),2))./std(squeeze(NINO_mn2(:,j,:)),1,2);
    end   
end

NINO34all = NINOmonth(:,:);    
% 
%     %% take ensemble averages
%     varmonth2 = inputs.varmonth;
%     isless = varmonth2 <= 2;
%     varmonth2 = (varmonth2 + isless*12);
% 
%     laglength = 3;    
%     
%     for i = 1:length(inputs.varmonthtomean)
%         monthstomean(i) = find(inputs.varmonth == inputs.varmonthtomean(i));
%     end
%     
%     %varmonths
%     datamonthall = permute(dataMonthArrange,[1,4,5,2,3]);
%     datamonthall = datamonthall(:,:,:,:);
% 
%     % Here I am regressing with the of the seasonal average
%     for j = 1:length(inputs.varmonth)
%         %together(:,:,:,:,j) = squeeze(datamonthall(:,:,:,varmonth2(j):12:end-inputs.varmonth(1)+inputs.varmonth(end)));
%         together(:,:,:,:,j) = squeeze(datamonthall(:,:,:,varmonth2(j):12:end));
%     end
%     temp = nanmean(together(:,:,:,:,inputs.varmonthtomean-2),5);
% 
%     if inputs.removeENSO        
%         % finding maximum regression lag time
%         for j = 1:size(datamonthall,1) % members
%             for k = 1:length(varmonth2)
%                 for l = 1:size(datamonthall,2) % latitudes
%                     for m = 1:size(datamonthall,3) % longitudes
%                         for lag = 1:laglength % lag
%                             regressors = [ones(length(squeeze(NINO34all(j,varmonth2(k)-lag+1:12:end-lag))),1),...
%                                 squeeze(NINO34all(j,varmonth2(k)-lag+1:12:end-lag))'];
%                             [b(j,k,l,m,lag,:)] = regress(squeeze(datamonthall(j,l,m,varmonth2(k):12:end-lag)),...
%                                 regressors);                        
%                             % finding largest lag correlation
%                         end
%                         [~,llc(j,k,l,m)] = max(abs(squeeze(b(j,k,l,m,:,2))));
%                         blag(j,k,l,m,:) = b(j,k,l,m,llc(j,k,l,m),:);
%                         dataVarMonth(j,:,k,l,m) = squeeze(datamonthall(j,l,m,varmonth2(k):12:end-laglength)) - ...
%                             squeeze(b(j,k,l,m,llc(j,k,l,m),2))*...
%                             squeeze(NINO34all(j,varmonth2(k)-llc(j,k,l,m)+1:12:end-laglength))';                        
%                     end
%                 end
%             end        
%         end
%         dataVarMonthAve = squeeze(nanmean(dataVarMonth(:,:,monthstomean,:,:),3));        
%     else        
%         dataVarMonthAve = permute(temp,[1,4,2,3]);
%         dataVarMonth = permute(together,[1,4,5,2,3]);
%     end        

end