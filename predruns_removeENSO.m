function [dataVarMonthAve,dataVarMonth,NINO34all] = predruns_removeENSO(dataMonthArrange,lats,lons,inputs,ClLevel)

%% calculate ENSO and remove it from the surface temperature data

latlimits = [-5 5];
lonlimits = [190 240];

latindex = lats >= latlimits(1) & lats <= latlimits(2);
lonindex = lons >= lonlimits(1) & lons <= lonlimits(2);
    
for j = 1:12
    NINO_mn(:,j,:,:) = squeeze(nanmean(dataMonthArrange(:,j,:,latindex,lonindex),5));
    NINO_mn2(:,j,:) = squeeze(nanmean(NINO_mn(:,j,:,:),4));
    NINOmonth(:,j,:) = (squeeze(NINO_mn2(:,j,:)) - nanmean(squeeze(NINO_mn2(:,j,:)),2))./std(squeeze(NINO_mn2(:,j,:)),1,2);
end   

NINO34all = NINOmonth(:,:);    
save(['/Volumes/MyBook/work/data/predruns/output/NINO34/',ClLevel,inputs.timeperiodvar(1),'-',inputs.timeperiodvar(1),'_',num2str(inputs.detrend)],'NINO34all','NINOmonth');

if inputs.removeENSO
    filename = ['/Volumes/MyBook/work/data/predruns/output/data/',ClLevel,'_TS_ninoremoved_',...
        monthnames(inputs.varmonthtomean,1,1),num2str(inputs.timeperiodvar(1)),'-',...
        num2str(inputs.timeperiodvar(2)),num2str(abs(inputs.lats(1))),'-',...
    num2str(abs(inputs.lats(2))),'_',num2str(inputs.detrend)];
else
    filename = ['/Volumes/MyBook/work/data/predruns/output/data/',ClLevel,'_TS_',...
        monthnames(inputs.varmonthtomean,1,1),num2str(inputs.timeperiodvar(1)),'-',...
        num2str(inputs.timeperiodvar(2)),num2str(abs(inputs.lats(1))),'-',...
    num2str(abs(inputs.lats(2))),'_',num2str(inputs.detrend)];
end
if ~exist([filename,'.mat'],'file')

    %% take ensemble averages
    varmonth2 = inputs.varmonth;
    isless = varmonth2 <= 2;
    varmonth2 = (varmonth2 + isless*12);

    laglength = 3;    
    
    for i = 1:length(inputs.varmonthtomean)
        monthstomean(i) = find(inputs.varmonth == inputs.varmonthtomean(i));
    end
    
    %varmonths
    datamonthall = permute(dataMonthArrange,[1,4,5,2,3]);
    datamonthall = datamonthall(:,:,:,:);

    % Here I am regressing with the of the seasonal average
    for j = 1:length(inputs.varmonth)
        %together(:,:,:,:,j) = squeeze(datamonthall(:,:,:,varmonth2(j):12:end-inputs.varmonth(1)+inputs.varmonth(end)));
        together(:,:,:,:,j) = squeeze(datamonthall(:,:,:,varmonth2(j):12:end));
    end
    temp = nanmean(together(:,:,:,:,inputs.varmonthtomean-2),5);

    if inputs.removeENSO        
        % finding maximum regression lag time
        for j = 1:size(datamonthall,1) % members
            for k = 1:length(varmonth2)
                for l = 1:size(datamonthall,2) % latitudes
                    for m = 1:size(datamonthall,3) % longitudes
                        for lag = 1:laglength % lag
                            regressors = [ones(length(squeeze(NINO34all(j,varmonth2(k)-lag+1:12:end-lag))),1),...
                                squeeze(NINO34all(j,varmonth2(k)-lag+1:12:end-lag))'];
                            [b(j,k,l,m,lag,:)] = regress(squeeze(datamonthall(j,l,m,varmonth2(k):12:end-lag)),...
                                regressors);                        
                            % finding largest lag correlation
                        end
                        [~,llc(j,k,l,m)] = max(abs(squeeze(b(j,k,l,m,:,2))));
                        blag(j,k,l,m,:) = b(j,k,l,m,llc(j,k,l,m),:);
                        dataVarMonth(j,:,k,l,m) = squeeze(datamonthall(j,l,m,varmonth2(k):12:end-laglength)) - ...
                            squeeze(b(j,k,l,m,llc(j,k,l,m),2))*...
                            squeeze(NINO34all(j,varmonth2(k)-llc(j,k,l,m)+1:12:end-laglength))';                        
                    end
                end
            end        
        end
        dataVarMonthAve = squeeze(nanmean(dataVarMonth(:,:,monthstomean,:,:),3));        
    else        
        dataVarMonthAve = permute(temp,[1,4,2,3]);
        dataVarMonth = permute(together,[1,4,5,2,3]);
    end        

    save(filename,'dataVarMonthAve','dataVarMonth');
else
    load(filename);
end

end