function [] = latPresTrends(data,solar,SPE,pressure,dates,intpres,SPEpressure)

years = CCMI_years(dates,1);
SPE = permute(SPE,[1,3,2]);
SPEpressure = permute(repmat(SPEpressure,[1,size(SPE,1),size(SPE,3)]),[2,1,3]);
%% put onto regular pressure and monthly arrange
for i = 1:size(data,1)
    [dataRegPres(i,:,:),~] = intRegPres(squeeze(data(i,:,:)),...
        squeeze(pressure(i,:,:))./100,intpres);
end

%% put SPE onto regular pressure
for i = 1:size(SPE,1)
    SPERegPres(i,:,:) = intRegPres(squeeze(SPE(i,:,:)),...
        squeeze(SPEpressure(i,:,:)),intpres);
end


%% rearranging so all have the same dimensions
solar = solar';
solar = solar(:);
SPERegPres = permute(SPERegPres,[2,3,1]);
SPERegPres = reshape(SPERegPres,[size(SPERegPres,1),size(SPERegPres,2)*size(SPERegPres,3)]);
dataRegPres = dataRegPres(:,:,1:end-11);
%% removing solar cycle and SPE from each month separately
regressyears = [2000,2014];

%% testing orthogonality
soltemp = (solar(12:12:end)-nanmean(solar(12:12:end)))./std(solar(12:12:end));
SPEtemp = (squeeze(SPERegPres(27,12:12:end))' - nanmean(squeeze(SPERegPres(27,12:12:end))))./std(squeeze(SPERegPres(27,12:12:end)));
[bsol,~,rsol] = regress(soltemp,[ones(length(soltemp),1),SPEtemp]);
[bSPE,~,rSPE] = regress(SPEtemp,[ones(length(soltemp),1),soltemp]);


% b = regress(totest_n,...
%     [ones(length(totest_n),1),...
%     SPE_totest_n,...
%     solar_totest_n,...
%     [1:length(totest_n)]']);
b = zeros(12,size(dataRegPres,1),size(dataRegPres,2),4);

for i = 1:12
    solartouse(:,i) = (solar(12+i:12:end) - nanmean(solar(12+i:12:end)))./std(solar(12+i:12:end));    
    lngth = size(solartouse,1);
    for j = 1:size(dataRegPres,1)
        % appropriate lags for SPE
        if i >= 3 && i <= 9
            if j > 48
                SPEmonth = nanmean(cat(3,SPERegPres(:,12+i-5:12:end-12+i),SPERegPres(:,12+i-4:12:end-12+i),...
                    SPERegPres(:,12+i-3:12:end-12+i),SPERegPres(:,12+i-2:12:end-12+i),SPERegPres(:,12+i-1:12:end-12+i),...
                    SPERegPres(:,12+i:12:end-12+i)),3);
            else
                SPEmonth = nanmean(cat(3,SPERegPres(:,12+i-3:12:end-12+i),SPERegPres(:,12+i-2:12:end-12+i),...
                   SPERegPres(:,12+i-1:12:end-12+i),SPERegPres(:,12+i:12:end-12+i)),3);
%                 SPEmonth = nanmean(cat(3,SPERegPres(:,12+i-5:12:end-12+i),SPERegPres(:,12+i-4:12:end-12+i),...
%                     SPERegPres(:,12+i-3:12:end-12+i),SPERegPres(:,12+i-2:12:end-12+i),SPERegPres(:,12+i-1:12:end-12+i),...
%                     SPERegPres(:,12+i:12:end-12+i)),3);
%                 SPEmonth = nanmean(cat(3,SPERegPres(:,12+i-2:12:end-12+i),...
%                     SPERegPres(:,12+i-1:12:end-12+i),SPERegPres(:,12+i:12:end-12+i)),3);
            end
        else
            if j < 49
                SPEmonth = nanmean(cat(3,SPERegPres(:,12+i-5:12:end-12+i),SPERegPres(:,12+i-4:12:end-12+i),...
                    SPERegPres(:,12+i-3:12:end-12+i),SPERegPres(:,12+i-2:12:end-12+i),SPERegPres(:,12+i-1:12:end-12+i),...
                    SPERegPres(:,12+i:12:end-12+i)),3);
            else
                SPEmonth = nanmean(cat(3,SPERegPres(:,12+i-3:12:end-12+i),SPERegPres(:,12+i-2:12:end-12+i),...
                    SPERegPres(:,12+i-1:12:end-12+i),SPERegPres(:,12+i:12:end-12+i)),3);
%                 SPEmonth = nanmean(cat(3,SPERegPres(:,12+i-5:12:end-12+i),SPERegPres(:,12+i-4:12:end-12+i),...
%                     SPERegPres(:,12+i-3:12:end-12+i),SPERegPres(:,12+i-2:12:end-12+i),SPERegPres(:,12+i-1:12:end-12+i),...
%                     SPERegPres(:,12+i:12:end-12+i)),3);
            end
        end
        for k = 1:size(dataRegPres,2) 
            if sum(SPEmonth(k,:)) == 0
                SPEtouse(:,i,j,k) = zeros(size(SPEmonth,2),1);
            else
                SPEtouse(:,i,j,k) = ((SPEmonth(k,:) - nanmean(SPEmonth(k,:)))./std(SPEmonth(k,:)))';
            end            
            datatouse(:,i,j,k) = squeeze((dataRegPres(j,k,12+i:12:end) - nanmean(dataRegPres(j,k,12+i:12:end))));
            lineartouse = (([1:lngth] - nanmean(1:lngth))./std(1:lngth))';
            
%             % first off, orthoganilise by removing residuals
%             [~,~,rsol] = regress(solartouse(:,i),[ones(lngth,1),SPEtouse(:,i,j,k),lineartouse]);
%             if sum(SPEmonth(k,:)) == 0
%                 rSPE = [];
%             else
%                 [~,~,rSPE] = regress(SPEtouse(:,i,j,k),[ones(lngth,1),solartouse(:,i),lineartouse]);
%             end
                        
            if sum(SPEmonth(k,:)) == 0                
                b(i,j,k,[1,3,4]) = regress(datatouse(:,i,j,k),[ones(lngth,1),solartouse(:,i),lineartouse]);
                
            else
                 b(i,j,k,:) = regress(datatouse(:,i,j,k),[ones(lngth,1),SPEtouse(:,i,j,k),solartouse(:,i),lineartouse]);                
            end
        end
    end
end
%%
mon = 2;
lat = 83;
pres = 27;
figure
plot(datatouse(:,mon,lat,pres),'k');
hold on
plot(squeeze(b(mon,lat,pres,1))+squeeze(b(mon,lat,pres,2))*SPEtouse(:,mon,lat,pres)+b(mon,lat,pres,3)*solartouse(:,mon)+b(mon,lat,pres,4)*lineartouse,'--');
plot(datatouse(:,mon,lat,pres) - squeeze(b(mon,lat,pres,1)) - squeeze(b(mon,lat,pres,2))*SPEtouse(:,mon,lat,pres) - b(mon,lat,pres,3)*solartouse(:,mon));

end
