function [b,ball,pvalue,uncertainty] = ozoneRegressionEnsAve(data,numyears,ARmodel)

%% 

%%
for i = 1:size(data,2) % pressure
    for j = 1:size(data,3) % latitudes        
        for k = 1:12            
            b(i,j,k,:) = regress(data(k:12:end,i,j),[ones(size(data(k:12:end,i,j),1),1),...
                (1:size(data(k:12:end,i,j),1))']);
        end     
        [ball(i,j,:),~,ball_residuals(:,i,j),~,ballstats] = regress(data(:,i,j),[ones(size(data(:,i,j),1),1),...
                (1:size(data(:,i,j),1))']);
        
        AC2(i,j,:,:) = corrcoef(ball_residuals(1:end-1,i,j),ball_residuals(2:end,i,j),'rows','pairwise');
        AC1(i,j) = AC2(i,j,1,2);                
        %AC1(i,j) = corr(ball_residuals(1:end-1,i,j),ball_residuals(2:end,i,j));        
        twosigma(i,j) =  nanstd(ball_residuals(:,i,j),1)*2;%*sqrt((1+AC1(i,j))./(1-AC1(i,j)));
        uncertainty.u(i,j) = nanstd(ball_residuals(:,i,j),1)./((length(ball_residuals(:,i,j))./120).^(3/2));
        
        if uncertainty.u(i,j)*2*sqrt((1+AC1(i,j))./(1-AC1(i,j))) >= abs(ball(i,j,2))*120
            uncertainty.sig(i,j) = -1;
        else
            uncertainty.sig(i,j) = 0;
        end
        
        pvalue(1,i,j) = ballstats(3)*sqrt((1+AC1(i,j))./(1-AC1(i,j)));      
%         linmodel(:,i,j) = (sum([ones(size(data(:,i,j),1),1),...
%                 (1:size(data(:,i,j),1))']'.*squeeze(ball(i,j,:))))';
    end    
end

pvalue (pvalue > .046) = 1;
pvalue (pvalue <= .046) = 0;

abc = ball(:,:,2)*204 >= twosigma./2;

% %% obtaining uncertainty
% 
% for i = 1:size(data,2) % pressure
%     i
%     for j = 1:size(data,3) % latitudes                        
%         for k = 1:10000
%             b_allresrand = ball_residuals(randsample(1:size(ball_residuals(:,i,j),1),size(ball_residuals(:,i,j),1)),i,j);
%             newdata = linmodel(:,i,j) + b_allresrand;    
%             ballnew(k,:) = regress(newdata,[ones(size(newdata,1),1),...
%                     (1:size(newdata,1))']);            
%         end
%         ball_sig(i,j) = std(ballnew(:,2),0,1)*sqrt((1+AC1(i,j))./(1-AC1(i,j)));
%     end
% end
% 
% if ARmodel == 1
% elseif ARmodel == 2
%     end
end


% find standard deviation of individual mon

% if ARmodel == 2
%         for i = 1:size(O3Anomaly.percent_residuals,2)
%             for j = 1:size(O3Anomaly.percent_residuals,3)
%                 ARcoeffs(:,i,j) = regress(O3Anomaly.percent_residuals(3:end,i,j),...
%                     [ones(length(O3Anomaly.percent_residuals(1:end-2,i,j)),1),...
%                     O3Anomaly.percent_residuals(1:end-2,i,j),...
%                     O3Anomaly.percent_residuals(2:end-1,i,j)]);
%                 %AR2coeffs(:,i,j) = aryule(O3Anomaly.percent_residuals(1:end-1,i,j),1);
% %                 O3Anomaly.percent_residuals_months_AR2(:,i,j) = ...
% %                     O3Anomaly.percent_residuals(:,i,j)*ARcoeffs(1,i,j)+...
% %                     O3Anomaly.percent_residuals(:,i,j)*ARcoeffs(2,i,j)+...
% %                     O3Anomaly.percent_residuals(:,i,j)*ARcoeffs(3,i,j);   
%                 
%                 O3Anomaly.percent_residuals_months_AR2(:,i,j) = ...
%                     [O3Anomaly.percent_residuals_months(1:end-2,i,j)*ARcoeffs(1,i,j)+...
%                     O3Anomaly.percent_residuals_months(1:end-2,i,j)*ARcoeffs(2,i,j)+...
%                     O3Anomaly.percent_residuals_months(2:end-1,i,j)*ARcoeffs(3,i,j);O3Anomaly.percent_residuals_months(end-1:end,i,j)];
%                 
%             end
%         end
%     elseif ARmodel == 1 %Y = bX(t-1) + e
%         for i = 1:size(O3Anomaly.percent_residuals,2)
%             for j = 1:size(O3Anomaly.percent_residuals,3)
%                 
%                 %dwp = dwtest(
%                 
%                 ARcoeffs(:,i,j) = regress(O3Anomaly.percent_residuals_months(2:end,i,j),...
%                     [ones(length(O3Anomaly.percent_residuals_months(1:end-1,i,j)),1),...
%                     O3Anomaly.percent_residuals_months(1:end-1,i,j)]);                    
%                 %AR2coeffs(:,i,j) = aryule(O3Anomaly.percent_residuals(1:end-1,i,j),1);
%                 O3Anomaly.percent_residuals_months_AR2(:,i,j) = ...
%                     [O3Anomaly.percent_residuals_months(1:end-1,i,j)*ARcoeffs(1,i,j)+...
%                     O3Anomaly.percent_residuals_months(1:end-1,i,j)*ARcoeffs(2,i,j);O3Anomaly.percent_residuals_months(end,i,j)];
%                 
%                 O3Anomaly.percent_residuals_months_ARtest(:,i,j) = ...
%                     O3Anomaly.percent_residuals_months(:,i,j)*ARcoeffs(1,i,j)+...
%                     O3Anomaly.percent_residuals_months(:,i,j)*ARcoeffs(2,i,j);
%                 
%                 O3Anomaly.percent_residuals_months_ARgone(:,i,j) = ...
%                     [O3Anomaly.percent_residuals_months(2:end,i,j)- ... 
%                     O3Anomaly.percent_residuals_months(1:end-1,i,j)*ARcoeffs(1,i,j)- ...
%                     O3Anomaly.percent_residuals_months(1:end-1,i,j)*ARcoeffs(2,i,j);O3Anomaly.percent_residuals_months(end,i,j)];
%                     
%             end
%         end
%     end
    