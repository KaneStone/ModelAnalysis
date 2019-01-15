function [b,predictorsout,O3Anomaly,uncertainty] = ozoneRegressionTrends(data,QBOdata,solardata,SPEdata,...
    ENSOdata,HFSdata,HFNdata,NO2data,aerosols,years,chemonly,linear,ARmodel,var)

% take regression of ozone anomalies using supplied predictors. 
% QBO data is broken down into 6 proxies

% all data is in the form [pres,lat,time] where the pressure has been
% interpolated to SWOOSH pressure levels

% this function uses two methods: 1. taking the regression on a monthly basis,
% then taking a linear trend of the residuals, 2. Using fourier pairs to
% account for seasonal cycle and including a linear term.

% need to create monthly normalised predictors as well.

numyears = length(years(1):years(2));
numofpred = 1;

%% creating annual cycle if used

ann(:,1) = sin(2.*pi.*(1:numyears*12)./12)';    
ann(:,2) = cos(2.*pi.*(1:numyears*12)./12)';    
ann(:,3) = sin(2.*pi.*2.*(1:numyears*12)./12)';    
ann(:,4) = cos(2.*pi.*2.*(1:numyears*12)./12)';    

%% creating 10 QBO predictors
artificialQBO = 0;
    
if sum(QBOdata(:)) ~= 0
    QBO(:,1) = 2.*(QBOdata(1,:) - min(QBOdata(1,:)))./(max(QBOdata(1,:))-min(QBOdata(1,:)))'-1;    
    
    if artificialQBO
        for i = 1:20
            QBOtemp(:,i) = circshift(QBO(:,1),i);
            dotproduct(i) = dot(QBO(:,1),ENSOtemp(:,i));
            [~,minind] = min(abs(dotproduct));    
        end
        QBO2temp = null(QBO(:,1).');
        QBO(:,2) = [QBO2temp(1,:),QBO2temp(1,size(data,3)-28)]';
        QBO(:,2) = circshift(QBO(:,1),3)';
        QBO(:,2) = (QBOdata(2,:) - nanmean(QBOdata(2,:)));
    end
    
    QBO(:,2) = 2.*(QBOdata(2,:) - min(QBOdata(2,:)))./(max(QBOdata(2,:))-min(QBOdata(2,:)))'-1;    
    
    % creating fourier pairs
    QBO(:,3) = QBO(:,1).*sin(2.*pi.*(1:numyears*12)./12)';    
    QBO(:,4) = QBO(:,2).*sin(2.*pi.*(1:numyears*12)./12)'; 
    QBO(:,5) = QBO(:,1).*cos(2.*pi.*(1:numyears*12)./12)';
    QBO(:,6) = QBO(:,2).*cos(2.*pi.*(1:numyears*12)./12)';               
    QBO(:,7) = QBO(:,1).*sin(2.*pi.*2.*(1:numyears*12)./12)';
    QBO(:,8) = QBO(:,2).*sin(2.*pi.*2.*(1:numyears*12)./12)';               
    QBO(:,9) = QBO(:,1).*cos(2.*pi.*2.*(1:numyears*12)./12)';
    QBO(:,10) = QBO(:,2).*cos(2.*pi.*2.*(1:numyears*12)./12)';               
    numofpred = numofpred+size(QBO,2);
    
    for i = 1:12
        QBOmonth(:,1,i) = 2.*(QBOdata(1,i:12:end) - min(QBOdata(1,i:12:end)))./...
            (max(QBOdata(1,i:12:end))-min(QBOdata(1,i:12:end)))'-1;    
        QBOmonth(:,2,i) = 2.*(QBOdata(2,i:12:end) - min(QBOdata(2,i:12:end)))./...
            (max(QBOdata(2,i:12:end))-min(QBOdata(2,i:12:end)))'-1;    
    end    
else
    QBO = [];
    QBOmonth = [];
end
%% creating heatflux predictor

artificialHF = 0;
if HFSdata
    HFS = 2.*(HFSdata' - min(HFSdata'))./(max(HFSdata') - min(HFSdata'))-1;
    HFN = 2.*(HFNdata' - min(HFNdata'))./(max(HFNdata') - min(HFNdata'))-1;    
    if artificialHF    
        for i = 1:20
            ENSOtemp(:,i) = circshift(HFS(:,1),i);
            ENSOtemp2(:,i) = circshift(HFN(:,1),i);
            dotproduct(i) = dot(HFS(:,1),ENSOtemp(:,i));
            dotproduct2(i) = dot(HFN(:,1),ENSOtemp2(:,i));
            [~,minind] = min(abs(dotproduct));    
            [~,minind2] = min(abs(dotproduct2));    
        end
        HFS(:,2) = circshift(HFS(:,1),2);
        HFN(:,2) = circshift(HFN(:,1),3);
    end        
    
    for i = 1:12
        HFSmonth(:,i) = 2.*(HFSdata(i:12:end)' - min(HFSdata(i:12:end)'))./(max(HFSdata(i:12:end)') - min(HFSdata(i:12:end)'))-1;
        HFNmonth(:,i) = 2.*(HFNdata(i:12:end)' - min(HFNdata(i:12:end)'))./(max(HFNdata(i:12:end)') - min(HFNdata(i:12:end)'))-1;
    end        
else
    HFS = [];
    HFN = [];
    HFSmonth = [];
    HFNmonth = [];
end

%%  aerosols
if aerosols
    aer = [2.*(aerosols' - min(aerosols'))./(max(aerosols') - min(aerosols'))-1]';    
    for i = 1:12        
        aermonth(:,i) = 2.*(aerosols(i:12:end)' - min(aerosols(i:12:end)'))./(max(aerosols(i:12:end)') - min(aerosols(i:12:end)'))-1;        
    end        
else
    aer = [];
    aermonth = [];
end

%% creating ENSO predictor

artificialENSO = 0;
if ENSOdata
    ENSO(:,1) = 2.*(ENSOdata-min(ENSOdata))./(max(ENSOdata)-min(ENSOdata))-1;
    
    if artificialENSO
        for i = 1:20
            ENSOtemp(:,i) = circshift(ENSO(:,1),i);
            dotproduct(i) = dot(ENSO,ENSOtemp(:,i));
            [~,minind] = min(abs(dotproduct));
        end
        ENSO(:,2) = circshift(ENSO(:,1),13)';
    end    
    numofpred = numofpred+2;
    
    for i = 1:12
        ENSOmonth(:,i) = 2.*(ENSOdata(i:12:end)-min(ENSOdata(i:12:end)))./(max(ENSOdata(i:12:end))-min(ENSOdata(i:12:end)))-1;
    end
    
else
    ENSO = [];
    ENSOmonth = [];
end

%% creating solar and SPE/NO2 predictors
artificialsolar = 0;
if solardata
    solar(:,1) = 2.*(solardata - min(solardata))./(max(solardata)-min(solardata))-1;
    if artificialsolar
        for i = 1:50
            solartemp(:,i) = circshift(solar(:,1),i);
            dotproduct(i) = dot(solar,solartemp(:,i));
            [~,minind] = min(abs(dotproduct));
        end
    end
    numofpred = numofpred+1;
    
    for i = 1:12
        solartemp1(:,i) = solardata(i:12:end);
        solarmonth1(:,i) = 2.*(solardata(i:12:end) - min(solardata(i:12:end)))./(max(solardata(i:12:end))-min(solardata(i:12:end)))-1;
    end
    solarmonth2 = nanmean(solartemp1,2); 
    solarmonth2 = repmat(2.*(solarmonth2 - min(solarmonth2))./(max(solarmonth2)-min(solarmonth2))-1,[1,12]);        
    solarmonth3 = 2.*(solartemp1 - min(solartemp1))./(max(solartemp1)-min(solartemp1))-1;        
    solarmonth = solarmonth3;
else
    solar = [];
    solarmonth = [];
end

%% creating SPE predictor
SPElevstouse = [22,27,29,29,29,29,29,29,29,29,29,29]; %= .23 hPa
if sum(SPEdata(:)) ~= 0
    numofpred = numofpred+1;
    % take 5 month weighted average
    for i = 1:size(data,3)        
        for j = 1:5
            SPEtemp(i,j) = SPEdata(SPElevstouse(j),12+i-j);
        end        
    end
    SPE = nanmean(SPEtemp,2);        
    SPE = (SPE - min(SPE))./(max(SPE) - min(SPE));
else
    %SPE = zeros(size(data));
    SPE = [];
end


%% creating NO2 predictor

if sum(NO2data(:)) ~= 0        
    
    % regress NO2 data by solar cycle and linear
    for i = 1:size(NO2data,1)
        for j = 1:size(NO2data,2)
            NO2data_regress(i,j,:) = regress(squeeze(NO2data(i,j,:)),[ones(length(squeeze(NO2data(i,j,:))),1),solar,(1:length(squeeze(NO2data(i,j,:))))']);
            % remove solar from NO2 proxy
            NO2_solarremoved(i,j,:) = squeeze(NO2data(i,j,:)) - NO2data_regress(i,j,2)*solar-NO2data_regress(i,j,3)*(1:length(squeeze(NO2data(i,j,:))))';
        end
    end
    
    
    NO2 = 2.*(NO2data - min(NO2data,[],3))./(max(NO2data,[],3)-min(NO2data,[],3))-1;
    %NO2 = 2.*(NO2_solarremoved - min(NO2_solarremoved,[],3))./(max(NO2_solarremoved,[],3)-min(NO2_solarremoved,[],3))-1;
    
    
    
    detrendNO2 = 0;
    if detrendNO2
        for i = 1:size(NO2data,1)
            for j = 1:size(NO2data,2)
                %NO2data2 = detrend(squeeze(NO2data(i,j,:))
            end
        end
                
        NO2data1 = 1;
    end
    
    for i = 1:12
        NO2month(:,:,i,:) = 2.*(NO2data(:,:,i:12:end) - min(NO2data(:,:,i:12:end),[],3))./(max(NO2data(:,:,i:12:end),[],3)-min(NO2data(:,:,i:12:end),[],3))-1;
        %NO2month(:,:,i,:) = 2.*(NO2_solarremoved(:,:,i:12:end) - min(NO2_solarremoved(:,:,i:12:end),[],3))./(max(NO2_solarremoved(:,:,i:12:end),[],3)-min(NO2_solarremoved(:,:,i:12:end),[],3))-1;
    end    
else
    NO2 = [];
    NO2month = [];
end

%NO2month(:,:,:,8:end) = 0;

%% constructing linear vector and vector of ones
vectorofones = ones(numyears*12,1);
if linear
    linearvector = (1:numyears*12)';    
else
    linearvector = [];    
end
monthvector = vectorofones(1:12:end);


%% constructing O3 anomalies as: (O3(month) - mean(O3(month)))./mean(O3(month)).

% percent anomalies
for i = 1:12
    O3Anomaly.percent(:,:,i,:) = (data(:,:,i:12:end) - nanmean(data(:,:,i:12:end),3))./nanmean(data(:,:,i:12:end),3)*100;
    if strcmp(var,'O3')
        O3Anomaly.ppmv(:,:,i,:) = (data(:,:,i:12:end) - nanmean(data(:,:,i:12:end),3))*1e6;
        O3Anomaly.ppmv2(:,:,i,:) = (data(:,:,i:12:end) - nanmean(data(:,:,:),3))*1e6;
    else
        O3Anomaly.ppmv(:,:,i,:) = (data(:,:,i:12:end) - nanmean(data(:,:,i:12:end),3));
        O3Anomaly.ppmv2(:,:,i,:) = (data(:,:,i:12:end) - nanmean(data(:,:,:),3));
    end
end
O3Anomaly.percent = permute(O3Anomaly.percent(:,:,:),[3,1,2]);
O3Anomaly.ppmv2 = permute(O3Anomaly.ppmv2(:,:,:),[3,1,2]);
O3Anomaly.ppmv = permute(O3Anomaly.ppmv(:,:,:),[3,1,2]);

%% take the regression
% Here I am taking both the full regression including linear trends, and
% also trends without the linear predictor so that I can remove the other
% predictors from the dataset.
for i = 1:size(O3Anomaly.ppmv,2) % pressure

    for j = 1:size(O3Anomaly.ppmv,3) % latitudes     
        
        if j <= 48
            HFtouse = HFS;
            HFtousemonth = HFSmonth;
        else
            HFtouse = HFN;
            HFtousemonth = HFNmonth;
        end
               
        if i >= 13 && sum(NO2data(:)) ~= 0    
            %predictorstouse = [vectorofones,linearvector,QBO,ENSO,HFtouse,aer,squeeze(NO2(i,j,:)),solar];                                            
            predictorstouse = [vectorofones,linearvector,QBO,solar,aer,ENSO,HFtouse,squeeze(NO2(i,j,:)),];                                            
        else
            %predictorstouse = [vectorofones,linearvector,QBO,ENSO,HFtouse,aer,solar];                                            
            predictorstouse = [vectorofones,linearvector,QBO,solar,aer,ENSO,HFtouse];                                            
        end
        predind = 1:size(predictorstouse,2);               
        % regression using fourier pair method
        [b.percent(i,j,predind),~,O3Anomaly.percent_residuals(:,i,j),~,~] = regress(O3Anomaly.percent(:,i,j),predictorstouse);        
        [b.ppmv(i,j,predind),~,O3Anomaly.ppmv_residuals(:,i,j)] = regress(O3Anomaly.ppmv(:,i,j),predictorstouse);      
        %p.percent(i,j,predind) = stats(3,:);        
        uncertainty.u(i,j) = nanstd(O3Anomaly.percent_residuals(:,i,j),0)./((length(O3Anomaly.percent_residuals(:,i,j))./120).^(3/2));        
        AC2(i,j,:,:) = corrcoef(O3Anomaly.percent_residuals(1:end-1,i,j),O3Anomaly.percent_residuals(2:end,i,j),'rows','pairwise');
        AC1(i,j) = AC2(i,j,1,2);
        %AC1(i,j) = corr(O3Anomaly.percent_residuals(1:end-1,i,j),O3Anomaly.percent_residuals(2:end,i,j));
         if i == 27 && j == 27
             a = 1;
         end
        if uncertainty.u(i,j)*2*sqrt((1+AC1(i,j))./(1-AC1(i,j))) >= abs(b.percent(i,j,2))*120 || isnan(uncertainty.u(i,j))
            uncertainty.sig(i,j) = -1;
        else
            uncertainty.sig(i,j) = 0;
        end
        % regression using monthly method
        for k = 1:12    
            
            if sum(QBOdata(:)) ~= 0
                QBOpred = QBOmonth(:,:,k);
            else
                QBOpred = [];
            end
            if ENSOdata
                ENSOpred = ENSOmonth(:,k);
            else
                ENSOpred = [];
            end
            if HFSdata
                HFpred = HFtousemonth(:,k);                       
            else
                HFpred = [];
            end    
            if aerosols
                aerpred = aermonth(:,k);                       
            else
                aerpred = [];
            end    
            if solardata
                solarpred = solarmonth(:,k);
            else
                solarpred = [];
            end
            if sum(NO2data(:)) ~= 0
                NO2pred = squeeze(NO2month(i,j,k,:));
            else
                NO2pred = [];
            end
            
            if chemonly                                                 
                if i >= 13 && sum(NO2data(:)) ~= 0                                        
                    predictors_month = [monthvector,NO2pred,solarpred];
                else
                    predictors_month = [monthvector,solarpred];
                end                                
            else
                if i >= 13 && sum(NO2data(:)) ~= 0                                        
                    predictors_month = [monthvector,QBOpred,ENSOpred,HFpred,aerpred,NO2pred,solarpred];
                else
                    predictors_month = [monthvector,QBOpred,ENSOpred,HFpred,aerpred,solarpred];
                end
            end       
            predindmonth = 1:size(predictors_month,2);                      
            
            [bmonth.ppmv(i,j,k,predindmonth),~,O3Anomaly.ppmv_residuals_months(i,j,k,:),~,stats(i,j,k,:)] = regress(O3Anomaly.ppmv(k:12:end,i,j),...
                predictors_month);
            
            [bmonth.percent(i,j,k,predindmonth),~,O3Anomaly.percent_residuals_months(i,j,k,:),~,stats_percent(i,j,k,:)] = regress(O3Anomaly.percent(k:12:end,i,j),...
                predictors_month);                       
            if i == 13
                a = 1;
            end
            if predindmonth == 1
                O3Anomaly.regmodel_months(i,j,k,:) = (predictors_month'.*squeeze(bmonth.ppmv(i,j,k,predindmonth)))';            
                O3Anomaly.percent_regmodel_months(i,j,k,:) = (predictors_month'.*squeeze(bmonth.percent(i,j,k,predindmonth)))';
            else
                O3Anomaly.regmodel_months(i,j,k,:) = (sum(predictors_month'.*squeeze(bmonth.ppmv(i,j,k,predindmonth))))';            
                O3Anomaly.percent_regmodel_months(i,j,k,:) = (sum(predictors_month'.*squeeze(bmonth.percent(i,j,k,predindmonth))))';
            end
            
        end
                        
        O3Anomaly.regmodel(:,i,j) = (sum(predictorstouse(:,:)'.*squeeze(b.ppmv(i,j,predind))))';
        O3Anomaly.percent_regmodel(:,i,j) = (sum(predictorstouse(:,:)'.*squeeze(b.percent(i,j,predind))))';
        
    end    
end

% create continuous time series
b.percent (b.percent == 0) = NaN;
O3Anomaly.regmodel_months = permute(O3Anomaly.regmodel_months(:,:,:),[3,1,2]);
O3Anomaly.percent_regmodel_months = permute(O3Anomaly.percent_regmodel_months(:,:,:),[3,1,2]);
O3Anomaly.ppmv_residuals_months = permute(O3Anomaly.ppmv_residuals_months(:,:,:),[3,1,2]);
O3Anomaly.percent_residuals_months = permute(O3Anomaly.percent_residuals_months(:,:,:),[3,1,2]);

%% AR2 model

if ARmodel
    if ARmodel == 2
        for i = 1:size(O3Anomaly.percent_residuals,2)
            for j = 1:size(O3Anomaly.percent_residuals,3)
                ARcoeffs(:,i,j) = regress(O3Anomaly.percent_residuals(3:end,i,j),...
                    [ones(length(O3Anomaly.percent_residuals(1:end-2,i,j)),1),...
                    O3Anomaly.percent_residuals(1:end-2,i,j),...
                    O3Anomaly.percent_residuals(2:end-1,i,j)]);
                %AR2coeffs(:,i,j) = aryule(O3Anomaly.percent_residuals(1:end-1,i,j),1);
%                 O3Anomaly.percent_residuals_months_AR2(:,i,j) = ...
%                     O3Anomaly.percent_residuals(:,i,j)*ARcoeffs(1,i,j)+...
%                     O3Anomaly.percent_residuals(:,i,j)*ARcoeffs(2,i,j)+...
%                     O3Anomaly.percent_residuals(:,i,j)*ARcoeffs(3,i,j);   
                
                O3Anomaly.percent_residuals_months_AR2(:,i,j) = ...
                    [O3Anomaly.percent_residuals_months(1:end-2,i,j)*ARcoeffs(1,i,j)+...
                    O3Anomaly.percent_residuals_months(1:end-2,i,j)*ARcoeffs(2,i,j)+...
                    O3Anomaly.percent_residuals_months(2:end-1,i,j)*ARcoeffs(3,i,j);O3Anomaly.percent_residuals_months(end-1:end,i,j)];
                
            end
        end
    elseif ARmodel == 1 %Y = bX(t-1) + e
        for i = 1:size(O3Anomaly.percent_residuals,2)
            for j = 1:size(O3Anomaly.percent_residuals,3)
                
                %dwp = dwtest(
                
                ARcoeffs(:,i,j) = regress(O3Anomaly.percent_residuals_months(2:end,i,j),...
                    [ones(length(O3Anomaly.percent_residuals_months(1:end-1,i,j)),1),...
                    O3Anomaly.percent_residuals_months(1:end-1,i,j)]);                    
                %AR2coeffs(:,i,j) = aryule(O3Anomaly.percent_residuals(1:end-1,i,j),1);
                O3Anomaly.percent_residuals_months_AR2(:,i,j) = ...
                    [O3Anomaly.percent_residuals_months(1:end-1,i,j)*ARcoeffs(1,i,j)+...
                    O3Anomaly.percent_residuals_months(1:end-1,i,j)*ARcoeffs(2,i,j);O3Anomaly.percent_residuals_months(end,i,j)];
                
                O3Anomaly.percent_residuals_months_ARtest(:,i,j) = ...
                    O3Anomaly.percent_residuals_months(:,i,j)*ARcoeffs(1,i,j)+...
                    O3Anomaly.percent_residuals_months(:,i,j)*ARcoeffs(2,i,j);
                
                O3Anomaly.percent_residuals_months_ARgone(:,i,j) = ...
                    [O3Anomaly.percent_residuals_months(2:end,i,j)- ... 
                    O3Anomaly.percent_residuals_months(1:end-1,i,j)*ARcoeffs(1,i,j)- ...
                    O3Anomaly.percent_residuals_months(1:end-1,i,j)*ARcoeffs(2,i,j);O3Anomaly.percent_residuals_months(end,i,j)];
                    
            end
        end
    end
   
end

if sum(SPEdata(:)) ~= 0 
    predictorsout = [predictorstouse];
else
    predictorsout = predictorstouse;
end

end
