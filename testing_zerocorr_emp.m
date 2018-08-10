% testing HSS using a leave-one out cross validation method on randomly
% sampled data
clear variables
%close all

% define length of predictor and predictant inputs
testlength = 30;
% define number of times to repeat method
notimes = 2000;

% define whether to normalize to median or not
normtomed = 'onlydata'; %'onlydata' 'no' 'all'

for j = 1:notimes
    % creating random data
    data = rand(testlength,1); predictor = (rand(testlength,1)-.5)./.5; r(j) = corr(data,predictor);
    predictorautocorr = predictor(2:end-1);
    %predictorautocorr = predictor;
    alldataupper = prctile(predictorautocorr,80);
    alldatalower = prctile(predictorautocorr,20);
    alldataupperind = find(predictorautocorr >= alldataupper);
    alldatalowerind = find(predictorautocorr <= alldatalower);
    upperlower = [alldataupperind',alldatalowerind'];
    %normalizing to mean
    %data = data - nanmean(data);
    %predictor = predictor - nanmean(predictor);
    count = 2;
    for i = 1:numel(data)-2
        sd12 = 1:length(data);
        %sd12 (sd12 == i ) = [];
        sd12 (sd12 == count | sd12 == count-1 | sd12 == count+1) = [];

        % extracing training dat
        data_train = data(sd12);%-nanmedian(data(sd12));
        data_train_med(i) = median(data(sd12));
        predictor_train = predictor(sd12);

        %construct linear regression model
        b(i,:) = regress(data_train,[predictor_train,ones(size(predictor_train))]);

        % predictdata
        upper_ind = find(predictor_train >= prctile(predictor_train,80));
        lower_ind = find(predictor_train <= prctile(predictor_train,20));
        %prediction(i) = b(i,1)*predictor(i) + b(i,2);
        %actualdata(i) = data(i)-nanmean(nanmean(data(sd12)));
        datadiff(i) = nanmean(data_train(upper_ind)) - nanmean(data_train(lower_ind));       
        %predictorleftsign(i) = sign(predictor(i));        
        predictorleftsign(i) = sign(predictor(count));        
        %actualdata(i) = data(i);
        actualdata(i) = data(count);
        count = count+1;
        
    end
    predictionsign = zeros(size(actualdata));
    predictionsign (predictorleftsign == 1 & sign(datadiff) == 1) = 1;
    predictionsign (predictorleftsign == -1 & sign(datadiff) == -1) = 1;
    predictionsign (predictorleftsign == -1 & sign(datadiff) == 1) = -1;
    predictionsign (predictorleftsign == 1 & sign(datadiff) == -1) = -1;
    
    % calculating Heidke skill score
    if strcmp(normtomed,'all')
        actualdatasign = sign(actualdata-data_train_med);
        %predictionsign =  sign(prediction-data_train_med);        
    elseif strcmp(normtomed,'onlydata')
        actualdatasign = sign(actualdata-data_train_med);
        %predictionsign =  sign(prediction);
    elseif strcmp(normtomed,'no')
        actualdatasign = sign(actualdata);
        %predictionsign =  sign(prediction);
    end
    
    
    forHSS = predictionsign;
    forHSS (forHSS ~= actualdatasign) = -2;
    forHSS (forHSS == actualdatasign) = 1;
    forHSS (forHSS == -2) = 0;

    HSStest(j) = (sum(forHSS,2) - size(forHSS,2)./2)./...
        (size(forHSS,2)-size(forHSS,2)./2)*100;  
    
    % extract upper and lower
    forHSSpct = forHSS(upperlower);
    
    HSStest_pct(j) = (sum(forHSSpct,2) - size(forHSSpct,2)./2)./...
        (size(forHSSpct,2)-size(forHSSpct,2)./2)*100;  
    
end

HSSmean = nanmean(HSStest);
HSSmean_pct = nanmean(HSStest_pct);

%% Sorting correlations by closest to zero
[sortvalues,sortindex] = sort(abs(r));

%% plotting 
fsize = 18;
lwidth = 2;
titles = {'HSS','HSS','Correlation'};
ylabels = {'HSS','HSS','r'};
MarkSize = 20:-1:11;

toplot = [HSStest;HSStest_pct;r]';

fig = figure;
set(fig,'color','white','position',[100 100 800 1000]);
for i = 1:3
    subplot(3,1,i)
    plot(toplot(:,i),'LineWidth',lwidth);
    hold on
    for j = 1:10
        ph(j) = plot(sortindex(j),toplot(sortindex(j),i),'o','MarkerSize',MarkSize(j),...
            'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor','k');
    end
    set(gca,'fontsize',fsize);    
    ylabel(ylabels{i},'fontsize',fsize+2);
    xlabel('Trial number','fontsize',fsize+2);
    title(titles{i},'fontsize',fsize+4);
    if i == 1        
        annotation('textbox',get(gca,'position'),'String',['Ave. HSS = ',num2str(HSSmean)],...
            'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left',...
            'fontsize',fsize+2,'EdgeColor','none','fontweight','bold');   
    elseif i == 2
        annotation('textbox',get(gca,'position'),'String',['Ave. HSS = ',num2str(HSSmean_pct)],...
            'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left',...
            'fontsize',fsize+2,'EdgeColor','none','fontweight','bold');   
    end
end

annotation('textbox',[0 1 1 0],'String',['Data test lenth = ',num2str(testlength),', Normalize to median = ',normtomed],...
            'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center',...
            'fontsize',fsize+6,'EdgeColor','none','fontweight','bold');   
        
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/HSStesting/'...
    'empHSStest_TestLength_',num2str(testlength),'_Norm2Med_',normtomed,'.pdf'];
export_fig(filename,'-pdf');