% testing HSS using a leave-one out cross validation method on randomly
% sampled data
clear variables
%close all

% define length of predictor and predictant inputs
testlength = 30;
% define number of times to repeat method
notimes = 500;

% define whether to normalize to median or not
normtomed = 'no'; %'onlydata' 'no' 'all'

for j = 1:notimes
    % creating random data
    data = rand(testlength,1); predictor = rand(testlength,1)-1; r(j) = corr(data,predictor);

    %normalizing to mean
    %data = data - nanmean(data);
    %predictor = predictor - nanmean(predictor);

    for i = 1:numel(data)
        sd12 = 1:length(data);
        sd12 (sd12 == i) = [];

        % extracing training dat
        data_train = data(sd12)-nanmean(data(sd12));
        data_train_med(i) = median(data_train);
        predictor_train = predictor(sd12);

        %construct linear regression model
        b(i,:) = regress(data_train,[predictor_train,ones(size(predictor_train))]);

        % predictdata

        prediction(i) = b(i,1)*predictor(i) + b(i,2);
        actualdata(i) = data(i)-nanmean(nanmean(data(sd12)));

    end
    % calculating Heidke skill score
    if strcmp(normtomed,'all')
        actualdatasign = sign(actualdata-data_train_med);
        predictionsign =  sign(prediction-data_train_med);        
    elseif strcmp(normtomed,'onlydata')
        actualdatasign = sign(actualdata-data_train_med);
        predictionsign =  sign(prediction);
    elseif strcmp(normtomed,'no')
        actualdatasign = sign(actualdata);
        predictionsign =  sign(prediction);
    end
    
    
    forHSS = predictionsign;
    forHSS (forHSS ~= actualdatasign) = -2;
    forHSS (forHSS == actualdatasign) = 1;
    forHSS (forHSS == -2) = 0;

    HSStest(j) = (sum(forHSS,2) - size(forHSS,2)./2)./...
        (size(forHSS,2)-size(forHSS,2)./2)*100;  
end

HSSmean = nanmean(HSStest);

%% Sorting correlations by closest to zero
[sortvalues,sortindex] = sort(abs(r));

%% plotting 
fsize = 18;
lwidth = 2;
titles = {'HSS','Correlation'};
ylabels = {'HSS','r'};
MarkSize = 20:-1:11;

toplot = [HSStest;r]';

fig = figure;
set(fig,'color','white','position',[100 100 800 1000]);
for i = 1:2
    subplot(2,1,i)
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
    end
end

annotation('textbox',[0 1 1 0],'String',['Data test lenth = ',num2str(testlength),', Normalize to median = ',normtomed],...
            'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center',...
            'fontsize',fsize+6,'EdgeColor','none','fontweight','bold');   
        
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/HSStesting/'...
    'HSStest_TestLength_',num2str(testlength),'_Norm2Med_',normtomed,'.pdf'];
export_fig(filename,'-pdf');