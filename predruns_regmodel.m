function [] = predruns_regmodel(Eachyear,TSdataout,toz_dataMonthArrange,tozmonth,lons,lats,hem,inputtime)

%for i = 1:size(Eachyear.lowindex,2)
for j = 1:size(Eachyear.lowindex,3)         
             varlower(:,j,:,:) = TSdataout.dataVarMonthAve.highcl(Eachyear.lowindex(:,tozmonth,j),j,:,:);
             varupper(:,j,:,:) = TSdataout.dataVarMonthAve.highcl(Eachyear.highindex(:,tozmonth,j),j,:,:);
             tozcomp(:,j) = [squeeze(toz_dataMonthArrange.highcl(Eachyear.lowindex(:,tozmonth,j),tozmonth,j));...
                 squeeze(toz_dataMonthArrange.highcl(Eachyear.highindex(:,tozmonth,j),tozmonth,j))];                  

end
varcomp = permute(cat(1,varlower,varupper),[3,4,1,2]);
varcomp = varcomp(:,:,:);
tozcomp = tozcomp(:);

if hem
    tozdates = [1995,2024];
    cn = 7;
else
    tozdates = [1995,2023];
    cn = 6;
end
%tozdates = [1955,1974];%[1955,1975]
directory = ['/Volumes/MyBook/work/data/predruns/','toz','/','highCl','/'];
tozfiles = dir([directory,'*.nc']);

if hem
    [~,~,~,~,toz_rest] = predruns_ReadInlayer_areaaverage(directory,tozfiles,'toz',tozdates,[63,90],1);
        toz_rest2 = toz_rest(:,:,end-cn:end);
else
    [~,~,~,~,toz_rest] = predruns_ReadInlayer_areaaverage(directory,tozfiles,'toz',tozdates,[-90,-63],1);
        toz_rest2 = toz_rest(:,:,end-cn:end);
end

%% load in the entire time series
if hem
    TS9524 = load('/Volumes/MyBook/work/data/predruns/output/data/TS_ninoremoved_1995-202463-90.mat');
else
    TS9524 = load('/Volumes/MyBook/work/data/predruns/output/data/TS_ninoremoved_1995-202490-63.mat');
end

%% find upper and lower percentiles in ozone
for i = 1:9
    pctupper(i) = prctile(squeeze(toz_dataMonthArrange.highcl(i,tozmonth,:)),80);
    pctlower(i) = prctile(squeeze(toz_dataMonthArrange.highcl(i,tozmonth,:)),20);
    upperind(i).u = find(squeeze(toz_dataMonthArrange.highcl(i,tozmonth,:)) >= pctupper(i));
    lowerind(i).u = find(squeeze(toz_dataMonthArrange.highcl(i,tozmonth,:)) <= pctlower(i));
    
    TSdifferences(i,:,:) = squeeze(nanmean(TSdataout.dataVarMonthAve.highcl(i,upperind(i).u,:,:),2) - nanmean(TSdataout.dataVarMonthAve.highcl(i,lowerind(i).u,:,:),2));
    for j = 1:size(TSdataout.dataVarMonthAve.highcl,3)
        for k = 1:size(TSdataout.dataVarMonthAve.highcl,4)
            pdiff(i,j,k) = ttest2(squeeze(TSdataout.dataVarMonthAve.highcl(i,lowerind(i).u,j,k)),squeeze(TSdataout.dataVarMonthAve.highcl(i,upperind(i).u,j,k)));
        end
    end
end

pdiff (pdiff == 0) = -1;
pdiff (pdiff == 1) = 0;

TSdifferences = permute(TSdifferences,[1,3,2]);
%% plot TS differences
plot_TSdiff = 0;
if plot_TSdiff
    for i = 1:9
        mtitle = {['No. ',num2str(i),'{, }',num2str(inputtime(1)),'-',num2str(inputtime(2))]};   

        %mtitle = ['Esemble mean correlations of ',num2str(abs(lats2(1))),'-',num2str(abs(lats2(2))),hemext,' toz and ',var];

        subplotmaps(TSdifferences(i,:,:),lons,lats,{'div','RdBu'},1,[],16,mtitle,'Longitude','Latitude','Temperature difference (K)','on',...
            [-5 5],22,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/',...
            'correlations/maps/',sprintf('%02d',i),'_difference_','NH',num2str(inputtime(1)),'_',num2str(inputtime(2))];

        export_fig(filename,'-png');
    end
end

%% take regression over inputtime (may not be the entire time series)

% data has already been detrended and anomalized (ozone may not be
% detrended, and may need to be)

for l = 1:size(TSdataout.dataVarMonthAve.highcl,1)
    for i = 1:size(TSdataout.dataVarMonthAve.highcl,3)
        for j = 1:size(TSdataout.dataVarMonthAve.highcl,4)
            b(l,j,i,:) = regress(squeeze(TSdataout.dataVarMonthAve.highcl(l,:,i,j))',...
                [ones(length(squeeze(TSdataout.dataVarMonthAve.highcl(l,:,i,j))),1),squeeze(toz_dataMonthArrange.highcl(l,tozmonth,:))]);
            modelpred(l,j,i,:) = b(l,j,i,1) + b(l,j,i,2)*squeeze(toz_dataMonthArrange.highcl(l,tozmonth,:));
            modelpredfuture(l,j,i,:) = b(l,j,i,1) + b(l,j,i,2)*squeeze(toz_rest2(l,tozmonth,:));
%             b(l,j,i,:) = regress(squeeze(TS9524.dataVarMonthAve.highcl(l,1:21,i,j))',...
%                 [ones(length(squeeze(TSdataout.dataVarMonthAve.highcl(l,1:21,i,j))),1),squeeze(toz_rest(l,tozmonth,1:21))]);
        end
    end
end

%% taking rolling correlations for selection criteria (if used)
rollwindow = 12;
for i = 1:size(TSdataout.dataVarMonthAve.highcl,1)
    %for j = 1:size(TSdataout.dataVarMonthAve.highcl,3)
    %    for k = 1:size(TSdataout.dataVarMonthAve.highcl,4)
    TS = permute(squeeze(TSdataout.dataVarMonthAve.highcl(i,:,:,:)),[3,2,1]);
            [correlations(i)] = rollingCorrelations(TS,...
                squeeze(toz_dataMonthArrange.highcl(i,tozmonth,:)),rollwindow);
            
            
    for j = 1:size(TSdataout.dataVarMonthAve.highcl,3)
        for k = 1:size(TSdataout.dataVarMonthAve.highcl,4)
            [rall(i,k,j),rallpval(i,k,j)] = corr(squeeze(TSdataout.dataVarMonthAve.highcl(i,:,j,k))',...
                squeeze(toz_dataMonthArrange.highcl(i,tozmonth,:)));
            palltoplot = rallpval(i,:,:);
            palltoplot (palltoplot <= .05) = 0;
            palltoplot (palltoplot > .05) = -1;
        end
    end
            
    %% correlation
    plot_corr = 0;
    if plot_corr

        mtitle = {['No. ',num2str(i),'{, }',num2str(inputtime(1)),'-',num2str(inputtime(2))]};   

        subplotmaps(rall(i,:,:),lons,lats,{'div','RdBu'},1,palltoplot,16,mtitle,'Longitude','Latitude','Correlation','on',...
            [-1 1],22,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/',...
            'correlations/maps/',sprintf('%02d',i),'_corr_','NH',num2str(inputtime(1)),'-',num2str(inputtime(2))];

        export_fig(filename,'-png');
    end

    %% t

    if hem
        %test(:,:,i,:) = (abs(correlations(i).r) - 2.*std(correlations(i).r,0,3)).*nanmean(correlations(i).r(:,:,1:end-1),3);
        test(:,:,i,:) = (abs(correlations(i).r) - 2.*std(correlations(i).r,0,3));
    else
        test(:,:,i,:) = (abs(correlations(i).r) - 2.*std(correlations(i).r,0,3));
        %test(:,:,i,:) = (abs(correlations(i).r)).*nanmean(correlations(i).r(:,:,1:end-1),3);
    end
    
end

%% plot std of rall correlations and temperature differences
plotstd = 0;
if plotstd
    stdofcd(1,:,:) = squeeze(std(rall,0,1));
    stdofcd2(1,:,:) = squeeze(std(TSdifferences,0,1));

    subplotmaps(stdofcd(1,:,:),lons,lats,{'seq','YlOrBr'},1,[],16,{['Standard deviation of correlations, ',num2str(inputtime(1)),'-',num2str(inputtime(2))]},'Longitude','latitude','Correlation','on',...
        [0 .5],22,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/',...
                'correlations/maps/',sprintf('%02d',i),'_stdcorr_','NH',num2str(inputtime(1)),'-',num2str(inputtime(2))];

    export_fig(filename,'-png');

    subplotmaps(stdofcd2(1,:,:),lons,lats,{'seq','YlOrBr'},1,[],16,{['Standard deviation of temperature differences, ',num2str(inputtime(1)),'-',num2str(inputtime(2))]},'Longitude','latitude','Correlation','on',...
        [0 4],22,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/',...
                'correlations/maps/',sprintf('%02d',i),'_stdtempdiff_','NH',num2str(inputtime(1)),'-',num2str(inputtime(2))];

    export_fig(filename,'-png');

end
%%
test2 = test(:,:,:,end);
%test2 = test(:,:,:,end).*nanmean(test,4);
test2 = permute(test2,[3,1,2]);
%% plotting rolling correlations
plotroll = 0;
if plotroll
    if hem == 1
        xlims = [0 90];
    else
        xlims = [-90 0];    
    end
    for i = 1:size(test2,1)
        subplotmaps(test2(i,:,:),lons,lats,{'seq','YlOrBr'},0,[],16,{['No.',num2str(i),' HighCl']},'Longitude','latitude','','on',...
            [0 1],11,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],'',1,[0 360],xlims,0,'none',1,'Miller Cylindrical');

        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/',...
            'correlations/Rolling/predmetric/Ind_','No11.',sprintf('%02d',i),monthnames(tozmonth,0,0),'toz','_','TS','_detrend'];

        export_fig(filename,'-png');    
    end
end

%% Extracting point areas follow the criteria: 
% of correlation above .5 and temperature differences above 4 over three
% areas

 Areasforanom(1,:,:) = [0 120;45,80]; % Russia
Areasforanom(2,:,:) = [230 290;30 70]; % USA
Areasforanom(3,:,:) = [30 120;15 45]; % Central Asia
Areasforanom(4,:,:) = [300 330;55 80]; % Greenland
tempfutureextract = [];
tozfutureextract = [];
pastmean = squeeze(nanmean(TSdataout.dataMonthArrangeMean.highcl(:,[3,4],:,:,:),2));
pastmean = cat(2,pastmean,pastmean(:,1:8,:,:));
temptouse = TS9524.dataVarMonthAve.highcl + squeeze(nanmean(TS9524.dataMonthArrangeMean.highcl(:,[3,4],:,:,:),2)) - ...
    pastmean;
for i = 1:size(b,1)
    for j = 1:size(b,2)
        for k = 1:size(b,3)
            tempfuture = squeeze(TS9524.dataVarMonthAve.highcl(i,end-cn:end,k,j));            
            tozfuture = squeeze(toz_rest2(i,tozmonth,:));
            futurepred = b(i,j,k,1)+b(i,j,k,2)*tozfuture;            
        end
    end
    %tozfutureextract = [tozfutureextract;futurepred(ind)];
    %tempfutureextract = [tempfutureextract,tempfuture(ind)];
    
    % find best fit of criteria within ach area.    
    for l = size(Areasforanom,1)
        latind = lats >= Areasforanom(l,2,1) & lats <= Areasforanom(l,2,2);
        lonind = lons >= Areasforanom(l,1,1) & lons <= Areasforanom(l,1,2);
        if abs(b(i,lonind,latind)) >= .5 && abs(TSdifferences(i,lonind,latind)) >= 4
            ind = [find(tozfuture >= pctupper(i));find(tozfuture <= pctlower(i))];                
            test = b(i,lonind,latind,2).*TSdifferences(i,lonind,latind);
            
        end                 
    end     
end
tempfutureextract = tempfutureextract';
figure;
plot(tozfutureextract)
hold on
plot(tempfutureextract)

a = sign(tempfutureextract);
c = sign(tozfutureextract);
percent = length(find(a==c))/length(a)*100;

%%
if hem
    areas = [1,75,22;3,275,35;8,170,67;5,260,50;7,220,63;4,150,82];
else
    areas = [5,20,-36;6,17,-34;1,147,-30;8,335,-38;4,290,-40;4,205,-60];%;3,275,35;3,40,55;5,260,50;7,220,63;4,150,82];
end
%areas = [4,150,82];

rollwindow = 12;
figure
for i = 1:size(areas,1)
    [~, lon(i)] = min(abs(lons - areas(i,2)));
    [~, lat(i)] = min(abs(lats - areas(i,3)));
    prediction(:,i) = b(areas(i,1),lon(i),lat(i),1) + b(areas(i,1),lon(i),lat(i),2)*squeeze(toz_dataMonthArrange.highcl(areas(i,1),tozmonth,:));
    temperature(:,i) = squeeze(TSdataout.dataVarMonthAve.highcl(areas(i,1),:,lat(i),lon(i)));
    toz(:,i) = toz_dataMonthArrange.highcl(areas(i,1),tozmonth,:);
    predrest(:,i) = b(areas(i,1),lon(i),lat(i),1) + b(areas(i,1),lon(i),lat(i),2)*squeeze(toz_rest2(areas(i,1),tozmonth,:));
    %predrest2(:,i) = squeeze(modelpredfuture(areas(i,1),lon(i),lat(i),:));
    temprest(:,i) = squeeze(TS9524.dataVarMonthAve.highcl(areas(i,1),:,lat(i),lon(i)));
    %temprest2(:,i) = squeeze(tempresttotest2(areas(i,1),lon(i),lat(i),:));
    tozall(:,i) = toz_rest(areas(i,1),tozmonth,:);
    hold on
    %plot(squeeze(test(lon(i),lat(i),areas(i,1),:)))    
    plot(squeeze(correlations(areas(i,1)).r(lon(i),lat(i),:)));
%     corrlength = length(temperature(:,i)) - rollwindow;
%     count = 1;
%     for j = 1:corrlength        
%         [correlations.r(i,j),correlations.p(i,j)] = corr(temperature(count:count+rollwindow-1,i),toz(count:count+rollwindow-1,i));        
%         count = count+1;
%     end        
end
legend
%%
num = 1;
% figure
% plot([prediction(:,num)-nanmean(prediction(:,num))./std(prediction(:,num)),temperature(:,num)-nanmean(temperature(:,num))./std(prediction(:,num))]);
figure
plot([(predrest(:,num) - nanmean(predrest(:,num)))./std(predrest(:,num)),(temprest(end-cn:end,num) - nanmean(temprest(end-cn:end,num)))./std(prediction(:,num))])
legend('pred','temp');
% temprest2 = temprest(end-7:end,:);
% plot(correlations.r')
% 
% figure
% plot(correlations.r')
% legend
% figure
% plot(correlations.p');
% 
% stand = std(correlations.r');
% 
% figure;
% test = correlations.r.*(1-correlations.p)./stand';
% plot(test');


%% plotting
createfig('medium','on');
for i = 1:size(areas,1)
    subplot(3,2,i)
    plot(1995:2016,temperature(:,i),'LineWidth',2);
    hold on
    temp = (predrest(:,i)-nanmean(predrest(:,i)))./std(predrest(:,i))*std(temprest(1:end-8,i))+nanmean(predrest(:,i));
    %temp = predrest(:,i);
    %plot(2017:2023,[temprest(end-6:end,i),predrest(:,i)],'LineWidth',3)
    plot(2017:2024,[temprest(end-cn:end,i),temp],'LineWidth',3)
    grid on
    [mx] = max([temperature(:,i);temprest(end-cn:end,i)]);
    [mn] = min([temperature(:,i);temprest(end-cn:end,i)]);
    set(gca,'fontsize',16)
    xlabel('Year','fontsize',18)
    ylabel('Temperature (K)','fontsize',18)
    ylim([mn-.25 mx+.25])
    title(['Member ',num2str(areas(i,1)),', ',num2str(areas(i,3)),'N, ',num2str(areas(i,2)),'E']);
    if i == size(areas,1)
        lh = legend('TS used in regression model','Future temperature','Predicted temperature');
        set(lh,'position',[0 0 1 .05],'orientation','horizontal','box','off','fontsize',16)
    end       
end

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/Predictions/','South3'];
        
    export_fig(filename,'-png');    

%% Read in the rest of the ozone

for i = 4%1:size(areas,1)
    figure;
    plot([(predrest(:,i) - nanmean(predrest(:,i)))./std(predrest(:,i)),(temprest2(:,i)-nanmean(temprest2(:,i)))./std(temprest2(:,i))]);%./std(temprest2(:,i))
end

for i = 4%1:size(areas,1)
    figure;
    plot([prediction(:,i),temperature(:,i)]);
    hold on
    plot(temprest(:,i));
end


%%


%% 
end