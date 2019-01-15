function [TSout,forcomposite,TScomposite_lowmean,TScomposite_highmean,...
    TScompositedifference,longitude,latitude,data,alldata_lowind,alldata_highind,TScompositelow,...
    TScompositehigh,TScompositedifference_individual] = ...
    predruns_composites(files,TSfiles,directory,TSdirectory,timeperiod,ifdetrend,percentile,...
    lats,ozonemonth,TSmonth)

zonalmean_reshape = [];
TSdata_reshape_shift = [];
TScomposite_lowmean = [];
TScomposite_highmean = [];

count = 1;
count2 = 1;
for i = 1:length(files)
    %read in weighted area average data [height,time]
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    [~,TSdata(i),~] = Read_in_netcdf([TSdirectory,TSfiles(i).name]);
    %zonalmean(i,:) = nanmean(data(i).toz,1);
    if i == 1
        for j = 1:length(lats)
            [~,latind(j)] = min(abs(lats(j)-data(i).lat));
        end
    end
     % construct year only vector
    years(i).y = CCMI_years(TSdata(i).date);
    dateindfirst = find(years(i).y == timeperiod(1),1);
    dateindlast = find(years(i).y == timeperiod(2),1,'last');
        
    % for composite
    alldata(:,:,count:count+dateindlast-dateindfirst) = data(i).toz(:,:,dateindfirst:dateindlast);            
    % creating individual low and high indicies
    for j = 1:12
        
        data_zonalmean = squeeze(nanmean(data(i).toz(:,latind(1):latind(2),dateindfirst+j-1:12:dateindlast),1));
        data_allmean(i,j,:) = weightedaverage(data_zonalmean,data(1).lat(latind(1):latind(2)));           
        lowpercentile(i,j) = prctile(data_allmean(i,j,:),percentile);
        highpercentile(i,j) = prctile(data_allmean(i,j,:),100-percentile);
        lowind(i,j).l = find(data_allmean(i,j,:) <= lowpercentile(i,j));
        highind(i,j).h = find(data_allmean(i,j,:) >= highpercentile(i,j));  
        for k = 1:length(TSdata(i).lat)
            TSmean1 = nanmean(TSdata(i).TS(:,k,j:12:end),3);
            TSdata_rearrange_detrend(i,j,k,:,:) = detrend(squeeze(TSdata(i).TS(:,k,j:12:end))','linear')...
                + repmat(TSmean1,[1,size(TSdata(i).TS(:,k,j:12:end),3)])';            
        end
        
    end
    TScompositelow_individual(i,:,:) = nanmean(squeeze(TSdata_rearrange_detrend(i,12,:,lowind(i,10).l,:)),2);
    TScompositehigh_individual(i,:,:) = nanmean(squeeze(TSdata_rearrange_detrend(i,12,:,highind(i,10).h,:)),2);    

    TSalldata(:,:,count:count+dateindlast-dateindfirst) = TSdata(i).TS(:,:,dateindfirst:dateindlast);
     
    for j = 1:size(TSdata(i).TS(:,:,dateindfirst:12:dateindlast),3)
        TSalldata_summer(:,:,count2) = nanmean(TSdata(i).TS(:,:,dateindfirst+11+(j-1)*12:dateindfirst+11+(j-1)*12+2),3);
        count2 = count2+1;
    end
    count = count+size(data(i).toz(:,:,dateindfirst:dateindlast),3);        
end
    
TScompositedifference_individual = TScompositelow_individual - TScompositehigh_individual;

%weightedaverage
alldata_weightmean = weightedaverage(squeeze(nanmean(alldata(:,latind(1):latind(2),:),1)),data(1).lat(latind(1):latind(2)));

for i = 1:12
    alldata_rearrange(i,:) = alldata_weightmean(:,i:12:end);
    TSalldata_rearrange(i,:,:,:) = TSalldata(:,:,i:12:end);
    alldata_lowpercentile(i) = prctile(alldata_rearrange(i,:),percentile);
    alldata_highpercentile(i) = prctile(alldata_rearrange(i,:),100-percentile);
    alldata_lowind(i).a = find(alldata_rearrange(i,:) <= alldata_lowpercentile(i));
    alldata_highind(i).a = find(alldata_rearrange(i,:) >= alldata_highpercentile(i));  
end


% detrend
longitude = TSdata(1).lon;
latitude = TSdata(1).lat;

%for i = 1:latitude
if ifdetrend
    TSmean = squeeze(nanmean(squeeze(TSalldata_rearrange),4));
    bp = 1:size(TSalldata_rearrange,4)/length(files):length(files)*size(TSalldata_rearrange,4)/length(files);
    for i = 1:12;
        for j = 1:length(latitude)
            TS_detrend(i,:,j,:) = detrend(squeeze(TSalldata_rearrange(i,:,j,:))','linear',bp)+...
                repmat(squeeze(TSmean(i,:,j))',[1,size(TSalldata_rearrange,4)])';
        end
    end 
    TS_detrend = permute(TS_detrend,[1,4,3,2]);
    forcomposite = TS_detrend;
else
    forcomposite = TSalldata_rearrange;
end

% if alldata_lowind(ozonemonth).a(end) == size(alldata_rearrange,2)
%     alldata_lowind(ozonemonth).a(end) = [];
% end
% if alldata_highind(ozonemonth).a(end) == size(alldata_rearrange,2)
%     alldata_highind(ozonemonth).a(end) = [];
% end
%     alldata_lowind(ozonemonth).a(end) = [];


if strcmp(TSmonth,'Summer')
%     TScompositelow = nanmean(cat(3,squeeze(forcomposite(12,:,:,alldata_lowind(ozonemonth).a)),...
%         squeeze(forcomposite(1,:,:,alldata_lowind(ozonemonth).a+1)),...
%         squeeze(forcomposite(2,:,:,alldata_lowind(ozonemonth).a+1))),3);
% 
%     TScompositehigh = nanmean(cat(3,squeeze(forcomposite(12,:,:,alldata_highind(ozonemonth).a)),...
%         squeeze(forcomposite(1,:,:,alldata_highind(ozonemonth).a+1)),...
%         squeeze(forcomposite(2,:,:,alldata_highind(ozonemonth).a+1))),3);

    TScompositelow = nanmean(TSalldata_summer(:,:,alldata_lowind(ozonemonth).a),3);
    TScompositehigh = nanmean(TSalldata_summer(:,:,alldata_highind(ozonemonth).a),3);
    TScompositedifference = TScompositelow - TScompositehigh;
    TSout = TSalldata_summer;
else
    TSout = squeeze(TSalldata_rearrange(TSmonth,:,:,:));
    TScompositelow = nanmean([squeeze(forcomposite(TSmonth,:,:,alldata_lowind(ozonemonth).a))],3);
    TScompositehigh = nanmean(squeeze(forcomposite(TSmonth,:,:,alldata_highind(ozonemonth).a)),3);
    TScompositedifference = TScompositelow - TScompositehigh;
end


%         % rearranging
%         zonalmean_reshape(i,j,:) = zonalmean(i,dateindfirst+j-1:12:dateindlast);
%         TSdata_reshape(i,:,:,j,:) = TSdata(i).TS(:,:,dateindfirst+j-1:12:dateindlast);
% 
%         % detrend
%         if ifdetrend            
%             for k = 1:length(TSdata(i).lat)
%                 %for l = 1:length(TSdata(i).lon)
%                     TSmean(i,k,j,:) = nanmean(squeeze(TSdata_reshape(i,:,k,j,:)),2);
%                     TS_detrend(i,k,j,:,:) = detrend(squeeze(TSdata_reshape(i,:,k,j,:))','linear')+...
%                         repmat(squeeze(TSmean(i,k,j,:))',[size(TSdata_reshape,5),1]);
%                 %end
%             end
%         end
% 
%         % finding 90th percentiles
%         lowpercentile(i,j) = prctile(zonalmean_reshape(i,j,:),percentile);
%         highpercentile(i,j) = prctile(zonalmean_reshape(i,j,:),100-percentile);
%         lowind(i,j).l = find(zonalmean_reshape(i,j,:) <= lowpercentile(i,j));
%         highind(i,j).h = find(zonalmean_reshape(i,j,:) >= highpercentile(i,j));                       
%     end    
% end
% if ifdetrend
%     TSdata_reshape_shift = circshift(permute(TS_detrend,[1,5,2,3,4]),[0,0,0,2,0]);
% else
%     TSdata_reshape_shift = circshift(TSdata_reshape,[0,0,0,2,0]);
% end

%% creating composites 
% This is wrong. I need to take the January correlation in regards to
% October.

% for j = 1:12
%     countlow = 1;
%     counthigh = 1;
%     for i = 1:length(files)
%         TScomposite_low(j,:,:,countlow:countlow+length(lowind(i,j).l)-1) = squeeze(TSdata_reshape_shift(i,:,:,j,lowind(i,j).l));
%         TScomposite_high(j,:,:,counthigh:counthigh+length(highind(i,j).h)-1) = squeeze(TSdata_reshape_shift(i,:,:,j,highind(i,j).h));                
%         countlow = countlow+length(lowind(i,j).l);
%         counthigh = counthigh+length(highind(i,j).h);
%     end    
% end
% TScomposite_low (TScomposite_low == 0) = NaN;
% TScomposite_high (TScomposite_high == 0) = NaN;
% 
% TScomposite_lowmean = nanmean(TScomposite_low,4);
% TScomposite_highmean = nanmean(TScomposite_high,4);
% TScompositedifference = TScomposite_lowmean - TScomposite_highmean;


end
