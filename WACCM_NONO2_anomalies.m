% NO+NO2 rate of decline anomaly

% produce the time series 2000-2015
% produce 2006 2007 2011 2015 anomalies.
% Look at 10 50 and 150 hPa levels
% Seasonal cycle anomalies

variable = 'NONO2';
pres = [150 50 10];
lats = [-90 -60];
runs = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};

%% Reading in WACCM data
[NumberDensity, MolConc, Pressure, Altitude, GMH, Latitudes, datedata] = ReadWACCMvertical(variable);

latindex = find(Latitudes > lats(1) & Latitudes < lats(2));

%% extracting latitudes and pressure levels

unipres = [1000:-50:150 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.05:.1];

count = 1;
for i = 1:length(runs) 
    for l = 1:size(NumberDensity.(runs{i}),3)
        for k = 1:size(NumberDensity.(runs{i}),1)
            vertical_uniform_ozone.(runs{i})(k,:,l) = interp1(log(squeeze(Pressure.(runs{i})(k,:,l))/100),...
                squeeze(NumberDensity.(runs{i})(k,:,l)),log(unipres));
        end        
        %rawdata(i,:,:) = squeeze(var20002014.(names{i})(rawlatindex,:,rawmonth:12:end)); 
    end            
    
    % vertical anolamy 
    if i > 1
        [~, yearmean(count,:,:)] = Verticalanomaly(vertical_uniform_ozone.(runs{i}),unipres,[-90 -60],...
            Latitudes,[2000 2014],[2000 2015],1:12,0);
        count = count+1;
    end
end

% autumn summer days (Mar-Aug)

%autsumdays = [60:243];

for i = 1:size(yearmean,1)
    NDlat(i).d = squeeze(nanmean(vertical_uniform_ozone.(runs{i})(latindex,:,:))); 
    for j = 1:length(pres)
        [~,presindex] = min(abs(pres(j)-unipres));
        NDlatpres(i,j,:) = yearmean(i,:,presindex);
        NDlatpresrate(i,j,:) = diff(yearmean(i,:,presindex));
    end
end

count = 2;
for i = 1:16
    NDlatpres_autsum(:,:,:,i) = NDlatpres(:,:,count:count+3);
    NDlatpres_autsum_rate(:,:,:,i) = NDlatpresrate(:,:,count:count+3);
    count = count+12;
end

%% plotting 
fsize = 18;
figure
fig = gcf;
set(fig,'color','white','position',[100 100 1000 1000]);
count = 1;
count1 = 1;
for i = 1:3;
        
%if i == 1 || i == 3 || i == 5
    subplot(3,1,i)
    ph = plot(squeeze(log(NDlatpres_autsum(1,count,:,[1:7,9:12,14:end]))),'color',[.6 .6 .6]);
    hold on 
    ph2 = plot(squeeze(log(NDlatpres_autsum(1,count,:,[8,13]))),'LineWidth',2);
    count = count+1;
% else
%     subplot(3,2,i)    
%     ph = plot(squeeze(NDlatpres_autsum_rate(1,count1,:,[1:7,9:12,14:end])),'color',[.6 .6 .6]);
%     hold on 
%     ph2 = plot(squeeze(NDlatpres_autsum_rate(1,count1,:,[8,13])),'LineWidth',2);
%     count1 = count1+1;
% end
set(gca,'xtick',1:1:4,'xticklabel',{'February','March','April','May'},'fontsize',fsize);
title([num2str(pres(i)),'hPa'],'fontsize',20);
end
lh = legend([ph(1);ph2],'All other years','2006','2011')%num2str([2006,2011]));
set(lh,'box','off');
mtit('NO+NO2','xoff',-.01,'yoff',.04,'fontsize',22);

export_fig('out.pdf','-pdf');
