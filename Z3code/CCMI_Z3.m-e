% Read in CESM-CCMI data and
clear all
variable = 'Z3';
hemisphere = 'south';
if strcmp(hemisphere,'south')
    ext = '6090S';
else
    ext = '6090N';
end
directory = ['/Volumes/ExternalOne/work/data/CESM-CCMI/',variable,'/zonalmean/',ext,'/'];
%histdirectory = ['/Volumes/My Book for Mac/work/data/CESM-CCMI/',variable,'/historical/zonalmean/'];
%histSDdirectory = ['/Volumes/My Book for Mac/work/data/CESM-CCMI/',variable,'/historicalSD/zonalmean/'];
write_netcdf = 0;
line_plots = 0;
contour_plots = 0;
files = dir([directory,'*.nc*']);
% histfiles = dir([histdirectory,'*.nc']);
% histSDfiles = dir([histSDdirectory,'*.nc']);
% ensfiles = dir([directory,'ensmean/','*.nc*']);
separate_lineplots = 1;
%% Reading in ERA-Interim
if strcmp(variable,'T');
    ERAdata = ncread(['/Volumes/ExternalOne/work/data/ERA-Interim/',variable,'/',variable,'_ERA-Interim.nc'],'t');
elseif strcmp(variable,'Z3');
    ERAdata = ncread(['/Volumes/ExternalOne/work/data/ERA-Interim/',variable,'/',variable,'_ERA-Interim.nc'],'z')./9.80665;
end
ERApressure = ncread('/Volumes/ExternalOne/work/data/ERA-Interim/Z3/Z3_ERA-Interim.nc','level');
ERAlatitude = ncread('/Volumes/ExternalOne/work/data/ERA-Interim/Z3/Z3_ERA-Interim.nc','latitude');
ERAzonalmean = squeeze(nanmean(ERAdata(:,101:end,:,:),1));

for j = 1:size(ERAzonalmean,2)
    Z3weighted_ERA(j,:) = weightedaverage(squeeze(ERAzonalmean(:,j,:)),ERAlatitude(101:end));        
end

for j = 1:12
    Z3_ERA_rearrange(j,:,:) = Z3weighted_ERA(:,j:12:240);
    Z3_ERA_rearrange2_1(j,:,:) = Z3weighted_ERA(:,j+240:12:276);
    Z3_ERA_rearrange2_2(j,:,:) = Z3weighted_ERA(:,j+288:12:456);
    Z3_ERA_rearrange2 = cat(3,Z3_ERA_rearrange2_1,Z3_ERA_rearrange2_2);
    %Z3_ERA_rearrange = circshift(Z3_ERA_rearrange,[6,0,0]);
    for k = 1:size(Z3_ERA_rearrange,2)
        [bERA(j,k,:),bintERA(j,k,:,:)] = ... %(month,pressure,b)
            regress(squeeze(Z3_ERA_rearrange(j,k,:)),[ones(1,size(Z3_ERA_rearrange(j,k,:),3));1:size(Z3_ERA_rearrange(j,:,:),3)]');                            
        [bERA2(j,k,:),bintERA2(j,k,:,:)] = ... %(month,pressure,b)
            regress(squeeze(Z3_ERA_rearrange2(j,k,:)),[ones(1,size(Z3_ERA_rearrange2(j,k,:),3));1:size(Z3_ERA_rearrange2(j,:,:),3)]');                            
    end        
end

%%

for i = 1:length(files);    
    [data(i).attributes, data(i).data,data(i).info] = Read_in_netcdf([directory,files(i).name]);
    data(i).data.(variable) = squeeze(data(i).data.(variable));
    data(i).data.date = squeeze(data(i).data.date(1,1,:));
    
    for j = 1:size(data(i).data.(variable),2)
        Z3weighted(i).wa(j,:) = weightedaverage(squeeze(data(i).data.(variable)(:,j,:)),data(i).data.lat);        
    end
    Pressureweighted(i).wa = weightedaverage(squeeze(data(i).data.PS),data(i).data.lat);
    Pressureweighted(i).wa = repmat(data(i).data.ap.*100000,1,size(Pressureweighted(i).wa,2)) + ...
        repmat(data(i).data.b,1,size(Pressureweighted(i).wa,2)) ...
        .* repmat(Pressureweighted(i).wa,length(data(i).data.lev),1);        
    Pressureweighted(i).wa = Pressureweighted(i).wa./100;
    years(i).y = num2str(data(i).data.date);
    for j = 1:length(years(i).y)
        yearsfinal(i).y(j,:) = str2double(years(i).y(j,1:4));
    end
    yearsfinal(i).y = circshift(yearsfinal(i).y,1);
    yearsfinal(i).y(1) = yearsfinal(i).y(2); 
end

plot(yearsfinal(4).y(12:12:end),squeeze(data(7).data.(variable)(5,49,12:12:end)))
hold on
plot(yearsfinal(4).y(12:12:end),squeeze(data(8).data.(variable)(5,49,12:12:end)))
plot(yearsfinal(4).y(12:12:end),squeeze(data(9).data.(variable)(5,49,12:12:end)))

for i = 1:length(histfiles)
    [histdata(i).attributes, histdata(i).histdata,histdata(i).info] = Read_in_netcdf([histdirectory,histfiles(i).name]);
    histdata(i).histdata.(variable) = squeeze(histdata(i).histdata.(variable));
    histdata(i).histdata.date = squeeze(histdata(i).histdata.date(1,1,:));
    
    for j = 1:size(histdata(i).histdata.(variable),2)
        histZ3weighted(i).wa(j,:) = weightedaverage(squeeze(histdata(i).histdata.(variable)(:,j,:)),histdata(i).histdata.lat);        
    end
    histPressureweighted(i).wa = weightedaverage(squeeze(histdata(i).histdata.PS),histdata(i).histdata.lat);
    histPressureweighted(i).wa = repmat(histdata(i).histdata.ap.*100000,1,size(histPressureweighted(i).wa,2)) + ...
        repmat(histdata(i).histdata.b,1,size(histPressureweighted(i).wa,2)) ...
        .* repmat(histPressureweighted(i).wa,length(histdata(i).histdata.lev),1);        
    histPressureweighted(i).wa = histPressureweighted(i).wa./100;
    histyears(i).y = num2str(histdata(i).histdata.date);
    for j = 1:length(histyears(i).y)
        histyearsfinal(i).y(j,:) = str2double(histyears(i).y(j,1:4));
    end
    histyearsfinal(i).y = circshift(histyearsfinal(i).y,1);
    histyearsfinal(i).y(1) = histyearsfinal(i).y(2); 
end

[SDhistdata(1).attributes, SDhistdata(1).histdata,SDhistdata(1).info] = Read_in_netcdf([histSDdirectory,histSDfiles.name]);
SDhistdata.histdata.(variable) = squeeze(SDhistdata.histdata.(variable));
%SDhistdata.histdata.date = squeeze(SDhistdata.histdata.date(1,1,:));
    
for j = 1:size(SDhistdata.histdata.(variable),2)
    SDhistZ3weighted.wa(j,:) = weightedaverage(squeeze(SDhistdata.histdata.(variable)(:,j,:)),SDhistdata.histdata.lat);        
end

[~,PSdata_temp,~,~] = ReadinCESMDateandPS(0,1);
SDhistdata.histdata.PS = squeeze(nanmean(PSdata_temp(18).p.PS(:,1:16,:))); 
SDhistPressureweighted.wa = weightedaverage(squeeze(SDhistdata.histdata.PS),SDhistdata.histdata.lat);
SDhistPressureweighted.wa = repmat(SDhistdata.histdata.ap.*100000,1,size(SDhistPressureweighted.wa,2)) + ...
    repmat(SDhistdata.histdata.b,1,size(SDhistPressureweighted.wa,2)) ...
    .* repmat(SDhistPressureweighted.wa,length(SDhistdata.histdata.lev),1);        
SDhistPressureweighted.wa = SDhistPressureweighted.wa./100;
% SDhistyears.y = num2str(SDhistdata.histdata.date);
% for j = 1:length(SDhistyears.y)
%     SDhistyearsfinal.y(j,:) = str2double(SDhistyears.y(j,1:4));
% end
%SDhistyearsfinal.y = circshift(SDhistyearsfinal.y,1);
SDhistyearsfinal.y = histyearsfinal(1).y(289:end);
%%
% for i = 1:length(ensfiles);
%     [ensdata(i).attributes,ensdata(i).data,ensdata(i).info] = Read_in_netcdf([directory,'ensmean/',ensfiles(i).name]);
%     ensdata(i).data.Z3 = squeeze(ensdata(i).data.Z3);
%     ensdata(i).data.date = squeeze(ensdata(i).data.date(1,1,:));
%     
%     for j = 1:size(ensdata(i).data.Z3,2)
%         ensZ3weighted(i).wa(j,:) = weightedaverage(squeeze(ensdata(i).data.Z3(:,j,:)),ensdata(i).data.lat);
%     end
%     ensPressureweighted(i).wa = weightedaverage(squeeze(data(i).data.PS),data(i).data.lat);
%     ensPressureweighted(i).wa = repmat(data(i).data.ap.*100000,1,size(ensPressureweighted(i).wa,2)) + ...
%         repmat(data(i).data.b,1,size(ensPressureweighted(i).wa,2)) ...
%         .* repmat(ensPressureweighted(i).wa,length(data(i).data.lev),1);        
%     ensPressureweighted(i).wa = ensPressureweighted(i).wa./100;
%     
% end

%% take mirror trends and plot

pres = [1000:-50:150 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.05:.15 ... 
    .1:-.005:.015 ...
    .01:-.0005:.0015 ...
    .001:-.00005:.00015];

% first put onto regular grid by interpolating pressure

for i = 1:length(Z3weighted)
    for j = 1:size(Z3weighted(i).wa,2)
        %Z3weightedregp(i).wa(:,j)  = interp1(log(Pressureweighted(i).wa(:,j)),Z3weighted(i).wa(:,j),log(data(i).data.lev),'linear','extrap');
        Z3weightedregp(i).wa(:,j)  = interp1(log(Pressureweighted(i).wa(:,j)),Z3weighted(i).wa(:,j),log(double(ERApressure)),'linear','extrap');
    end
   
end

for i = 1:length(histZ3weighted)
    for j = 1:size(histZ3weighted(i).wa,2)
        %Z3weightedregp(i).wa(:,j)  = interp1(log(Pressureweighted(i).wa(:,j)),Z3weighted(i).wa(:,j),log(data(i).data.lev),'linear','extrap');
        histZ3weightedregp(i).wa(:,j)  = interp1(log(histPressureweighted(i).wa(:,j)),histZ3weighted(i).wa(:,j),log(double(ERApressure)),'linear','extrap');
    end
   
end

% for i = 1:length(ensZ3weighted)
%     for j = 1:size(ensZ3weighted(i).wa,2)
%         ensZ3weightedregp(i).wa(:,j)  = interp1(log(ensPressureweighted(i).wa(:,j)),ensZ3weighted(i).wa(:,j),log(double(ERApressure)),'linear','extrap');
%     end
% end



%% take linear trends over apropriate times
for i = 1:length(Z3weighted)-3
    dateind_1(i,:) = find(data(i).data.date >= 19790201 & data(i).data.date <= 19990101);    
    dateind_2_temp = find(data(i).data.date >= 19990201 & data(i).data.date <= 20150101);    
    dateind_2_temp (data(i).data.date(dateind_2_temp) >= 20020201 & data(i).data.date(dateind_2_temp) <= 20030101) = [];
    dateind_2(i,:) = dateind_2_temp;
    Z3_dateextract_1(i).wa = Z3weightedregp(i).wa(:,dateind_1(i,:));
    Z3_dateextract_2(i).wa = Z3weightedregp(i).wa(:,dateind_2(i,:));
            
    for j = 1:12
        Z3_dateextract_1_rearrange(i).wa(j,:,:) = Z3_dateextract_1(i).wa(:,j:12:end);
        Z3_dateextract_2_rearrange(i).wa(j,:,:) = Z3_dateextract_2(i).wa(:,j:12:end);
        
        for k = 1:size(Z3_dateextract_1_rearrange(i).wa,2)
            [bWAC1(i,j,k,:),bintWAC1(i,j,k,:,:)] = ... %(run,month,pressure,b)
                regress(squeeze(Z3_dateextract_1_rearrange(i).wa(j,k,:)),[ones(1,size(Z3_dateextract_1_rearrange(i).wa(j,k,:),3));1:size(Z3_dateextract_1_rearrange(i).wa(j,:,:),3)]');                
            [bWAC2(i,j,k,:),bintWAC2(i,j,k,:,:)] = ... %(run,month,pressure,b)
                regress(squeeze(Z3_dateextract_2_rearrange(i).wa(j,k,:)),[ones(1,size(Z3_dateextract_2_rearrange(i).wa(j,k,:),3));1:size(Z3_dateextract_2_rearrange(i).wa(j,:,:),3)]');                
        end        
    end
       
end

for i = 1:length(histZ3weighted)
    histdateind_1(i,:) = find(histdata(i).histdata.date >= 19790201 & histdata(i).histdata.date <= 19990101);    
    histdateind_2_alltemp(i,:) = find(histdata(i).histdata.date >= 19990201 & histdata(i).histdata.date <= 20150101);    
    histdateind_2_temp = histdateind_2_alltemp(i,:);
    histdateind_2_temp (histdata(i).histdata.date(histdateind_2_temp) >= 20020201 & histdata(i).histdata.date(histdateind_2_temp) <= 20030101) = [];
    histdateind_2(i,:) = histdateind_2_temp;
    histZ3_dateextract_1(i).wa = histZ3weightedregp(i).wa(:,histdateind_1(i,:));
    histZ3_dateextract_2(i).wa = histZ3weightedregp(i).wa(:,histdateind_2(i,:));
            
    for j = 1:12
        histZ3_dateextract_1_rearrange(i).wa(j,:,:) = histZ3_dateextract_1(i).wa(:,j:12:end);
        histZ3_dateextract_2_rearrange(i).wa(j,:,:) = histZ3_dateextract_2(i).wa(:,j:12:end);
        
        for k = 1:size(histZ3_dateextract_1_rearrange(i).wa,2)
            [histbWAC1(i,j,k,:),histbintWAC1(i,j,k,:,:)] = ... %(run,month,pressure,b)
                regress(squeeze(histZ3_dateextract_1_rearrange(i).wa(j,k,:)),[ones(1,size(histZ3_dateextract_1_rearrange(i).wa(j,k,:),3));1:size(histZ3_dateextract_1_rearrange(i).wa(j,:,:),3)]');                
            [histbWAC2(i,j,k,:),histbintWAC2(i,j,k,:,:)] = ... %(run,month,pressure,b)
                regress(squeeze(histZ3_dateextract_2_rearrange(i).wa(j,k,:)),[ones(1,size(histZ3_dateextract_2_rearrange(i).wa(j,k,:),3));1:size(histZ3_dateextract_2_rearrange(i).wa(j,:,:),3)]');                
        end        
    end
       
end


% for i = 1:length(ensZ3weighted)-1
%     dateind_1(i,:) = find(data(i).data.date >= 19790201 & data(i).data.date <= 19990101);    
%     dateind_2_temp = find(data(i).data.date >= 19990201 & data(i).data.date <= 20150101);    
%     dateind_2_temp (data(i).data.date(dateind_2_temp) >= 20020201 & data(i).data.date(dateind_2_temp) <= 20030101) = [];
%     dateind_2(i,:) = dateind_2_temp;
%     ensZ3_dateextract_1(i).wa = ensZ3weightedregp(i).wa(:,dateind_1(i,:));
%     ensZ3_dateextract_2(i).wa = ensZ3weightedregp(i).wa(:,dateind_2(i,:));
%             
%     for j = 1:12
%         ensZ3_dateextract_1_rearrange(i).wa(j,:,:) = ensZ3_dateextract_1(i).wa(:,j:12:end);
%         ensZ3_dateextract_2_rearrange(i).wa(j,:,:) = ensZ3_dateextract_2(i).wa(:,j:12:end);
%         
%         for k = 1:size(Z3_dateextract_1_rearrange(i).wa,2)
%             [ensbWAC1(i,j,k,:),ensbintWAC1(i,j,k,:,:)] = ... %(run,month,pressure,b)
%                 regress(squeeze(ensZ3_dateextract_1_rearrange(i).wa(j,k,:)),[ones(1,size(ensZ3_dateextract_1_rearrange(i).wa(j,k,:),3));1:size(ensZ3_dateextract_1_rearrange(i).wa(j,:,:),3)]');                
%             [ensbWAC2(i,j,k,:),ensbintWAC2(i,j,k,:,:)] = ... %(run,month,pressure,b)
%                 regress(squeeze(ensZ3_dateextract_2_rearrange(i).wa(j,k,:)),[ones(1,size(ensZ3_dateextract_2_rearrange(i).wa(j,k,:),3));1:size(ensZ3_dateextract_2_rearrange(i).wa(j,:,:),3)]');                
%         end        
%     end
%        
% end

bERAtoplot = circshift(flip(bERA,2),[7,0,0])*10;
bERAtoplot2 = circshift(flip(bERA2,2),[7,0,0])*10;
% plot(squeeze(Z3_dateextract_1_rearrange(1).wa(12,10,:)),'k')
% hold on
% plot(squeeze(Z3_dateextract_1_rearrange(2).wa(12,10,:)),'k')
% plot(squeeze(Z3_dateextract_1_rearrange(3).wa(12,10,:)),'k')
% 
% plot(squeeze(Z3_ERA_rearrange(12,10,:))./9.80665)

%% saving output
lattosave = histdata(1).histdata.lat;
REFC1 = histZ3weightedregp;
yearsREFC1 = histyearsfinal.y;
REFC2(1).wa = Z3weightedregp(1).wa;
REFC2(2).wa = Z3weightedregp(2).wa;
REFC2(3).wa = Z3weightedregp(3).wa;
yearsREFC2 = yearsfinal(1:3).y;
GHG(1).wa = Z3weightedregp(4).wa;
GHG(2).wa = Z3weightedregp(5).wa;
GHG(3).wa = Z3weightedregp(6).wa;
yearsGHG = yearsfinal(4:6).y;
ODS(1).wa = Z3weightedregp(7).wa;
ODS(2).wa = Z3weightedregp(8).wa;
ODS(3).wa = Z3weightedregp(9).wa;
yearsODS = yearsfinal(7:9).y;
%ODS2000(1).wa = Z3weightedregp(10).wa;
%ODS2000(2).wa = Z3weightedregp(11).wa;
%ODS2000(3).wa = Z3weightedregp(12).wa;
%yearsODS2000 = yearsfinal(10:12).y;
ERAInterim = Z3weighted_ERA;
count = 1;
yearcount = 1979;
for i = 1:floor(size(Z3weighted_ERA,2)/12)
    yearsERA(count:count+11) = repmat(yearcount,12,1);
    yearcount = yearcount+1;
    count = count+12;
end
yearsERA = [yearsERA,yearsERA(end)+1];
save(['/Volumes/My Book for Mac/work/data/CESM-CCMI/',variable,'/output/REFC1_6090S_ERApressure.mat'],'REFC1','ERApressure','yearsREFC1');
save(['/Volumes/My Book for Mac/work/data/CESM-CCMI/',variable,'/output/REFC2_6090S_ERApressure.mat'],'REFC2','ERApressure','yearsREFC2');
save(['/Volumes/My Book for Mac/work/data/CESM-CCMI/',variable,'/output/SENC2fGHG_6090S_ERApressure.mat'],'GHG','ERApressure','yearsGHG');
save(['/Volumes/My Book for Mac/work/data/CESM-CCMI/',variable,'/output/SENC2fODS_6090S_ERApressure.mat'],'ODS','ERApressure','yearsODS');
%save(['/Volumes/My Book for Mac/work/data/CESM-CCMI/',variable,'/output/SENC2fODS2000_6090S_ERApressure.mat'],'ODS2000','ERApressure','yearsODS2000');
save(['/Volumes/My Book for Mac/work/data/CESM-CCMI/',variable,'/output/ERA-Interim_6090S_ERApressure.mat'],'ERAInterim','ERApressure','yearsERA');
%% plotting history
if contour_plots
    
    titles = {'REF-C2 - no.1', 'REF-C2 - no.2', 'REF-C2 - no.3','fGHG-1960 - no.1',...
        'fGHG-1960 - no.2', 'fGHG-1960 - no.3',...
        'fODS-1960 - no.1','fODS-1960 - no.2', 'fODS-1960 - no.3'};
    plotmtitle = 'GPH trends, 1979-1998';
    bWAC1toplot = circshift(flip(bWAC1(:,:,:,2),3),[0,7,0])*10;

    pressure = flip(double(ERApressure));
    prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];
        presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],5,[],[],2,1};
        logprestick = log(prestick);

    subplotmaps(bWAC1toplot,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','m/decade','on',...
        [-400 400],22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);

    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/GPH_allmembers_1979-1998.pdf'],'-pdf');

    %% plotting history REF-C1
    titles = {'REF-C1 - no.1', 'REF-C1 - no.2', 'REF-C1 - no.3','REF-C1 - no.4','REF-C1 - no.5','ERA-Interim'};
    plotmtitle = 'GPH trends, 1999-2014';
    histbWAC1toplot = circshift(flip(histbWAC1(:,:,:,2),3),[0,7,0])*10;

    histbWAC1toplot = permute(cat(3,squeeze(histbWAC1toplot(1,:,:)),squeeze(histbWAC1toplot(2,:,:))...
        ,squeeze(histbWAC1toplot(3,:,:)),squeeze(histbWAC1toplot(4,:,:)),squeeze(histbWAC1toplot(5,:,:)),bERAtoplot(:,:,2)),[3,1,2]);

    subplotmaps(histbWAC1toplot,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','m/decade','on',...
        [-400 400],22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);

    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/GPH_REF-C1_allmembers_1979-1998.pdf'],'-pdf');

    %% plotting future REF-C1
    titles = {'REF-C1 - no.1', 'REF-C1 - no.2', 'REF-C1 - no.3','REF-C1 - no.4','REF-C1 - no.5','ERA-Interim'};
    plotmtitle = 'GPH trends, 1979-1998';
    histbWAC1toplot = circshift(flip(histbWAC2(:,:,:,2),3),[0,7,0])*10;

    histbWAC1toplot = permute(cat(3,squeeze(histbWAC1toplot(1,:,:)),squeeze(histbWAC1toplot(2,:,:))...
        ,squeeze(histbWAC1toplot(3,:,:)),squeeze(histbWAC1toplot(4,:,:)),squeeze(histbWAC1toplot(5,:,:)),bERAtoplot2(:,:,2)),[3,1,2]);

    subplotmaps(histbWAC1toplot,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','m/decade','on',...
        [-400 400],22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);

    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/GPH_REF-C1_allmembers_1999-2014.pdf'],'-pdf');

    %% plotting future
    bWAC2toplot = circshift(flip(bWAC2(:,:,:,2),3),[0,7,0])*10;
    plotmtitle = 'GPH trends, 1999-2014';
    subplotmaps(bWAC2toplot,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','m/decade','on',...
        [-400 400],22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);
    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/GPH_allmembers_1999-2016.pdf'],'-pdf');
    %% hist ensemble
    titles = {'ERA-Interim', 'REF-C2', 'SEN-C2-fGHG','SEN-C2-fODS'};    
    plotmtitle = 'Ensemble average GPH trends, 1979-1998';
    bWAC1toplot = circshift(flip(bWAC1(:,:,:,2),3),[0,7,0])*10;
    bWAC1toplot_average(1,:,:) = nanmean(bWAC1toplot(1:3,:,:),1);  
    bWAC1toplot_average(2,:,:) = nanmean(bWAC1toplot(4:6,:,:),1);  
    bWAC1toplot_average(3,:,:) = nanmean(bWAC1toplot(7:9,:,:),1);  

    ensbWAC1toplot = permute(cat(3,squeeze(bWAC1toplot_average(1,:,:)),squeeze(bWAC1toplot_average(2,:,:))...
        ,squeeze(bWAC1toplot_average(3,:,:)),bERAtoplot(:,:,2)),[3,1,2]);

    ensbWAC1toplot = circshift(ensbWAC1toplot,[1,0,0]); 

    pressure = flip(double(ERApressure));
    prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];
        presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],5,[],[],2,1};
        logprestick = log(prestick);

    subplotmaps(ensbWAC1toplot,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','m/decade','on',...
        [-400 400],22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);
    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/GPH_ensemble_1979-1998.pdf'],'-pdf');
    %% future ensemble
    titles = {'ERA-Interim', 'REF-C2', 'SEN-C2-fGHG','SEN-C2-fODS'};    

    plotmtitle = 'Ensemble average GPH trends, 1999-2014';
    bWAC2toplot = circshift(flip(bWAC2(:,:,:,2),3),[0,7,0])*10;
    bWAC2toplot_average(1,:,:) = nanmean(bWAC2toplot(1:3,:,:),1);  
    bWAC2toplot_average(2,:,:) = nanmean(bWAC2toplot(4:6,:,:),1);  
    bWAC2toplot_average(3,:,:) = nanmean(bWAC2toplot(7:9,:,:),1);  

    ensbWAC1toplot2 = permute(cat(3,squeeze(bWAC2toplot_average(1,:,:)),squeeze(bWAC2toplot_average(2,:,:))...
        ,squeeze(bWAC2toplot_average(3,:,:)),bERAtoplot2(:,:,2)),[3,1,2]);
    ensbWAC1toplot2 = circshift(ensbWAC1toplot2,[1,0,0]); 

    pressure = flip(double(ERApressure));
    prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];
        presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],5,[],[],2,1};
        logprestick = log(prestick);

    subplotmaps(ensbWAC1toplot2,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','m/decade','on',...
        [-200 200],22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);
    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/GPH_ensemble_1999-2016.pdf'],'-pdf');

    %% producing line plots
    
end
monthnames = {'January','February','March','April','May','June','July','August','September',...
        'October','November','December'};
%%
fsize = 18;
if separate_lineplots
    % historical ensembles (REF-C1
    colors2 = cbrewer('qual','Set1',5);
    mon = 10;    
    count = 1;
    for i = 1:6
        fig1 = figure;
        set(fig1,'color','white','position',[100 100 700 800]);   
        pressure = 2;    
        [~,presind] = min(abs(ERApressure - pressure));
        if i <= 5
            for j = 1:2
                s(i,j) = subplot(2,1,j);
                if j == 1
                    plot(histyearsfinal(i).y(mon:12:end),histZ3weightedregp(i).wa(presind,mon:12:end),'LineWidth',2)
                    hold on
                    plot(histyearsfinal(i).y(histdateind_1(i,mon:12:end)),squeeze(histbWAC1(i,mon,presind,2)*[1:length(histyearsfinal(i).y(histdateind_1(i,mon:12:end)))]+histbWAC1(i,mon,presind,1)),'LineWidth',3)
                    plot(histyearsfinal(i).y(histdateind_2_alltemp(i,mon:12:end)),squeeze(histbWAC2(i,mon,presind,2)*[1:length(histyearsfinal(i).y(histdateind_2_alltemp(i,mon:12:end)))]+histbWAC2(i,mon,presind,1)),'LineWidth',3)
                    xlabel('year','fontsize',fsize+2)
                    ylabel('Geopotential height (m)','fontsize',fsize+2);                    
                    set(gca,'fontsize',fsize);                    
                    title([monthnames{mon},' zonal mean at ',num2str(pressure),' hPa'],'fontsize',fsize);
                end                
                if j == 2;
                    pressure = 100;    
                    [~,presind] = min(abs(ERApressure - pressure));                    
                    plot(histyearsfinal(i).y(mon:12:end),histZ3weightedregp(i).wa(presind,mon:12:end),'LineWidth',2)
                    hold on
                    plot(histyearsfinal(i).y(histdateind_1(i,mon:12:end)),squeeze(histbWAC1(i,mon,presind,2)*[1:length(histyearsfinal(i).y(histdateind_1(i,mon:12:end)))]+histbWAC1(i,mon,presind,1)),'LineWidth',3)
                    plot(histyearsfinal(i).y(histdateind_2_alltemp(i,mon:12:end)),squeeze(histbWAC2(i,mon,presind,2)*[1:length(histyearsfinal(i).y(histdateind_2_alltemp(i,mon:12:end)))]+histbWAC2(i,mon,presind,1)),'LineWidth',3)
                    xlabel('year','fontsize',fsize+2)
                    ylabel('Geopotential height (m)','fontsize',fsize+2);                    
                    set(gca,'fontsize',fsize);                    
                    title([monthnames{mon},' zonal mean at ',num2str(pressure),' hPa'],'fontsize',fsize);
                    
                end                                
            end 
            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/Z3_',num2str(i),'_',num2str(mon),'_','2and100hPa','hPa.pdf'];
            export_fig(filename,'-pdf');
        else
            plot(histyearsfinal(1).y(mon:12:end),nanmean([squeeze(histZ3weightedregp(1).wa(presind,mon:12:end));squeeze(histZ3weightedregp(2).wa(presind,mon:12:end));...
                squeeze(histZ3weightedregp(3).wa(presind,mon:12:end));...
                squeeze(histZ3weightedregp(4).wa(presind,mon:12:end));squeeze(histZ3weightedregp(5).wa(presind,mon:12:end))],1),'LineWidth',2,'LineStyle','-')
            hold on
            
            histbWAC1_average = squeeze(nanmean(histbWAC1(:,mon,presind,:),1));
            histbWAC2_average = squeeze(nanmean(histbWAC2(:,mon,presind,:),1));
            
            plot(histyearsfinal(i-1).y(histdateind_1(i-1,mon:12:end)),squeeze(histbWAC1_average(2)*[1:length(histyearsfinal(i-1).y(histdateind_1(i-1,mon:12:end)))]+histbWAC1_average(1)),'LineWidth',3)
            plot(histyearsfinal(i-1).y(histdateind_2_alltemp(i-1,mon:12:end)),squeeze(histbWAC2_average(2)*[1:length(histyearsfinal(i-1).y(histdateind_2_alltemp(i-1,mon:12:end)))]+histbWAC2_average(1)),'LineWidth',3)
            
            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/Z3_','ensavg','_',num2str(mon),'_','2and100hPa','hPa.pdf'];
            export_fig(filename,'-pdf');
            
        end
        
    end
ax.ColorOrderIndex = 1;
    fig1 = figure;
    set(fig1,'color','white','position',[100 100 700 800]);   
    pressure = 100;    
    [~,presind] = min(abs(ERApressure - pressure));
    for i = 1:5
        plot(histyearsfinal(i).y(mon:12:end),histZ3weightedregp(i).wa(presind,mon:12:end),'LineWidth',1,'color',[.7 .7 .7])
        hold on                
    end
    
    plot(histyearsfinal(1).y(mon:12:end),nanmean([squeeze(histZ3weightedregp(1).wa(presind,mon:12:end));squeeze(histZ3weightedregp(2).wa(presind,mon:12:end));...
        squeeze(histZ3weightedregp(3).wa(presind,mon:12:end));...
        squeeze(histZ3weightedregp(4).wa(presind,mon:12:end));squeeze(histZ3weightedregp(5).wa(presind,mon:12:end))],1),'LineWidth',2,'LineStyle','-','color',colors2(1,:))

    histbWAC1_average = squeeze(nanmean(histbWAC1(:,mon,presind,:),1));
    histbWAC2_average = squeeze(nanmean(histbWAC2(:,mon,presind,:),1));

    plot(histyearsfinal(i).y(histdateind_1(i,mon:12:end)),squeeze(histbWAC1_average(2)*[1:length(histyearsfinal(i).y(histdateind_1(i,mon:12:end)))]+histbWAC1_average(1)),'LineWidth',4,'color',colors2(2,:))
    plot(histyearsfinal(i).y(histdateind_2_alltemp(i,mon:12:end)),squeeze(histbWAC2_average(2)*[1:length(histyearsfinal(i).y(histdateind_2_alltemp(i,mon:12:end)))]+histbWAC2_average(1)),'LineWidth',4,'color',colors2(3,:))
    
    plot(1979:2016,Z3weighted_ERA(presind,mon:12:end),'color','k','LineWidth',3,'LineStyle','--')
    hold on
    plot(1979:1998,bERA(mon,presind,2)*[1:length(1979:1998)]+bERA(mon,presind,1),'color',colors2(2,:),'LineWidth',3)
    plot(1999:2016,bERA2(mon,presind,2)*[1:length(1999:2016)]+bERA2(mon,presind,1),'color',colors2(3,:),'LineWidth',3)
    
    xlabel('year','fontsize',fsize+2)
    ylabel('Geopotential height (m)','fontsize',fsize+2);                    
    set(gca,'fontsize',fsize);                    
    title([monthnames{mon},' zonal mean at ',num2str(pressure),' hPa'],'fontsize',fsize);

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/Z3_','combine','_',num2str(mon),'_',num2str(pressure),'hPa.pdf'];
    export_fig(filename,'-pdf');
    
end

%%
if line_plots
    fsize = 18;
    month = 12;
    pressure = 2;
    latitude = -70;
    [~,presind] = min(abs(data(1).data.lev - pressure));
    [~,SDpresind] = min(abs(SDhistdata(1).histdata.lev - pressure));
    [~,latind] = min(abs(data(1).data.lat - latitude));
    
    [~,ERApresind] = min(abs(ERApressure - pressure));
    [~,ERAlatind] = min(abs(ERAlatitude - latitude));
    
%     colors = {'k','k','k','b','b','b','r','r','r'};     
%     for i = 1:9;        
%         hold on
%         plot(squeeze(data(i).data.Z3(6,30,10:12:end)),colors{i});
%     end
    colors2 = cbrewer('qual','Set1',5);
    fig = figure;
    set(fig,'color','white','position',[100 100 1000 700]);
    hold on        
    count = 1;
    for i = [1,4,7];
        plot(yearsfinal(i).y(month:12:end),nanmean([Z3weighted(i).wa(presind,month:12:end);Z3weighted(i+1).wa(presind,month:12:end);...
            Z3weighted(i+2).wa(presind,month:12:end)],1),'color',colors2(count+1,:),'LineWidth',3)
        count = count+1;
    end    
    
    plot(histyearsfinal(1).y(month:12:end),nanmean([squeeze(histZ3weighted(1).wa(presind,month:12:end));squeeze(histZ3weighted(2).wa(presind,month:12:end));...
        squeeze(histZ3weighted(3).wa(presind,month:12:end));...
        squeeze(histZ3weighted(4).wa(presind,month:12:end));squeeze(histZ3weighted(5).wa(presind,month:12:end))],1),'color',colors2(1,:),'LineWidth',3,'LineStyle','-')
    
    plot(1979:2016,Z3weighted_ERA(ERApresind,month:12:end),'color',[.7 .7 .7],'LineWidth',3,'LineStyle','--')
    plot(1979:2014,SDhistZ3weighted.wa(SDpresind,month:12:end),'color','k','LineWidth',3,'LineStyle','-')
    
    set(gca,'fontsize',fsize);
    set(gca,'xtick',1950:10:2100,'xticklabel',1950:10:2100);
    lh = legend('REF-C2','SEN-C2-fGHG','SEN-C2-fODS','REF-C1 (prescribed ocean)','ERA-Interim','REF-C1SD');
    set(lh,'location','SouthEast','box','off','fontsize',fsize+2);
    set(gca,'box','on');
    xlabel('year','fontsize',fsize+2);
    ylabel('Geopotential height (m)','fontsize',fsize+2);
    title([monthnames{month},' ensemble average between 90 and 60?S at ',num2str(pressure),' hPa'],'fontsize',fsize+2);
    
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/Z3_',num2str(month),'_',num2str(pressure),'hPa.pdf'];
    export_fig(filename,'-pdf');
end
