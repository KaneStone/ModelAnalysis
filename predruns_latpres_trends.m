%predruns solar cycle stuff

clear all 
%close all

%% read in predruns

var = 'O3';
month = 10;
ClLevel = 'highCl';
timeperiodhigh = [2005,2015];%[1955,1975]
timeperiod = [2005,2014];%[1955,1975]

vardirectory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/zonalmean/'];
varfiles = dir([vardirectory,'*.nc']);

[data.(ClLevel),years.(ClLevel),composite.(ClLevel),dataMonthArrange.(ClLevel)]...
    = predruns_ReadInlayer(vardirectory,varfiles,var,timeperiodhigh,[-90 90],0);

% calculate pressure
for i = 1:length(varfiles) 
    data.(ClLevel)(i).Pressure = repmat(data.(ClLevel)(i).hyam,[1,size(data.(ClLevel)(i).PS,1),size(data.(ClLevel)(i).PS,2)]).*100000 ...
        + repmat(data.(ClLevel)(i).hybm,[1,size(data.(ClLevel)(i).PS,1),size(data.(ClLevel)(i).PS,2)]).*...
        permute(repmat(data.(ClLevel)(i).PS,[1,1,size(data.(ClLevel)(i).hyam,1)]),[3,1,2]);
    
    data.(ClLevel)(i).Pressure = permute(data.(ClLevel)(i).Pressure,[2,1,3]);
end



%% extracting times and taking regular pressure

yearsindex(1) = find(years.(ClLevel)(1).y == timeperiod(1),1,'first');
yearsindex(2) = find(years.(ClLevel)(1).y == timeperiod(2),1,'last');

for i = 1:length(varfiles) 
    data.(ClLevel)(i).timeextract = data.(ClLevel)(i).O3(:,:,yearsindex(1):yearsindex(2)); % starts on December for seasons
    data.(ClLevel)(i).pressuretimeextract = data.(ClLevel)(i).Pressure(:,:,yearsindex(1):yearsindex(2)); % starts on December for seasons
    for j = 1:size(data.(ClLevel)(i).O3,1)
        [regpres(i,j,:,:),pres] = intRegPres(squeeze(data.(ClLevel)(i).timeextract(j,:,:)),squeeze(data.(ClLevel)(i).pressuretimeextract(j,:,:))./100);
    end
end

regpres_ensmean = squeeze(nanmean(regpres,1));

%% take weighted average
lats = [-90,-80;-80,-70;-70,-60;-60,-50;-50,-40;-40,-30;-30,-20;-20,-10;-10,0;...
    0,10;10,20;20,30;30,40;40,50;50,60;60,70;70,80;80,90];

for i = 1:size(lats,1)
    latindex = find(data.(ClLevel)(1).lat >= lats(i,1) & data.(ClLevel)(1).lat < lats(i,2));
    for j = 1:size(regpres_ensmean,2)
        regpres_wa(i,j,:) = weightedaverage(squeeze(regpres_ensmean(latindex,j,:)),data.(ClLevel)(1).lat(latindex));
    end
end


%% extract month array
for i = 1:12
    regpres_wa_month(i).r = regpres_wa(:,:,i:12:end);
end

%% taking regression
for k = 1:length(regpres_wa_month)
    for i = 1:size(regpres_wa_month(1).r,1)
        for j = 1:size(regpres_wa_month(1).r,2)
            [b(k,i,j,:),~] = regress(squeeze(regpres_wa_month(k).r(i,j,:)),...
                [ones(length(squeeze(regpres_wa_month(k).r(i,j,:))),1),(1:length(squeeze(regpres_wa_month(k).r(i,j,:))))']);       
        end
    end
end


%% plotting


for i = 1:size(regpres_wa_month(1).r,1)/2
    bsouth = circshift(permute(squeeze(b(:,1:9,:,2)),[2,1,3]),[0,7,0])*10*1e6;
    bnorth = circshift(permute(squeeze(b(:,10:18,:,2)),[2,1,3]),[0,7,0])*10*1e6;
end

    prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
        5,[],[],2,1,[],[],[],[],.5,[],[],.2,.1};
    logprestick = log(prestick);

    titles = {'1','2','3','4','5','6','7','8','9'};
    plotmtitle =  'trends';
    
    subplotmaps(bnorth,1:12,log(pres),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','ppmv/decade','on',...
        [-.6 .6],22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(.1) log(300)],1,'-',0,[]);



%% read in swoosh

[~,Swoosh,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/combinedo3qeq_swoosh-v02.1-198401-201501-latpress-2.5deg-L31.nc');
temp = repmat(1984:2014,12,1);
Swoosh.years = [temp(:);2015];

yearsindexsw(1) = find(Swoosh.years == timeperiod(1),1,'first');
yearsindexsw(2) = find(Swoosh.years == timeperiod(2),1,'last');

sw_textract = Swoosh.combinedo3qeq(:,:,yearsindexsw(1):yearsindexsw(2));

%% take weighted average
lats = [-90,-80;-80,-70;-70,-60;-60,-50;-50,-40;-40,-30;-30,-20;-20,-10;-10,0;...
    0,10;10,20;20,30;30,40;40,50;50,60;60,70;70,80;80,90];

for i = 1:size(lats,1)
    latindex = find(Swoosh.lat >= lats(i,1) & Swoosh.lat < lats(i,2));
    for j = 1:size(sw_textract,2)
        swoosh_wa(i,j,:) = weightedaverage(squeeze(sw_textract(latindex,j,:)),Swoosh.lat(latindex));
    end
end

%% extracting months
for i = 1:12
   sw_textract_months(i).s = swoosh_wa(:,:,i:12:end);
end

%% taking regression
for k = 1:length(sw_textract_months)
    for i = 1:size(sw_textract_months(1).s,1)
        for j = 1:size(sw_textract_months(1).s,2)
            [bs(k,i,j,:),~] = regress(squeeze(sw_textract_months(k).s(i,j,:)),...
                [ones(length(squeeze(sw_textract_months(k).s(i,j,:))),1),(1:length(squeeze(sw_textract_months(k).s(i,j,:))))']);       
        end
    end
end

%% taking trends
%% plotting


for i = 1:size(regpres_wa_month(1).r,1)/2
    bsouth = circshift(permute(squeeze(bs(:,1:9,:,2)),[2,1,3]),[0,7,0])*10;
    bnorth = circshift(permute(squeeze(bs(:,10:18,:,2)),[2,1,3]),[0,7,0])*10;
end

    prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
        5,[],[],2,1,[],[],[],[],.5,[],[],.2,.1};
    logprestick = log(prestick);

    titles = {'test','test','test','test','test','test','test','test','test'};
    plotmtitle =  'test';
    
    subplotmaps(bnorth,1:12,log(Swoosh.level),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','ppmv/decade','on',...
        [-.6 .6],22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(.1) log(300)],1,'-',0,[]);
