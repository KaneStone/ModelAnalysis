clear all
var = 'O3';
area = '7590S';
ClLevel = 'highCl';
dates = [1995,2005];
directory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/',area,'/'];
files = dir([directory,'*.nc']);

% Read in future
[data.future,pressure.future,years.future,dataRegPres.future,regPres.future,dataMonthArrange.future] = ...
    predruns_ReadIn(directory,files,var,dates);

ClLevel = 'lowCl';
pastdates = [1965,1975];
directory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/',area,'/'];
pastfiles = dir([directory,'*.nc']);

% Read in past
[data.past,pressure.past,years.past,dataRegPres.past,regPres.past,dataMonthArrange.past] = ...
    predruns_ReadIn(directory,pastfiles,var,pastdates);

%% Taking differences of future period from past period
pastAverage = nanmean(dataMonthArrange.past,4);
futureAverage = nanmean(dataMonthArrange.future,4);
for i = 1:size(futureAverage,1)
    futureDifference(i,:,:) = futureAverage(i,:,:) - pastAverage(end,:,:);    
end

OctPastAve = squeeze(pastAverage(end,12,:));
OctFutureAve = squeeze(futureAverage(end,12,:));

plot(squeeze(pastAverage(1:end-1,10,:)),log(regPres.future(1,:)),'y');
hold on
plot(squeeze(futureAverage(1:end-1,10,:)),log(regPres.future(1,:)),'c');
plot(OctPastAve,log(regPres.future(1,:)),'r','LineWidth',2);
plot(OctFutureAve,log(regPres.future(1,:)),'b','LineWidth',2);
figure
plot(OctPastAve-OctFutureAve,log(regPres.future(1,:)));

%% plotting

if strcmp(var,'T')
    titlevar = 'Temperature';
elseif strcmp(var,'Z3')
    titlevar = 'Geopotential height';
elseif strcmp(var,'O3')
    titlevar = 'Ozone';
end

btoplot = circshift(futureDifference,[0,7,0]);
if mod(size(btoplot,1),2)
    btoplot(end+1,:,:,:) = zeros(size(btoplot,2),size(btoplot,3),size(btoplot,4));
    btoplot (btoplot == 0) = NaN;
end

plotmtitle = [area,' ',num2str(dates(1)),'-',num2str(dates(2)),', ',titlevar,' difference from ',num2str(pastdates(1)),'-',num2str(pastdates(2))];

titles = {'No. 1','No. 2','No. 3','No. 4','No. 5','No. 6','No. 7','No. 8','No. 9','Ensemble average'}; 

if strcmp(var,'T');
     contourlims = [-15 15];
     maxheight = 1;
     ctitle = 'Kelvin';
elseif strcmp(var,'Z3');
    contourlims = [-1000 1000];
    maxheight = .1;
    ctitle = 'm';
elseif strcmp(var,'O3');
    contourlims = [-2e-6 2e-6];
    maxheight = .1;
    ctitle = 'mol/mol';
end

prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
if maxheight > 10
    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
        5,[],[],2,1,[],[],[],[],.5,[],[],.2,.1};
else
    presticklabel = {1000,[],[],[],[],[],[],[],[],100,[],[],[],[],[],[],[],[],10,[],[],[],[],...
        [],[],[],[],1,[],[],[],[],[],[],[],[],.1};
end
logprestick = log(prestick);

subplotmaps(btoplot,1:12,log(regPres.future(1,:)),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
    contourlims,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(maxheight) log(700)],1,'-',0);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/',var,'/differences/','Individualfuture_differences',...
    num2str(dates(1)),'-',num2str(dates(2)),'_from_',num2str(pastdates(1)),'-',num2str(pastdates(2)),'_',area,'.png'];
export_fig(filename,'-png');
