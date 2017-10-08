
% Need to add a year onto variable that is using the percentiles

clear all
var = 'O3';
tozvar = 'toz';
tozmonth = 10;
varmonth = 'summer';
percentile = 10;
%% Read in highcl variable
area = '7590S';
ClLevel = 'highCl';
dates = [1995,2016];
directory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/',area,'/'];
files = dir([directory,'*.nc']);

[data.highcl,pressure.highcl,years.highcl,dataRegPres.highcl,regPres.highcl,...
    dataMonthArrange.highcl,dataRegPres_composite.highcl] = predruns_ReadIn(directory,files,var,dates);

% dataMonthArrange = [no. members + ens avg, months, pressures, years]

%% Read in lowcl variable
ClLevel = 'lowCl';
pastdates = [1955,1976];
directory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/',area,'/'];
pastfiles = dir([directory,'*.nc']);

[data.lowcl,pressure.lowcl,years.lowcl,dataRegPres.lowcl,regPres.lowcl,...
    dataMonthArrange.lowcl,dataRegPres_composite.lowcl] = predruns_ReadIn(directory,pastfiles,var,pastdates);

%% Read in TOZ highcl and take percentiles
ClLevel = 'highCl';
lats = [-90,-75];
tozdates = [1995,2015];
directory = ['/Volumes/MyBook/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,~] = ...
    predruns_ReadInlayer(directory,tozfiles,tozvar,tozdates,lats);

[pct_highcl] = predruns_varPercentiles(toz_composite.highcl.montharrange,tozmonth,percentile,length(tozfiles));

%% Read in TOZ lowcl and take percentiles
ClLevel = 'lowCl';
tozpastdates = [1955,1975];
directory = ['/Volumes/MyBook/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfilespast = dir([directory,'*.nc']);
[toz_data.lowcl,toz_years.lowcl,toz_varweighted.lowcl,toz_composite.lowcl,~] = ...
    predruns_ReadInlayer(directory,tozfilespast,tozvar,tozpastdates,lats);

[pct_lowcl] = predruns_varPercentiles(toz_composite.lowcl.montharrange,tozmonth,percentile,length(tozfilespast));

%% Extract high and low variable years and average

for i = 1:length(tozfiles)
    for j = 1:12
        if j <= 4
            varextract.highcl.lowind(i,j).a = squeeze(dataMonthArrange.highcl(i,j,:,pct_highcl.lowindrestruct(i).a+1));                
            varextract.highcl.highind(i,j).a = squeeze(dataMonthArrange.highcl(i,j,:,pct_highcl.highindrestruct(i).a+1));                
        else
            varextract.highcl.lowind(i,j).a = squeeze(dataMonthArrange.highcl(i,j,:,pct_highcl.lowindrestruct(i).a));
            varextract.highcl.highind(i,j).a = squeeze(dataMonthArrange.highcl(i,j,:,pct_highcl.highindrestruct(i).a));
        end
    end
    
end
    
for i = 1:length(tozfilespast)
    for j = 1:12
        if j <= 4
            varextract.lowcl.lowind(i,j).a = squeeze(dataMonthArrange.lowcl(i,j,:,pct_lowcl.lowindrestruct(i).a+1));                
            varextract.lowcl.highind(i,j).a = squeeze(dataMonthArrange.lowcl(i,j,:,pct_lowcl.highindrestruct(i).a+1));                
        else
            varextract.lowcl.lowind(i,j).a = squeeze(dataMonthArrange.lowcl(i,j,:,pct_lowcl.lowindrestruct(i).a));
            varextract.lowcl.highind(i,j).a = squeeze(dataMonthArrange.lowcl(i,j,:,pct_lowcl.highindrestruct(i).a));
        end
    end
    
end

for j = 1:12
    varextractmean.highcl.lowind(j,:) = nanmean([varextract.highcl.lowind(:,j).a],2);
    varextractmean.highcl.highind(j,:) = nanmean([varextract.highcl.highind(:,j).a],2);
    varextractmean.lowcl.lowind(j,:) = nanmean([varextract.lowcl.lowind(:,j).a],2);
    varextractmean.lowcl.highind(j,:) = nanmean([varextract.lowcl.highind(:,j).a],2);
end

vardifference.highcl = varextractmean.highcl.lowind - varextractmean.highcl.highind;
vardifference.lowcl = varextractmean.lowcl.lowind - varextractmean.lowcl.highind;

%% creating plotting array
btoplot = permute(cat(3,vardifference.lowcl,vardifference.highcl),[3,1,2]);
btoplot = circshift(btoplot,[0,8,0]);
%% plotting
if mod(size(btoplot,1),2)
    btoplot(end+1,:,:,:) = zeros(size(btoplot,2),size(btoplot,3),size(btoplot,4));
    btoplot (btoplot == 0) = NaN;
end

plotmtitle = [area,' ',var,' composite differences'];

titles = {['LowCl ',num2str(tozpastdates(1)),'-',num2str(tozpastdates(2))],['HighCl ',...
    num2str(tozdates(1)),'-',num2str(tozdates(2))]}; 

if strcmp(var,'T');
     contourlims = [-15,15];
     maxheight = .1;
     ctitle = 'K';
elseif strcmp(var,'Z3');
    contourlims = [-1500 1500];
    maxheight = .001;
    ctitle = 'm';
elseif strcmp(var,'O3');
    contourlims = [-1.5e-6 1.5e-6];
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
subplotmaps(btoplot,1:12,log(regPres.highcl(1,:)),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
    contourlims,22,1:12,{'M','J','J','A','S','O','N','D','J','F','M','A'},...
    fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(maxheight) log(700)],1,'-',0);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/',var,'/compositedifferences/',...
   var,'_',num2str(tozpastdates(1)),'-',num2str(tozpastdates(2)),'and',num2str(tozdates(1)),'-',num2str(tozdates(2)),'_',area,'.png'];
export_fig(filename,'-png');
