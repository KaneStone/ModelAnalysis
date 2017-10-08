% Predictability runs trend maps


clear all
var = 'T';
area = '7590S';
ClLevel = 'highCl';
trenddates = [2000,2016];
directory = ['/Volumes/MyBook/work/data/predruns/',var,'/',ClLevel,'/',area,'/'];
files = dir([directory,'*.nc']);

DepletionEraDifference = 1;
%% read in data and interpolate onto regular pressure levels.

for i = 1:length(files)
    %read in weighted area average data [height,time]
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    
    %calculate pressure from hybrid_height
    pressure(i).p = permute(repmat(data(i).hyam.*data(i).P0,1,size(data(i).(var),1),size(data(i).(var),3)) + ...
       repmat(data(i).hybm,1,size(data(i).(var),1),size(data(i).(var),3)) ...
       .* permute(repmat(data(i).PS,1,1,length(data(i).lev)),[3,1,2]),[2,1,3])./100; 
    
    %weighted average
    for j = 1:size(data(i).(var),2)
        varweighted(i).wa(j,:) = weightedaverage(squeeze(data(i).(var)(:,j,:)),data(i).lat);        
        presweighted(i).wa(j,:) = weightedaverage(squeeze(pressure(i).p(:,j,:)),data(i).lat);        
    end     
        
    % construct year only vector
    years(i).y = CCMI_years(data(i).date);
    
    % interpolate onto regular pressure
    [dataRegPres(i,:,:),regPres(i,:)] = intRegPres(varweighted(i).wa,presweighted(i).wa);
    
end

dataRegPres(length(files)+1,:,:) = nanmean(dataRegPres,1);
regPres(length(files)+1,:) = nanmean(regPres,1);
years(length(files)+1).y = years(1).y;
for i = 1:size(dataRegPres,1)
    % taking trends
    dateind = find(years(i).y >= trenddates(1) & years(i).y <= trenddates(2));  
    for j = 1:12        
        dataMonthArrange(i,j,:,:) = dataRegPres(i,:,dateind(1)+j-1:12:dateind(end));
        for k = 1:size(dataMonthArrange,3)
            [b(i,j,k,:),bint(i,j,k,:,:)] = regress(squeeze(dataMonthArrange(i,j,k,:)),...
                [ones(1,size(dataMonthArrange(i,j,:,:),4));...
                1:size(dataMonthArrange(i,j,:,:),4)]');
        end
    end    
    
end

% if DepletionEraDifference == 1
%     
%     [] = predruns_depEraDiff();
%     
% end

%% subplotmaps

btoplot = circshift(b(:,:,:,2),[0,7,0,0])*10;

if mod(size(btoplot,1),2)
    btoplot(end+1,:,:,:) = zeros(size(btoplot,2),size(btoplot,3),size(btoplot,4));
    btoplot (btoplot == 0) = NaN;
end

plotmtitle = [area,' ',var,' trends (',num2str(trenddates(1)),'-',num2str(trenddates(2)),')'];
if strcmp(ClLevel,'highCl');
    titles = {'No. 1','No. 2','No. 3','No. 4','No. 5','No. 6','No. 7','No. 8','No. 9','Ensemble average'}; 
else
    titles = {'No. 1','No. 2','No. 3','No. 4','No. 5','No. 6','No. 7','No. 8','No. 9',...
        'No. 10','Ensemble average','{}'}; 
end

if strcmp(var,'T');
     contourlims = [-3 3];
     maxheight = .1;
     ctitle = 'K/decade';
else
    contourlims = [-100 100];
    maxheight = 5;
    ctitle = 'm/decade';
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
%btoplot (btoplot < contourlimsfuture(1)) = contourlimsfuture(1)+1;
%btoplot (btoplot > contourlimsfuture(2)) = contourlimsfuture(2)-1;

% subplotmaps(btoplot,1:12,log(regPres(1).d),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
%     contourlims,14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
%     fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(700)],1,'-',0);

subplotmaps(btoplot,1:12,log(regPres(1,:)),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
    contourlims,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(maxheight) log(700)],1,'-',0);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/predruns/',var,'/trends/','Individualfuture_',...
    num2str(trenddates(1)),'-',num2str(trenddates(2)),'_',area,'.png'];
export_fig(filename,'-png');

