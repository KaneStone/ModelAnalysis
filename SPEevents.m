% plot SPE events
clear all

[~,data,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/WACCM/SPE/spes_1963-2014_c150717.nc');

%% getting months

temp = num2str(data.date);
for j = 1:length(temp)
    mons(j) = str2double(temp(j,5:6));
    years(j) = str2double(temp(j,1:4));
end
clearvars temp

%% extract time
timeindex = years >= 2001 & years <= 2001 & mons >= 11;
toplot = data.Prod(:,timeindex);
dates = data.date(timeindex);
toplot(toplot == 0) = NaN;

%% plotting everything
prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
    5,[],[],2,1,[],[],[],[],.5,[],[],.2,.1};
logprestick = log(prestick);

limits = [0,1,2.^(.7*(1:30))];

cbrew = cbrewer('seq','YlOrBr',32);
cbrew = flipud(cbrew);
cbrew(1,:) = [0 0 0];
contourf(dates,log(data.pressure),toplot,'LineStyle','none');
set(gca,'ytick',fliplr(logprestick),'yticklabel',fliplr(presticklabel));
set(gca,'YDir','reverse');
%colormap(cbrew);
colorbar
