% TCO percent difference lines

% Read in and construct trend data for maps

clear all;
clc;

runtype = 'monthly';
TCOdirectory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/TCO/'];

%addpath(genpath('/home/stonek/code'))

TCOfiles = dir([TCOdirectory,'TOZ*']);

yearmin = 20000201;
yearmax = 20150201;
remove2002 = 'none'; % 'all_lats' or 'Southern' or 'none'

cmap = flipud(cbrewer('div','RdBu',10));

TCOname = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
TCOnamel = {'CCMI','MAM','VC-MAM','Chem-only','Chem-only-fSSTs','Chem-only-noleap'};

fsize = 18;
lwidth = 2;
lats = [-90 -60];
anomyear = [2000 2014];
months = [1:12];
%lats = [[-90 -60];[-60 -30];[-30 0];[0 30];[30 60];[60 90]];

mon = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};

szlat = size(lats);
TCOzonallat = [];

for i = 1:length(TCOfiles)
    [info.(TCOname{i}), data.(TCOname{i}), attributes.(TCOname{i})] = ...
        Read_in_netcdf([TCOdirectory,TCOfiles(i).name]);
    
    %zonal mean data
    TCOzonal.(TCOname{i}) = squeeze(nanmean(data.(TCOname{i}).toz));
    
    for p = 1:szlat(1)
        latindex(p).p = find(data.(TCOname{i}).lat >= lats(p,1) & ...
            data.(TCOname{i}).lat <= lats(p,2));
    end
    
    %removing leap years
    dates = num2str(data.(TCOname{i}).date);
    monthday = str2num(dates(:,5:8));
    years =  str2num(dates(:,1:4));
    monthdayindex = find(monthday == 0229);
    if strcmp(runtype,'daily')
        TCOzonal.(TCOname{i})(:,monthdayindex) = []; 
        monthday(monthdayindex) = [];
        data.(TCOname{i}).date(monthdayindex) = [];
        interval = 365;
        %extracting time period
        monthday = monthday(data.(TCOname{i}).date >= 20000101 & data.(TCOname{i}).date < 20150101);
        TCOzonal20002014.(TCOname{i}) = TCOzonal.(TCOname{i})(:,...
            data.(TCOname{i}).date >= 20000101 & data.(TCOname{i}).date < 20150101);
        TCOzonal20002014.(TCOname{i})(:,973:1338) = NaN; 
    else
        interval = 12;
        %extracting time period        
        
        TCOzonal20002014.(TCOname{i}) = TCOzonal.(TCOname{i})(:,...
        data.(TCOname{i}).date >= yearmin & data.(TCOname{i}).date < yearmax);
    
        datetemp = data.(TCOname{i}).date(data.(TCOname{i}).date >= yearmin & data.(TCOname{i}).date < yearmax);
        %removing 2002 from September
    
        switch remove2002
            case 'all_lats'
                TCOzonal20002014.(TCOname{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            case 'Southern'        
                TCOzonal20002014.(TCOname{i})(data.(TCOname{i}).lat < 0,...
                datetemp >= 20021001 & datetemp < 20031001) = NaN; 
                
                Latitudes.(TCOname{i})(data.(TCOname{i}).lat < 0,...
                datetemp >= 20021001 & datetemp < 20031001) = NaN;                
        end
                        
    end
    
    if strcmp(runtype,'daily')
       time_xaxis = 365;
    else
       time_xaxis = 12; 
    end
    
    [WACCMyearmeananomaly(i,:),WACCMyearmean(i,:)] = TCOanomaly(TCOzonal20002014.(TCOname{i}),...
        lats,data.(TCOname{i}).lat,anomyear,[2000, 2015],months);
    
end


SBUVfilename = ['/Users/kanestone/work/projects/WACCM/netcdffiles/SBUV/',...
    'sbuv_v86_mod.int_lyr.70-15.za.r5fromMcPeters.txt'];
[SBUVdata, SBUVmonthyear, SBUVTCO, SBUVyearmean,SBUVyearmeananomaly] = ...
    ReadSBUV(SBUVfilename,lats,anomyear,months);

WACCMyearmean = [zeros(6,30), WACCMyearmean];
WACCMyearmean (WACCMyearmean == 0) = NaN;
WACCMyearmeananomaly = [zeros(6,30), WACCMyearmeananomaly];
WACCMyearmeananomaly (WACCMyearmeananomaly == 0) = NaN;

combined = [WACCMyearmean; SBUVyearmean];
combined(:,1:30) = NaN;

combinedanom = [WACCMyearmeananomaly; SBUVyearmeananomaly];
combinedanom(:,1:30) = NaN;

volcanoes = WACCMyearmeananomaly(2,:) - WACCMyearmeananomaly(3,:);
dynamics = WACCMyearmeananomaly(3,:) - WACCMyearmeananomaly(4,:);
chemistry = WACCMyearmeananomaly(5,:);

if length(months) == 1
    montitle = mon{months};
elseif length(months) == 12
    montitle = 'Yearly';
elseif length(months) ~= 1 || length(months) ~= 12
    montitle = sprintf('%02d',months);
end

% plotting
%%
% Year mean anomaly
figure;
fig = gcf;
set(fig,'color','white','position',[100 100 840 560]);
fh = plot(WACCMyearmeananomaly([2:3,5],:)','Linewidth',lwidth);
hold on
sbuvh = plot(SBUVyearmeananomaly,'--k','Linewidth',lwidth);
plot(1:50,zeros(1,50),'-.','color',[.3 .3 .3]);
ylim([floor(min(combinedanom(:))) ceil(max(combinedanom(:)))]);
set(gca,'ytick',-20:1:30,'yticklabel',-20:1:30,'xtick',1:1:50,'xticklabel',1970:1:2030,...
    'fontsize',fsize-2);
xlim([30 46])
xlabel('Year','fontsize',fsize);
ylabel('Percent TCO','fontsize',fsize);
title([montitle,' difference from the ',num2str(anomyear(1)),'-',num2str(anomyear(2)),' mean, ',...
    num2str(lats(1)),' to ',num2str(lats(2)),'{\circ}N'],'fontsize',fsize+2);
lh = legend([fh',sbuvh],'MAM','VC-MAM','Chem-only-fSSTs','SBUV','Location','NorthWest');
set(lh,'box','off','fontsize',fsize-2);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/TimeSeries/Difference','_',montitle,'_',...
    num2str(anomyear(1)),'-',num2str(anomyear(2)),'_',num2str(lats(1)),'Sto',num2str(lats(2)),'N.pdf'];
export_fig(fig,filename,'-nofontswap');

%Contributions
forlimits = [volcanoes,dynamics,chemistry];

figure;
fig2 = gcf;
set(fig2,'color','white','position',[100 100 840 560]);
dh(1) = plot(volcanoes,'LineWidth',lwidth);
hold on
dh(2) = plot(dynamics,'LineWidth',lwidth);
dh(3) = plot(chemistry,'LineWidth',lwidth);
plot(1:50,zeros(1,50),'-.','color',[.3 .3 .3]);

ylim([floor(min(forlimits(:))) ceil(max(forlimits(:)))]);
set(gca,'ytick',-20:1:30,'yticklabel',-20:1:30,'xtick',1:1:50,'xticklabel',1970:1:2030,'fontsize',fsize-2);

%ylim([-1 2]);
%set(gca,'ytick',-20:1:30,'yticklabel',-20:1:30,'xtick',1:1:50,'xticklabel',1970:1:2030,'fontsize',fsize-2);

xlim([30 46])
xlabel('Year','fontsize',fsize);
ylabel('Percent TCO','fontsize',fsize);
title([montitle,' Contributions to the difference, ',...
    num2str(lats(1)),' to ',num2str(lats(2)),'{\circ}N'],'fontsize',fsize+2);
lh = legend(dh,'Volcanoes','Dynamics','Chemistry','Location','NorthWest');
set(lh,'box','off','fontsize',fsize-2);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/TimeSeries/Contributions','_',montitle,'_',...
    num2str(anomyear(1)),'-',num2str(anomyear(2)),'_',num2str(lats(1)),'Sto',num2str(lats(2)),'N.pdf'];
export_fig(fig2,filename,'-nofontswap');

%Everything
figure;
fig3 = gcf;
set(fig3,'color','white','position',[100 100 840 560]);
eh = plot(WACCMyearmean','LineWidth',lwidth);
hold on
sbuvh2 = plot(SBUVyearmean,'--k','Linewidth',2);
plot(1:50,zeros(1,50));
ylim([floor(min(combined(:))/2)*2 ceil(max(combined(:))/2)*2]);
set(gca,'ytick',100:5:500,'yticklabel',100:5:500,'xtick',1:1:50,'xticklabel',1970:1:2030,'fontsize',fsize-2);
xlim([30 46])
xlabel('Year','fontsize',fsize);
ylabel('Dobson Units','fontsize',fsize);
title([montitle,' TCO zonal mean, ',...
    num2str(lats(1)),' to ',num2str(lats(2)),'{\circ}N'],'fontsize',fsize+2);
lh = legend([eh', sbuvh2],[TCOname,'SBUV'],'Location','NorthWest');
set(lh,'box','off','fontsize',fsize-2);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/TimeSeries/TimeSeries','_',montitle,'_',...
    num2str(anomyear(1)),'-',num2str(anomyear(2)),'_',num2str(lats(1)),'Sto',num2str(lats(2)),'N.pdf'];
export_fig(fig3,filename,'-nofontswap');

% figure;
% fig2 = gcf;
% set(fig2,'color','white','position',[100 100 840 560]);
% plot(volcanoes);
% hold on
% plot(dynamics);
% plot(chemistry);
% plot(1:50,zeros(1,50));
% ylim([-50 50]);
% 
% set(gca,'ytick',-5:1:20,'yticklabel',-5:1:20,'xtick',1:1:50,'xticklabel',1970:1:2030);
% xlim([30 46])