clear variables

%% user inputs

inputs = predruns_inputs;

%% Read in model (user inputs located in function)
[surfacedata,tozdata] = predruns_ReadinModel(inputs,0);

% obtaining number and names of fields 
fields = fieldnames(surfacedata);
nofields = length(fields);

%% Read in NDISC data
if inputs.compareERA
    predruns_readinNDISC(inputs);
end


%% Calculating the upper and lower percentiles of toz based on user inputs defined below
% user inputs

for i = 1:nofields
    nofiles = length(tozdata.(fields{i}).data);
    [pct.(fields{i}),tozextract.(fields{i}),Eachyear.(fields{i})] = predruns_varPercentiles(...
        tozdata.(fields{i}).toz_composite.montharrange,tozdata.(fields{i}).dataMonthArrange,...
        inputs.tozmonth,inputs.percentile,nofiles);
end

longitude = surfacedata.(fields{1}).data(1).lon;
latitude = surfacedata.(fields{1}).data(1).lat;

%% find ozone extremes, mean, and take difference

for i = 1:nofiles
    extremeslow(:,:,:,i) = permute(squeeze(surfacedata.highCl.dataMonthArrange(i,inputs.varmonthtomean,pct.highCl.ind.lowind(i,:),:,:)),[2,3,1]);
    extremeshigh(:,:,:,i) = permute(squeeze(surfacedata.highCl.dataMonthArrange(i,inputs.varmonthtomean,pct.highCl.ind.highind(i,:),:,:)),[2,3,1]);        
end

lowmean = nanmean(extremeslow(:,:,:),3);
highmean = nanmean(extremeshigh(:,:,:),3);
%difference = nanmean(extremeshigh(:,:,:),3) - nanmean(extremeslow(:,:,:),3);
if strcmp(inputs.var,'SNOWHICE')
    %difference = permute(reshape((highmean - lowmean)./lowmean*100,[1,size(highmean)]),[1,3,2]);
    difference = permute(reshape((highmean - lowmean),[1,size(highmean)]),[1,3,2]);
    difference (isnan(difference)) = 0;
elseif strcmp(inputs.var,'ICEFRAC')
    difference = permute(reshape((highmean - lowmean),[1,size(highmean)]),[1,3,2]);
end

combinesnow = 0;
if combinesnow
    inputs.var = 'SNOWHLND';
    %% Read in model (user inputs located in function)
    [surfacedata2,~] = predruns_ReadinModel(inputs,0);

    % obtaining number and names of fields 
    fields2 = fieldnames(surfacedata2);
    nofields2 = length(fields2);
    
    %% find ozone extremes, mean, and take difference

    for i = 1:nofiles
        extremeslow2(:,:,:,i) = permute(squeeze(surfacedata2.highCl.dataMonthArrange(i,inputs.varmonthtomean,pct.highCl.ind.lowind(i,:),:,:)),[2,3,1]);
        extremeshigh2(:,:,:,i) = permute(squeeze(surfacedata2.highCl.dataMonthArrange(i,inputs.varmonthtomean,pct.highCl.ind.highind(i,:),:,:)),[2,3,1]);        
    end

    extremeslowcomp = extremeslow2(:,:,:)+extremeslow(:,:,:);
    extremeshighcomp = extremeshigh2(:,:,:)+extremeshigh(:,:,:);
    
    p = ttest2(permute(extremeslowcomp,[3,2,1]),permute(extremeshighcomp,[3,2,1]));
    
    p (p == 1) = -2;
    p (p == 0) = -1;
    p (p == -2) = 0;
    
    lowmean2 = nanmean(extremeslow2(:,:,:),3);
    highmean2 = nanmean(extremeshigh2(:,:,:),3);
    %difference = nanmean(extremeshigh(:,:,:),3) - nanmean(extremeslow(:,:,:),3);
    %difference2 = permute(reshape((highmean2 - lowmean2)./lowmean2*100,[1,size(highmean2)]),[1,3,2]);
    difference2 = permute(reshape((highmean2 - lowmean2),[1,size(highmean2)]),[1,3,2]);
    difference2 (isnan(difference2)) = 0;
    differencecombine = difference+difference2;
    
end

%% plot data
if strcmp(inputs.var,'SNOWHLND')
    title = {['Change in ',monthnames(inputs.varmonthtomean,0,0),' snow depth due to ',monthnames(inputs.tozmonth,0,0), ' ozone extremes']};
    ctitle = {'m'};
    clims = [-.2 .2];
    titext = 'SNOWLAND';
    toplotc = differencecombine;
    toplotp = p;
elseif strcmp(inputs.var,'SNOWHICE')
    title = {['Change in ',monthnames(inputs.varmonthtomean,0,0),' snow depth due to ',monthnames(inputs.tozmonth,0,0), ' ozone extremes']};
    ctitle = {'m'};
    clims = [-.2 .2];
    titext = 'SNOWLICE';
    if combinesnow
        toplotc = differencecombine;
        toplotp = p;
    else
        toplotc = difference;
        toplotp = [];
    end                    
    
elseif strcmp(inputs.var,'ICEFRAC')
    title = {['Change in ',monthnames(inputs.varmonthtomean,0,0),' sea ice fraction due to ',monthnames(inputs.tozmonth,0,0), ' ozone extremes']};
    ctitle = {'ice fraction'};
    clims = [-.1 .1];
    titext = 'ICEFRAC';
    toplotc = difference;
    toplotp = [];
end

filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
subplotmaps(toplotc(1,:,:),longitude,latitude,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',ctitle,'on',...
    clims,22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

filename = [filedir,monthnames(inputs.varmonthtomean,1,1),'_',titext,'_',...
    monthnames(inputs.tozmonth,1,1),'_Arcticozoneextremes_over',...
    num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
%print(filename,'-depsc');
export_fig(filename,'-png');

clearvars title

%% plot contour lines
matdefcols = get(groot,'DefaultAxesColorOrder');

fraction = .15;

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewer('qual','Set2',8);
createfig('medium','on')

m_proj('Stereographic','lon',-180,'lat',90,'rad',abs(50))           
[~,h] = m_contour(longitude,latitude,squeeze(lowmean(:,:)),[0,fraction]);
set(h,'LineWidth',3,'color',cbrewqual(1,:));
hold on
[~,h2] = m_contour(longitude,latitude,squeeze(highmean(:,:)),[0,fraction]);
set(h2,'LineWidth',3,'color',cbrewqual(2,:));
m_coast('patch',cbrewqual2(7,:));
m_grid('ytick',0:15:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',20,'LineWidth',2);    
th = title('Change in sea ice extent relative to Arctic ozone extemes','fontsize',20);
set(th,'units','normalized')
set(th,'position',[.5 1.045 0])
filename = [filedir,monthnames(inputs.varmonthtomean,1,1),'_iceextenddiff_from_',...
    monthnames(inputs.tozmonth,1,1),'_Arcticozoneextremes_over',...
    num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
export_fig(filename,'-png');

%% zoom in on barents kara sea area (75N 40E)
createfig('medium','on')
x = longitude;
x(size(x,1)/2+1:size(x,1)) = x(size(x,1)/2+1:size(x,1)) - 360;
x = [x(size(x,1)/2+1:size(x,1));x(1:size(x,1)/2)];
datalow = circshift(lowmean,size(lowmean,2)/2,2);
datahigh = circshift(highmean,size(highmean,2)/2,2);
m_proj('lambert','lon',[0 120],'lat',[65 85]);
[~,h] = m_contour(x,latitude,squeeze(datalow(:,:)),[0,fraction]);
hold on
set(h,'LineWidth',3,'color',cbrewqual(1,:));
hl = plot(1,1,'color',cbrewqual(1,:),'LineWidth',3);

[~,h2] = m_contour(x,latitude,squeeze(datahigh(:,:)),[0,fraction]);
set(h2,'LineWidth',3,'color',cbrewqual(2,:));
hl2 = plot(1,1,'color',cbrewqual(2,:),'LineWidth',3);
m_coast('patch',cbrewqual2(7,:));
m_grid('box','fancy','tickdir','in');
th = title(['Ensemble average ',monthnames(inputs.varmonthtomean,0,0),' sea ice extent relative to ',...
    monthnames(inputs.tozmonth,0,0),' Arctic ozone extemes'],'fontsize',20);
set(th,'units','normalized')
set(th,'position',[.5 1.045 0])
filename = [filedir,monthnames(inputs.varmonthtomean,1,1),'_iceextenddiff_from_',...
    monthnames(inputs.tozmonth,1,1),'_Arcticozoneextremes_over_barentskara',...
    num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
lh = legend([hl,hl2],'relative to TCO lower 20th percentile','relative to TCO upper 20th percentile');
set(lh,'box','off','fontsize',20,'location','southoutside')
export_fig(filename,'-png');
