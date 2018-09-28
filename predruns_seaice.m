clear variables

%% user inputs

inputs = predruns_inputs;

%% Read in model (user inputs located in function)
[surfacedata,tozdata] = predruns_ReadinModel(inputs,0);

% obtaining number and names of fields 
fields = fieldnames(surfacedata);
nofields = length(fields);


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
difference = permute(reshape(highmean - lowmean,[1,size(highmean)]),[1,3,2]);

%% plot data
if strcmp(inputs.var,'SNOWHLND')
    title = {['Change in ',monthnames(inputs.varmonthtomean,0,0),' land snow depth due to ',monthnames(inputs.tozmonth,0,0), ' ozone extremes']};
    ctitle = {'meters'};
    clims = [-.05 .05];
    titext = 'SNOWLAND';
elseif strcmp(inputs.var,'SNOWHICE')
    title = {['Change in ',monthnames(inputs.varmonthtomean,0,0),' ice snow depth due to ',monthnames(inputs.tozmonth,0,0), ' ozone extremes']};
    ctitle = {'meters'};
    clims = [-.05 .05];
    titext = 'SNOWLICE';
elseif strcmp(inputs.var,'ICEFRAC')
    title = {['Change in ',monthnames(inputs.varmonthtomean,0,0),' sea ice fraction due to ',monthnames(inputs.tozmonth,0,0), ' ozone extremes']};
    ctitle = {'ice fraction'};
    clims = [-.2 .2];
    titext = 'ICEFRAC';
end

filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
subplotmaps(difference(1,:,:),longitude,latitude,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',ctitle,'on',...
    [-.05,.05],22,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');

filename = [filedir,monthnames(inputs.varmonthtomean,1,1),'_',titext,'_',monthnames(inputs.tozmonth,1,1),'_Arcticozoneextremes_over',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
%print(filename,'-depsc');
export_fig(filename,'-png');

%% plot contour lines
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
filename = [filedir,monthnames(inputs.varmonthtomean,1,1),'_iceextenddiff_from_',monthnames(inputs.tozmonth,1,1),'_Arcticozoneextremes_over',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
%print(filename,'-depsc');
export_fig(filename,'-png');
