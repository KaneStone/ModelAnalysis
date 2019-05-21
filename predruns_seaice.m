clear variables

%% user inputs

inputs = predruns_inputs;
inputs.compareSSW = 0;
%% Read in model (user inputs located in function)
if strcmp(inputs.var,'ICEFRAC')
    seaice = 1;
else
    seaice = 0;
end
[surfacedata,tozdata] = predruns_ReadinModel(inputs,seaice);

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
        inputs.tozmonth,inputs.percentile,nofiles,inputs.varmonth);
end

longitude = surfacedata.(fields{1}).data(1).lon;
latitude = surfacedata.(fields{1}).data(1).lat;

%% if compare to SSWs
if inputs.compareSSW
    SSWin.lat = 60;
    SSWin.pres = 10;
    SSWin.plot = 0;
    [SSWSPVcomposite] = predruns_SSWSPV(SSWin);
end

%% find ozone extremes, mean, and take difference
seasons = 0;
if seasons 
    mons = [4,5,6;7,8,9;10,11,12];
else
    mons = [1:12]';
end

%% Read in SSW and ozone coincidence extreme
coincidence = 0;
if coincidence
    coincidencenoSPV = load('/Volumes/ExternalOne/work/data/predruns/output/TCOwindextremes/coincidence_wind40and0.mat');
    coincidencenoSSW = load('/Volumes/ExternalOne/work/data/predruns/output/TCOwindextremes/coincidence_wind48and-2.mat');
    
end

%%

for l = 1:length(mons)
    inputs.varmonthtomean = mons(l,:);
for i = 1:nofiles
    
    
    if inputs.compareSSW
        extremeslow_wind(i).w(:,:,:) = permute(squeeze(nanmean(surfacedata.highCl.dataMonthArrange(i,inputs.varmonthtomean,SSWSPVcomposite.model.SPV(i+1,:),:,:),2)),[2,3,1]); %+1 to account for only 9 toz files
        extremeshigh_wind(i).w(:,:,:) = permute(squeeze(nanmean(surfacedata.highCl.dataMonthArrange(i,inputs.varmonthtomean,SSWSPVcomposite.model.SSW(i+1,:),:,:),2)),[2,3,1]);        
    end
    
    if ~coincidence
        extremeslow(:,:,:,i) = permute(squeeze(nanmean(surfacedata.highCl.dataMonthArrange(i,inputs.varmonthtomean,pct.highCl.ind.lowind(i,:),:,:),2)),[2,3,1]);
        extremeshigh(:,:,:,i) = permute(squeeze(nanmean(surfacedata.highCl.dataMonthArrange(i,inputs.varmonthtomean,pct.highCl.ind.highind(i,:),:,:),2)),[2,3,1]);        
    else
        if ~isempty(coincidencenoSPV.ozoneind.lower.noSPV(i+1).a)
            extremeslowtemp(i).a = permute(squeeze(nanmean(surfacedata.highCl.dataMonthArrange(i,inputs.varmonthtomean,coincidencenoSPV.ozoneind.lower.noSPV(i+1).a,:,:),2)),[2,3,1]);
            if size(extremeslowtemp(i).a,2) == 1
                extremeslowtemp(i).a = permute(extremeslowtemp(i).a,[2,3,1]);
            else
                extremeslowtemp(i).a = permute(extremeslowtemp(i).a,[3,1,2]);
            end
        else
            extremeslowtemp(i).a = zeros(1,96,144);
            extremeslowtemp(i).a (extremeslowtemp(i).a == 0) = NaN;
        end
        if ~isempty(coincidencenoSSW.ozoneind.upper.noSSW(i+1).a)        
            extremeshightemp(i).a = permute(squeeze(nanmean(surfacedata.highCl.dataMonthArrange(i,inputs.varmonthtomean,coincidencenoSSW.ozoneind.upper.noSSW(i+1).a,:,:),2)),[2,3,1]);        
            if size(extremeshightemp(i).a,2) == 1
                extremeshightemp(i).a = permute(extremeshightemp(i).a,[2,3,1]);
            else
                extremeshightemp(i).a = permute(extremeshightemp(i).a,[3,1,2]);
            end
        else
            extremeshightemp(i).a = zeros(1,96,144);
            extremeshightemp(i).a (extremeshightemp(i).a == 0) = NaN;
        end
    end
    
end

%%
if coincidence
    extremeslow = permute(cat(1,extremeslowtemp(:).a),[2,3,1]);
    extremeshigh = permute(cat(1,extremeshightemp(:).a),[2,3,1]);
end

%% plotting individuals
plotind = 0;
if plotind
    for i = 1:size(extremeslow,4)
        toplotdifference_temp = permute(repmat(nanmean(extremeshigh(:,:,:,i),3) - nanmean(extremeslow(:,:,:,i),3),1,1,1),[3,2,1]);
        toplotdifferencefinal = cat(2,toplotdifference_temp,toplotdifference_temp(:,1,:));
        longtouse = [longitude;longitude(1)+360];
        
        title = {['No. ',sprintf('%02d',i),', Change in ',monthnames(inputs.varmonthtomean,1,'single'),' sea ice fraction due to ',monthnames(inputs.tozmonth,0,'single'), ' ozone extremes']};
        ctitle = {'Ice fraction'};
        clims = [-.4 .4];
        titext = 'ICEFRAC';
        
        
        filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
        subplotmaps(toplotdifferencefinal,longtouse,latitude,{'div','RdBu'},1,p,16,title,'Longitude','Latitude',ctitle,'on',...
            clims,18,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic');

        filename = [filedir,'Ind',' No. ',sprintf('%02d',i),'_',monthnames(inputs.varmonthtomean,1,'single'),'_',titext,'_',...
            monthnames(inputs.tozmonth,1,'long'),'_Arcticozoneextremes_over',...
            num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
        %print(filename,'-depsc');
        export_fig(filename,'-png');
    end
end

%%
windmetric = 0;
if windmetric
    lowmean = nanmean(cat(3,extremeslow_wind(:).w),3);
    highmean = nanmean(cat(3,extremeshigh_wind(:).w),3);
    windname_ext = 'wind';
else
    lowmean = nanmean(extremeslow(:,:,:),3);
    highmean = nanmean(extremeshigh(:,:,:),3);
    windname_ext = 'Arcticozoneextremes';
end

lowmean (lowmean < 0) = 0;
highmean (highmean < 0) = 0;

%difference = nanmean(extremeshigh(:,:,:),3) - nanmean(extremeslow(:,:,:),3);
if strcmp(inputs.var,'SNOWHICE')
    %difference = permute(reshape((highmean - lowmean)./lowmean*100,[1,size(highmean)]),[1,3,2]);
    difference = permute(reshape((highmean - lowmean),[1,size(highmean)]),[1,3,2]);
    difference (isnan(difference)) = 0;
elseif strcmp(inputs.var,'ICEFRAC')
    difference = permute(reshape((highmean - lowmean),[1,size(highmean)]),[1,3,2]);
elseif strcmp(inputs.var,'SNOWHLND')
    difference = permute(reshape((highmean - lowmean),[1,size(highmean)]),[1,3,2]);
    difference (isnan(difference)) = 0;
end

differencefinal = difference;

%% calculate significance

extremeslowcomp = extremeslow(:,:,:);
extremeshighcomp = extremeshigh(:,:,:);

p = ttest2(permute(extremeslowcomp,[3,2,1]),permute(extremeshighcomp,[3,2,1]),'Alpha',.1);
    
p (p == 1) = -2;
p (p == 0) = 1;
p (p == -2) = 0.9;
p (isnan(p)) = 1;

%%
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
        extremeslow2(:,:,:,i) = permute(squeeze(nanmean(surfacedata2.highCl.dataMonthArrange(i,inputs.varmonthtomean,pct.highCl.ind.lowind(i,:),:,:),2)),[2,3,1]);
        extremeshigh2(:,:,:,i) = permute(squeeze(nanmean(surfacedata2.highCl.dataMonthArrange(i,inputs.varmonthtomean,pct.highCl.ind.highind(i,:),:,:),2)),[2,3,1]);        
    end

    extremeslowcomp = extremeslow2(:,:,:)+extremeslow(:,:,:);
    extremeshighcomp = extremeshigh2(:,:,:)+extremeshigh(:,:,:);
    
    p = ttest2(permute(extremeslowcomp,[3,2,1]),permute(extremeshighcomp,[3,2,1]),'Alpha',.1);
    
    p (p == 1) = -2;
    p (p == 0) = 1;
    p (p == -2) = 0.9;
    p (isnan(p)) = 1;
    
    lowmean2 = nanmean(extremeslow2(:,:,:),3);
    highmean2 = nanmean(extremeshigh2(:,:,:),3);
    %difference = nanmean(extremeshigh(:,:,:),3) - nanmean(extremeslow(:,:,:),3);
    %difference2 = permute(reshape((highmean2 - lowmean2)./lowmean2*100,[1,size(highmean2)]),[1,3,2]);
    difference2 = permute(reshape((highmean2 - lowmean2),[1,size(highmean2)]),[1,3,2]);
    difference2 (isnan(difference2)) = 0;
    differencefinal = differencefinal+difference2;
    
end
alldifferencelandorice(l,:,:) = difference;
alldifference(l,:,:) = differencefinal;
pfinal(l,:,:) = p;
lowmeanfinal(l,:,:) = lowmean;
highmeanfinal(l,:,:) = highmean;
lowmeanfinal_ind(l,:,:,:,:) = extremeslow;
highmeanfinal_ind(l,:,:,:,:) = extremeshigh;

end
%% plot data
if coincidence
    coincidence_ext = 'windTCOcoincidence';
else
    coincidence_ext = '';
end
plotsig = 0;
plot_concdiff = 0;
if plot_concdiff
    if strcmp(inputs.var,'SNOWHLND')
        if windmetric
            titles = {['Change in ',monthnames(inputs.varmonthtomean,1,'single'),' snow depth due to ',monthnames(inputs.tozmonth,1,'long'), ' ozone extremes']};
        else
            titles = {['(a): ', monthnames(mons(1,:),1,'single')],['(b): ', monthnames(mons(2,:),1,'single')],['(c): ', monthnames(mons(3,:),1,'single')]};
        end
        ctitle = {'m'};
        clims = [-.1 .1];
        titext = 'SNOWLAND';
        patchcoast = 0;
        if combinesnow
            toplotc = alldifference;
        else
            toplotc = alldifferencelandorice;
        end
        toplotp = pfinal;
    elseif strcmp(inputs.var,'SNOWHICE')
        if windmetric
            titles = {['Change in ',monthnames(inputs.varmonthtomean,1,'single'),' snow depth due to ',monthnames(inputs.tozmonth,1,'long'), ' ozone extremes']};
        else
            titles = {['(a): ', monthnames(mons(1,:),1,'single')],['(b): ', monthnames(mons(2,:),1,'single')],['(c): ', monthnames(mons(3,:),1,'single')]};
        end
        ctitle = {'m'};
        clims = [-.2 .2];
        titext = 'SNOWLICE';
        patchcoast = 0;
        if combinesnow
            toplotc = alldifference;
            toplotp = p;
        else
            toplotc = alldifferencelandorice;
            toplotp = [];
        end                    

    elseif strcmp(inputs.var,'ICEFRAC')
        if windmetric
            titles = {['Change in ',monthnames(inputs.varmonthtomean,1,'single'),' sea ice fraction due to winter SSWs and SPVs']};
        else
            %title = {['Change in ',monthnames(inputs.varmonthtomean,1,'single'),' sea ice fraction due to ',monthnames(inputs.tozmonth,1,'long'), ' ozone extremes']};
            titles = {['(a): ', monthnames(mons(1,:),1,'single')],['(b): ', monthnames(mons(2,:),1,'single')],['(c): ', monthnames(mons(3,:),1,'single')]};
        end
        ctitle = {'Ice fraction'};
        clims = [-.3 .3];
        titext = 'ICEFRAC';
        toplotc = alldifference;
        toplotp = pfinal;
        patchcoast = 1;
    end

    toplotc2 = cat(2,toplotc,toplotc(:,1,:));
    longtouse = [longitude;longitude(1)];
    toplotpext = cat(2,toplotp,toplotp(:,1,:));
    filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
%     [fig,sh] = subplotmaps(toplotc2,longtouse,latitude,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',ctitle,'on',...
%         clims,22,[longtouse(1:24:end-1)]-180,[longtouse(1:24:end-1)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic',patchcoast);
    
    ytick = 0:15:90;
    xtick = longtouse(1:24:end-1);
    fsize = 16;
    nlevels = 22;
    cstep = (clims(2)-clims(1))/(nlevels-2);
    contourlevels = [clims(1)-cstep:cstep:clims(2)+cstep];
    cbrew = cbrewernowhite({'div','RdBu'},nlevels);
    
    cbrew = [cbrew];
    
    ylimit = [45 90];
    xlimit = [0 360];
    fig = createfig('largelandscape','on');
    ppos = get(gca,'position');
    for i = 1:size(toplotc,1)
        sp(i) = subplot(1,3,i);
        spos(i,:) = get(sp(i),'position');
        
        if i == 1
            set(sp(i),'position',[spos(i,1)-spos(i,1)./1.2,spos(i,2)-spos(i,2),spos(i,[3,4])*1.3]);
            spos(i,:) = get(sp(i),'position');                
        else            
            set(sp(i),'position',[spos(i-1,1)+spos(i-1,3)+spos(1,1),spos(i,2)-spos(i,2),spos(i,[3,4])*1.3]);
            spos(i,:) = get(sp(i),'position');
        end
        
        
        m_proj('Stereographic','lon',0,'lat',ylimit(2),'rad',abs(ylimit(2)-ylimit(1)))
        [~,h] = m_contourf(longtouse,latitude,squeeze(toplotc2(i,:,:))',contourlevels,'LineStyle','none');                     
        
        m_grid('ytick',double(ytick) ,'xtick',6,'XaxisLocation','bottom','fontsize',fsize);                    
        if i == size(toplotc,1)
            ch = colorbar;
            set(sp(i),'position',spos(i,:));
            set(get(ch,'ylabel'),'string',ctitle,'fontsize',fsize+2)
            cbaxloc = get(ch,'Position');
            set(ch,'YTick',[clims(1):cstep*2:clims(2)],'fontsize',fsize)
            cbarrow;
        end
        colormap(cbrew);          
        caxis([clims(1)-cstep clims(2)+cstep]);
        
        th(i) = title(titles{i},'fontsize',fsize+2);
        thpos = get(th(i),'position');
        set(th(i),'Position',[thpos(1),thpos(2)+.1,thpos(3)]);
        
        if patchcoast
            m_coast('patch',[.5 .5 .5],'edgecolor','none','LineWidth',1);
        else                
            m_coast('color','k','LineWidth',1);
        end
        
        hold on

        % create highlighted significance
        
%         ax2(i).h = axes;        
%         m_proj('Stereographic','lon',0,'lat',ylimit(2),'rad',abs(ylimit(2)-ylimit(1)))        
%         set(ax2(i).h,'position',get(sp(i),'position'));
%         ax2(i).h.Visible = 'off';
%         %m_grid('ytick',double(ytick) ,'xtick',6,'XaxisLocation','bottom','fontsize',fsize);   
%         set(ax2(i).h,'color','none');        
        

        if plotsig
            [c2,h2] = m_contourf(longtouse,latitude,squeeze(toplotpext(i,:,:))',[.9,1],'LineStyle','none','color',[1 1 1]);    
            [c3,h3] = m_contour(longtouse,latitude,squeeze(toplotpext(i,:,:))',[.9,.9],'LineStyle','-','color',[.6 .6 .6],'LineWidth',2);            
            %colormap([[.7 .7 .7];[.7 .7 .7]]);
            drawnow;
            hFills = h2.FacePrims;  % array of TriangleStrip objects

            [hFills(1).ColorType] = deal('truecoloralpha');  % default = 'truecolor'
            [hFills(2).ColorType] = deal('truecoloralpha');  % default = 'truecolor'

            for idx = 1:2 
             hFills(idx).ColorData(4) = 125;   % default=255
             hFills(idx).ColorData(1) = 210;   % default=255
             hFills(idx).ColorData(2) = 210;   % default=255
             hFills(idx).ColorData(3) = 210;   % default=255

            end
        end
        %[c3,h3] = m_contour(longtouse,latitude,squeeze(toplotpext(i,:,:))',[.9,.9],'LineStyle','-','color',[.7 .7 .7],'LineWidth',2);            
        
        
        
    end
    
    
    %set(get(sh,'Title'),'Position',[.5 1.03 0])

    filename = [filedir,sprintf('%02d',inputs.varmonthtomean),'_',titext,'_',...
        monthnames(inputs.tozmonth,1,'single'),'_',windname_ext,'_',...
        num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',coincidence_ext];
    %print(filename,'-dsvg','-painters');
    export_fig(filename,'-png','-r200');
    %print(gcf,filename,'epsc')
    clearvars title
end

%% plot contour lines

titlab = {'a: ','b: ','c: '};
cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewer('qual','Set2',8);
fraction = .15;
plot_iceextent = 0;
filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
if plot_iceextent && seaice
    
    
    createfig('largelandscape','on')
    for i = 1:size(lowmeanfinal,1)
        subplot(1,size(lowmeanfinal,1),i)
        matdefcols = get(groot,'DefaultAxesColorOrder');
               
        m_proj('Stereographic','lon',0,'lat',90,'rad',abs(90-45))
        [~,h] = m_contour(longitude,latitude,squeeze(lowmeanfinal(i,:,:)),[0,fraction]);
        hold on
        h12 = plot(1,1,'LineWidth',3,'color',cbrewqual(1,:));
        set(h,'LineWidth',3,'color',cbrewqual(1,:));    
        [~,h2] = m_contour(longitude,latitude,squeeze(highmeanfinal(i,:,:)),[0,fraction]);
        h22 = plot(1,1,'LineWidth',3,'color',cbrewqual(2,:));
        set(h2,'LineWidth',3,'color',cbrewqual(2,:));
        m_coast('patch',[.5 .5 .5],'edgecolor','none','LineWidth',1);    
        m_grid('ytick',0:15:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',20,'LineWidth',2);    

        tit = {[titlab{i},monthnames(mons(i,:),1,'single')]};     
        th = title(tit);
        set(th,'units','Normalize');
        set(th,'Position',[.5 1.10 0],'fontsize',20);
        if i == 2
            lh = legend([h12,h22],'Low ozone','High ozone');
            set(lh,'box','off','fontsize',20,'orientation','horizontal','location','south');
            lhpos = get(lh,'position');            
            set(lh,'position',[lhpos(1) lhpos(2)-.1 lhpos(3) lhpos(4)]);
        end
        
        
    end
    filename = [filedir,sprintf('%02d',inputs.varmonthtomean),'_iceextenddiff_from_',...
        monthnames(inputs.tozmonth,1,'single'),'_Arcticozoneextremes_over',...
        num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
    export_fig(filename,'-png');
end

%% Zoom in on different seas
Barents = 0;
Seasplot.lon = [25,80;90,235;315,378;100,200];
Seasplot.lat = [65,80;66,80;55,85;40,65];
seanames = {'BKS','LEC','GS','OB'};
fraction  = .15;
if Barents
    
    %x(size(x,1)/2+1:size(x,1)) = x(size(x,1)/2+1:size(x,1)) - 360;
    %x = [x(size(x,1)/2+1:size(x,1));x(1:size(x,1)/2)];
    for i = 1:size(Seasplot.lon,1)
        for j = 1:size(lowmeanfinal,1)
            createfig('medium','on') 
            if i == 3
                datalow = cat(3,lowmeanfinal(j,:,:),lowmeanfinal(j,:,1:11));%circshift(lowmeanfinal(j,:,:),size(lowmeanfinal(j,:,:),3)/2,2);
                datahigh = cat(3,highmeanfinal(j,:,:),highmeanfinal(j,:,1:11));%circshift(highmeanfinal(j,:,:),size(highmeanfinal(j,:,:),3)/2,2);                
                x = [longitude;longitude(1:11)+360];
            else
                datalow = lowmeanfinal(j,:,:);%circshift(lowmeanfinal(j,:,:),size(lowmeanfinal(j,:,:),3)/2,2);
                datahigh = highmeanfinal(j,:,:);%circshift(highmeanfinal(j,:,:),size(highmeanfinal(j,:,:),3)/2,2);
                x = longitude;
            end
            m_proj('lambert','lon',[Seasplot.lon(i,1) Seasplot.lon(i,2)],'lat',[Seasplot.lat(i,1) Seasplot.lat(i,2)]);
            [~,h] = m_contour(x,latitude,squeeze(datalow(1,:,:)),[0,fraction]);
            hold on
            set(h,'LineWidth',3,'color',cbrewqual(1,:));
            hl = plot(1,1,'color',cbrewqual(1,:),'LineWidth',3);

            [~,h2] = m_contour(x,latitude,squeeze(datahigh(1,:,:)),[0,fraction]);
            set(h2,'LineWidth',3,'color',cbrewqual(2,:));
            hl2 = plot(1,1,'color',cbrewqual(2,:),'LineWidth',3);
            m_coast('patch',cbrewqual2(7,:));
            m_grid('box','fancy','tickdir','in');
            th = title(['Ensemble average ',monthnames(mons(j,:),1,'single'),' sea ice extent relative to ',...
                monthnames(inputs.tozmonth,1,'long'),' Arctic ozone extemes'],'fontsize',20);
            set(th,'units','normalized')
            set(th,'position',[.5 1.045 0])
            filename = [filedir,sprintf('%02d',j),'_iceextenddiff_from_',...
                monthnames(inputs.tozmonth,1,'long'),'_Arcticozoneextremes_over_',seanames{i},...
                num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
            lh = legend([hl,hl2],'relative to TCO lower 20th percentile','relative to TCO upper 20th percentile');
            set(lh,'box','off','fontsize',20,'location','southoutside')
            export_fig(filename,'-png');
            close 1
        end
    end
end    
    
%% seas for accumulation metric

% Barents Kara sea; Laptet, East Siberian, and Chukchi seas, Greenland
% Seas; Sea of Okhotsk and Berling Sea

% Calculating the mean
toplot = 1;

% flag locations below threshold
%% Read in area grid cells
areacell = gridcellarea(latitude,longitude);

%% flag ensemble and individuals
% ensemble
below.low = lowmeanfinal;
below.high = highmeanfinal;
below.low (lowmeanfinal < .15) = 0;
below.low (lowmeanfinal >= .15) = 1;
below.high (highmeanfinal < .15) = 0;
below.high (highmeanfinal >= .15) = 1;

% individual
lowmeanind = squeeze(nanmean(lowmeanfinal_ind,4)); %[months,lat,lon,fileno.]
highmeanind = squeeze(nanmean(highmeanfinal_ind,4));

indbelow.low = lowmeanind;
indbelow.high = highmeanind;
indbelow.low (lowmeanind < .15) = 0;
indbelow.low (lowmeanind >= .15) = 1;
indbelow.high (highmeanind < .15) = 0;
indbelow.high (highmeanind >= .15) = 1;

%% difference in sea ice area
below.area.low = lowmeanfinal;
below.area.high = highmeanfinal;
below.area.low (lowmeanfinal < .15) = 0;
below.area.high (highmeanfinal < .15) = 0;

indbelow.area.low = lowmeanind;
indbelow.area.high = highmeanind;
indbelow.area.low (lowmeanind < .15) = 0;
indbelow.area.high (highmeanind < .15) = 0;

%% calculating sea ice extent for different regions
% Barents Kara sea; Laptet, East Siberian, and Chukchi seas, Greenland
% Seas; Sea of Okhotsk and Berling Sea
Seasplot.lon = [25,80;90,235;310,18;100,200]; %310,378
Seasplot.lat = [65,80;66,80;55,85;40,65];
fraction = .15;

for i = 1:size(Seasplot.lon,1)
    latind = latitude >= Seasplot.lat(i,1) & latitude <= Seasplot.lat(i,2);
    if i == 3
        lonind = longitude >= Seasplot.lon(i,1) | longitude < Seasplot.lon(i,2);
    else
        lonind = longitude >= Seasplot.lon(i,1) & longitude < Seasplot.lon(i,2);
    end
    seaiceextent.ensemble.low(i).a = permute(areacell(latind,lonind).*permute(below.low(:,latind,lonind),[2,3,1]),[3,1,2]);
    seaiceextent.ensemble.low(i).a = sum(seaiceextent.ensemble.low(i).a(:,:),2);
    seaiceextent.ensemble.high(i).a = permute(areacell(latind,lonind).*permute(below.high(:,latind,lonind),[2,3,1]),[3,1,2]);
    seaiceextent.ensemble.high(i).a = sum(seaiceextent.ensemble.high(i).a(:,:),2);
    seaiceextent.ensemble.difference(:,i) = seaiceextent.ensemble.high(i).a - seaiceextent.ensemble.low(i).a;
    
    seaiceextent.ind.low(i).a = permute(areacell(latind,lonind).*permute(indbelow.low(:,latind,lonind,:),[2,3,4,1]),[4,3,1,2]);
    seaiceextent.ind.low(i).a = sum(seaiceextent.ind.low(i).a(:,:,:),3);
    seaiceextent.ind.high(i).a = permute(areacell(latind,lonind).*permute(indbelow.high(:,latind,lonind,:),[2,3,4,1]),[4,3,1,2]);
    seaiceextent.ind.high(i).a = sum(seaiceextent.ind.high(i).a(:,:,:),3);
    seaiceextent.ind.difference(:,:,i) = seaiceextent.ind.high(i).a - seaiceextent.ind.low(i).a;
    
    seaicearea.ensemble.low(i).a = permute(areacell(latind,lonind).*permute(below.area.low(:,latind,lonind),[2,3,1]),[3,1,2]);
    seaicearea.ensemble.low(i).a = sum(seaicearea.ensemble.low(i).a(:,:),2);
    seaicearea.ensemble.high(i).a = permute(areacell(latind,lonind).*permute(below.area.high(:,latind,lonind),[2,3,1]),[3,1,2]);
    seaicearea.ensemble.high(i).a = sum(seaicearea.ensemble.high(i).a(:,:),2);
    seaicearea.ensemble.difference(:,i) = seaicearea.ensemble.high(i).a - seaicearea.ensemble.low(i).a;
    
    seaicearea.ind.low(i).a = permute(areacell(latind,lonind).*permute(indbelow.area.low(:,latind,lonind,:),[2,3,4,1]),[4,3,1,2]);
    seaicearea.ind.low(i).a = sum(seaicearea.ind.low(i).a(:,:,:),3);
    seaicearea.ind.high(i).a = permute(areacell(latind,lonind).*permute(indbelow.area.high(:,latind,lonind,:),[2,3,4,1]),[4,3,1,2]);
    seaicearea.ind.high(i).a = sum(seaicearea.ind.high(i).a(:,:,:),3);
    seaicearea.ind.difference(:,:,i) = seaicearea.ind.high(i).a - seaicearea.ind.low(i).a;
    
end

seaicearea.ind.std = squeeze(std(seaicearea.ind.difference,0,2));
seaiceextent.ind.std = squeeze(std(seaiceextent.ind.difference,0,2));
%% plot
createfig('medium','on')
lwidth = 3;
fsize = 18;
colors = cbrewer('qual','Set1',10);
Seanames = {'Barents-Kara sea','Laptet, East Siberian, and Chukchi seas','Greenland Sea','Sea of Okhotsk and Berling Sea'};
lims = [[-2e5 2e5];[-2e5 4e5];[-2.5e5 1e5];[-2e5 8.5e5]];
yticks = [.5e5,1e5,.5e5,1e5];
% sea ice extent (sum up grid point area for all grid points above .15)
for i = 1:size(Seasplot.lon,1)
    sp(i) = subplot(2,2,i);
    sppos(i,:) = get(sp(i),'position');
    
    %MOVE PLOTS TO MAKE LOOK NICE
    plot([-1 13],[0 0],'--k');
    hold on
    eh(i) = errorbar(mons,seaiceextent.ensemble.difference(:,i),seaiceextent.ind.std(:,i),'color',colors(1,:),'LineWidth',lwidth-1);    
    ph(i) = plot(mons,seaiceextent.ensemble.difference(:,i),'LineWidth',lwidth,'Color',...
        colors(1,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
    
    set(gca,'xtick',mons,'xticklabel',mons,'ytick',-8e5:yticks(i):8e5,'yticklabel',-8:yticks(i)/1e5:8,'fontsize',fsize);
    xlim([.5 12.5])
    ylim(lims(i,:));
    if i == 1 || i == 3
        ylabel(['Sea ice extent( ',char(215),'10^5 km^2)'],'fontsize',fsize+2)
    end
    if i == 3 || i == 4
        xlabel('Month','fontsize',fsize+2);        
    end
    title(Seanames{i},'fontsize',fsize);
end

% sea ice area (sum up fraction of sea ice for all grid point above .15)
createfig('medium','on')
lims = [[-1.25e5 1.75e5];[-.5e5 4e5];[-1.5e5 1e5];[-1e5 5e5]];
yticks = [.5e5,.5e5,.5e5,1e5];
for i = 1:size(Seasplot.lon,1)
    subplot(2,2,i)
    plot([-1 13],[0 0],'--k');
    hold on
    eh(i) = errorbar(mons,seaicearea.ensemble.difference(:,i),seaicearea.ind.std(:,i),'color',colors(1,:),'LineWidth',lwidth-1);    
    ph(i) = plot(mons,seaicearea.ensemble.difference(:,i),'LineWidth',lwidth,'Color',...
        colors(1,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
    set(gca,'xtick',mons,'xticklabel',mons,'ytick',-8e5:yticks(i):8e5,'yticklabel',-8:yticks(i)/1e5:8,'fontsize',fsize);
    xlim([.5 12.5])
    ylim(lims(i,:));
    if i == 1 || i == 3
        ylabel(['Sea ice area( ',char(215),'10^5 km^2)'],'fontsize',fsize+2)
    end
    if i == 3 || i == 4
        xlabel('Month','fontsize',fsize+2);        
    end
    title(Seanames{i},'fontsize',fsize);
    
end

% both on same graph
createfig('medium','on')
lims = [[-2e5 2e5];[-2e5 4e5];[-2.5e5 1e5];[-2e5 8.5e5]];
yticks = [.5e5,1e5,.5e5,1e5];
for i = 1:size(Seasplot.lon,1)
    subplot(2,2,i)
    plot([-1 13],[0 0],'--k');
    hold on
    eh(i) = errorbar(mons,seaiceextent.ensemble.difference(:,i),seaiceextent.ind.std(:,i),'color',colors(1,:),'LineWidth',lwidth-1);    
    ph(i) = plot(mons,seaiceextent.ensemble.difference(:,i),'LineWidth',lwidth,'Color',...
        colors(1,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
    eh(i) = errorbar(mons,seaicearea.ensemble.difference(:,i),seaicearea.ind.std(:,i),'color',colors(2,:),'LineWidth',lwidth-1);    
    ph(i) = plot(mons,seaicearea.ensemble.difference(:,i),'LineWidth',lwidth,'Color',...
        colors(2,:),'Marker','o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
    
    
    set(gca,'xtick',mons,'xticklabel',mons,'ytick',-8e5:yticks(i):8e5,'yticklabel',-8:yticks(i)/1e5:8,'fontsize',fsize);
    xlim([.5 12.5])
    ylim(lims(i,:));
    if i == 1 || i == 3
        ylabel(['Sea ice area( ',char(215),'10^5 km^2)'],'fontsize',fsize+2)
    end
    if i == 3 || i == 4
        xlabel('Month','fontsize',fsize+2);        
    end
    title(Seanames{i},'fontsize',fsize);
    
end

%[areadifferencecorrect] = calculatearea(lowmeanfinal,highmeanfinal,longitude,latitude,toplot,[2:4:36],[2:4:18,26:4:36]);

%% Calculating individuals
toplot = 1;
tactogetherind = [2:4:36];
tactogetherindlow = [2:4:36];
for i = 1:nofiles    
    lowmeanind = squeeze(nanmean(lowmeanfinal_ind(:,:,:,:,i),4));
    highmeanind = squeeze(nanmean(highmeanfinal_ind(:,:,:,:,i),4));
    [areadifferencecorrectind(i,:,:)] = calculatearea(lowmeanind,highmeanind,longitude,latitude,toplot,tactogetherind,tactogetherindlow);
    pause;
end

%% plot differences

createfig('medium','on');

% areadifferencecorrectlog = areadifferencecorrect;
% areadifferencecorrectlog (areadifferencecorrectlog < 0) = -log10(abs(areadifferencecorrect(areadifferencecorrectlog < 0)));
% areadifferencecorrectlog (areadifferencecorrect > 0) = log10(areadifferencecorrect(areadifferencecorrect > 0));
for i = 1:size(areadifferencecorrect,2)
    subplot(2,2,i)
    plot(mons,areadifferencecorrect(:,i),'LineWidth',3);
end
%set(gca, 'YScale', 'log')
