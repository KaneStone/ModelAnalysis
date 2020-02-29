clear variables

%% user inputs

inputs = predruns_inputs;
inputs.compareSSW = 0;
inputs.fraction = .15;

%% Read in model (user inputs located in function)

% Barents Kara sea; Laptet, East Siberian, and Chukchi seas, Greenland
% Seas; Sea of Okhotsk and Berling Sea

if strcmp(inputs.var,'ICEFRAC')
    seaice = 1;
%     Seasplot.lon = [20,94;95,235;310,380;130,200];
%     Seasplot.lat = [65,80;64,80;55,85;45,65];

    Seasplot.lon = [20,100;100,235;310,380]; %Barents was 94
    Seasplot.lat = [65,85;45,85;55,85];
    
%     Seasplot.lon = [50,100;100,235;310,410]; %Barents was 94
%     Seasplot.lat = [65,85;45,85;55,85];
    
% including Canada
    %Seasplot.lon = [20,94;95,235;310,380;235,310];
    %Seasplot.lat = [65,80;45,80;55,85;45,80];
    
% inlcuding Bering pass    
%     Seasplot.lon = [20,94;95,235;310,380;180,200];
%     Seasplot.lat = [65,80;45,80;55,85;55,75];
    
else
    seaice = 0;
%     Seasplot.obs.lon = [93,163,205];
%     Seasplot.obs.lat = [65,63,60];
%     Seasplot.mod.lon = [93,163,205];
%     Seasplot.mod.lat = [65,62,60];
%     Seasplot.lon = [75,100;155,170;205,235]; % Northern Siberia, Eastern Siberia, Alaska
%     Seasplot.lat = [64,75;60,65;48,64];

    Seasplot.lon = [75,100;75,115;205,235]; % Northern Siberia, Southern Siberia, Alaska
    Seasplot.lat = [64,75;52,60;50,64];

    Seasplot.obs.lon = Seasplot.lon;
    Seasplot.obs.lat = Seasplot.lat;
    Seasplot.mod.lon = Seasplot.lon;
    Seasplot.mod.lat = Seasplot.lat;
    
    
end

%% read in model data
[surfacedata,tozdata] = predruns_ReadinModel(inputs);

% read in surface temperature calculating ENSO
inputemp = inputs.var;
inputs.var = 'TS';
[modeltemp,~] = predruns_ReadinModel(inputs);
inputs.var = inputemp;

surfacedata.highCl.dataMonthArrangeTS = modeltemp.highCl.dataMonthArrange;

modeldim.lat = surfacedata.highCl.data(1).lat ; 
modeldim.lon = surfacedata.highCl.data(1).lon;

% obtaining number and names of fields 
fields = fieldnames(surfacedata);
nofields = length(fields);

%% remove enso
if inputs.removeENSO
    [~,surfacedata.highCl.dataMonthArrange_re,surfacedata.highCl.ENSO] = predruns_removeENSO(surfacedata.highCl.dataMonthArrange,surfacedata.highCl.dataMonthArrangeTS,modeldim.lat,modeldim.lon,inputs,'highCl');
    surfacedata.highCl.dataMonthArrange_re = permute(surfacedata.highCl.dataMonthArrange_re,[1,3,2,4,5]);
    surfacedata.highCl.dataMonthArrange_re = cat(2,surfacedata.highCl.dataMonthArrange_re(:,1:2,:,:,:),surfacedata.highCl.dataMonthArrange_re); 
end

%% Read in observed arctic data and Bodeker scientific total column ozone
inputs.obstouse = 'MERRA'; %'MERRA', 'Goddard', or 'ERA'
if inputs.compareERA    
    %first input is either 'seaiceindex','Goddard', or 'ERA'
    
    % read in temperature
    if inputs.removeENSO
        ERAtempall = ReadinERA('/Volumes/ExternalOne/work/data/ERA-Interim/TS/TS_ERA-Interim.nc');
        observedTS = ERAtempall.t2m(:,:,13:end-2);
%         for i = 1:12            
%             ERAtemp(i,:,:,:) = ERAtempall.t2m(:,:,i:12:end);
%         end
%         ERAtemp = permute(ERAtemp,[1,4,2,3]);
    else
        ERAtemp = [];
    end
    [observations] = predruns_readinObservedIceandSnow(inputs,modeldim,observedTS);
    observations.TS = observedTS;
else
    observations.latitude = [90:-1.5:-90]';
    observations.longitude = [0:1.5:360-1.5]';
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
inputs.seasons = 1;
if inputs.seasons 
    mons = [4,5,6;7,8,9;10,11,12];
    %mons = [4;9;12];
else
    mons = inputs.varmonth';
    %mons = [1:12]';
end

% Read in SSW and ozone coincidence extreme
inputs.coincidence = 0;
inputs.nocoincidence = 0;

% calculate difference in surface metric based on input extremes (either ozone or wind)
differences = Seaice_TakeDifference(inputs,mons,surfacedata,nofiles,pct,observations,[],.05);

%% plotting individuals
plotind = 0;
combine = 0;
if plotind
    Seaice_PlotSeaIceIndividual(inputs,differences.difference.ind,latitude,longitude,mons,seasons)
end

%% plotting ensemble
plotens = 1;
combine = 1;
if inputs.seasons && plotens    
    Seaice_PlotSeaIceEnsemble(inputs,differences,latitude,longitude,observations.latitude,...
        observations.longitude,mons,combine,Seasplot)
end

%% plot sea ice extent
lineplots = 1;
if lineplots
    if ~inputs.seasons && strcmp(inputs.var,'ICEFRAC') % do not detrend sea ice in inputs, but do detrend ozone.
        inputs.clines = 0;
        Seaice_Extent(differences,Seasplot,latitude,longitude,observations.latitude,...
            observations.longitude,inputs,mons,surfacedata,observations,pct);
    elseif ~inputs.seasons
        Seaice_SnowAverage(differences,Seasplot,latitude,longitude,observations.latitude,...
            observations.longitude,inputs,mons,surfacedata,observations,pct);
    end
end

%% construct regression model line plots
regline = 1;
if regline
    predruns_cryregression(differences,Seasplot,latitude,longitude,observations.latitude,...
            observations.longitude,inputs,mons,surfacedata,observations,tozdata.highCl.dataMonthArrange,observations.toz.zm)
end

%% plot comparison line plots for all ozone and sea ice extent

seaice_comparisonLinePlots(surfacedata.highCl.dataMonthArrange,observations,...
    tozdata,observations.toz.zm,Seasplot,latitude,longitude,observations.latitude,observations.longitude,inputs,pct);

%% prediction

seaice_leaveout_prediction(surfacedata,tozdata,inputs,latitude,longitude,differences,observations,Seasplot)

%dataVarMonthAve,dataVarMonth,tozdata.(inputs.ClLevel{1}).dataMonthArrange,inputs,latitude,longitude,pct,differences.indmonths.individual,surfacedata.(fields{i}).dataMonthArrange
%% plot data

% I am leaving this in for now as I am going to need the significance stuff
% later.

% if coincidence
%     coincidence_ext = 'windTCOcoincidence';
% else
%     coincidence_ext = '';
% end
% plotsig = 0;
% plot_concdiff = 0;
% if plot_concdiff
%     if strcmp(inputs.var,'SNOWHLND')
%         if windmetric
%             titles = {['Change in ',monthnames(inputs.varmonthtomean,1,'single'),' snow depth due to ',monthnames(inputs.tozmonth,1,'long'), ' ozone extremes']};
%         else
%             titles = {['(a): ', monthnames(mons(1,:),1,'single')],['(b): ', monthnames(mons(2,:),1,'single')],['(c): ', monthnames(mons(3,:),1,'single')]};
%         end
%         ctitle = {'m'};
%         clims = [-.1 .1];
%         titext = 'SNOWLAND';
%         patchcoast = 0;
%         if combinesnow
%             toplotc = alldifference;
%         else
%             toplotc = alldifferencelandorice;
%         end
%         toplotp = pfinal;
%     elseif strcmp(inputs.var,'SNOWHICE')
%         if windmetric
%             titles = {['Change in ',monthnames(inputs.varmonthtomean,1,'single'),' snow depth due to ',monthnames(inputs.tozmonth,1,'long'), ' ozone extremes']};
%         else
%             titles = {['(a): ', monthnames(mons(1,:),1,'single')],['(b): ', monthnames(mons(2,:),1,'single')],['(c): ', monthnames(mons(3,:),1,'single')]};
%         end
%         ctitle = {'m'};
%         clims = [-.2 .2];
%         titext = 'SNOWLICE';
%         patchcoast = 0;
%         if combinesnow
%             toplotc = alldifference;
%             toplotp = p;
%         else
%             toplotc = alldifferencelandorice;
%             toplotp = [];
%         end                    
% 
%     elseif strcmp(inputs.var,'ICEFRAC')
%         if windmetric
%             titles = {['Change in ',monthnames(inputs.varmonthtomean,1,'single'),' sea ice inputs.fraction due to winter SSWs and SPVs']};
%         else
%             %title = {['Change in ',monthnames(inputs.varmonthtomean,1,'single'),' sea ice inputs.fraction due to ',monthnames(inputs.tozmonth,1,'long'), ' ozone extremes']};
%             titles = {['(a): ', monthnames(mons(1,:),1,'single')],['(b): ', monthnames(mons(2,:),1,'single')],['(c): ', monthnames(mons(3,:),1,'single')]};
%         end
%         ctitle = {'Ice inputs.fraction'};
%         clims = [-.3 .3];
%         titext = 'ICEFRAC';
%         toplotc = alldifference;
%         toplotp = pfinal;
%         patchcoast = 1;
%     end
% 
%     toplotc2 = cat(2,toplotc,toplotc(:,1,:));
%     longtouse = [longitude;longitude(1)];
%     toplotpext = cat(2,toplotp,toplotp(:,1,:));
%     filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
% %     [fig,sh] = subplotmaps(toplotc2,longtouse,latitude,{'div','RdBu'},1,[],16,title,'Longitude','Latitude',ctitle,'on',...
% %         clims,22,[longtouse(1:24:end-1)]-180,[longtouse(1:24:end-1)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic',patchcoast);
%     
%     ytick = 0:15:90;
%     xtick = longtouse(1:24:end-1);
%     fsize = 16;
%     nlevels = 22;
%     cstep = (clims(2)-clims(1))/(nlevels-2);
%     contourlevels = [clims(1)-cstep:cstep:clims(2)+cstep];
%     cbrew = cbrewernowhite({'div','RdBu'},nlevels);
%     
%     cbrew = [cbrew];
%     
%     ylimit = [45 90];
%     xlimit = [0 360];
%     fig = createfig('largelandscape','on');
%     ppos = get(gca,'position');
%     for i = 1:size(toplotc,1)
%         sp(i) = subplot(1,3,i);
%         spos(i,:) = get(sp(i),'position');
%         
%         if i == 1
%             set(sp(i),'position',[spos(i,1)-spos(i,1)./1.2,spos(i,2)-spos(i,2),spos(i,[3,4])*1.3]);
%             spos(i,:) = get(sp(i),'position');                
%         else            
%             set(sp(i),'position',[spos(i-1,1)+spos(i-1,3)+spos(1,1),spos(i,2)-spos(i,2),spos(i,[3,4])*1.3]);
%             spos(i,:) = get(sp(i),'position');
%         end
%         
%         
%         m_proj('Stereographic','lon',0,'lat',ylimit(2),'rad',abs(ylimit(2)-ylimit(1)))
%         [~,h] = m_contourf(longtouse,latitude,squeeze(toplotc2(i,:,:))',contourlevels,'LineStyle','none');                     
%         
%         m_grid('ytick',double(ytick) ,'xtick',6,'XaxisLocation','bottom','fontsize',fsize);                    
%         if i == size(toplotc,1)
%             ch = colorbar;
%             set(sp(i),'position',spos(i,:));
%             set(get(ch,'ylabel'),'string',ctitle,'fontsize',fsize+2)
%             cbaxloc = get(ch,'Position');
%             set(ch,'YTick',[clims(1):cstep*2:clims(2)],'fontsize',fsize)
%             cbarrow;
%         end
%         colormap(cbrew);          
%         caxis([clims(1)-cstep clims(2)+cstep]);
%         
%         th(i) = title(titles{i},'fontsize',fsize+2);
%         thpos = get(th(i),'position');
%         set(th(i),'Position',[thpos(1),thpos(2)+.1,thpos(3)]);
%         
%         if patchcoast
%             m_coast('patch',[.5 .5 .5],'edgecolor','none','LineWidth',1);
%         else                
%             m_coast('color','k','LineWidth',1);
%         end
%         
%         hold on
% 
%         % create highlighted significance
%         
% %         ax2(i).h = axes;        
% %         m_proj('Stereographic','lon',0,'lat',ylimit(2),'rad',abs(ylimit(2)-ylimit(1)))        
% %         set(ax2(i).h,'position',get(sp(i),'position'));
% %         ax2(i).h.Visible = 'off';
% %         %m_grid('ytick',double(ytick) ,'xtick',6,'XaxisLocation','bottom','fontsize',fsize);   
% %         set(ax2(i).h,'color','none');        
%         
% 
%         if plotsig
%             [c2,h2] = m_contourf(longtouse,latitude,squeeze(toplotpext(i,:,:))',[.9,1],'LineStyle','none','color',[1 1 1]);    
%             [c3,h3] = m_contour(longtouse,latitude,squeeze(toplotpext(i,:,:))',[.9,.9],'LineStyle','-','color',[.6 .6 .6],'LineWidth',2);            
%             %colormap([[.7 .7 .7];[.7 .7 .7]]);
%             drawnow;
%             hFills = h2.FacePrims;  % array of TriangleStrip objects
% 
%             [hFills(1).ColorType] = deal('truecoloralpha');  % default = 'truecolor'
%             [hFills(2).ColorType] = deal('truecoloralpha');  % default = 'truecolor'
% 
%             for idx = 1:2 
%              hFills(idx).ColorData(4) = 125;   % default=255
%              hFills(idx).ColorData(1) = 210;   % default=255
%              hFills(idx).ColorData(2) = 210;   % default=255
%              hFills(idx).ColorData(3) = 210;   % default=255
% 
%             end
%         end
%         %[c3,h3] = m_contour(longtouse,latitude,squeeze(toplotpext(i,:,:))',[.9,.9],'LineStyle','-','color',[.7 .7 .7],'LineWidth',2);            
%         
%         
%         
%     end
%     
%     
%     %set(get(sh,'Title'),'Position',[.5 1.03 0])
% 
%     filename = [filedir,sprintf('%02d',inputs.varmonthtomean),'_',titext,'_',...
%         monthnames(inputs.tozmonth,1,'single'),'_',windname_ext,'_',...
%         num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',coincidence_ext];
%     %print(filename,'-dsvg','-painters');
%     export_fig(filename,'-png','-r200');
%     %print(gcf,filename,'epsc')
%     clearvars title
% end




