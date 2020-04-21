function [] = Seaice_PlotSeaIceIndividual(inputs,difference,latitude,longitude,mons,season,combine)

if season
    montype = 'single';
else
    montype = 'long';
end

if strcmp(inputs.var,'ICEFRAC')
    varname = 'sea ice fraction'; 
    clims = [-.25 .25];
    ctitle = {'Sea ice fraction'};
    patchcoast = 1;
else
    varname = 'snow cover';
    clims = [-.1 .1];
    ctitle = {'Meters'};   
    patchcoast = 0;
end

%% plotting all individuals on one
differencemean = squeeze(difference(:,:,2,:));
longtouse = [longitude;longitude(1)+360];
differencemean = cat(2,differencemean,differencemean(:,1,:));

filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/Individual/';
        
ytick = 60:15:90;
%xtick = longtouse(1:24:end-1);
fsize = 16;
nlevels = 22;
cstep = (clims(2)-clims(1))/(nlevels-2);
contourlevels = [clims(1)-cstep:cstep:clims(2)+cstep];
cbrew = cbrewernowhite({'div','RdBu'},nlevels);

ylimit = [44 90];
xlimit = [0 360];
createfig('largesquare','on');
for j = 1:9
    % plot observations
    sp(j) = subplot(3,3,j);
    sppos(j,:) = get(sp(j),'position');
    if j == 2 || j == 5 || j == 8
        set(sp(j),'position',[sppos(j,1)-.075,sppos(j,2),sppos(j,3:4)]);
    elseif j == 3 || j == 6 || j == 9
        set(sp(j),'position',[sppos(j,1)-.15,sppos(j,2),sppos(j,3:4)]);
    end
    sppos(j,:) = get(sp(j),'position');
    m_proj('Stereographic','lon',0,'lat',ylimit(2),'rad',abs(ylimit(2)-ylimit(1)));    
    axm = axesm('MapProjection','stereo','MapLatLimit',[45 90],'Grid','On','MLineLocation',60,...
        'PlineLocation',15,'MLabelLocation',60,'PLabelLocation',[60,75],'MeridianLabel','on','ParallelLabel','on','MLabelParallel',45,'LabelRotation','on','fontsize',fsize-2); 
    %m_grid('ytick',double(ytick) ,'xtick',6,'XaxisLocation','bottom','fontsize',fsize-2,'box','off');    
    framem; axis off; tightmap;
    
    contourfm(latitude,longtouse,differencemean(:,:,j),contourlevels,'LineStyle','none');
    colormap(cbrew);      
    caxis([clims(1)-cstep clims(2)+cstep]);        
    th = title(['No. ',sprintf('%02d',j)],'fontsize',fsize+2); 
    thpos = get(th,'position');
    set(th,'position',[thpos(1),thpos(2)+.1,thpos(3)]);
    load coastlines
    geoshow(coastlat, coastlon, 'Color', 'black')
    if patchcoast
        patchm(coastlat,coastlon,[.5 .5 .5])                                 
    else
        linem(coastlat,coastlon,'Color','k')  
    end

end

ch = colorbar;
set(sp(end),'position',sppos(end,:));
set(get(ch,'ylabel'),'string',ctitle,'fontsize',fsize+2)
cbaxloc = get(ch,'Position');
set(ch,'Position',[cbaxloc(1)+.01,cbaxloc(2),cbaxloc(3),.79]);
set(ch,'YTick',[clims(1):cstep*2:clims(2)],'fontsize',fsize)      
set(ch,'YaxisLocation','right','fontsize',fsize);          
cbarrow;
%axpos = tightplots(gcf,10,25,.05,.15);

minx = min(sppos(:,1));
maxx = max(sppos(:,1))+sppos(1,3);
maxy = max(sppos(:,2))+sppos(1,4);

annotation('textbox',[minx 1 maxx-minx 0],'String',...
    ['Change in JAS ','sea ice fraction due to ',monthnames(inputs.tozmonth,1,'long'),' ozone extremes'],...
    'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
    'EdgeColor','none','fontweight','bold'); 

% ch = subplotmaps(differencemean,longtouse,latitude,{'div','RdBu'},1,[],16,title,[],[],ctitle,'on',...
%     clims,18,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic',patchcoast);
% gh = get(gca,'title');
% ghpos = get(gh,'position');
% set(gh,'Position',[ghpos(1) ghpos(2)+ghpos(2)./50  ghpos(3)])     
filename = [filedir,'All_',monthnames([7,8,9],1,'single'),'_',inputs.var,'_',...
    monthnames(inputs.tozmonth,1,'long'),'_Arcticozoneextremes_over',...
    num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
%print(filename,'-depsc');
print(gcf,filename,'-depsc','-painters');
%export_fig(filename,'-png');

%%
for i = 1:size(difference,3)
    for j = 1:size(difference,4)
        longtouse = [longitude;longitude(1)+360];

        titles = {['No. ',sprintf('%02d',j),', Change in ',...
            monthnames(mons(i),1,montype),' ',...
            varname,' due to ',monthnames(inputs.tozmonth,1,'long'), ' ozone extremes']};        
        if combine
             toplotdifferencefinal = permute(difference(:,:,:,j),[3,2,1]);
        else
            toplotdifferencefinal = permute(difference(:,:,i,j),[3,2,1]);
        end
        toplotdifferencefinal = cat(2,toplotdifferencefinal,toplotdifferencefinal(:,1,:));

        filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/Individual/';
        
        ch = subplotmaps(toplotdifferencefinal,longtouse,latitude,{'div','RdBu'},1,[],16,title,[],[],ctitle,'on',...
            clims,18,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic',patchcoast);
        gh = get(gca,'title');
        ghpos = get(gh,'position');
        set(gh,'Position',[ghpos(1) ghpos(2)+ghpos(2)./50  ghpos(3)])     
        filename = [filedir,'Ind_','No.',sprintf('%02d',j),'_',sprintf('%02d',mons(i)),'_',inputs.var,'_',...
            monthnames(inputs.tozmonth,1,'long'),'_Arcticozoneextremes_over',...
            num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
        %print(filename,'-depsc');
        export_fig(filename,'-png');
        close 1
    end
    if combine
        break;
    end
end
end