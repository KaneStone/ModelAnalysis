function [] = Seaice_PlotSeaIceEnsemble(inputs,difference,latitude,longitude,obslat,obslon,mons,combine,Seasplot)

% convert to model dimensions
for i = 1:size(difference.observations.difference,1)
    for k = 1:size(difference.observations.difference,2)
        templatinterp(i,k,:) = interp1(obslat,squeeze(difference.observations.difference(i,k,:)),latitude);
    end
end
for i = 1:size(difference.observations.difference,1)
    for k = 1:size(templatinterp,3)
        difference.observations.difference_modeldim(i,:,k) = interp1(obslon,squeeze(templatinterp(i,:,k)),longitude);
    end
end

combineddifference_modeldim = cat(1,permute(difference.difference.ens,[3,2,1]),difference.observations.difference_modeldim);

if inputs.seasons
    montype = 'single';
else
    montype = 'long';
end

if strcmp(inputs.var,'ICEFRAC')
    varname = 'sea ice fraction'; 
    clims = [-.25 .25];
    ctitle = {'Sea ice fraction'};
    patchcoast = 1;
    titext = 'sea ice fraction';
elseif strcmp(inputs.var,'SNOWHICE')
    varname = 'snow cover';
    clims = [-.1 .1];
    ctitle = {'Meters'};   
    patchcoast = 0;
    titext = 'snow depth';
elseif strcmp(inputs.var,'PRECSC')
    varname = 'convective snow';
    clims = [-5e-10 5e-10];
    ctitle = {'m/s'};   
    patchcoast = 0;
    titext = 'convectivesnow';
elseif strcmp(inputs.var,'PRECSL')
    varname = 'large scale snow';
    clims = [-5e-9 5e-9];
    ctitle = {'m/s'};   
    patchcoast = 0;
    titext = 'largescalesnow';
end

ytick = 60:15:90;
%xtick = longtouse(1:24:end-1);
fsize = 16;
nlevels = 22;
cstep = (clims(2)-clims(1))/(nlevels-2);
contourlevels = [clims(1)-cstep:cstep:clims(2)+cstep];
cbrew = cbrewernowhite({'div','RdBu'},nlevels);

ylimit = [44 90];
xlimit = [0 360];
spind = {'d','e','a','f','b','c'};
% Barents Kara sea; Laptet, East Siberian, and Chukchi seas, Greenland
% Seas; Sea of Okhotsk and Berling Sea

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([4,1,3,8],:);


%difference.observations.ttest = cat(2,difference.observations.ttest,difference.observations.ttest(:,1,:));
difference.difference.ttest = permute(difference.difference.ttest,[3,2,1]);

for i = 1:size(combineddifference_modeldim,1)
    longtouse = [longitude;longitude(1)+360];   
    fig = figure;
    set(fig,'position',[100 100 1200 900],'color','white');
    if combine
        %toplotdifferencefinal = permute(combineddifference,[3,2,1]); 
        obsfinal = permute(cat(2,difference.observations.difference,difference.observations.difference(:,1,:)),[2,3,1]);
        
        obsfinal (obsfinal >= clims(2)) = clims(2);
        obsfinal (obsfinal < clims(1)) = clims(1);
        
        modelfinal = permute(cat(2,difference.difference.ens,difference.difference.ens(:,1,:)),[2,1,3]);
        
        modelfinal (modelfinal >= clims(2)) = clims(2);
        modelfinal (modelfinal < clims(1)) = clims(1);
        
        modellontoplot = [longitude;longitude(1)+360]; 
        obslontoplot =  [obslon;obslon(1)+360]; 
        slen = size(modelfinal,3);
        for j = 1:size(mons,1)            
            titles{j} = {['Observations: ',monthnames(mons(j,:),1,montype)]};                     
            titles{j+slen} = {['Model: ',monthnames(mons(j,:),1,montype)]};                     
        end          
    else
        %toplotdifferencefinal = permute(combineddifference(:,:,i),[3,2,1]);        
        titles = {['Ensemble ','Change in ',monthnames(mons(i,:),1,montype),...
            titext,' due to ',monthnames(inputs.tozmonth,1,'long'), ' ozone extremes']};    
    end
    %toplotdifferencefinal = cat(2,combineddifference_modeldim,combineddifference_modeldim(:,1,:));
    
    filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/Ensemble/';
                        
    if combine
        for j = 1:slen
            % plot observations
            sp(j) = subplot(2,3,j);
            sppos(j,:) = get(sp(j),'position');
            m_proj('Stereographic','lon',0,'lat',ylimit(2),'rad',abs(ylimit(2)-ylimit(1)));    
            axm = axesm('MapProjection','stereo','MapLatLimit',[45 90],'Grid','On','MLineLocation',60,...
                'PlineLocation',15,'MLabelLocation',60,'PLabelLocation',[60,75],'MeridianLabel','on','ParallelLabel','on','MLabelParallel',45,'LabelRotation','on','fontsize',fsize-2); 
            %m_grid('ytick',double(ytick) ,'xtick',6,'XaxisLocation','bottom','fontsize',fsize-2,'box','off');    
            framem; axis off; tightmap;
            if inputs.compareERA
                contourfm(obslat,obslontoplot,squeeze(obsfinal(:,:,j))',contourlevels,'LineStyle','none');
                colormap(cbrew);      
                caxis([clims(1)-cstep clims(2)+cstep]);        
                th = title(titles{j},'fontsize',fsize+2);            
                load coastlines
                geoshow(coastlat, coastlon, 'Color', 'black')
                if patchcoast
                    patchm(coastlat,coastlon,[.5 .5 .5])                                 
                else
                    linem(coastlat,coastlon,'Color','k')  
                end
                % plot significance
                hold on
                difference.observations.ttest (isnan(difference.observations.ttest)) = 0;
                f = figure;
                [H, ch] = contour(squeeze(difference.observations.ttest(j,:,:)),[0 1],'k');            
                close(f);
                figure(fig);
                contourm(obslat,obslon,squeeze(difference.observations.ttest(j,:,:))',[1 1],'k');
                oneloc = find(H(1,:) == 1);

                for pl = 1:length(oneloc)-1

    %                 tohatchcombinedlat = difference.observations.interplat(H(1,oneloc(pl)+1:oneloc(pl+1)-1));
    %                 tohatchcombinedlon = difference.observations.interplon(H(2,oneloc(pl)+1:oneloc(pl+1)-1));

                    tohatchcombinedlat = obslat(H(1,oneloc(pl)+1:oneloc(pl+1)-1));
                    tohatchcombinedlon = obslon(H(2,oneloc(pl)+1:oneloc(pl+1)-1));

                    m_hatch(tohatchcombinedlon,...
                        tohatchcombinedlat,'single',45,4,'color','k');                                
                end
            
                lathatch = repmat(obslat,[1,size(difference.observations.ttest(j,:,:),2)])';
                lonhatch = repmat(obslon,[1,size(difference.observations.ttest(j,:,:),3)]);
                %for k = size(difference.observations.ttest,3)                
                    latvec = lathatch(squeeze(difference.observations.ttest(j,:,:)) == 1);
                    lonvec = lonhatch(squeeze(difference.observations.ttest(j,:,:)) == 1);                                                
                %end      
                
                m_line(lonvec(1:2:end),latvec(1:2:end),'LineStyle','none','marker','.','color','k');
                                    
                hold on
                if j == 1
                    for k = 1:length(Seasplot.lon)
                        if strcmp(inputs.var,'ICEFRAC')
                            if k == 3
                                obslon2 = [obslon;obslon(1:20)+360];
                                lonstoplot = obslon2(obslon2 >= Seasplot.lon(k,1) & obslon2 <= Seasplot.lon(k,2));                                    
                            else
                                lonstoplot = obslon(obslon >= Seasplot.lon(k,1) & obslon <= Seasplot.lon(k,2));                
                                latstoplot = obslat(obslat >= Seasplot.lat(k,1) & obslat <= Seasplot.lat(k,2));
                            end
                        else
                            lonstoplot = obslon(obslon >= Seasplot.obs.lon(k,1) & obslon <= Seasplot.obs.lon(k,2));                
                            latstoplot = obslat(obslat >= Seasplot.obs.lat(k,1) & obslat <= Seasplot.obs.lat(k,2));
                        end                    

                        if strcmp(inputs.var,'ICEFRAC')
                            [x1,y1] = meshgrid(lonstoplot,latstoplot);
                            x1 = double(x1);
                            y1 = double(y1);
                            bound = boundary(x1(:),y1(:));

                            m_plot(x1(bound),y1(bound),...
                                'LineStyle','-','LineWidth',2,'color',cbrewqual2(k,:));%,'LineStyle',lstyle{k});   
                        else
    %                         m_plot(Seasplot.obs.lon(k),Seasplot.obs.lat(k),...
    %                             'LineStyle','-','Marker','p','MarkerFaceColor',...
    %                             cbrewqual2(k,:),'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',1);%,'LineStyle',lstyle{k});   
                            [x1,y1] = meshgrid(lonstoplot,latstoplot);
                            x1 = double(x1);
                            y1 = double(y1);
                            bound = boundary(x1(:),y1(:));                    
                            m_plot(x1(bound),y1(bound),...
                                'LineStyle','-','LineWidth',2,'color',cbrewqual2(k,:));%,'LineStyle',lstyle{k});   
                        end

                    end
                end
            end
            
            % plot model
            sp(j+slen) = subplot(2,3,j+slen);
            sppos(j+slen,:) = get(sp(j+slen),'position');                      
            
            m_proj('Stereographic','lon',0,'lat',ylimit(2),'rad',abs(ylimit(2)-ylimit(1)))    
            axm = axesm('MapProjection','stereo','MapLatLimit',[45 90],'Grid','On','MLineLocation',60,...
                'PlineLocation',15,'MLabelLocation',60,'PLabelLocation',[60,75],'MeridianLabel','on','ParallelLabel','on','MLabelParallel',45,'LabelRotation','on','fontsize',fsize-2); 
            %m_grid('ytick',double(ytick) ,'xtick',6,'XaxisLocation','bottom','fontsize',fsize-2);    
            framem; axis off; tightmap

            contourfm(latitude,modellontoplot,double(squeeze(modelfinal(:,:,j))'),contourlevels,'LineStyle','none');
            colormap(cbrew);      
            caxis([clims(1)-cstep clims(2)+cstep]);        
            th = title(titles{j+slen},'fontsize',fsize+2);            
            load coastlines
            geoshow(coastlat, coastlon, 'Color', 'black')
            if patchcoast
                patchm(coastlat,coastlon,[.5 .5 .5])                                 
            else
                linem(coastlat,coastlon,'Color','k')  
            end                        
            
%             m_proj('Stereographic','lon',0,'lat',ylimit(2),'rad',abs(ylimit(2)-ylimit(1)))                
% 
%             [~,h] = m_contourf(modellontoplot,latitude,squeeze(modelfinal(:,:,j))',contourlevels,'LineStyle','none');        
%             m_grid('ytick',double(ytick) ,'xtick',6,'XaxisLocation','bottom','fontsize',fsize);    
%             if patchcoast
%                 m_coast('patch',[.5 .5 .5],'edgecolor','none','LineWidth',1);
%             else                
%                 m_coast('color','k','LineWidth',1);
%             end
%             
%             colormap(cbrew);          
%             caxis([clims(1)-cstep clims(2)+cstep]);
%             th = title(titles{j+slen},'fontsize',fsize+2);            
                                              
            
            hold on
            
            % plot significance
            difference.difference.ttest (isnan(difference.difference.ttest)) = 0;
            f = figure;
            [H, ch] = contour(squeeze(difference.difference.ttest(j,:,:)),[0 1]);            
            close(f);
            figure(fig);
            
            oneloc = find(H(1,:) == 1);
         
            for pl = 1:length(oneloc)-1
                
                tohatchcombinedlat = latitude(H(1,oneloc(pl)+1:oneloc(pl+1)-1));
                tohatchcombinedlon = longitude(H(2,oneloc(pl)+1:oneloc(pl+1)-1));                
                if ~isempty(tohatchcombinedlat)
                    m_hatch(tohatchcombinedlon,...
                        tohatchcombinedlat,'single',45,4,'color','k');                                
                end                   
            end         
            
            lathatch = repmat(latitude,[1,size(difference.difference.ttest(j,:,:),2)])';
                lonhatch = repmat(longitude,[1,size(difference.difference.ttest(j,:,:),3)]);
                %for k = size(difference.observations.ttest,3)                
                    latvec = lathatch(squeeze(difference.difference.ttest(j,:,:)) == 1);
                    lonvec = lonhatch(squeeze(difference.difference.ttest(j,:,:)) == 1);                                                
                %end      
                
                %m_line(lonvec(1:2:end),latvec(1:2:end),'LineStyle','none','marker','.','color','k');
            
            
            if j == 1
                for k = 1:length(Seasplot.lon)
                    if strcmp(inputs.var,'ICEFRAC')
                        latstoplot = obslat(obslat >= Seasplot.lat(k,1) & obslat <= Seasplot.lat(k,2));
                        if k == 3 
                            obslon2 = [obslon;obslon(1:20)+360];
                            lonstoplot = obslon2(obslon2 >= Seasplot.lon(k,1) & obslon2 <= Seasplot.lon(k,2));                
                        else
                            lonstoplot = obslon(obslon >= Seasplot.lon(k,1) & obslon <= Seasplot.lon(k,2));                
                            
                        end
                        
                    else
                        lonstoplot = obslon(obslon >= Seasplot.obs.lon(k,1) & obslon <= Seasplot.obs.lon(k,2));                
                        latstoplot = obslat(obslat >= Seasplot.obs.lat(k,1) & obslat <= Seasplot.obs.lat(k,2));
                    end                                        

                    if strcmp(inputs.var,'ICEFRAC')
                        [x1,y1] = meshgrid(lonstoplot,latstoplot);
                        x1 = double(x1);
                        y1 = double(y1);
                        bound = boundary(x1(:),y1(:));                    
                        m_plot(x1(bound),y1(bound),...
                            'LineStyle','-','LineWidth',2,'color',cbrewqual2(k,:));%,'LineStyle',lstyle{k});   
                    else 
%                         m_plot(Seasplot.mod.lon(k),Seasplot.mod.lat(k),...
%                             'LineStyle','-','Marker','p','MarkerFaceColor',...
%                             cbrewqual2(k,:),'MarkerEdgeColor','k','MarkerSize',10,'LineWidth',1);%,'LineStyle',lstyle{k});   
                       [x1,y1] = meshgrid(lonstoplot,latstoplot);
                        x1 = double(x1);
                        y1 = double(y1);
                        bound = boundary(x1(:),y1(:));                    
                        m_plot(x1(bound),y1(bound),...
                            'LineStyle','-','LineWidth',2,'color',cbrewqual2(k,:));%,'LineStyle',lstyle{k});   
                    end

                end
            end
            
            
        end
        
        ch = colorbar;
        set(sp(end),'position',sppos(end,:));
        set(get(ch,'ylabel'),'string',ctitle,'fontsize',fsize+2)
        cbaxloc = get(ch,'Position');
        set(ch,'YTick',[clims(1):cstep*2:clims(2)],'fontsize',fsize)      
        set(ch,'YaxisLocation','right','fontsize',fsize);          
        
    end
%     ch = subplotmaps(toplotdifferencefinal,longtouse,latitude,{'div','RdBu'},1,[],fsize,title,[],[],ctitle,'on',...
%         clims,18,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic',patchcoast);
    
    % This is where I put on significance somehow
    %
    %
    %
    
    axpos = tightplots(gcf,10,25,.05,.15);
   
%     gh = get(gca,'title');
%     ghpos = get(gh,'position');
%     set(gh,'Position',[ghpos(1) ghpos(2)+ghpos(2)./50  ghpos(3)])     
    if inputs.nocoincidence
        coinext = 'nocoinc';
    elseif inputs.coincidence        
        coinext = 'coinc';
    else
        coinext = [];
    end
    
    if inputs.removeENSO
        ensoext = 're';
    else
        ensoext = [];
    end
    
    filename = [filedir,'Ensemble_',inputs.obstouse,'_',sprintf('%02d',mons(i,:)),'_',inputs.var,'_',...
        monthnames(inputs.tozmonth,1,'long'),'_Arcticozoneextremes_over',...
        num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',coinext,'_',ensoext];
    %print(filename,'-depsc');
    
    
    if combine        
        minx = min(axpos(:,1));
        maxx = max(axpos(:,1))+axpos(1,3);
        maxy = max(axpos(:,2))+axpos(1,4);
        annotation('textbox',[minx 1-(1-maxy)./2 maxx-minx 0],'String',...
            ['Change in ',titext,' due to ',monthnames(inputs.tozmonth,1,'long'),' ozone extremes'],...
            'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');     
    
        for g = 1:size(axpos,1)            
            annotation('textbox',[axpos(g,1),axpos(g,2),axpos(g,3:4)],'String',spind{g},...
                'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+2,... % 
                'EdgeColor','none','fontweight','bold');    
        end        
    
        print(filename,'-dpng');
        %print(filename,'-depsc');
        print(gcf,filename,'-depsc','-painters');
%         export_fig(filename,'-pdf');
        close 1
        break;
    else
        export_fig(filename,'-png');
        close 1
    end    
end
end