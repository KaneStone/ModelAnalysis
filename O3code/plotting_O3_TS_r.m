% plotting correlations
clear all
contour_plot_ensemble = 1;
contour_plot_individual = 0;
contour_plot_ensemble_hist = 0;
line_plot = 0;

constantyears = 1;
pastonly = 1;
poly = 1;
    
if pastonly
    pastext = 'past';
else
    pastext = [];
end
if constantyears
    file_ext = '_constantyears';
else
    file_ext = [];
end
hemisphere = 'south';
if strcmp(hemisphere,'north')
    hemisphere_ext = {'March-April','March'};
elseif strcmp(hemisphere,'south')
    hemisphere_ext = {'December-February','November'};
end
    

if strcmp(hemisphere,'south')            
    sites = [100,-5;75,-35;-130,-65;-110,42;-30,15;20,50;-160,66];            
elseif strcmp(hemisphere,'south')
    sites = [100,-5;75,-35;-130,-65;-110,42;-30,15;20,50;-160,66];            
end

ext = '200hPa';
load(['/Volumes/My Book for Mac/work/data/CESM-CCMI/O3/output/',hemisphere,'_',ext,'_O3_correlations',file_ext,'_',pastext]);
%load(['/Volumes/My Book for Mac/work/data/CESM-CCMI/O3/output/',hemisphere,'_',ext,'_O3_correlations_poly',file_ext]);
if poly
    r = rpoly;
end
names = fieldnames(TSdata_MAM);
names = names{1};

%runnames = {'REF-C2 no. 1 (1955-2010)','REF-C2 no. 2 (1955-2099)','REF-C2 no. 3 (1955-2099)',...
%    'SEN-C2-fGHG no. 1 (1960-2099)','SEN-C2-fGHG no. 2 (1960-2099)','SEN-C2-fGHG no. 3 (1960-2099)',...
%    'SEN-C2-fODS no. 1 (1960-2099)','SEN-C2-fODS no. 2 (1960-2099)','SEN-C2-fODS no. 3 (1960-2099)',...
%    'REF-C1 no. 1 (1955-2015)','REF-C1 no. 2 (1955-2015)','REF-C1 no. 3 (1955-2015)' ,...
%    'REF-C1 no. 4 (1955-2015)','REF-C1 no. 5 (1955-2015)','REF-C1SD no. 1 (1979-2015)',...
%    'REF-C2 ens (1955-2099)','SEN-C2-fGHG ens (1960-2099)','SEN-C2-fODS ens (1960-2099)',...
%    'REF-C1 ens (1955-2015)'};

runnames = {'REF-C2 no. 1','REF-C2 no. 2','REF-C2 no. 3',...
    'SEN-C2-fGHG no. 1','SEN-C2-fGHG no. 3','SEN-C2-fGHG no. 3',...
    'SEN-C2-fODS no. 1','SEN-C2-fODS no. 2','SEN-C2-fODS no. 3',...
    'REF-C1 no. 1','REF-C1 no. 2','REF-C1 no. 3','REF-C1 no. 4','REF-C1 no. 5',...
    'REF-C1SD no .1','REF-C2 ens','SEN-C2-fGHG ens','SEN-C2-fODS ens','REF-C1 ens'};

if pastonly
 runnames(15) = [];   
end

%% ensemble average contour plot
%ens_rtoplot = cat(3,squeeze(r(15).r(4,:,:)),squeeze(r(19).r(4,:,:)),squeeze(r(16).r(4,:,:)),...
%    squeeze(r(17).r(4,:,:)),squeeze(r(18).r(4,:,:)));
cbrew = flipud(cbrewer('div','RdBu',13));
fsize = 18;

if contour_plot_ensemble
    
    if pastonly
        noruns = 4;
        sx = 2;
        sy = 2;
        cbind = 4;
        createfig('medium','on')    
    elseif futureonly
        noruns = 3;
        sx = 1;
        sy = 3;
        cbind = 3;
        createfig('medium','on')    
    else noruns = 6;
        sx = 2;
        sy = 3;
        cbind = 5;
        createfig('large','on')    
    end
    count = 15;
    for i = 1:noruns
        if i < 6
        sp(i) = subplot(sy,sx,i); 
        sp_pos(i,:) = get(sp(i),'position');
        subpos(i,:) = get(sp(i),'outerposition');
        contourfm(latitude,[longitude;360],double([squeeze(r(count).r(:,:))',squeeze(r(count).r(1,:))']),14,'LineStyle','none','Color',[.7 .7 .7]); 
        colormap(cbrew);
        %contourfm(latitude,longitude,double(squeeze(r(count).r(4,:,:))')); 
        caxis([-.65 .65]);
        xlim([-181 181]);
        ylim([-91 91]);
        if i == cbind
            ch = colorbar('eastoutside');
            cbaxloc = get(ch,'position');            
            set(ch,'ytick',-.65:.1:.65);
            set(get(ch,'ylabel'),'string','r','fontsize',fsize+2)
        end
        set(gca,'fontsize',fsize-4,'box','on','xtick',-180:40:180,'xticklabel',-180:40:180,...
            'ytick',-90:20:90,'yticklabel',-90:20:90);
        load coastlines
        geoshow(coastlat, coastlon, 'Color', 'black','LineWidth',2)
        if i == cbind || i == cbind-1
            xlabel('Longitude','fontsize',fsize);
        end
        if mod(i,2)
            ylabel('Latitude','fontsize',fsize);
        else            
        end
        
        if i == 1 || i == 2
            set(sp(i),'position',[sp_pos(i,1)-.05,sp_pos(i,2)-.02,sp_pos(i,3)+.025,sp_pos(i,4)+.025]);
            
        elseif i == 3 || i == 4
            set(sp(i),'position',[sp_pos(i,1)-.05,sp_pos(i,2)+.02,sp_pos(i,3)+.025,sp_pos(i,4)+.025]);
        elseif i == 5 || i == 6
            set(sp(i),'position',[sp_pos(i,1)-.05,sp_pos(i,2)+.06,sp_pos(i,3)+.025,sp_pos(i,4)+.025]);
        end
        
        sp_pos(i,:) = get(sp(i),'position');
        subpos(i,:) = get(sp(i),'outerposition');
        title(runnames{count},'fontsize',fsize);               
        count = count+1;
            cbind = noruns;
        else
            
            %sp(i) = subplot(3,2,i);            
        end  
        hold on        
        h = plot(sites(:,1),sites(:,2),'Marker','pentagram','linestyle','none',...
                'MarkerEdgeColor',[77,175,74]./255,'MarkerFaceColor',[77,175,74]./255,'MarkerSize',15);                
        set(gca,'LineWidth',1,'box','on');
    end

%moveplots(sp,.08,.08);
set(ch,'Position',[sp_pos(2,1)+sp_pos(2,3)+.02, sp_pos(cbind,2),...
    cbaxloc(3), sp_pos(2,2)+sp_pos(2,4)-sp_pos(cbind,2)],...
    'Box','on','YAxisLocation','right','fontsize',fsize);
if strcmp(ext,'200hPa')    
        mtit([hemisphere_ext{1},' 200 hPa temperature - ',hemisphere_ext{2},' 70 hPa ozone correlations'],'xoff',-.01,'yoff',.02,'fontsize',fsize+4);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/200hPaT_70hPaO3_r',file_ext,'_',pastext,'_',hemisphere,'.pdf'];                
else
    mtit([hemisphere_ext{1},' surface temperature - ',hemisphere_ext{2},' 70 hPa ozone correlations'],'xoff',-.01,'yoff',.02,'fontsize',fsize+4);    
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/TS_70hPaO3_r',file_ext,'_',hemisphere,'.pdf'];                
end
export_fig(filename,'-pdf');
end
%%
if contour_plot_ensemble_hist
    createfig('large')
    count = 10;
    for i = 1:6     
        sp(i) = subplot(3,2,i); 
        sp_pos(i,:) = get(sp(i),'position');
        subpos(i,:) = get(sp(i),'outerposition');
        contourfm(latitude,[longitude;360],double([squeeze(r(count).r(4,:,:))',squeeze(r(count).r(4,1,:))]),14,'LineStyle','none','Color',[.7 .7 .7]); 
        colormap(cbrew);
        %contourfm(latitude,longitude,double(squeeze(r(count).r(4,:,:))')); 
        caxis([-.65 .65]);
        xlim([-181 181]);
        ylim([-91 91]);
        if i == 5
            ch = colorbar('eastoutside');
            cbaxloc = get(ch,'position');            
            set(ch,'ytick',-.65:.1:.65);
        end
        set(gca,'fontsize',fsize-4,'box','on','xtick',-180:40:180,'xticklabel',-180:40:180,...
            'ytick',-90:20:90,'yticklabel',-90:20:90);
        load coastlines
        geoshow(coastlat, coastlon, 'Color', 'black','LineWidth',2)
        if i == 5 || i == 6
            xlabel('Longitude','fontsize',fsize);
        end
        if mod(i,2)
            ylabel('Latitude','fontsize',fsize);
        else            
        end
        
        if i == 1 || i == 2
            set(sp(i),'position',[sp_pos(i,1)-.05,sp_pos(i,2)-.02,sp_pos(i,3)+.025,sp_pos(i,4)+.025]);
            
        elseif i == 3 || i == 4
            set(sp(i),'position',[sp_pos(i,1)-.05,sp_pos(i,2)+.02,sp_pos(i,3)+.025,sp_pos(i,4)+.025]);
        elseif i == 5 || i == 6
            set(sp(i),'position',[sp_pos(i,1)-.05,sp_pos(i,2)+.06,sp_pos(i,3)+.025,sp_pos(i,4)+.025]);
        end
        
        sp_pos(i,:) = get(sp(i),'position');
        subpos(i,:) = get(sp(i),'outerposition');
        title(runnames{count},'fontsize',fsize);               
        count = count+1;
                      
            %sp(i) = subplot(3,2,i);                    
        set(gca,'LineWidth',1,'box','on');
    end

%moveplots(sp,.08,.08);
set(ch,'Position',[sp_pos(2,1)+sp_pos(2,3)+.02, sp_pos(5,2),...
    cbaxloc(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(5,2)-.03],...
    'Box','on','YAxisLocation','right','fontsize',fsize);
if strcmp(ext,'200hPa')    
        mtit([hemisphere_ext{1},' 200 hPa temperature - ',hemisphere_ext{2},' 70 hPa ozone correlations'],'xoff',-.01,'yoff',.02,'fontsize',fsize+4);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/HistInd_200hPaT_70hPaO3_r',file_ext,'_',hemisphere,'.pdf'];                
else
    mtit([hemisphere_ext{1},' surface temperature - ',hemisphere_ext{2},' 70 hPa ozone correlations'],'xoff',-.01,'yoff',.02,'fontsize',fsize+4);    
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/HistInd_TS_70hPaO3_r',file_ext,'_',hemisphere,'.pdf'];                
end
export_fig(filename,'-pdf');
end

%% 
if contour_plot_individual
    createfig('large')
    count = 1;
    for i = 1:9        
        sp(i) = subplot(3,3,i); 
        sp_pos(i,:) = get(sp(i),'position');
        subpos(i,:) = get(sp(i),'outerposition');
        contourfm(latitude,[longitude;360],double([squeeze(r(count).r(4,:,:))',squeeze(r(count).r(4,1,:))]),14,'LineStyle','none','Color',[.7 .7 .7]); 
        colormap(cbrew);
        %contourfm(latitude,longitude,double(squeeze(r(count).r(4,:,:))')); 
        caxis([-.65 .65]);
        xlim([-181 181]);
        ylim([-91 91]);
        if i == 9
            ch = colorbar('eastoutside');
            cbaxloc = get(ch,'position');            
            set(ch,'ytick',-.65:.1:.65);
        end
        set(gca,'fontsize',fsize-4,'box','on','xtick',-180:40:180,'xticklabel',-180:40:180,...
            'ytick',-90:20:90,'yticklabel',-90:20:90);
        load coastlines
        geoshow(coastlat, coastlon, 'Color', 'black','LineWidth',2)
        if i == 7 || i == 8 || i == 9
            xlabel('Longitude','fontsize',fsize);
        end
        if i == 1 || i == 4 || i == 7
            ylabel('Latitude','fontsize',fsize);
        else
        
            
        %end
        end        
        if i == 1 || i == 2 || i == 3
            set(sp(i),'position',[sp_pos(i,1)-.04,sp_pos(i,2)-.04,sp_pos(i,3)+.025,sp_pos(i,4)+.025]);
            
        elseif i == 4 || i == 5 || i == 6
            set(sp(i),'position',[sp_pos(i,1)-.04,sp_pos(i,2)+.04,sp_pos(i,3)+.025,sp_pos(i,4)+.025]);
        elseif i == 7 || i == 8 || i == 9
            set(sp(i),'position',[sp_pos(i,1)-.04,sp_pos(i,2)+.12,sp_pos(i,3)+.025,sp_pos(i,4)+.025]);
        end
        sp_pos(i,:) = get(sp(i),'position');
        subpos(i,:) = get(sp(i),'outerposition');
        title(runnames{count},'fontsize',fsize);               
        count = count+1;                      
        set(gca,'LineWidth',1,'box','on');    
    end
    %moveplots(sp,.08,.2);
    set(ch,'Position',[sp_pos(3,1)+sp_pos(3,3)+.02, sp_pos(9,2),...
        cbaxloc(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(9,2)-.015],...
        'Box','on','YAxisLocation','right','fontsize',fsize);
    if strcmp(ext,'200hPa')    
        mtit([hemisphere_ext{1},' 200 hPa temperature - ',hemisphere_ext{2},' 70 hPa ozone correlations'],'xoff',-.01,'yoff',-.02,'fontsize',fsize+4);
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/Ind_200hPaT_70hPaO3_r',file_ext,'_',hemisphere,'.pdf'];                
    else
        mtit([hemisphere_ext{1},' surface temperature - ',hemisphere_ext{2},' 70 hPa ozone correlations'],'xoff',-.01,'yoff',-.02,'fontsize',fsize+4);    
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/Ind_TS_70hPaO3_r',file_ext,'_',hemisphere,'.pdf'];                
    end
    export_fig(filename,'-pdf');
end


%% testing line plot
if line_plot
    constantyears = 0;
    corryears = [1979 2099];
    if strcmp(hemisphere,'north')
        offset = 0;
    else
        offset = 12;
    end
    %lon = 165;
    %lat = 52;
    detr = 1;
    for j = 16;
    fig = figure;
    set(fig,'color','white','position',[50 50 1500 900]);
    iadd = 14;
    mon = 11;
    for i = 1:size(sites,1)
        [~,lonindex] = min(abs(longitude-sites(i,1)-180));
        [~,latindex] = min(abs(latitude-sites(i,2)));
        subplot(ceil(size(sites,1)./2),2,i)
        bp(j) = min(find(years(j).y == 1995)-1)./12;
        if constantyears
            yearstart(i) = max(find(years(iadd+i).y == corryears(1)))./12;
            yearend(i) = max(find(years(iadd+i).y == corryears(2)))./12;   
            bp(i) = bp(i) - yearstart(i);    
        end
        
        yyaxis left
        if constantyears
            plot(years(iadd+i).y(yearstart(i)*12:12:yearend(i)*12-offset),detrend(squeeze(TSdata_MAM(14+i).(names)(lonindex,latindex,yearstart(i):yearend(i)))')+...
                nanmean(squeeze(TSdata_MAM(iadd+i).(names)(lonindex,latindex,yearstart(i):yearend(i)))'));
            set(gca,'Ydir','reverse');
            yyaxis right    
            plot(years(iadd+i).y(yearstart(i)*12:12:yearend(i)*12-offset),...
                detrend(squeeze(O3weighted(14+i).wa(1,yearstart(i)*12-12+mon:12:yearend(i)*12-offset))','linear',bp(i)));
        else
            %plot(years(iadd+i).y(1:12:end-offset),detrend(squeeze(TSdata_MAM(14+i).(names)(lonindex,latindex,:))')+...
            %    nanmean(squeeze(TSdata_MAM(iadd+i).(names)(lonindex,latindex,:))'));
            if detr
               p1 = polyfit(years(j).y(1+offset:12:bp(j)*12),squeeze(TSdata_MAM(j).(names)(lonindex,latindex,1:bp(j)-offset/12)),3);
                x1 = years(j).y(1+offset:12:bp(j)*12);
                y1 = polyval(p1,x1);                                                

                p2 = polyfit(years(j).y(bp(j)*12+12:12:end),squeeze(TSdata_MAM(j).(names)(lonindex,latindex,bp(j):end)),2);                        
                x2 = years(j).y(bp(j)*12+12:12:end);
                y2 = polyval(p2,x2);
                Temp_detrend(j).a(i,:) = [squeeze(TSdata_MAM(j).(names)(lonindex,latindex,1:bp(j)-offset/12))-y1;squeeze(TSdata_MAM(j).(names)(lonindex,latindex,bp(j):end))-y2];
                    
                %ozone
                po1 = polyfit(years(j).y(1:12:bp(j)*12-12),squeeze(O3weighted(j).wa(4,mon:12:bp(j)*12-12))',3);
                xo1 = years(j).y(1:12:bp(j)*12-12);
                yo1 = polyval(po1,xo1);                                                

                po2 = polyfit(years(j).y(bp(j)*12:12:end-12),squeeze(O3weighted(j).wa(4,bp(j)*12-12+mon:12:end-12))',2);                        
                xo2 = years(j).y(bp(j)*12:12:end-12);
                yo2 = polyval(po2,xo2);                          
                O3_detrend = [squeeze(O3weighted(j).wa(4,mon:12:bp(j)*12-12))'-yo1;squeeze(O3weighted(j).wa(4,bp(j)*12-12+mon:12:end-12))'-yo2];
                
                plot(years(j).y(1:12:end-offset),Temp_detrend(j).a(i,:))
                ylim([-4 4]);
                hold on
                yyaxis right
                plot(years(j).y(1:12:end-offset),O3_detrend)
                ylim([-6e-7 6e-7]);
                
                
                else
                  plot(years(iadd+j).y(1:12:end-offset),squeeze(TSdata_MAM(14+j).(names)(lonindex,latindex,:))')                
                set(gca,'Ydir','reverse')
                
                yyaxis right    
                plot(years(iadd+j).y(1:12:end-offset),squeeze(O3weighted(14+j).wa(4,mon:12:end-offset))');
                
                
                
            end
        end

    %     figure
    %     yyaxis left
    %     plot(years(14+i).y(1:12:end-offset),squeeze(TSdata_MAM(14+i).(names)(lonindex,latindex,:))');
    %     set(gca,'Ydir','reverse');
    %     yyaxis right    
    %     plot(years(14+i).y(1:12:end-offset),squeeze(O3weighted(14+i).wa(1,3:12:end-offset))');

    end
    end
end

%%
mon = 11;
plot(years(1).y(1:12:end)',squeeze(O3weighted(1).wa(4,mon:12:end))');
hold on
%plot(years(1).y(469+mon:12:end)',squeeze(O3weighted(1).wa(4,470+mon:12:end))')
%% 
hold on
abc = squeeze(O3weighted(1).wa(1,mon:12:469))';
p = polyfit(years(1).y(1:12:468)',abc',3);
x1 = years(1).y(1:12:468)';
y1 = polyval(p,x1);
plot(years(1).y(1:12:468)',abc);
hold on
plot(x1,y1);



abc2 = squeeze(O3weighted(1).wa(1,468+mon:12:end))';
p = polyfit(years(1).y(468+mon:12:end)',abc2',2);
x1 = years(1).y(468+mon:12:end)';
y2 = polyval(p,x1);
plot(years(1).y(468+mon:12:end)',abc2);
hold on
plot(x1,y2);

%%
final = [squeeze(O3weighted(1).wa(4,mon:12:469))'-y1';squeeze(O3weighted(1).wa(4,468+mon:12:end))'-y2'];
figure;
plot(final);

%%
