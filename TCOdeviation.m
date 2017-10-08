% TCO deviation

TCOdirectory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/TCO/'];

%addpath(genpath('/home/stonek/code'))

TCOfiles = dir([TCOdirectory,'TOZ*']);

yearmin = 20040201;
yearmax = 20160201;
remove2002 = 'Southern'; % 'all_lats' or 'Southern' or 'none'

TCOname = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};

halon = zeros(6,144);
halat = zeros(6,96);
halon (halon == 0) = NaN;
halat (halat == 0) = NaN;

lats = [[30 60];[60 90]];
%lats = [[-90 -60];[-60 -30];[-30 0];[0 30];[30 60];[60 90]];

szlat = size(lats);
TCOzonallat = [];

yearstoexclude = 2015;
deviation_year = 2015;
vert = 0;
allyears = 2014:2015;
for i = 1:length(TCOfiles)
    
    [info.(TCOname{i}), data.(TCOname{i}), attributes.(TCOname{i})] = ...
        Read_in_netcdf([TCOdirectory,TCOfiles(i).name]);
    
    %time extract
    TCO.(TCOname{i}) = data.(TCOname{i}).toz(:,:,data.(TCOname{i}).date >= yearmin & data.(TCOname{i}).date <= yearmax);
    date.(TCOname{i}) = data.(TCOname{i}).date(data.(TCOname{i}).date >= yearmin & data.(TCOname{i}).date <= yearmax);
    %calculate deviation
    TCO.(TCOname{i}) = cat(3,TCO.(TCOname{i}),zeros(size(TCO.(TCOname{i}),1),size(TCO.(TCOname{i}),2),1));
    TCO.(TCOname{i}) (TCO.(TCOname{i}) == 0) = NaN;
    
    yearindex = find(date.(TCOname{i}) == str2double([num2str(deviation_year),'0201']));
    
    for j = 1:12
        Mean(:,:,j) = squeeze(nanmean(TCO.(TCOname{i})(:,:,j:12:end-12),3));
    end

    deviation(i,:,:,:) = squeeze(TCO.(TCOname{i})(:,:,yearindex:yearindex+11))-Mean;
%     deviation2 = squeeze(data(:,deviationindex+1,:))-Mean2;
%     deviation_percent = (squeeze(data(:,deviationindex,:))-Mean1)./Mean1*100;
%     deviation2_percent = (squeeze(data(:,deviationindex+1,:))-Mean2)./Mean2*100;
%     deviation_final = [deviation(1:12,:); deviation2(1:12,:)];
%     deviation_final_percent = [deviation_percent(1:12,:); deviation2_percent(1:12,:)];
        
end

%% Reading in SAD_SULFC

SADdirectory = ['/Volumes/My Book for Mac/work/data/WACCM/SAD_SULFC/'];

%addpath(genpath('/home/stonek/code'))

SADfiles = dir([SADdirectory,'SAD*']);

SADname = {'MAM','VCMAM'};

for i = 1:length(SADfiles)
    
    [SADinfo.(SADname{i}), SADdata.(SADname{i}), SADattributes.(SADname{i})] = ...
        Read_in_netcdf([SADdirectory,SADfiles(i).name]);
    
end

SADlevel = sum(SADdata.MAM.SAD_SULFC(:,:,SADdata.MAM.lev <= 200,:),3);

%% plotting
fsize = 16;
createfig('large');

cmap = flipud(cbrewer('div','RdBu',18));
cmapfinal = cmap(1:10,:);
intervals = -170:20:30;

SADcmap = cbrewer('seq','YlOrBr',12);
SADintervals = 3e-7:7e-8:1.0e-6;

daytitle = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};
montoplot = 9;
montoplot = montoplot-1;
count = 1;
for i = 1:6
    sp = subplot(3,2,i);
    sp_pos(i,:) = get(sp,'position');
    if i == 1 || i == 3 || i == 5 
        sep = contourf(data.CCMI.lat,data.CCMI.lon,squeeze(deviation(2,:,:,count+montoplot)),intervals);
        ylabel('Longitude','fontsize',fsize+2);
        set(sp,'position',[sp_pos(i,1)+.025,sp_pos(i,2:4)])                
        
    else sep = contourf(data.CCMI.lat,data.CCMI.lon,squeeze(deviation(3,:,:,count+montoplot)),intervals);
        set(sp,'position',[sp_pos(i,1)-.025,sp_pos(i,2:4)])
    end
    colormap(cmapfinal);           

    set(gca,'fontsize',fsize-2,'xtick',-90:10:90);
    sp_pos(i,:) = get(sp,'position');
    caxis([min(intervals),max(intervals)]);
    if i == 5 || i == 6
        xlabel('Latitude','fontsize',fsize+2);
    end
    
    %moving suplots
    
    
    if i == 6;
        ch = colorbar;       
        set(ch,'Ticks',intervals)
        set(get(ch,'ylabel'),'string','Dobson Units','fontsize',fsize+2)
        ch_pos = get(ch,'position');
            set(ch,'Position',[sp_pos(i,1)+sp_pos(i,3)+sp_pos(i,3)/15, ch_pos(2),...
            ch_pos(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(6,2)],...
            'Box','on','YAxisLocation','right','fontsize',fsize);
    end 
    
    if i == 1 || i == 3 || i == 5
        haxes1 = gca;
        haxes1_pos(i,:) = get(haxes1,'Position');
        haxes2 = axes('Position',haxes1_pos(i,:),...
          'XAxisLocation','bottom',...
          'YAxisLocation','left',...
          'Color','none');
        hold on        
        sep = contour(data.CCMI.lat,data.CCMI.lon,squeeze(SADlevel(:,:,count+montoplot)),SADintervals,'LineWidth',1.5);        
        colormap(haxes2,SADcmap);           
        set(gca,'fontsize',fsize-2);
        sp_pos(i,:) = get(sp,'position');
        caxis([min(SADintervals),max(SADintervals)]);
        
    end
    
    if i == 2 || i == 4 || i == 6
        title(['VC-MAM ',daytitle{count+montoplot},' ',num2str(deviation_year),' deviation'],'fontsize',fsize+4);
        count = count+1;
    else title(['MAM ',daytitle{count+montoplot},' 2015 deviation'],'fontsize',fsize+4);
    end
end
mtit([num2str(deviation_year),' anomaly from 2004-2014 mean'],'xoff',-.01,'yoff',.04,'fontsize',fsize+8);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
    num2str(deviation_year),'_WACCM_TCO_anomaly_','_',daytitle{montoplot+1},'-',daytitle{montoplot+3},'.png'];
    export_fig(filename,'-png')
