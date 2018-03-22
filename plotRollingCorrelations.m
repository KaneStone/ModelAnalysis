function [] = plotRollingCorrelations(rollcorr,ozonetoplot,loc,numsubplots,subplotorientation,...
    titles,fsize,rollwindow,rollyears,latitudes,longitudes,modelext)

%% plotting rolling correlations
cbrewqual = cbrewer('qual','Set1',10);
ozonesmooth = smooth(ozonetoplot,rollwindow,'moving');

createfig('largeportrait','on')

%-----------midlat-positive----------------------
% latstoplot = [-26,-30,-22];
% lonstoplot = [146,295,25];


for i = 1:numsubplots
    subplot(subplotorientation(1),subplotorientation(2),i);
    
    set(gca,'ycolor','k')
    if mod(rollyears,2)
        plot(rollyears,ozonesmooth(ceil(rollwindow/2):end-ceil(rollwindow/2)),'LineWidth',3,'color','k');
    else
        plot(rollyears,ozonesmooth(rollwindow/2+1:end-rollwindow/2),'LineWidth',3,'color','k');
    end
    forleg{1} = 'TCO';    
    set(gca,'fontsize',fsize);
    ylabel('DU','fontsize',fsize+2);       
    if i == 1 || i == 3
        set(gca,'ycolor','k')
    else        
        set(gca,'ycolor','k','ydir','reverse')
    end
    
    yyaxis right
    set(gca,'YDir','reverse','ycolor','k')

    for j = 1:length(loc(i).lats)
        [~,rolllatind] = min(abs(latitudes - loc(i).lats(j))); 
        [~,rolllonind] = min(abs(longitudes - loc(i).lons(j)));

%         rolllatind = latitudes >= loc(i).lats(j,1) & latitudes <= loc(i).lats(j,2); 
%         rolllonind = longitudes >= loc(i).lons(j,1) & longitudes <= loc(i).lons(j,2); 
%         toplot = permute(rollcorr(rolllonind,rolllatind,:),[3,1,2]);
%         toplot = nanmean(toplot(:,:),2);
        toplot = squeeze(rollcorr(rolllonind,rolllatind,:));
        ph(j) = plot(rollyears,toplot,'LineWidth',2,'color',cbrewqual(j,:));
        hold on
        forleg{j+1} = [num2str(loc(i).lats(j)),'N, ',num2str(loc(i).lons(j)),'E'];
    end
    ylabel('Rolling correlation','fontsize',fsize+2);
    xlabel('Year','fontsize',fsize+2);    
    title(titles{i},'fontsize',fsize+2);
    lh = legend(forleg);    
    if i == 1 || i == 3
        set(lh,'box','off','location','SouthEast','fontsize',fsize);
    else
        set(lh,'box','off','location','SouthWest','fontsize',fsize);
    end
        
    clearvars forleg
end


if loc(1).lats(1) > 0
    hemext = 'North';
else
    hemext = 'South';
end

annotation('textbox',[0 .905 1 .1],'String',[hemext,'ern Hemisphere - selected rolling correlations ', num2str(rollwindow),' year window, ', num2str(rollyears(1)-rollwindow/2-1),'-',num2str(rollyears(end)+rollwindow/2)],...
    'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,'EdgeColor','none','fontweight','bold')

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/correlations/Rolling/',hemext,'_',modelext,'_Rolling_',num2str(rollwindow),'year_',num2str(rollyears(1)),'-',num2str(rollyears(end))];
export_fig(filename,'-pdf');

% ------------midlat-negative---------------------
% latstoplotneg = [-30,-20];
% lonstoplotneg = [210,170];
% subplot(2,2,2);
% for i = 1:length(latstoplotneg)
%     [~,rolllatind] = min(abs(ERAdata.latitude - latstoplotneg(i))); 
%     [~,rolllonind] = min(abs(ERAdata.longitude - lonstoplotneg(i)));
%             
%     ph(i) = plot(rollyears,squeeze(rollcorr.r(rolllonind,rolllatind,:)),'LineWidth',2,'color',cbrewqual(i,:));
%     hold on
%     forleg{i} = [num2str(latstoplotneg(i)),'N, ',num2str(lonstoplotneg(i)),'E'];
% end
% yyaxis right
% set(gca,'YDir','reverse')
% plot(rollyears,ozonesmooth(ceil(rollwindow/2):end-ceil(rollwindow/2)),'LineWidth',3,'color','k');
% forleg{i+1} = 'Normalized ozone';
% lh = legend(forleg);    
% set(lh,'box','off','location','NorthEast','fontsize',fsize)
% 
% ylabel('Normalized TOZ','fontsize',fsize+2);
% title('Negative correlations - mid-latitudes','fontsize',fsize+2);
% ------------polar-positive---------------------
% latstoplot = [-60,-55,-60];
% lonstoplot = [60,140,240];
% 
% subplot(2,2,3);
% for i = 1:length(latstoplot)
%     [~,rolllatind] = min(abs(ERAdata.latitude - latstoplot(i))); 
%     [~,rolllonind] = min(abs(ERAdata.longitude - lonstoplot(i)));
%         
%     yyaxis left
%     ph(i) = plot(rollyears,squeeze(rollcorr.r(rolllonind,rolllatind,:)),'LineWidth',2,'color',cbrewqual(i,:));
%     hold on
%     forleg{i} = [num2str(latstoplot(i)),'N, ',num2str(lonstoplot(i)),'E'];
% end
% ylabel('Rolling correlation','fontsize',fsize+2);
% yyaxis right
% set(gca,'YDir','reverse')
% plot(rollyears,ozonesmooth(ceil(rollwindow/2):end-ceil(rollwindow/2)),'LineWidth',3,'color','k');
% forleg{i+1} = 'Normalized ozone';
% lh = legend(forleg);    
% set(lh,'box','off','location','SouthEast','fontsize',fsize)
% xlabel('Year','fontsize',fsize+2);
% title('Positive correlations - polar-latitudes','fontsize',fsize+2);
% 
% ------------polar-negative---------------------
% latstoplotneg = [-50,-47];
% lonstoplotneg = [320,57];
% subplot(2,2,4);
% for i = 1:length(latstoplotneg)
%     [~,rolllatind] = min(abs(ERAdata.latitude - latstoplotneg(i))); 
%     [~,rolllonind] = min(abs(ERAdata.longitude - lonstoplotneg(i)));
%         
%     yyaxis left
%     ph(i) = plot(rollyears,squeeze(rollcorr.r(rolllonind,rolllatind,:)),'LineWidth',2,'color',cbrewqual(i,:));
%     hold on
%     forleg{i} = [num2str(latstoplotneg(i)),'N, ',num2str(lonstoplotneg(i)),'E'];
% end
% yyaxis right
% plot(rollyears,ozonesmooth(ceil(rollwindow/2):end-ceil(rollwindow/2)),'LineWidth',3,'color','k');
% forleg{i+1} = 'Normalized ozone';
% lh = legend(forleg);    
% set(lh,'box','off','location','NorthEast','fontsize',fsize)
% 
% ylabel('Normalized TOZ','fontsize',fsize+2);
% xlabel('Year','fontsize',fsize+2);
% title('Negative correlations - polar-latitudes','fontsize',fsize+2);




end