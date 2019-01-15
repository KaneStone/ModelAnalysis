function [] = predruns_plotdiffmonths(differences,lons,lats,moninlength)
%%
    
    complats = [66,70,52,35];
    complons = [310,60,255,60];
    
    [~,lats_ind] = min(abs(lats - complats));
    [~,lons_ind] = min(abs(lons - complons));
    
    for i = 1:moninlength
        for j = 1:size(lons_ind,2)
            diffarea = squeeze(differences(i,lons_ind(j),...
                lats_ind(j)));
            meandiff(i,j) = nanmean(diffarea(:));
        end
    end

    %%    
    cbrewqual = cbrewer('qual','Set1',10);
    cbrewqual2 = cbrewqual([3,4,10,8],:);
    lnames = {['Greenland, (',num2str(abs(complons(1))),'{\circ}','W, ',num2str(abs(complats(1))),'{\circ}','N)'],...
            ['Russia, (',num2str(abs(complons(2))),'{\circ}','W, ',num2str(abs(complats(2))),'{\circ}','N)'],...
            ['North America, (',num2str(abs(complons(3))),'{\circ}','W, ',num2str(abs(complats(3))),'{\circ}','N)'],...
            ['Himalayas, (',num2str(abs(complons(4))),'{\circ}','W, ',num2str(abs(complats(4))),'{\circ}','N)']};
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/compositedifferences/diffAreaLines/North_',num2str(20),'th_percentile_lines'];
    xticklab = {'Mar','Apr','May','Jun'};    
    
    createfig('medium','on')
    lstyle = {':','-','-.','--'};        
    for i = 1:size(meandiff,2)
        plot(meandiff(:,i),'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
        hold on
    end
    hold on
    plot([0 5],[0,0],'--k','lineWidth',2);
    xlim([.9,4.1]);
    set(gca,'fontsize',20,'xtick',1:5,'xticklabel',xticklab)
    xlabel('Month','fontsize',22);
    ylabel('Temperature difference (K)','fontsize',22);
    title('Model ensemble composite surface temperature difference','fontsize',24)
    lh = legend(lnames);
    set(lh,'box','off','fontsize',24,'location','SouthEast')
    
    export_fig(filename,'-pdf');
end
