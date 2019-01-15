function [] = predruns_plotSkillLines(GSS,waclon,waclat);

%% plot line plots
fsize = 16;
% obslats = [70,67,50,30];%Greenland,Russia,USA,India
% obslons = [310,70,260,78];
complats = [68,67,52,36];
complons = [310,70,255,75];
% for i = 1:length(obslats)
%     [~,latind(i)] = min(abs(obslats(i)-ERAdata.latitude));
%     [~,lonind(i)] = min(abs(obslons(i)-ERAdata.longitude));
% end
for i = 1:length(complats)
    [~,modlatind(i)] = min(abs(complats(i)-waclat));
    [~,modlonind(i)] = min(abs(complons(i)-waclon));
end

lnames = {['Greenland, (',num2str(abs(obslons(1))),'{\circ}','W, ',num2str(abs(obslats(1))),'{\circ}','N)'],...
    ['Russia, (',num2str(abs(obslons(2))),'{\circ}','W, ',num2str(abs(obslats(2))),'{\circ}','N)'],...
    ['North America, (',num2str(abs(obslons(3))),'{\circ}','W, ',num2str(abs(obslats(3))),'{\circ}','N)'],...
    ['Himalayas, (',num2str(abs(obslons(4))),'{\circ}','W, ',num2str(abs(obslats(4))),'{\circ}','N)']};
lnamesmod = {['Greenland, (',num2str(abs(complons(1))),'{\circ}','W, ',num2str(abs(complats(1))),'{\circ}','N)'],...
    ['Russia, (',num2str(abs(complons(2))),'{\circ}','W, ',num2str(abs(complats(2))),'{\circ}','N)'],...
    ['North America, (',num2str(abs(complons(3))),'{\circ}','W, ',num2str(abs(complats(3))),'{\circ}','N)'],...
    ['Himalayas, (',num2str(abs(complons(4))),'{\circ}','W, ',num2str(abs(complats(4))),'{\circ}','N)']};

xticklab = {'March','April','May','June','July'};

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([3,4,10,8],:);

fig = figure;
set(fig,'position',[100 100 1200 400],'color','white');
lstyle = {':','-','-.','--'};        
sp(1) = subplot(1,2,1);
sp_pos(1,:) = get(sp(1),'position');
for i = 1:length(obslats)    
    pho(i) = plot(differences_ind(:,lonind(i),latind(i)),'o','linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
    hold on
end
plot([0 5],[0,0],'--k','lineWidth',2);
xlim([.5,5.5]);
ylim([-7 7]);
set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
xlabel('Month','fontsize',fsize+2);
ylabel('Temperature difference (K)','fontsize',fsize+2);
title('Observed temperature difference','fontsize',fsize+4)
lh = legend(pho,lnames);
set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
sp(2) = subplot(1,2,2);
totimes = [-1.5,-.5,.5,1.5];
for i = 1:length(complats)      
    %plot(squeeze(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:))','linewidth',1,'LineStyle',lstyle{i},'color',cbrewqual2(i,:)./1.5);
    ph(i) = plot([1:5]+(totimes(i)./20),squeeze(nanmean(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:),1)),'o',...
        'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
    hold on
    errorbar([1:5]+(totimes(i)./20),squeeze(nanmean(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:),1)),squeeze(std(modeldifferences.indmonths.individual(:,modlonind(i),modlatind(i),:),0,1)),...
        'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
end
sp_pos(2,:) = get(sp(2),'position');
set(sp(2),'position',[sp_pos(2,1)-.05,sp_pos(2,2),sp_pos(2,3),sp_pos(1,4)]);
set(sp(1),'position',[sp_pos(1,1),sp_pos(1,2),sp_pos(1,3),sp_pos(1,4)]);
%plot(meandiff,'linewidth',3);
hold on
plot([0 5],[0,0],'--k','lineWidth',2);
xlim([.5,5.5]);
ylim([-7 7]);
set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
xlabel('Month','fontsize',fsize+2);

title('Ensemble mean temperature difference','fontsize',fsize+4)
lh = legend(ph,lnamesmod);
set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/North_20th_percentile_lines'];
export_fig(filename,'-pdf');
end
