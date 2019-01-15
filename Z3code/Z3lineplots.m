function [] = Z3lineplots(month,pres,data,data_b,data_b2,years,ERApressure,trendyears,...
    SDsave,SDsaveyears,sim,linetitle,area,xlimits)

monthnames = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};

fsize = 18;
[~,presind] = min(abs(double(ERApressure) - pres));
colors2 = cbrewer('qual','Set1',5);
fig = figure;
set(fig,'color','white','position',[100 100 1000 700]);
hold on        
for i = 1:length(data)   
    
    hold on
    if i == length(data)
        
        ph(5) = plot(SDsaveyears.y(month:12:end),SDsave.wa(presind,month:12:end)-nanmean(SDsave.wa(presind,month:12:end))+...
            nanmean(data(i).wa(presind,find(years(1).y == SDsaveyears.y(1),1)-1+month:12:find(years(1).y == SDsaveyears.y(end),1)-1)),'color','k','LineWidth',2);
        ph(3) = plot(years(i).y(month:12:end),data(i).wa(presind,month:12:end),'color',[177,89,40]/255,'LineWidth',2); 
        ph(4) = plot(trendyears(1,1):trendyears(1,2),data_b(i,month,presind,2)*[1:length(trendyears(1,1):trendyears(1,2))]+data_b(i,month,presind,1),'LineWidth',3,'color',[228,26,28]/255);
        plot(trendyears(2,1):trendyears(2,2),data_b2(i,month,presind,2)*[1:length(trendyears(2,1):trendyears(2,2))]+data_b2(i,month,presind,1),'LineWidth',3,'color',[228,26,28]/255);       
    else
        ph(1) = plot(years(i).y(month:12:end),data(i).wa(presind,month:12:end),'color',[.7 .7 .7],'LineWidth',1);
        ph(2) = plot(trendyears(1,1):trendyears(1,2),data_b(i,month,presind,2)*[1:length(trendyears(1,1):trendyears(1,2))]+data_b(i,month,presind,1),'LineWidth',3,'color',[.5 .5 .5]);
        plot(trendyears(2,1):trendyears(2,2),data_b2(i,month,presind,2)*[1:length(trendyears(2,1):trendyears(2,2))]+data_b2(i,month,presind,1),'LineWidth',3,'color',[.5 .5 .5])    ;
    end
    
end    

set(gca,'fontsize',fsize);
set(gca,'xtick',1950:10:2100,'xticklabel',1950:10:2100);
lh = legend(ph,'Ensemble members','Ensemble trends','Ensemble average','Ensemble average trend','Specified Dynamics (scaled)');
set(lh,'location','NorthEast','box','off','fontsize',fsize+2);
set(gca,'box','on');
xlabel('Year','fontsize',fsize+2);
ylabel(linetitle,'fontsize',fsize+2);
title([monthnames{month},'{ }',sim,', ',area,', ',num2str(pres),' hPa'],'fontsize',fsize+2);
xlim(xlimits);
end
