% Produce trend lines for zonal mean WACCM data.

TCOdirectory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/TCO/'];

plotzonaltrends = 0;
barplottrends = 1;

latitudes = [[-90 -60];[-60 -30];[30 60];[60 90]];
%latitudes = [[-60 60];[-15 15];[-90 0];[0 90]];

szlat = size(latitudes);

yearmin = 20000201;
yearmax = 20150201;
yearmax1 = num2str(yearmax);
yearmax1 = str2double(yearmax1(1:4));

remove2002 = 'Southern'; % 'all_lats' or 'Southern' or 'none'
%addpath(genpath('/home/stonek/code'))

cbrew = cbrewer('qual','Set1',10);
cbrew([5,7],:) = [];
cbrewfade = cbrewer('qual','Pastel1',10);

TCOfiles = dir([TCOdirectory,'TOZ*']);

cmap = flipud(cbrewer('div','RdBu',10));

TCOname = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
TCOzonallat = zeros(6,12,6,16);
TCOzonallat (TCOzonallat == 0) = NaN;

TCOzonallatdiff = zeros(6,12,6,16);
TCOzonallatdiff (TCOzonallatdiff == 0) = NaN;

for i = 1:length(TCOfiles)
    [info.(TCOname{i}), data.(TCOname{i}), attributes.(TCOname{i})] = ...
        Read_in_netcdf([TCOdirectory,TCOfiles(i).name]);
    
    %finding latitudes
    for p = 1:szlat(1)
        latindex(p).p = find(data.(TCOname{i}).lat >= latitudes(p,1) & ...
            data.(TCOname{i}).lat <= latitudes(p,2));
    end
    
    %zonal mean data
    TCOzonal.(TCOname{i}) = squeeze(nanmean(data.(TCOname{i}).toz));
    
    %removing leap years
    dates = num2str(data.(TCOname{i}).date);
    monthday = str2num(dates(:,5:8));
    years =  str2num(dates(:,1:4));
    monthdayindex = find(monthday == 0229);
    
    interval = 12;
    %extracting time period

    TCOzonal20002014.(TCOname{i}) = TCOzonal.(TCOname{i})(:,...
        data.(TCOname{i}).date >= yearmin & data.(TCOname{i}).date < yearmax);
    datetemp = data.(TCOname{i}).date(data.(TCOname{i}).date >= yearmin & data.(TCOname{i}).date < yearmax);
    %removing 2002 from September
    
    switch remove2002
        case 'all_lats'
            TCOzonal20002014.(TCOname{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        case 'Southern'        
            TCOzonal20002014.(TCOname{i})(data.(TCOname{i}).lat < 0,...
                datetemp >= 20021001 & datetemp < 20031001) = NaN; 
    end
    
    if yearmax > 20151101
        %adding in NaNs for missing December
        TCOzonal20002014.(TCOname{i}) = cat(2,TCOzonal20002014.(TCOname{i}),TCOzonal20002014.(TCOname{i})(:,end));
        TCOzonal20002014.(TCOname{i})(:,end) = NaN; 
    end
    % computing differences
    
    for k = 1:12
        temp = length(nanmean(TCOzonal20002014.(TCOname{i})(latindex(1).p,k:12:end))');
        for p = 1:szlat(1)
            [~, TCOzonallat(p,k,i,1:temp)] = TCOanomaly(TCOzonal20002014.(TCOname{i}),...
                latitudes(p,:),data.(TCOname{i}).lat,[2000 yearmax1-1],[2000 yearmax1],k);
            %TCOzonallat1(p,k,i,1:temp) = nanmean(TCOzonal20002014.(TCOname{i})(latindex(p).p,k:12:end))';          
            stats(p,k,i) = regstats(squeeze(TCOzonallat(p,k,i,:)),...
                1:length(squeeze(TCOzonallat(p,k,i,:))),'linear',{'tstat','rsquare','yhat','r'});
            [b(p,k,i,:),bint(p,k,i,:,:),~,~] = regress(squeeze(TCOzonallat(p,k,i,:)),...
                [1:length(squeeze(TCOzonallat(p,k,i,:)));ones(1,length(squeeze(TCOzonallat(p,k,i,:))))]',.2);
            conf(p,k,i,:) = bint(p,k,i,:,1) - b(p,k,i,:);
        end        
    end        
    
end

%% Computing differences
diffnames = {'Volcanoes','Dynamics','Chemistry','ssts'};
TCOzonaldiff.(diffnames{1}) = TCOzonal20002014.(TCOname{2})-TCOzonal20002014.(TCOname{3});
TCOzonaldiff.(diffnames{2}) = TCOzonal20002014.(TCOname{3})-TCOzonal20002014.(TCOname{4});
TCOzonaldiff.(diffnames{3}) = TCOzonal20002014.(TCOname{5});
TCOzonaldiff.(diffnames{4}) = TCOzonal20002014.(TCOname{4})-TCOzonal20002014.(TCOname{5});

for i = 1:length(diffnames);
    for k = 1:12
        temp = length(nanmean(TCOzonaldiff.(diffnames{i})(latindex(1).p,k:12:end))');
        for p = 1:szlat(1)
            [~, TCOzonallatdiff(p,k,i,1:temp)] = TCOanomaly(TCOzonaldiff.(diffnames{i}),...
                latitudes(p,:),data.(TCOname{i}).lat,[2000 yearmax1-1],[2000 yearmax1],k);
            %TCOzonallat1(p,k,i,1:temp) = nanmean(TCOzonal20002014.(TCOname{i})(latindex(p).p,k:12:end))';          
            diffstats(p,k,i) = regstats(squeeze(TCOzonallatdiff(p,k,i,:)),...
                1:length(squeeze(TCOzonallatdiff(p,k,i,:))),'linear',{'tstat','rsquare','yhat','r'});
            pdiff(p,k,i) = diffstats(p,k,i).tstat.pval(2);    
            [bdiff(p,k,i,:),bdiffint(p,k,i,:,:),~,~] = regress(squeeze(TCOzonallatdiff(p,k,i,:)),...
                [1:length(squeeze(TCOzonallatdiff(p,k,i,:)));ones(1,length(squeeze(TCOzonallatdiff(p,k,i,:)))),]',.2);
            diffconf(p,k,i,:) = bdiffint(p,k,i,:,2) - bdiff(p,k,i,:);
        end        
    end   
end

%% plot line trends
%universal label names and graph constants
mon = {'January','February','March','April','May','June','July','August','September','October','November','December'};
mons = {'J','F','M','A','M','J','J','A','S','O','N','D'};
lnames = {'CCMI','MAM','VC-MAM','Chem-only','Chem-only-fixedSSTs','Chem-only-noleap'};
diff_namelist = {'Volcanoes','SSTs','Dynamics','Chemistry'};
fsize = 20;
lsize = 3;

if plotzonaltrends
    
    titles = {'90-60S','60-30S','30-0S','0-30N','30-60N','60-90N'};    
    for j = 1:6;
        plottitle = titles{j};
        for k = 1:12;
            figure;
            fig = gcf;
            set(fig,'position',[100 100 840 560],'Visible','off');
            %ap = get(gca,'position');
            %set(gca,'position',[ap(1), ap(2), ap(3)-.2, ap(4)]);
            for i = 1:6
                ph2(i) = plot(squeeze(TCOzonallat(j,k,i,:)),'--o','LineWidth',lsize-1);
                set(ph2(i),'color',cbrew(i,:));            
                hold on
            end
            for i = 1:6
                ph(i) = plot((1:length(squeeze(TCOzonallat(j,k,i,:))))*stats(j,k,i).tstat.beta(2)...
                    +stats(j,k,i).tstat.beta(1),'LineWidth',lsize+1);
                set(ph(i),'color',cbrew(i,:));
            end
            set(gca,'xtick',1:2:20,'xticklabel',2000:2:2020,'fontsize',fsize-2);        
            xlabel('Year','fontsize',fsize);
            ylabel('Dobson Units','fontsize',fsize);
            title(['TCO','{ }',plottitle,'{ }', mon{k}],'fontsize',fsize+2);
            xlim([0 17]);
            lh = legend(ph,lnames,'location','NorthEastOutside','fontsize',fsize-4,'box','off');

            annotation('textbox',[.69 .60 .25 .05],...
                    'String',{['{\beta}' '{ }' setstr(177) '{ }' '80% conf']},...
                    'fontsize',fsize-2,...
                    'Color','k',...
                    'EdgeColor','none');

            for l = 1:6
                annotation('textbox',[.69 .60-(l/20) .25 .05],...
                    'String',{[sprintf('%02.2f',(stats(j,k,l).tstat.beta(2))) '{ }' setstr(177) '{ }' sprintf('%02.2f',conf(j,k,l,1))]},...
                    'fontsize',fsize-2,...
                    'Color',cbrew(l,:),...
                    'EdgeColor','none');
            end

            hold off
            export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/linetrends/',...
                sprintf('%02d',k),'WACCM_trends_',plottitle,'.pdf']);
            close(fig);
        end
    end
end

% %calculating differences
% diff(:,:,1,:) = squeeze(b(:,:,2,:) - b(:,:,3,:));
% diff(:,:,2,:) = squeeze(b(:,:,4,:) - b(:,:,5,:));
% diff(:,:,3,:) = squeeze(b(:,:,3,:) - b(:,:,4,:));
% diff(:,:,4,:) = b(:,:,5,:);
% 
% %ttest2(squeeze(b(:,:,2,:),b(:,:,3,:)));
% 
% diffconf(:,:,1,:) = squeeze(conf(:,:,2,:) - conf(:,:,3,:));
% diffconf(:,:,2,:) = squeeze(conf(:,:,4,:) - conf(:,:,5,:));
% diffconf(:,:,3,:) = squeeze(conf(:,:,3,:) - conf(:,:,4,:));
% diffconf(:,:,4,:) = conf(:,:,5,:);

%% plotting
difference = 0;
if difference
   b = bdiff;
   conf = diffconf;
end

if barplottrends
    if ~difference
        %btitles = {'SouthMidandPole_30degreebin.pdf','NorthMidandPole_30degreebin.pdf'};
        btitles = {'SouthMidandPole_30degreebin1.pdf','NorthMidandPole_30degreebin1.pdf'};
    else
        %btitles = {'DiffSouthMidandPole_30degreebin.pdf','DiffNorthMidandPole_30degreebin.pdf'};
        btitles = {'DiffSouthMidandPole_30degreebin.pdf','DiffNorthMidandPole_30degreebin1.pdf'};
    end
    count = 0;
    count1 = 1;
    for j = 1:2
        figure;
        fig = gcf;
        set(fig,'color','white','position',[100 100 840 840],'visible','off');     
        %axespos = get(gca,'position');
        %set(gca,'position',[axespos(1) axespos(2)+.2 axespos(3) axespos(4)-.2]);
        if szlat(1) == 1
            iend = 1;
        else iend = szlat(1)/2;
        end
        for i = 1:iend;
        
            sp(i) = subplot(length(1:iend),1,i);
            sppos = get(sp(i),'position');
            set(sp(i),'position',[sppos(1),sppos(2)+.2./szlat(1)/2,sppos(3),sppos(4)-.2./szlat(1)/2]);
            regcoef = squeeze(b(i+count,:,:,1));
            regconf = squeeze(conf(i+count,:,:,1)); 
            
            if ~difference
                regcoef(:,[1,6]) = [];
                regconf(:,[1,6]) = [];
            end
            
            [phbar, herrbar] = barwitherr(regconf,regcoef,1);
            set(phbar(i),'BarWidth',1);
            bxd = get(herrbar(i),'xdata');
            for k = 1:length(phbar)
                set(herrbar(k),'color',cbrew(k,:));
                set(phbar(k),'EdgeColor',cbrew(k,:),'LineWidth',1.5);
            end            
            colormap(cbrewfade(1:size(regcoef,2),:));

            hold on
            plot(repmat(1.5:1:size(regcoef,1)-.5,11,1),repmat(-4:1:6,size(regcoef,1)-1,1)','--k','LineWidth',.5);   
            if ~difference
                ylim([-3 5]);
                title([num2str(latitudes(count1,1)),' to ',num2str(latitudes(count1,2)),'{\circ}N{ }',...
                    'TCO trends (2000-2014)'],'fontsize',fsize+2);
            else
                ylim([-2 2]);
                title([num2str(latitudes(count1,1)),' to ',num2str(latitudes(count1,2)),'{\circ}N{ }',...
                    'TCO trend contributions (2000-2014)'],'fontsize',fsize+2);
            end
            xlim([.5 size(regcoef,1)+.5]);
            set(gca,'xticklabel',mons,'ytick',-5:1:5,'fontsize',fsize-2)
            ylabel('DU/year','fontsize',fsize);
            xlabel('month','fontsize',fsize);            
            count1 = count1 + 1;            
        
        end        
        count = count+iend;
         %creating axes for legend
        if ~difference
            TCOnameleg = {lnames{1,2},lnames{1,3},lnames{1,4},lnames{1,5}};
        else
            TCOnameleg = diffnames;
        end
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
            'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        lh = legend(phbar,TCOnameleg);        
        set(lh,'fontsize',fsize,'box','off','Orientation','horizontal','position',[0.5, .03, .01, .01]);
        if ~difference
            export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/barplots/',...
                btitles{j}]);
        else
            export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/barplots/',...
                btitles{j}]);
        end
    end           
end

close all



