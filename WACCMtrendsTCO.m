% Read in and construct trend data for maps

%clear all;
clc;
clear all;
plotlinediff = 0;
meanvoldiff = 0;
plot_maps = 1;

runtype = 'monthly';
TCOdirectory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/TCO/'];

%addpath(genpath('/home/stonek/code'))

TCOfiles = dir([TCOdirectory,'TOZ*']);

yearmin = 20000201;
yearmax = 20150201;
remove2002 = 'Southern'; % 'all_lats' or 'Southern' or 'none'

cmap = flipud(cbrewer('div','RdBu',10));

TCOname = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};

halon = zeros(6,144);
halat = zeros(6,96);
halon (halon == 0) = NaN;
halat (halat == 0) = NaN;

lats = [[30 60];[60 90]];
%lats = [[-90 -60];[-60 -30];[-30 0];[0 30];[30 60];[60 90]];

szlat = size(lats);
TCOzonallat = [];

for i = 1:length(TCOfiles)
    [info.(TCOname{i}), data.(TCOname{i}), attributes.(TCOname{i})] = ...
        Read_in_netcdf([TCOdirectory,TCOfiles(i).name]);
    
    %zonal mean data
    TCOzonal.(TCOname{i}) = squeeze(nanmean(data.(TCOname{i}).toz));
    
    for p1 = 1:szlat(1)
        latindex(p1).p = find(data.(TCOname{i}).lat >= lats(p1,1) & ...
            data.(TCOname{i}).lat <= lats(p1,2));
    end
    
    %removing leap years
    dates = num2str(data.(TCOname{i}).date);
    monthday = str2num(dates(:,5:8));
    years =  str2num(dates(:,1:4));
    monthdayindex = find(monthday == 0229);
    if strcmp(runtype,'daily')
        TCOzonal.(TCOname{i})(:,monthdayindex) = []; 
        monthday(monthdayindex) = [];
        data.(TCOname{i}).date(monthdayindex) = [];
        interval = 365;
        %extracting time period
        monthday = monthday(data.(TCOname{i}).date >= 20000101 & data.(TCOname{i}).date < 20150101);
        TCOzonal20002014.(TCOname{i}) = TCOzonal.(TCOname{i})(:,...
            data.(TCOname{i}).date >= 20000101 & data.(TCOname{i}).date < 20150101);
        TCOzonal20002014.(TCOname{i})(:,973:1338) = NaN; 
    else
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
                
                Latitudes.(TCOname{i})(data.(TCOname{i}).lat < 0,...
                datetemp >= 20021001 & datetemp < 20031001) = NaN;                
        end
                        
    end
    
    if strcmp(runtype,'daily')
       time_xaxis = 365;
    else
       time_xaxis = 12; 
    end

    %calculate yearly trends through least square
    for j = 1:time_xaxis
        toztrends.(TCOname{i})(:,j) = trend(TCOzonal20002014.(TCOname{i})(:,j:time_xaxis:end),[],2);
        for k = 1:size(TCOzonal20002014.(TCOname{i}),1)
            x = TCOzonal20002014.(TCOname{i})(k,j:interval:end)';
            %[b(i,j,k,:), bint(i,j,k,:,:), r(i,j,k,:), rint(i,j,k,:,:),...
            %     stats(i,j,k,:)] = regress(TCOzonal20002014.(TCOname{i})(k,j:12:end)',regfun); 
            stats = regstats(x,1:length(x),'linear',{'tstat','rsquare','yhat','r'});
            stats2 = regress(x,[ones(length(x),1),(1:length(x))']);
            p(i,j,k) = stats.tstat.pval(2);
            h(i,j,k) = stats.tstat.pval(1);
            b1(i,j,k) = stats.tstat.beta(1);
            b(i,j,k) = stats.tstat.beta(2);
            rsquare(i,j,k) = stats.rsquare;
            residuals(i,j,k,:) = stats.r;  
            yhat(i,j,k,:) = stats.yhat;
            [h_osttest(i,j,k), p_osttest(i,j,k)] = ttest2(squeeze(yhat(i,j,k,:)),...
                TCOzonal20002014.(TCOname{i})(k,j:12:end)',.05,'right');
        end
    end
    toztrends.(TCOname{i}) = circshift(toztrends.(TCOname{i}),[0 -1]);
    test = 1;
   
    for k = 1:12
        for m = 1:szlat(1)
            %TCOzonallat(m,k,i,:) = nanmean(TCOzonal20002014.(TCOname{i})(latindex(m).p,k:12:end))';            
            temp = TCOzonal20002014.(TCOname{i})(latindex(m).p,k:12:end);
            WLats = repmat(data.(TCOname{i}).lat,[1,size(temp,2)]);
            WLats = WLats(latindex(m).p,:);
            switch remove2002
            case 'all_lats'            
            case 'Southern'                        
                                
                WLats (isnan(temp)) = NaN;
                
            end
                                    
            W = cosd(WLats);
            TCOzonallat(m,k,i,:) = nansum(temp.*W)./nansum(W);
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
    for k = 1:time_xaxis
        temp = length(nanmean(TCOzonaldiff.(diffnames{i})(latindex(1).p,k:12:end))');
        for j = 1:size(TCOzonal20002014.(TCOname{i}),1)
%             [~, TCOzonallatdiff(j,k,i,1:temp)] = TCOanomaly(TCOzonaldiff.(diffnames{i}),...
%                 latitudes(j,:),data.(TCOname{i}).lat,[2000 2014],[2000 2015],k);
            %TCOzonallat1(p,k,i,1:temp) = nanmean(TCOzonal20002014.(TCOname{i})(latindex(p).p,k:12:end))';          
            diffstats(j,k,i) = regstats(squeeze(TCOzonaldiff.(diffnames{i})(j,k:12:end)),...
                1:length(squeeze(TCOzonaldiff.(diffnames{i})(j,k:12:end))),'linear',{'tstat','rsquare','yhat','r'});
            pdiff(i,k,j) = diffstats(j,k,i).tstat.pval(2);            
            [bdiff(j,k,i,:),bdiffint(j,k,i,:,:),~,~] = regress(squeeze(TCOzonaldiff.(diffnames{i})(j,k:12:end))',...
                [ones(1,length(squeeze(TCOzonaldiff.(diffnames{i})(j,k:12:end))));1:length(squeeze(TCOzonaldiff.(diffnames{i})(j,k:12:end)))]',.8);
            bdiff1(i,k,j) = diffstats(j,k,i).tstat.beta(2);
            diffconf(j,k,i,:) = bdiffint(j,k,i,:,2) - bdiff(j,k,i,:);
        end        
    end   
end

%% 
for i = 1:szlat(1)  
    if lats(i,1) < 0 && lats(i,1) > 0
        titles2{i} = {[num2str(abs(lats(i,1))),'{\circ}S-',num2str(abs(lats(i,2))),'{\circ}N']};    
    elseif lats(i,1) < 0 && lats(i,1) < 0
        titles2{i} = {[num2str(abs(lats(i,1))),'-',num2str(abs(lats(i,2))),'{\circ}S']};
    elseif lats(i,1) >= 0 && lats(i,1) >= 0
        titles2{i} = {[num2str(abs(lats(i,1))),'-',num2str(abs(lats(i,2))),'{\circ}N']};
    end
end


if meanvoldiff
    count = 1;
    for k = 1:szlat
        for i = 1:12
            TCOmeanvoldiff(count,i,:) = squeeze((TCOzonallat(k,i,2,:) - TCOzonallat(k,i,3,:))./TCOzonallat(k,i,3,:)*100);
        end
        count = count + 1;
    end


    TCOmeanvoldiffplot(:,:,:) = TCOmeanvoldiff(:,[3,6,8,10],:);

    legendtitles3 = {'March','June','August','October'};    

    linesfilename2 = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/','percentVolcanoes3090N.pdf'];                

    linefig2 = lineplots(TCOmeanvoldiffplot,[-3.1 .1],'Year','% of TCO',titles2,14,2,legendtitles3,'% TCO difference (MAM - VC-MAM)',meanvoldiff);
    set(linefig2,'position',[100 100 840 560]);
    export_fig(linefig2,linesfilename2,'-pdf');
    close(linefig2);
end

%%
%absolute differences
volcabsdiff(1,:,:) = TCOzonal20002014.(TCOname{2})(:,:)' - TCOzonal20002014.(TCOname{3})(:,:)';

titles = {'CCMI','MAM','VC-MAM','Chem-only','Chem-only fixed SSTs','Chem-only noleap'};

%plotting latitude linies
mon = {'January','February','March','April','May','June','July','August','September','October','November','December'};

if plotlinediff
    for i = 1:12;
        linefig = lineplots(squeeze(TCOzonallat(:,i,:,:)),[],'Year','Dobson Units',titles2,18,2,titles,mon{i},meanvoldiff);
        linesfilename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/',sprintf('%02d',i),'WACCM_ZonalMeanTCO_',mon{i},'.pdf'];                
        set(linefig,'position',[100 100 840 840]);
        export_fig(linefig,linesfilename,'-pdf');
        close(linefig);
    end
end

if plot_maps
    %calculating differences
    diff(1,:,:) = squeeze(b(2,:,:) - b(3,:,:));
    diff(2,:,:) = squeeze(b(4,:,:) - b(5,:,:));
    diff(3,:,:) = squeeze(b(3,:,:) - b(4,:,:));
    diff(4,:,:) = b(5,:,:);

    longitudes = data.(TCOname{1}).lon;
    latitudes = data.(TCOname{1}).lat;

    ptest = p;
    for i = 1:numel(p)
        if p(i) > .2
            p(i) = p(i)-1.01;
        end
    end
    for i = 1:numel(pdiff)
        if pdiff(i) > .2
            pdiff(i) = pdiff(i)-1.01;
        end
    end

    if strcmp(runtype,'daily');
        x_intervals = 1:365;
    %    b = cat(2,b,b(:,1:31,:));
    %    p = cat(2,p,p(:,1:31,:));
    %    diff = cat(2,diff,diff(:,1:31,:));
    else
        x_intervals = 1:13;
        b = cat(2,b,b(:,1,:));
        p = cat(2,p,p(:,1,:));
        bdiff1 = cat(2,bdiff1,bdiff1(:,1,:));
        pdiff = cat(2,pdiff,pdiff(:,1,:));
    end

    titlediff = {'MAM - VC-MAM','Chem only - Chem only fSSTs','VC-MAM - Chem-only','Chem-only'};
    diff_namelist = {'Volcanoes (MAM - VC-MAM)','SSTs (Chem-only - Chem-only-fSSTs)','Dynamics (VC-MAM - Chem-only)','Chemistry (Chem-only-fSSTs)'};
    %set(gca,'xtick',[1,32,60,92,123,155,186,208,239,271,302,334],'xticklabel',['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'],'ytick',-90:30:90,'yticklabel',-90:30:90,'fontsize',fsize)
    if strcmp(runtype,'daily')
        xt = [1,32,60,92,123,155,186,208,239,271,302,334];
    else
        xt = 1:1:12;
    end
    months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];

    fig = subplotmaps(double(b),x_intervals,double(latitudes),{'div','RdBu'},1,p,14,titles,'Months',...
        'Latitude ({\circ}N)','DU/year','on',[-2,2],22,xt,months,-90:30:90,-90:30:90,'WACCM 2000-2015 trends',...
        0,[1 12]);
    %rgb2cm
    
    set(fig,'position',[100 100 840 840]);

    bdiff1 (bdiff1 < -1) = -1;
    fig2 = subplotmaps(double(bdiff1),x_intervals,double(latitudes),{'div','RdBu'},1,pdiff,14,diffnames,'Months',...
        'Latitude ({\circ}N)','DU/year','on',[-1,1],22,xt,months,-90:30:90,-90:30:90,'WACCM 2000-2015 trend contributions',...
        0,[1 12]);

    set(fig2,'position',[100 100 840 560]);

    absdiffig = subplotmaps(double(volcabsdiff),1:180,double(latitudes),...
        {'seq','YlOrBr'},0,[],20,{'TCO difference due to volcanoes (MAM - VC-MAM)'},'Year','Latitude ({\circ}N)','Dobson Units',...
        'on',[-12,0],13,1:24:180,2000:2:2015,-90:30:90,-90:30:90);

    set(absdiffig,'position',[100 100 840 560]);

    if strcmp(runtype,'daily')
         filename = '/home/stonek/projects/ozoneVolcanoes/figures/WACCM_TCOtrends19992014_daily.eps';

         filename2 = '/home/stonek/projects/ozoneVolcanoes/figures/WACCM_TCOtrends19992014_daily_diff.eps'; 
    else filename = '/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/WACCM_TCOtrends20002014.png';
        filename2 = '/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/WACCM_TCOtrends20002014_diff.png';
        filename3 = '/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/VolcanicDiff20002014.png';

    end
    %rgb2cm

    export_fig(fig,filename,'-png');

    export_fig(fig2,filename2,'-png');

    export_fig(absdiffig,filename3,'-png');
end

