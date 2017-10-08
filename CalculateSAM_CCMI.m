% Calculate CCMI SAM
clear all

constantyears = 1;

directory = ['/Volumes/My Book for Mac/work/data/CESM-CCMI/PSL/'];
files = dir([directory,'*.nc*']);
files([10,11,12,19]) = [];
    
[years,~,yearsfiles,~] = ReadinCESMDateandPS(1,0);

%% matching existing files
x = zeros(1,length(yearsfiles));
for i = 1:length(yearsfiles)
    for j = 1:length(files)
        if regexp(files(j).name,yearsfiles(i).name)
            x(i) = 1;        
        end
    end
end
years(~x) = [];

%% read in O3 data
for i = 1:length(files);    
    [data(i).attributes, data(i).data,data(i).info] = Read_in_netcdf([directory,files(i).name]);
end
data(length(files)+1).data.PSL = nanmean(cat(4,data(1).data.PSL,data(2).data.PSL,data(3).data.PSL),4);
data(length(files)+2).data.PSL = nanmean(cat(4,data(4).data.PSL,data(5).data.PSL,data(6).data.PSL),4);
data(length(files)+3).data.PSL = nanmean(cat(4,data(7).data.PSL,data(8).data.PSL,data(9).data.PSL),4);
data(length(files)+4).data.PSL = nanmean(cat(4,data(10).data.PSL,data(11).data.PSL,data(12).data.PSL,...
    data(13).data.PSL,data(14).data.PSL),4);

years(length(files)+1).y = years(1).y;
years(length(files)+2).y = years(4).y;
years(length(files)+3).y = years(7).y;
years(length(files)+4).y = years(10).y;

%% calculating SAM
basedates = [1979,2014];
lats = [-40,-65];
[~,latindex(1)] = min(abs(data(1).data.lat - lats(1)));    
[~,latindex(2)] = min(abs(data(2).data.lat - lats(2)));    
%% testing
for i = 1:12
    test1(i,:) = squeeze(nanmean(data(15).data.PSL(:,latindex(1),i:12:end)));
    test1_mean(i) = nanmean(test1(i,:),2);
    test1_std(i) = nanstd(test1(i,:),0,2);
    test2(i,:) = squeeze(nanmean(data(15).data.PSL(:,latindex(2),i:12:end)));
    test2_mean(i) = nanmean(test2(i,:),2);
    test2_std(i) = nanstd(test2(i,:),0,2);
    test1_norm(i,:) = (test1(i,:) - test1_mean(i))./test1_std(i);
    test2_norm(i,:) = (test2(i,:) - test2_mean(i))./test2_std(i);
end

% SAM = test1_norm - test2_norm;
% plot(SAM(2,:))
% hold on
%plot(SAM(1,:))
%plot(SAM(2,:))

% test1_norm = (test1-nanmean(test1))./nanstd(test1);
% test2_norm = (test2-nanmean(test2))./nanstd(test2);
% test_sam = test1_norm - test2_norm;
% plot(test_sam(2:12:end))

%% monthly
for j = 1:length(data)
    for k = 1:length(basedates)
        dateindex(k,j) = find(years(j).y == basedates(k),1,'last')./12;            
    end            
    for l = 1:12
        for i = 1:length(lats)
            [~,latindex(i)] = min(abs(data(1).data.lat - lats(i)));    
            PSL_zm(i,j).a(l,:) = squeeze(nanmean(data(j).data.PSL(:,latindex(i),l:12:end)));        
            
            PSL_base(i,j).a(l) = nanmean(PSL_zm(i,j).a(l,dateindex(1,j):dateindex(2,j)),2);
            PSL_base_std(i,j).a(l) = nanstd(PSL_zm(i,j).a(l,dateindex(1,j):dateindex(2,j)),0,2);
            PSL_anomaly(i,j).a(l,:) = (PSL_zm(i,j).a(l,:) - PSL_base(i,j).a(l))./PSL_base_std(i,j).a(l);
        end        
        PSL_SAM(j).a(l,:) = PSL_anomaly(1,j).a(l,:) - PSL_anomaly(2,j).a(l,:);
    end
    for l = 1:4
        for i = 1:length(lats)
            [~,latindex(i)] = min(abs(data(1).data.lat - lats(i)));               
            if l == 1
                PSL_zm_seasons(i,j).a(l,:) = [zeros(1,1),nanmean([PSL_zm(i,j).a(1,2:end);PSL_zm(i,j).a(2,2:end);PSL_zm(i,j).a(12,1:end-1)])];                
                PSL_zm_seasons(i,j).a (PSL_zm_seasons(i,j).a == 0) = NaN;
            elseif l == 2
                PSL_zm_seasons(i,j).a(l,:) = nanmean([PSL_zm(i,j).a(3,1:end);PSL_zm(i,j).a(4,1:end);PSL_zm(i,j).a(5,1:end)]);                
            elseif l == 3 
                PSL_zm_seasons(i,j).a(l,:) = nanmean([PSL_zm(i,j).a(6,1:end);PSL_zm(i,j).a(7,1:end);PSL_zm(i,j).a(8,1:end)]);                
            elseif l == 4
                PSL_zm_seasons(i,j).a(l,:) = nanmean([PSL_zm(i,j).a(9,1:end);PSL_zm(i,j).a(10,1:end);PSL_zm(i,j).a(11,1:end)]);                
            end
            PSL_base_seasons(i,j).a(l) = nanmean(PSL_zm_seasons(i,j).a(l,dateindex(1,j):dateindex(2,j)),2);
            PSL_base_std_seasons(i,j).a(l) = nanstd(PSL_zm_seasons(i,j).a(l,dateindex(1,j):dateindex(2,j)),0,2);
            PSL_anomaly_seasons(i,j).a(l,:) = (PSL_zm_seasons(i,j).a(l,:) - PSL_base_seasons(i,j).a(l))./PSL_base_std_seasons(i,j).a(l);
        end        
        PSL_SAM_seasons(j).a(l,:) = PSL_anomaly_seasons(1,j).a(l,:) - PSL_anomaly_seasons(2,j).a(l,:);
    end
end

%% smoothing data
%smoothtype = 'moving';
smoothtype = '';
if strcmp(smoothtype,'rlowess');
    smoothlength = .25;
    smoothlength2 = 10;
else
    smoothlength = 20;
    smoothlength2 = 20;
end
for i = 1:length(data)
    for j = 1:size(PSL_SAM(1).a,1)
        smooth_PSL_SAM(i).a(j,:) = smooth(PSL_SAM(i).a(j,:),smoothlength,smoothtype);        
    end
    smooth_PSL_SAM(i).a = [zeros(size(PSL_SAM(1).a,1),smoothlength2/2),...
        smooth_PSL_SAM(i).a(:,smoothlength2/2+1:end-smoothlength2/2),zeros(size(PSL_SAM(1).a,1),smoothlength2/2)];
    smooth_PSL_SAM(i).a (smooth_PSL_SAM(i).a == 0) = NaN;
    for j = 1:size(PSL_SAM_seasons(1).a,1)
        smooth_PSL_SAM_seasons(i).a(j,:) = smooth(PSL_SAM_seasons(i).a(j,:),smoothlength,smoothtype);                
    end
    smooth_PSL_SAM_seasons(i).a = [zeros(size(PSL_SAM_seasons(1).a,1),smoothlength2/2),...
        smooth_PSL_SAM_seasons(i).a(:,smoothlength2/2+1:end-smoothlength2/2),zeros(size(PSL_SAM_seasons(1).a,1),smoothlength2/2)];
    smooth_PSL_SAM_seasons(i).a (smooth_PSL_SAM_seasons(i).a == 0) = NaN;
end

%% plotting
runnames = {'REF-C2 no. 1','REF-C2 no. 2','REF-C2 no. 3',...
    'SEN-C2-fGHG no. 1','SEN-C2-fGHG no. 3','SEN-C2-fGHG no. 3',...
    'SEN-C2-fODS no. 1','SEN-C2-fODS no. 2','SEN-C2-fODS no. 3',...
    'REF-C1 no. 1','REF-C1 no. 2','REF-C1 no. 3','REF-C1 no. 4','REF-C1 no. 5',...
    'REF-C1SD no .1','REF-C2 ens','SEN-C2-fGHG ens','SEN-C2-fODS ens','REF-C1 ens'};

plotnames = {'REF-C2','SEN-C2-fGHGs','SEN-C2-fODSs',...
    'REF-C1','REF-C1 - specified dynamics','Ensemble averages'};

fsize = 18;

% creating plotting array
season = 1;
toplot(1).a = [smooth_PSL_SAM_seasons(1).a(season,:);smooth_PSL_SAM_seasons(2).a(season,:);...
    smooth_PSL_SAM_seasons(3).a(season,:);nanmean([smooth_PSL_SAM_seasons(1).a(season,:);smooth_PSL_SAM_seasons(2).a(season,:);...
    smooth_PSL_SAM_seasons(3).a(season,:)])];
toplotyears(1).a = years(1).y;
toplot(2).a = [smooth_PSL_SAM_seasons(4).a(season,:);smooth_PSL_SAM_seasons(5).a(season,:);...
    smooth_PSL_SAM_seasons(6).a(season,:);nanmean([smooth_PSL_SAM_seasons(4).a(season,:);smooth_PSL_SAM_seasons(5).a(season,:);...
    smooth_PSL_SAM_seasons(6).a(season,:)])];
toplotyears(2).a = years(4).y;
toplot(3).a = [smooth_PSL_SAM_seasons(7).a(season,:);smooth_PSL_SAM_seasons(8).a(season,:);...
    smooth_PSL_SAM_seasons(9).a(season,:);nanmean([smooth_PSL_SAM_seasons(7).a(season,:);smooth_PSL_SAM_seasons(8).a(season,:);...
    smooth_PSL_SAM_seasons(9).a(season,:)])];
toplotyears(3).a = years(7).y;
toplot(4).a = [smooth_PSL_SAM_seasons(10).a(season,:);smooth_PSL_SAM_seasons(11).a(season,:);...
    smooth_PSL_SAM_seasons(12).a(season,:);smooth_PSL_SAM_seasons(13).a(season,:);...
    smooth_PSL_SAM_seasons(14).a(season,:);nanmean([smooth_PSL_SAM_seasons(10).a(season,:);smooth_PSL_SAM_seasons(11).a(season,:);...
    smooth_PSL_SAM_seasons(12).a(season,:);smooth_PSL_SAM_seasons(13).a(season,:);...
    smooth_PSL_SAM_seasons(14).a(season,:)])];
toplotyears(4).a = years(10).y;
toplot(5).a = smooth_PSL_SAM_seasons(15).a(season,:);    
toplotyears(5).a = years(15).y;

createfig('large');
colors = cbrewer('qual','Set1',10);
for i = 1:6

    sp(i) = subplot(3,2,i);
    sp_pos(i,:) = get(sp(i),'position');

    if i < 6
        for j = 1:size(toplot(i).a,1)
            if j == size(toplot(i).a,1)
                plot(toplotyears(i).a(1:12:end),toplot(i).a(j,:),'LineWidth',3,'color',colors(2,:))
            else
                plot(toplotyears(i).a(1:12:end),toplot(i).a(j,:),'color',[.7 .7 .7],'LineWidth',2)            
            end
            hold on
        end
    else    
        for j = 1:4;
            plot(toplotyears(j).a(1:12:end),toplot(j).a(end,:),'LineWidth',3,'color',colors(j,:))
            hold on            
        end
        lh = legend('REF-C2','SEN-C2-fGHGs','SEN-C2-fODSs','REF-C1');
        set(lh,'fontsize',fsize-6,'location','NorthWest','box','off');
       
    end
    xlim([1954 2100])
    ylim([-2.5 2.5])
    title(plotnames{i},'fontsize',fsize+2);
    set(gca,'fontsize',fsize-2,'xtick',1960:20:2100,'xticklabel',1960:20:2100);
    if mod(i,2) 
        ylabel('SAM','fontsize',fsize+2)
    end
    if i == 5 || i == 6
        xlabel('Year','fontsize',fsize+2);
    end    



end

moveplots(sp,.08,.08);
mtit('Southern annular mode index','xoff',-.01,'yoff',.04,'fontsize',fsize+4);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/SAM/SAM',smoothtype,'.pdf'];                
export_fig(filename,'-pdf');
