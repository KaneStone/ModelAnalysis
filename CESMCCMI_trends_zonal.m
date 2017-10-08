% Calculate Temperature trends for CESM CCMI simulations.
clear all

read_in = 0;
remove_SD = 1;
line_plots = 0;
plot_all = 1;
plot_ens = 1;
hindcast = 1;
directory = '/Volumes/MyBook/work/data/CESM-CCMI/T/50hPa/';
files = dir([directory,'9030*']);

if remove_SD
    for i = 1:length(files)
        if regexp(files(i).name,'REFC1SD');
            a(i) = 1;
        else a(i) = 0;
        end
    end
    files = files(~a);
end

trendyears = [1960,2000];

if read_in
    %% Read in data

    for i = 1:length(files);
        [~,data(i).d,~] = Read_in_netcdf([directory,files(i).name]);
        years_temp = num2str(data(i).d.date);
        for j = 1:length(years_temp)
            years(i).y(j) = str2double(years_temp(j,1:4));
        end    
        years(i).y = circshift(years(i).y,[0,1]);
        years(i).y(1) = years(i).y(2);
        if years(i).y(end) >= trendyears(2)
            dateind(i,:) = find(years(i).y >= trendyears(1) & years(i).y <= trendyears(2));
            dataatdate(i,:,:,:) = squeeze(data(i).d.T(:,:,:,dateind(i,:)));
        else
            hindcast = 0;
            continue
        end

        for l = 1:length(data(1).d.lon)
            for j = 1:length(data(1).d.lat)
                for k = 1:12
                    [b(i,l,j,k,:),bint(i,l,j,k,:,:)] = regress(squeeze(dataatdate(i,l,j,k:12:end)),...
                        [ones(1,size(dataatdate(i,l,j,k:12:end),4));...
                        1:size(dataatdate(i,l,j,k:12:end),4)]');      
                end
            end   
        end
    end

    %% taking ensemble averages

    lon = data(1).d.lon;
    lat = data(1).d.lat;

    data_ens(1).d = nanmean(cat(5,data(1).d.T,data(2).d.T,data(3).d.T),5);
    data_ens(2).d = nanmean(cat(5,data(4).d.T,data(5).d.T,data(6).d.T),5);
    data_ens(3).d = nanmean(cat(5,data(7).d.T,data(8).d.T,data(9).d.T),5);
    if hindcast
        data_ens(4).d = nanmean(cat(5,data(10).d.T,data(11).d.T,data(12).d.T,data(13).d.T,data(14).d.T),5);
    end
    years_ens(1).y = years(1).y;
    years_ens(2).y = years(1).y;
    years_ens(3).y = years(1).y;
    if hindcast
        years_ens(4).y = years(end).y;
    end

    for i = 1:length(data_ens);

        dateind(i,:) = find(years_ens(i).y >= trendyears(1) & years_ens(i).y <= trendyears(2));
        dataatdate_ens(i,:,:,:) = squeeze(data_ens(i).d(:,:,:,dateind(i,:)));

        for l = 1:length(lon)
            for j = 1:length(lat)
                for k = 1:12
                [bens(i,l,j,k,:),bintens(i,l,j,k,:,:)] = regress(squeeze(dataatdate_ens(i,l,j,k:12:end)),...
                    [ones(1,size(dataatdate_ens(i,l,j,k:12:end),4));...
                    1:size(dataatdate_ens(i,l,j,k:12:end),4)]');      
                end
            end   
        end
    end


    %% constructing eddy components
    b_latmean = nanmean(b,2);
    b_latmean_ext = repmat(b_latmean,[1,size(b,2),1,1,1]);
    b_eddy = b - b_latmean_ext;
    b_eddy = cat(2,b_eddy(:,end,:,:,:),b_eddy(:,:,:,:,:));
    b = cat(2,b(:,end,:,:,:),b(:,:,:,:,:));

    bens_latmean = nanmean(bens,2);
    bens_latmean_ext = repmat(bens_latmean,[1,size(bens,2),1,1,1]);
    bens_eddy = bens - bens_latmean_ext;
    bens_eddy = cat(2,bens_eddy(:,end,:,:,:),bens_eddy(:,:,:,:,:));
    bens = cat(2,bens(:,end,:,:,:),bens(:,:,:,:,:));

    if hindcast
        nameprint = {'Ref-C2_no.1','Ref-C2_no.2','Ref-C2_no.3','Sen-C2-fGHG_no.1',...
            'Sen-C2-fGHG_no.2','Sen-C2-fGHG_no.3','Sen-C2-fODS_no.1','Sen-C2-fODS_no.2',...
            'Sen-C2-fODS_no.3'};
        nameprint_ens = {'Ref-C2','Sen-C2-fGHG','Sen-C2-fODS'};
    else
        nameprint = {'Ref-C2_no.1','Ref-C2_no.2','Ref-C2_no.3','Sen-C2-fGHG_no.1',...
            'Sen-C2-fGHG_no.2','Sen-C2-fGHG_no.3','Sen-C2-fODS_no.1','Sen-C2-fODS_no.2',...
            'Sen-C2-fODS_no.3'};
        nameprint_ens = {'Ref-C2','Sen-C2-fGHG','Sen-C2-fODS'};
    end

    save(['/Volumes/MyBook/work/data/CESM-CCMI/output/T_trendsandeddytrends',num2str(trendyears(1)),'-',num2str(trendyears(2)),'.mat'],'b','b_eddy','bens','bens_eddy','lat','lon','nameprint','nameprint_ens');
else
    if line_plots
        for i = 1:length(files);
            [~,data(i).d,~] = Read_in_netcdf([directory,files(i).name]);
        end
    end
    load(['/Volumes/MyBook/work/data/CESM-CCMI/output/T_trendsandeddytrends',num2str(trendyears(1)),'-',num2str(trendyears(2)),'.mat']);
end

%% plotting trends and eddy component 
lontoplot = [lon(end);lon];
mon = {'January','February','March','April','May','June','July','August','September','October',...
    'November','December'};
titles = {'trend','eddy component'};
if plot_all
    contourlims = [-2.5 2.5];
    nameprint2 = nameprint;
    for i = 1:length(nameprint2)
        nameprint2{i} (nameprint2{i} == '_') = ' ';
    end

    for i = 1:size(b,1);
        for j = 10;        
            btoplot = permute(cat(3,squeeze(b(i,:,:,j,2))*10,squeeze(b_eddy(i,:,:,j,2))*10),[3,1,2]);

            fig = subplotmaps(btoplot,lontoplot,lat,{'div','RdBu'},1,[],18,titles,'Longitude {\circ}N','Latitude {\circ}N','T/decade','off',...
                contourlims,22,0:20:360,0:20:360,...
                -90:10:-30,-90:10:30,[num2str(trendyears(1)),'-',num2str(trendyears(2)),', ', mon{j},' - ',nameprint2{i}],1,[0 360],[-90 -30],0,'-',1);        

            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/ZonalAsymmetry/TemperatureTrends/',mon{j},'_',num2str(trendyears(1)),'-',num2str(trendyears(2)),'_',nameprint{i}];            
            export_fig(filename,'-png')
            close(fig);
        end
    end
end
%% plot ensemble
if plot_ens
    contourlims = [-2.5 2.5];
    nameprint2 = nameprint_ens;
    for i = 1:length(nameprint2)
        nameprint2{i} (nameprint2{i} == '_') = ' ';
    end

    for i = 1:length(nameprint_ens)
        for j = 10;        
            benstoplot = permute(cat(3,squeeze(bens(i,:,:,j,2))*10,squeeze(bens_eddy(i,:,:,j,2))*10),[3,1,2]);

            fig = subplotmaps(benstoplot,lontoplot,lat,{'div','RdBu'},1,[],18,titles,'Longitude {\circ}N','Latitude {\circ}N','T/decade','off',...
                contourlims,22,0:20:360,0:20:360,...
                -90:10:-30,-90:10:30,[num2str(trendyears(1)),'-',num2str(trendyears(2)),', ', mon{j},' - ',nameprint2{i}],1,[0 360],[-90 -30],0,'-',1);        
            %tight_subplot(2,1,[.01 .03],[.1 .01],[.01 .01]) 
            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/ZonalAsymmetry/TemperatureTrends/',mon{j},'_',num2str(trendyears(1)),'-',num2str(trendyears(2)),'_',nameprint_ens{i}];            
            export_fig(filename,'-png')
            close(fig);
        end
    end
end

%% line plots
if line_plots
    datestart = [1955,1960,1960,1955];
    data_ens(1).d = squeeze(nanmean(cat(5,data(1).d.T,data(2).d.T,data(3).d.T),5));
    data_ens(2).d = squeeze(nanmean(cat(5,data(4).d.T,data(5).d.T,data(6).d.T),5));
    data_ens(3).d = squeeze(nanmean(cat(5,data(7).d.T,data(8).d.T,data(9).d.T),5));
    data_ens(4).d = squeeze(nanmean(cat(5,data(10).d.T,data(11).d.T,data(12).d.T,data(13).d.T,data(14).d.T),5));
    
    ensmem = 4;
    nameprint_ensall = {'Ref-C2','Sen-C2-fGHG','Sen-C2-fODS','Ref-C1'};
    latlines = [-50,-60,-70,-90];
    lonlines = [0 120 240];

    cbrew = cbrewer('qual','Set1',12);

    count = 1;
    for i = 1:length(lonlines)
        for j = 1:length(latlines)
            legheader{count,1} = {[num2str(latlines(j)),'S, ',num2str(lonlines(i)),'E']};
            count = count+1;
        end
    end

    for i = 1:length(latlines)
        [~, latind(i)] = min(abs(lat - latlines(i)));
    end

    for i = 1:length(lonlines)
        [~, lonind(i)] = min(abs(lon - lonlines(i)));
    end

    linestoplot(1,:,:) = [data_ens(ensmem).d(lonind(1),latind(1),:);data_ens(ensmem).d(lonind(1),latind(2),:);data_ens(ensmem).d(lonind(1),latind(3),:);data_ens(ensmem).d(lonind(1),latind(4),:);...
        data_ens(ensmem).d(lonind(2),latind(1),:);data_ens(ensmem).d(lonind(2),latind(2),:);data_ens(ensmem).d(lonind(2),latind(3),:);data_ens(ensmem).d(lonind(2),latind(4),:);...
        data_ens(ensmem).d(lonind(3),latind(1),:);data_ens(ensmem).d(lonind(3),latind(2),:);data_ens(ensmem).d(lonind(3),latind(3),:);data_ens(ensmem).d(lonind(3),latind(4),:)];

    createfig('medium','on')
    for i = 1:12
        ph(i) = plot(squeeze(linestoplot(:,i,10:12:end))','LineWidth',2,'color',cbrew(i,:));
        hold on
    end
    xlabel('year','fontsize',20);
    ylabel('Kelvin','fontsize',20);
    set(gca,'fontsize',18,'xtick',1:20:200,'xticklabel',datestart(ensmem):20:2150);
    title([nameprint_ensall{ensmem},' - October temperature time series'],'fontsize',22);

    lh = legend(ph,{'-50S, 0E';'-60S, 0E';'-70S, 0E';'-90S, 0E';'-50S, 120E';'-60S, 120E';'-70S, 120E';'-90S, 120E';'-50S, 240E';'-60S, 240E';'-70S, 240E';'-90S, 240E'},'location','EastOutside');
    set(lh,'fontsize',18);

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/ZonalAsymmetry/TemperatureTrends/LinePlots_',mon{10},'_',nameprint_ensall{ensmem},num2str(trendyears(1)),'-',num2str(trendyears(2))];            
    export_fig(filename,'-png')
    clearvars linestoplot
end
