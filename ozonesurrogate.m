% Attempting to use ozone as a surrogate for itself to look at dynamical influences

%% Read in data

[~,data,~] = Read_in_netcdf('/Volumes/MyBook/work/data/CESM-CCMI/O3/O3_f.e11.FWTREFC1.f19_f19.ccmi30.002.nc');

%% calculate pressure
Pressure = permute(repmat(data.ap.*100000,1,size(data.PS)),[2,3,1,4]) + ...
        permute(repmat(data.b,1,size(data.PS)),[2,3,1,4]) ...
        .* permute(repmat(data.PS,1,1,1,length(data.lev)),[1,2,4,3]); 
Pressure_zonmean = squeeze(nanmean(Pressure,1));
%%
preslev = 50;
mon = 7;
lat = -30; 
[~,presind] = min(abs(preslev-data.lev));


data_at_pres = squeeze(data.O3(:,1:48,presind,mon:12:end));
%data_at_lat = squeeze(data.O3(:,latind,:,mon:12:end));
data_at_pres_zonmean = squeeze(nanmean(data_at_pres));
%data_at_lat_zonmean = squeeze(nanmean(data_at_lat));

data_zonmean = squeeze(nanmean(data.O3(:,1:48,:,mon:12:end)));

%% normalized anomalies
for j = 1:size(data_zonmean,2)
    for i = 1:size(data_zonmean,1)
        toplot(i,j,:) = (squeeze(data_zonmean(i,j,:)) - squeeze(nanmean(data_zonmean(i,j,:))))./squeeze(std(data_zonmean(i,j,:)));    
        if i >= 2
            toplotdiff(i-1,j,:) = toplot(i,j,:) - toplot(i-1,j,:);       
        end
    end
end
% for i = 1:size(data_at_lat_zonmean,1)
%     toplot_lat(i,:) = (data_at_lat_zonmean(i,:) - nanmean(data_at_lat_zonmean(i,:)))./std(data_at_lat_zonmean(i,:));
%     if i >= 2
%         toplotdiff_lat(i-1,:) = toplot_lat(i,:) - toplot_lat(i-1,:);
%     end
%     
% end

%% plotting
lat = -20;
[~,latind] = min(abs(lat-data.lat));
% figure;
% contourf(squeeze(toplotdiff(:,presind,:)))
% caxis([-.5 .5]);
% colorbar

 
 for i = size(data_zonmean,1)-1:-1:1
    
    
prestoplot = squeeze(Pressure_zonmean(i,:,:));
prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];
    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],5,[],[],2,1};
    logprestick = log(prestick);

%monthnames = {'January','February','March','April','May','June','July','August','September',...
%    'October','November','December'};

titles = {'placeholder'};    

plotmtitle = 'placeholder';
ctitle = 'placeholder'

subplotmaps(squeeze(toplotdiff(i,:,:)),1955:1:2015,log(prestoplot),{'div','RdBu'},1,[],16,titles,'Year','Pressure (hPa)',ctitle,'on',...
    [-.5 .5],14,1:5,60,1955:5:2015,...
    fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(.1) log(1000)],1,'-',0);

%contourf(1:60,log(data.lev),squeeze(toplotdiff(i,:,:)));
%caxis([-.5 .5]);
%colorbar
%set(gca,'YDir','reverse');
%set(gca,'ytick',log(data.lev),'yticklabel',data.lev)
export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/Ozonesurrogates/',num2str(data.lat(i)),'_','prestime.png'],'-png');
close(fig);
 end
 
 
 