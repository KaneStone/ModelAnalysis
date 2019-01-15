% Attempting to use ozone as a surrogate for itself to look at dynamical influences
clear all
%% Read in data
if ~exist('/Volumes/ExternalOne/work/data/WACCM/output/REFC1_01.mat','file')
    [~,data,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/CESM-CCMI/O3/O3_f.e11.FWTREFC1.f19_f19.ccmi30.002.nc');

    %% calculate pressure
    Pressure = permute(repmat(data.ap.*100000,1,size(data.PS)),[2,3,1,4]) + ...
            permute(repmat(data.b,1,size(data.PS)),[2,3,1,4]) ...
            .* permute(repmat(data.PS,1,1,1,length(data.lev)),[1,2,4,3]); 
    Pressure_zonmean = squeeze(nanmean(Pressure,1));

    %% Interpolate onto regular pressure levels
 regpres.O3 = zeros(size(data.O3,1),size(data.O3,2),52,size(data.O3,4));
    for i = 1:size(data.O3,1)
        i
        for j = 1:size(data.O3,2)
            [regpres.O3(i,j,:,:),regpres.pres] = intRegPres(squeeze(data.O3(i,j,:,:)),squeeze(Pressure(i,j,:,:)));
        end
    end
regpres.lat = data.lat;
regpres.lon = data.lon;
regpres.date = data.date;
    save('/Volumes/ExternalOne/work/data/WACCM/output/REFC1_01.mat','regpres','-v7.3');
else
    load('/Volumes/ExternalOne/work/data/WACCM/output/REFC1_01.mat');
end

%%
preslev = 100;
mon = 7;
lat = -30; 
[~,presind] = min(abs(regpres.pres-preslev));

for i = 1:12
    data_at_pres(:,:,:,i) = squeeze(regpres.O3(:,1:48,presind,i:12:end));
end
%data_at_lat = squeeze(data.O3(:,latind,:,mon:12:end));
% data_at_pres_zonmean = squeeze(nanmean(data_at_pres));
% %data_at_lat_zonmean = squeeze(nanmean(data_at_lat));
% 
% data_zonmean = squeeze(nanmean(data.O3(:,1:48,:,mon:12:end)));

%% normalized anomalies
clearvars toplotdiff toplot
for j = 1:size(data_at_pres,1) % longitude
    for i = 1:size(data_at_pres,2) % latitude
        for k = 1:12
            toplot(j,i,:,k) = (squeeze(data_at_pres(j,i,:,k)) - squeeze(nanmean(data_at_pres(j,i,:,k))))./squeeze(std(data_at_pres(j,i,:,k)));    
            if i >= 2
                toplotdiff(j,i-1,:,k) = toplot(j,i,:,k) - toplot(j,i-1,:,k);       
                r(j,i-1,k) = corr(squeeze(toplot(j,i,:,k)),squeeze(toplot(j,i-1,:,k)));       
            end
        end
    end
end

%% correlations of longitude lines 
for j = 1:size(data_at_pres,3) % year
    for i = 1:size(data_at_pres,2) % latitude
        for k = 1:12
            toplotlonlines(j,i,:,k) = (squeeze(data_at_pres(:,i,j,k)) - squeeze(nanmean(data_at_pres(:,i,j,k))))./squeeze(std(data_at_pres(:,i,j,k)));    
            if i >= 2
                toplotdifflonlines(j,i-1,:,k) = toplotlonlines(j,i,:,k) - toplotlonlines(j,i-1,:,k);       
                rlonlines(j,i-1,k) = corr(squeeze(toplotlonlines(j,i,:,k)),squeeze(toplotlonlines(j,i-1,:,k)));       
            end
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
%%
contourf(1955:2014,regpres.lat(1:47),rlonlines(:,:,12)',20,'LineStyle','none')
caxis([.5 1])
colorbar
%%
close all
for i = 1:60
    figure;
    subplot(2,1,1)
    contourf(regpres.lon,regpres.lat(1:47),toplotdiff(:,:,i,7)');
    caxis([-.5 .5]);
    subplot(2,1,2)
    contourf(regpres.lon,regpres.lat(1:47),r(:,:,i,7)');
    colorbar
end

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
ctitle = 'placeholder';

subplotmaps(toplotdiff(i,:,:),1955:2014,log(prestoplot),{'div','RdBu'},1,[],16,titles,'Year','Pressure (hPa)',ctitle,'on',...
    [-.5 .5],18,1:5:60,1955:5:2015,fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(.1) log(1000)],1,'-',0,'none');

%subplotmaps(differences2,data.LENS.longitude,log(LENSverens.lev),{'div','RdBu'},1,h,18,contourtitle2,'Longitude','Pressure (hPa)','Kelvin','on',...
%    [-2,2],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),contourtitle,1,[0 360],[log(1) log(1000)],1,'-',0,'none');

%contourf(1:60,log(data.lev),squeeze(toplotdiff(i,:,:)));
%caxis([-.5 .5]);
%colorbar
%set(gca,'YDir','reverse');
%set(gca,'ytick',log(data.lev),'yticklabel',data.lev)
export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/Ozonesurrogates/',num2str(data.lat(i)),'_','prestime.png'],'-png');
close(fig);
 end
 
 
 
