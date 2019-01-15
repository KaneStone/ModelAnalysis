clear all

lats = [50 60; 60 70];
compareyears = [2000,2015];
highclcompareyears = [1998,2024];

%% import SWOOSH
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201611-latpress-10deg-L31.nc');
% SWOOSHlatindex = find(SWOOSH.lat > lats(1) & SWOOSH.lat < lats(2));
% SWOOSH.combinedo3q = SWOOSH.combinedo3q(:,:,1:end-11);
% SWOOSHyears = 1984:2015;
% SWOOSHtimeindex = find(SWOOSHyears >= compareyears(1) & SWOOSHyears <= compareyears(2));

%% Read in highcl runs
directory = '/Volumes/ExternalOne/work/data/predruns/O3/highCl/zonalmean/';
files = dir([directory,'*.nc']);
for k = 1:size(lats,2)
    
    for i = 1:length(files)

        [~,highcldata(i),~] = Read_in_netcdf([directory,files(i).name]);      
        
        if i == 1
            latindex = find(highcldata(i).lat >= lats(1,k) & highcldata(i).lat <= lats(2,k));
        end
        
        highcldataPressure(i).p =  permute(repmat(highcldata(i).hyam*100000,[1,size(highcldata(i).PS)]),[2,1,3]) + ...
            permute(repmat(highcldata(i).hybm,[1,size(highcldata(i).PS)]),[2,1,3]) .* ...
            double(permute(repmat(highcldata(i).PS,[1,1,length(highcldata(i).lev)]),[1,3,2]));       
        if i == 1
            highclyears = 1995:2024;
        end
        for j = 1:size(highcldata(i).O3,1)
            [higcl_dataRegPres(i).h(:,j,:),~] = intRegPres(squeeze(highcldata(i).O3(j,:,:)),...
                squeeze(highcldataPressure(i).p(j,:,:))./100,SWOOSH.level);
        end
        for j = 1:12
            higcl_ra(i).h(j,:,:,:) = higcl_dataRegPres(i).h(:,:,j:12:end);
        end

        hcldateindex = find(highclyears >= highclcompareyears(1) & highclyears <= highclcompareyears(2));    
        higcl_ra_te(i).h(:,:,:,k) = permute(circshift(squeeze(nanmean(higcl_ra(i).h(:,:,latindex,hcldateindex),3)),[7,0,0]),[3,1,2]);
    end
    
end
higcl_ra_te(i+1).h = nanmean(cat(5,higcl_ra_te(:).h),5);

%% extracting locations

lev = 2;
monind = 12;
monind2 = monind - 5;
if monind > 12
    monind = monind -12;
end
[~,levind] = min(abs(SWOOSH.level - lev));

highcllines = cat(5,higcl_ra_te(:).h);
highcllines = squeeze(highcllines(:,monind2,levind,:,:))*1e6;

%%

for i = 1:size(highcllines,3)
    
    means = squeeze(nanmean(highcllines,1));
    stds = squeeze(std(highcllines,1));
    
end

meansall = repmat(reshape(means,[1,size(means)]),[size(highcllines,1),1,1]);
stdsall = repmat(reshape(stds,[1,size(stds)]),[size(highcllines,1),1,1]);

stdsall = circshift(stdsall,[0,1,0]);

anomaly = (highcllines - meansall)./stdsall;

anomaly2(:,1,:) = anomaly(:,1,:) - anomaly(:,2,:);% + meansall(:,1,:);
%anomaly2(:,2,:) = anomaly(:,1,:) - anomaly(:,2,:);% + meansall(:,2,:);

plot(squeeze(anomaly2(:,:,10)));


