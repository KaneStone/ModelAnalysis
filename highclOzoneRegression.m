% highcl runs regression
clear all

timperiod = [2000 2024];

%% import SWOOSH
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201611-latpress-10deg-L31.nc');

%% Read in highcl runs
directory = '/Volumes/MyBook/work/data/predruns/O3/highCl/zonalmean/';
files = dir([directory,'*.nc']);
for i = 1:length(files)
    [~,highcldata(i),~] = Read_in_netcdf([directory,files(i).name]);
    
    if i == 1
        highclyears = 1995:2024;
        highclyears = repmat(highclyears,[12,1]);
        highclyears = highclyears(:);        
    end
    highcldata(i).O3 = highcldata(i).O3(:,:,highclyears >= timperiod(1) & highclyears <= timperiod(2));
    highcldata(i).PS = highcldata(i).PS(:,highclyears >= timperiod(1) & highclyears <= timperiod(2));
    
    highcldataPressure(i).p =  permute(repmat(highcldata(i).hyam*100000,[1,size(highcldata(i).PS)]),[2,1,3]) + ...
        permute(repmat(highcldata(i).hybm,[1,size(highcldata(i).PS)]),[2,1,3]) .* ...
        double(permute(repmat(highcldata(i).PS,[1,1,length(highcldata(i).lev)]),[1,3,2]));           
    for j = 1:size(highcldata(i).O3,1)
        [higcl_dataRegPres(i).h(:,j,:),~] = intRegPres(squeeze(highcldata(i).O3(j,:,:)),...
            squeeze(highcldataPressure(i).p(j,:,:))./100,SWOOSH.level);
    end
end

%% Read in Uwind
directory = '/Volumes/MyBook/work/data/predruns/U/highCl/';
filesU = dir([directory,'*.nc']);
for i = 1:length(filesU)
    [~,highclU(i),~] = Read_in_netcdf([directory,filesU(i).name]);
    highclUzonalmean(i).U = squeeze(nanmean(highclU(i).U(:,:,[2,7],:),1));
    highclUzonalmean(i).U = highclUzonalmean(i).U(:,highclyears >= timperiod(1) & highclyears <= timperiod(2));
end

%% taking ensemble averages
higcl_dataRegPres(10).h = nanmean(cat(4,higcl_dataRegPres(:).h),4);
highclUzonalmean(10).U = nanmean(cat(3,highclUzonalmean(:).U),3);

%% remove QBO from each month indipendently
for l = 1:10
    QBOyearly(:,l) = (highclUzonalmean(l).U(2,:) - nanmean(highclUzonalmean(l).U(2,:),2))./std(highclUzonalmean(l).U(2,:),0,2);
    temp = null(QBOyearly(:,l).','r');
    QBOyearly2(:,l) = [(temp(1,:)' - nanmean(temp(1,:),2))./std(temp(1,:),0,2);-QBOyearly(end,l)];
    
    for j = 1:size(higcl_dataRegPres(1).h,1)
        for k = 1:size(higcl_dataRegPres(1).h,2)            
            for i = 1:12
                O3monthmean(j,k,l,i) = squeeze(nanmean(higcl_dataRegPres(l).h(j,k,i:12:end),3));
                O3anomtemp(j,k,l,i,:) = (squeeze(higcl_dataRegPres(l).h(j,k,i:12:end)) - ...
                    O3monthmean(j,k,l,i))./O3monthmean(j,k,l,i)*100;
            end
            numyears = size(O3anomtemp,5);
            O3anomyearly(:,j,k,l) = squeeze(O3anomtemp(j,k,l,:));
            
            % Predictors
            onesvectors = ones(length(QBOyearly(:,l)),1);            
            linearvector = (1:length(QBOyearly(:,l)))';
            QBOsin(:,l) = QBOyearly(:,l).*sin(2.*pi.*(1:numyears*12)./12)';
            QBOcos(:,l) = QBOyearly(:,l).*cos(2.*pi.*(1:numyears*12)./12)';                        
            QBOsin2(:,l) = QBOyearly2(:,l).*sin(2.*pi.*(1:numyears*12)./12)';            
            QBOcos2(:,l) = QBOyearly2(:,l).*cos(2.*pi.*(1:numyears*12)./12)';                        
            byearly(j,k,l,:) = regress(O3anomyearly(:,j,k,l),[onesvectors,QBOyearly(:,l),QBOsin(:,l),...
                QBOcos(:,l),QBOyearly2(:,l),QBOsin2(:,l),QBOcos2(:,l),linearvector]);
            
            %testing
%             if j == 26
%                 figure
%                 plot(O3anomyearly(:,j,k,l))
%                 hold on
%                 plot(byearly(1)*onesvectors + QBOyearly(:,l)*byearly(j,k,l,2) + QBOsin(:,l)*byearly(j,k,l,3) + ...
%                     QBOcos(:,l)*byearly(j,k,l,4) + QBOyearly2(:,l)*byearly(j,k,l,5) + QBOsin2(:,l)*byearly(j,k,l,6) + ...
%                     QBOcos2(:,l)*byearly(j,k,l,7) + linearvector*byearly(j,k,l,8));
%             end
            
        end
    end
    
    for i = 1:12
        QBO(:,l,i) = (highclUzonalmean(l).U(2,i:12:end) - nanmean(highclUzonalmean(l).U(2,i:12:end),2))./std(highclUzonalmean(l).U(2,i:12:end),0,2);
        temp = null(QBO(:,l,i).','r');
        QBO2(:,l,i) = [(temp(1,:)' - nanmean(temp(1,:),2))./std(temp(1,:),0,2);-QBO(end,l,i)];
        
        for j = 1:size(higcl_dataRegPres(1).h,1)
            for k = 1:size(higcl_dataRegPres(1).h,2)            
                O3anom(:,i,j,k,l) = squeeze(higcl_dataRegPres(l).h(j,k,i:12:end) - ...
                    nanmean(higcl_dataRegPres(l).h(j,k,i:12:end),3));
%                 b(i,j,k,l,:) = regress(O3anom(:,i,j,k,l),[ones(length(QBO(:,l,i)),1),...
%                     QBO(:,l,i),QBO2(:,l,i)]);%,(1:length(QBO2(:,l,i)))'
                b(i,j,k,l,:) = regress(O3anom(:,i,j,k,l),[ones(length(QBO(:,l,i)),1),...
                    QBO(:,l,i)]);%,QBO2(:,l,i)]);%,(1:length(QBO2(:,l,i)))'
                O3anomQBOrm(:,i,j,k,l) = O3anom(:,i,j,k,l) - ...
                    squeeze(b(i,j,k,l,2))*QBO(:,l,i);% - ...
                    %squeeze(b(i,j,k,l,3))*QBO2(:,l,i);                
                blinearQBOrm(i,j,k,l,:) = regress(O3anomQBOrm(:,i,j,k,l),[ones(length(QBO(:,l,i)),1),(1:length(QBO(:,l,i)))']);
                blinear(i,j,k,l,:) = regress(O3anom(:,i,j,k,l),[ones(length(QBO(:,l,i)),1),(1:length(QBO(:,l,i)))']);
            end
        end
    end
end


%%
createfig('medium','on');
mon = [1,8,2,9];
lev = [15,26,26,27];
lat = [49,18,96-13,96-13];
for i = 1:4
    subplot(2,2,i)
    %plot((O3anom(:,mon(i),lev,lat(i),10) - squeeze(b(mon(i),lev,lat(i),10,2))*QBO(:,10,mon(i)) - squeeze(b(mon(i),lev,lat(i),10,3))*QBO2(:,10,mon(i)))*1e6);
    hold on
    plot(nanmean(O3anomQBOrm(:,mon(i),lev(i),lat(i),1:9),5)*1e6);
    plot(O3anom(:,mon(i),lev(i),lat(i),10)*1e6);
end
%plot(nanmean(O3anom(:,12,27,80,:),5))

%% 
plot(O3anomyearly(:,lev(2),96,10));

hold on
plot(byearly(lev(2),96,10,8)*(1:numyears*12))


%% plot

prestick = [300,200,100,90:-10:10,9:-1:1];
presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
    5,[],3,2,1};
logprestick = log(prestick);

titles = {'Ensemble average 2000-2014 yearly ozone trends','2'};

%bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
% bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
%bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
subplotmaps(bforplot,highcldata(1).lat,log(SWOOSH.level),{'div','RdBu'},1,[],22,titles,'Latitude','Pressure (hPa)','ppmv/decade','on',...
    [-5 5],22,-90:10:90,-90:10:90,...
    fliplr(logprestick),fliplr(presticklabel),{''} ,1,[-90 90],[log(1.1) log(300)],1,'none',0,'');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','Highcl_ensave_LatPres_QBOremoved'];

export_fig(filename,'-pdf');

%% plot
latind = 83;
bforplot2 = permute(reshape(squeeze(b(:,:,latind,10,2))',[1,31,12]),[1,3,2])*1e6*10;

latfortit = num2str(highcldata(1).lat(latind));

titles = {['Ensemble average 2000-2014 ozone trends at ',latfortit,'N']};

subplotmaps(bforplot2,1:12,log(SWOOSH.level),{'div','RdBu'},1,[],22,titles,'Month','Pressure (hPa)','ppmv/decade','on',...
    [-1 1],22,1:12,{'J','F','M','A','M','J','J','A','S','O','N','D'},...
    fliplr(logprestick),fliplr(presticklabel),{''} ,1,[1 12],[log(1.1) log(200)],1,'none',0,'');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','Highcl_ensave_MonthPres_QBOremoved_',latfortit];

export_fig(filename,'-pdf');
%%
createfig('medium','on');
subplot(2,1,1);
contourf(highcldata(1).lat,log(SWOOSH.level),squeeze(nanmean(b(:,:,:,10,4),1)));
set(gca,'Ydir','reverse','ytick',flipud(log(SWOOSH.level)),'yticklabel',flipud(SWOOSH.level));
colorbar
caxis([-2e-8 2e-8])

% subplot(3,1,2);
% contourf(highcldata(1).lat,log(SWOOSH.level),squeeze(b(2,:,:,10,4)));
% set(gca,'Ydir','reverse','ytick',flipud(log(SWOOSH.level)),'yticklabel',flipud(SWOOSH.level));
% colorbar
% caxis([-2e-8 2e-8])

subplot(2,1,2);
contourf(1:12,log(SWOOSH.level),squeeze(b(:,:,13,10,4))');
set(gca,'Ydir','reverse','ytick',flipud(log(SWOOSH.level)),'yticklabel',flipud(SWOOSH.level));
colorbar
caxis([-2e-8 2e-8])