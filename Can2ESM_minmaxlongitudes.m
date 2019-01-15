%% Can2ESM minimum and maximum locations

% Read in Can2ESM Temperature data
clear all
preslev = '50hPa';
directory = ['/Volumes/ExternalOne/work/data/CanESM2/',preslev,'/3090S/'];
files = dir([directory,'*.nc']);
read_in = 1;
plot_ens = 0;
plot_all = 0;
plot_ensall = 1;
line_plots = 0;
timeperiod = [1960,1980];
sims = {'historical-Natural','Historical','Anthropogenic aerosols','Stratospheric ozone'};
structnames = {'Natural','Historical','Aerosols','StratOzone'};
structnames2 = {'N1','N2','N3','N4','N5','N6','N7','N8','N9','N10','N11','N12','N13','N14','N15','N16','N17','N18','N19','N20',...
    'N21','N22','N23','N24','N25','N26','N27','N28','N29','N30','N31','N32','N33','N34','N35','N36','N37','N38','N39','N40',...
    'N41','N42','N43','N44','N45','N46','N47','N48','N49','N50',...
    'H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12','H13','H14','H15','H16','H17','H18','H19','H20',...
    'H21','H22','H23','H24','H25','H26','H27','H28','H29','H30','H31','H32','H33','H34','H35','H36','H37','H38','H39','H40',...
    'H41','H42','H43','H44','H45','H46','H47','H48','H49','H50',...
    'A1','S1','A2','S2','A3','S3','A4','S4','A5','S5','A6','S6','A7','S7','A8','S8','A9','S9','A10','S10',...
    'A11','S11','A12','S12','A13','S13','A14','S14','A15','S15','A16','S16','A17','S17','A18','S18','A19','S19','A20','S20',...
    'A21','S21','A22','S22','A23','S23','A24','S24','A25','S25','A26','S26','A27','S27','A28','S28','A29','S29','A30','S30',...
    'A31','S31','A32','S32','A33','S33','A34','S34','A35','S35','A36','S36','A37','S37','A38','S38','A39','S39','A40','S40',...
    'A41','S41','A42','S42','A43','S43','A44','S44','A45','S45','A46','S46','A47','S47','A48','S48','A49','S49','A50','S50'};
if read_in
    [ta,years,ta_ens,ensyears,ta_ensall,ensallyears,lon,lat] = Read_in_Can2ESMfortrends(directory,files,'ta');
end

%% restructurig temperature ensemble
for i = 1:size(ta_ensall,1)
    ta_ensalltogether(i,:,:,:) = permute(cat(3,ta_ensall(i,:).w),[3,1,2]);
end

for i = 1:size(ta,1)
    ta_alltogether(i,:,:,:) = permute(cat(3,ta(i,:).w),[3,1,2]);
end

%% plotting longitude lines
years2 = ensallyears(1).y;
dates = [1995,2024];
dateind = find(years2 >= dates(1) & years2 <= dates(2));
dates2 = [1955,1979];
dateind2 = find(years2 >= dates2(1) & years2 <= dates2(2));

dates3 = [2050,2079];
dateind3 = find(years2 >= dates3(1) & years2 <= dates3(2));

fsize = 20;
createfig('medium','on');
cbrew2 = cbrewer('qual','Paired',10);
lonlines = permute(nanmean(cat(4,squeeze(cat(4,ta_alltogether(51:100,:,9,dateind(1)-1+9:12:dateind(length(dateind))))),...
squeeze(cat(4,ta_alltogether(51:100,:,9,dateind(1)-1+10:12:dateind(length(dateind))))),...
squeeze(cat(4,ta_alltogether(51:100,:,9,dateind(1)-1+11:12:dateind(length(dateind)))))),4),[2,3,1]);




lonlines2 = permute(nanmean(cat(4,squeeze(cat(4,ta_alltogether(51:100,:,9,dateind2(1)-1+9:12:dateind2(length(dateind2))))),...
    squeeze(cat(4,ta_alltogether(51:100,:,9,dateind2(1)-1+10:12:dateind2(length(dateind2))))),...
    squeeze(cat(4,ta_alltogether(51:100,:,9,dateind2(1)-1+11:12:dateind2(length(dateind2)))))),4),[2,3,1]);

lonlines3 = permute(nanmean(cat(4,squeeze(cat(4,ta_alltogether(51:100,:,9,dateind3(1)-1+9:12:dateind3(length(dateind3))))),...
    squeeze(cat(4,ta_alltogether(51:100,:,9,dateind3(1)-1+10:12:dateind3(length(dateind3))))),...
    squeeze(cat(4,ta_alltogether(51:100,:,9,dateind3(1)-1+11:12:dateind3(length(dateind3)))))),4),[2,3,1]);

lonlinesstrat = permute(nanmean(cat(4,squeeze(cat(4,ta_alltogether(102:2:end,:,9,dateind(1)-1+9:12:dateind(length(dateind))))),...
    squeeze(cat(4,ta_alltogether(102:2:end,:,9,dateind(1)-1+10:12:dateind(length(dateind))))),...
    squeeze(cat(4,ta_alltogether(102:2:end,:,9,dateind(1)-1+11:12:dateind(length(dateind)))))),4),[2,3,1]);

lonlinesstrat2 = permute(nanmean(cat(4,squeeze(cat(4,ta_alltogether(102:2:end,:,9,dateind2(1)-1+9:12:dateind2(length(dateind2))))),...
    squeeze(cat(4,ta_alltogether(102:2:end,:,9,dateind2(1)-1+10:12:dateind2(length(dateind2))))),...
    squeeze(cat(4,ta_alltogether(102:2:end,:,9,dateind2(1)-1+11:12:dateind2(length(dateind2)))))),4),[2,3,1]);

lonlinesstrat3 = permute(nanmean(cat(4,squeeze(cat(4,ta_alltogether(102:2:end,:,9,dateind3(1)-1+9:12:dateind3(length(dateind3))))),...
    squeeze(cat(4,ta_alltogether(102:2:end,:,9,dateind3(1)-1+10:12:dateind3(length(dateind3))))),...
    squeeze(cat(4,ta_alltogether(102:2:end,:,9,dateind3(1)-1+11:12:dateind3(length(dateind3)))))),4),[2,3,1]);

test = reshape(lonlines,size(lonlines,1),size(lonlines,2)*size(lonlines,3));
test2 = reshape(lonlines2,size(lonlines2,1),size(lonlines2,2)*size(lonlines2,3));
test3 = reshape(lonlines3,size(lonlines3,1),size(lonlines3,2)*size(lonlines3,3));

amplitude = max(test)-min(test);
ampstd = std(amplitude);
ampmean = nanmean(amplitude);

plot(lon,test,'color',cbrew2(1,:))
hold on
plot(lon,test3,'color',cbrew2(3,:))
plot(lon,test2,'color',cbrew2(5,:))

ph(1) = plot(lon,nanmean(test,2),'LineWidth',5,'color',cbrew2(2,:));
ph(2) = plot(lon,nanmean(test2,2),'--','LineWidth',5,'color',cbrew2(6,:));
ph(3) = plot(lon,nanmean(test3,2),':','LineWidth',5,'color',cbrew2(4,:));
amplitude2 = max(test2)-min(test2);
ampstd2 = std(amplitude2);
ampmean2 = nanmean(amplitude2);

amplitude3 = max(test3)-min(test3);
ampstd3 = std(amplitude3);
ampmean3 = nanmean(amplitude3);

set(gca,'fontsize',fsize+2);
xlim([-5 365]);
xlabel('Longitude ({\circ}E)','fontsize',fsize+4);
ylabel('Temperature (K)','fontsize',fsize+4);
title('Austral spring temperatures at 65{\circ}S and 50hPa','fontsize',fsize+6);
lh = legend(ph,'1995-2025 (high chlorine)','1955-1979 (low chlorine)','2050-2079 (high GHGs - low chlorine)');
set(lh,'fontsize',fsize+4,'box','off','location','south')

filename = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TalkFigures/CanESM2Amplitudes_65S';
export_fig(filename,'-png');

%% strat
teststrat = reshape(lonlinesstrat,size(lonlinesstrat,1),size(lonlinesstrat,2)*size(lonlinesstrat,3));
teststrat2 = reshape(lonlinesstrat2,size(lonlinesstrat2,1),size(lonlinesstrat2,2)*size(lonlinesstrat2,3));
teststrat3 = reshape(lonlinesstrat3,size(lonlinesstrat3,1),size(lonlinesstrat3,2)*size(lonlinesstrat3,3));

plot(lon,teststrat,'color',cbrew2(1,:))
hold on
plot(lon,teststrat3,'color',cbrew2(3,:))
plot(lon,teststrat2,'color',cbrew2(5,:))

ph(1) = plot(lon,nanmean(teststrat,2),'LineWidth',5,'color',cbrew2(2,:));
ph(2) = plot(lon,nanmean(teststrat2,2),'--','LineWidth',5,'color',cbrew2(6,:));
ph(3) = plot(lon,nanmean(teststrat3,2),':','LineWidth',5,'color',cbrew2(4,:));
amplitude2 = max(teststrat2)-min(teststrat2);
ampstd2 = std(amplitude2);
ampmean2 = nanmean(amplitude2);

amplitude3 = max(teststrat3)-min(teststrat3);
ampstd3 = std(amplitude3);
ampmean3 = nanmean(amplitude3);

set(gca,'fontsize',fsize+2);
xlim([-5 365]);
xlabel('Longitude ({\circ}E)','fontsize',fsize+4);
ylabel('Temperature (K)','fontsize',fsize+4);
title('Austral spring temperatures at 65{\circ}S and 50 hPa','fontsize',fsize+6);
lh = legend(ph,'1995-2025 (high chlorine)','1955-1979 (low chlorine)','2050-2079 (low chlorine)');
set(lh,'fontsize',fsize+4,'box','off','location','south')

filename = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TalkFigures/CanESM2stratAmplitudes_65S.png';
export_fig(filename,'r200','-png');


%% finding maximum and minimum for all members
lats = [-50,-55,-60,-65,-70,-75,-80,-85];
for i = 1:size(ta_alltogether)
    for k = 1:length(lats)
        for j = 1:12        
            [~,latind] = min(abs(lats(k)-lat));   
            [tempminvalue,tempminind] = min(squeeze(ta_alltogether(i,:,latind,j:12:end)));        
            [tempmaxvalue,tempmaxind] = max(squeeze(ta_alltogether(i,:,latind,j:12:end)));       
            Can2ESM.minvalue_m(i).m(j,k,:) = tempminvalue;
            Can2ESM.minind_m(i).m(j,k,:) = tempminind;
            Can2ESM.maxvalue_m(i).m(j,k,:) = tempmaxvalue;
            Can2ESM.maxind_m(i).m(j,k,:) = tempmaxind;
        end
        
        [tempminvalue,tempminind] = min(nanmean(cat(3,squeeze(ta_alltogether(i,:,latind,9:12:end)),...
            squeeze(ta_alltogether(i,:,latind,10:12:end)),...
            squeeze(ta_alltogether(i,:,latind,11:12:end))),3),[],1);        
        [tempmaxvalue,tempmaxind] = max(nanmean(cat(3,squeeze(ta_alltogether(i,:,latind,9:12:end)),...
            squeeze(ta_alltogether(i,:,latind,10:12:end)),...
            squeeze(ta_alltogether(i,:,latind,11:12:end))),3),[],1);        
        Can2ESM.minvalue_m_spr(i).m(k,:) = tempminvalue;
        Can2ESM.minind_m_spr(i).m(k,:) = tempminind;                                    
        Can2ESM.maxvalue_m_spr(i).m(k,:) = tempmaxvalue;
        Can2ESM.maxind_m_spr(i).m(k,:) = tempmaxind;
        
        
    end
    Can2ESM.minlongitude_m(i).m = lon(Can2ESM.minind_m(i).m);
    Can2ESM.maxlongitude_m(i).m = lon(Can2ESM.maxind_m(i).m);
    Can2ESM.amplitude_m(i).m = Can2ESM.maxvalue_m(i).m - Can2ESM.minvalue_m(i).m;
    
    Can2ESM.minlongitude_m_spr(i).m = lon(Can2ESM.minind_m_spr(i).m);
    Can2ESM.maxlongitude_m_spr(i).m = lon(Can2ESM.maxind_m_spr(i).m);
    Can2ESM.amplitude_m_spr(i).m = Can2ESM.maxvalue_m_spr(i).m - Can2ESM.minvalue_m_spr(i).m;
end

for i = 1:size(ta_alltogether,1)    
    for k = 1:length(lats)   
        for j = 1:12
            for l = 1:size(Can2ESM.maxlongitude_m(i).m(j,k,:),3)
                if Can2ESM.maxlongitude_m(i).m(j,k,l) > 300
                    Can2ESM.maxlongitude_m(i).m(j,k,l) = Can2ESM.maxlongitude_m(i).m(j,k,l)-360;
                end
                if Can2ESM.minlongitude_m(i).m(j,k,l) < 100
                    Can2ESM.minlongitude_m(i).m(j,k,l) = Can2ESM.minlongitude_m(i).m(j,k,l)+360;
                end                                                
            end
        end
        for l = 1:size(Can2ESM.maxlongitude_m_spr(i).m(k,:),2)
            if Can2ESM.maxlongitude_m_spr(i).m(k,l) > 300
                Can2ESM.maxlongitude_m_spr(i).m(k,l) = Can2ESM.maxlongitude_m_spr(i).m(k,l)-360;
            end
            if Can2ESM.minlongitude_m_spr(i).m(k,l) < 100
                Can2ESM.minlongitude_m_spr(i).m(k,l) = Can2ESM.minlongitude_m_spr(i).m(k,l)+360;                           
            end
        end
    end
end

%%
ens = [1:50;51:100;101:2:200;102:2:200];
for i = 1:length(structnames)
    Can2ESM.maxlongitude2.(structnames{i}) = nanmean(cat(4,Can2ESM.maxlongitude_m(ens(i,:)).m),4);
    Can2ESM.minlongitude2.(structnames{i}) = nanmean(cat(4,Can2ESM.minlongitude_m(ens(i,:)).m),4);
    Can2ESM.maxlongitude_spr2.(structnames{i}) = nanmean(cat(3,Can2ESM.maxlongitude_m_spr(ens(i,:)).m),3);
    Can2ESM.minlongitude_spr2.(structnames{i}) = nanmean(cat(3,Can2ESM.minlongitude_m_spr(ens(i,:)).m),3);
    Can2ESM.amplitude2.(structnames{i}) = nanmean(cat(4,Can2ESM.amplitude_m(ens(i,:)).m),4);
    Can2ESM.amplitude_spr2.(structnames{i}) = nanmean(cat(3,Can2ESM.amplitude_m_spr(ens(i,:)).m),3);
end

Can2ESM.longitude = lon;
Can2ESM.latitude = lats;

%%
for i = 1:size(ta_ensalltogether,1)    
    for k = 1:length(lats)
        [~,latind] = min(abs(lats(k)-lat));            
        for j = 1:12                    
            [tempminvalue,tempminind] = min(squeeze(ta_ensalltogether(i,:,latind,j:12:end)));        
            [tempmaxvalue,tempmaxind] = max(squeeze(ta_ensalltogether(i,:,latind,j:12:end)));                   
            Can2ESM.minvalue.(structnames{i})(j,k,:) = tempminvalue;
            Can2ESM.minind.(structnames{i})(j,k,:) = tempminind;                                    
            Can2ESM.maxvalue.(structnames{i})(j,k,:) = tempmaxvalue;
            Can2ESM.maxind.(structnames{i})(j,k,:) = tempmaxind;            
        end
        
        [tempminvalue,tempminind] = min(nanmean(cat(3,squeeze(ta_ensalltogether(i,:,latind,9:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,10:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,11:12:end))),3),[],1);        
        [tempmaxvalue,tempmaxind] = max(nanmean(cat(3,squeeze(ta_ensalltogether(i,:,latind,9:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,10:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,11:12:end))),3),[],1);               
        Can2ESM.minvalue_spr.(structnames{i})(k,:) = tempminvalue;
        Can2ESM.minind_spr.(structnames{i})(k,:) = tempminind;                                    
        Can2ESM.maxvalue_spr.(structnames{i})(k,:) = tempmaxvalue;
        Can2ESM.maxind_spr.(structnames{i})(k,:) = tempmaxind;
    end    
    
    Can2ESM.minlongitude.(structnames{i}) = lon(Can2ESM.minind.(structnames{i}));
    Can2ESM.maxlongitude.(structnames{i}) = lon(Can2ESM.maxind.(structnames{i}));
    Can2ESM.amplitude.(structnames{i}) = Can2ESM.maxvalue.(structnames{i}) - Can2ESM.minvalue.(structnames{i});
    
    Can2ESM.minlongitude_spr.(structnames{i}) = lon(Can2ESM.minind_spr.(structnames{i}));
    Can2ESM.maxlongitude_spr.(structnames{i}) = lon(Can2ESM.maxind_spr.(structnames{i}));
    Can2ESM.amplitude_spr.(structnames{i}) = Can2ESM.maxvalue_spr.(structnames{i}) - Can2ESM.minvalue_spr.(structnames{i});
    
end
for i = 1:size(ta_ensalltogether,1)    
    for k = 1:length(lats)     
        for j = 1:12
            for l = 1:size(Can2ESM.maxlongitude.(structnames{i})(j,k,:),3)
                if Can2ESM.maxlongitude.(structnames{i})(j,k,l) > 300
                    Can2ESM.maxlongitude.(structnames{i})(j,k,l) = Can2ESM.maxlongitude.(structnames{i})(j,k,l)-360;
                end
                if Can2ESM.minlongitude.(structnames{i})(j,k,l) < 100
                    Can2ESM.minlongitude.(structnames{i})(j,k,l) = Can2ESM.minlongitude.(structnames{i})(j,k,l)+360;
                end            
            end
        end

        for l = 1:size(Can2ESM.maxlongitude_spr.(structnames{i})(k,:),2)
            if Can2ESM.maxlongitude_spr.(structnames{i})(k,l) > 300
                Can2ESM.maxlongitude_spr.(structnames{i})(k,l) = Can2ESM.maxlongitude_spr.(structnames{i})(k,l)-360;
            end
            if Can2ESM.minlongitude_spr.(structnames{i})(k,l) < 100
                Can2ESM.minlongitude_spr.(structnames{i})(k,l) = Can2ESM.minlongitude_spr.(structnames{i})(k,l)+360;                           
            end
        end        
    
    end
end
Can2ESM.years = years;
Can2ESM.ensyears = ensallyears;
Can2ESM.ensnames = structnames;
Can2ESM.names = structnames2;
save('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/Can2ESM.mat','Can2ESM');

%% plotting amplitude
createfig('medium','on')
contourf(amplitude.H,14)
%caxis([0,22]);
colorbar
title('Historical');

%% plotting
fsize = 18;
cbrew = cbrewer('qual','Set1',10);
createfig('large','on');
%subplot(2,1,1)
m_proj('Stereographic','lon',-180,'lat',-90,'rad',50)
m_coast('color','k','LineWidth',3);
m_grid('ytick',[-90 -80 -70 -60 -50 -40 -30 -20],'XaxisLocation','top','fontsize',fsize);
for i = 1:length(structnames)
    mh(i) = m_line(squeeze(minlongitude.(structnames{i})(10,:))',lats','color',cbrew(i,:),'LineWidth',3);    
    m_line(squeeze(maxlongitude.(structnames{i})(10,:))',lats','color',cbrew(i,:),'LineWidth',3);
    
    %m_line(squeeze(minlong_mean.(fieldnames{i})(:,10)),lats','color',cbrew(i,:),'LineWidth',3);    
    %m_line(squeeze(maxlong_mean.(fieldnames{i})(:,10)),lats','color',cbrew(i,:),'LineWidth',3);
    
end
legend(mh,'HistoricalNat','Historical','Anthropogenic aerosols','Stratospheric ozone');
title([num2str(timeperiod(1)),'-',num2str(timeperiod(2))],'position',[0,.98,],'fontsize',20);

