%% LENS minimum and maximum locations

% Read in Can2ESM Temperature data
clear all
preslev = '50hPa';
directory = ['/Volumes/MyBook/work/data/LENS/',preslev,'/3090S/'];
files = dir([directory,'*.nc']);
read_in = 1;
plot_ens = 0;
plot_all = 0;
plot_ensall = 1;
line_plots = 0;
timeperiod = [1995,2014];
sims = {'historical-Natural','Historical','Anthropogenic aerosols','Stratospheric ozone'};
structnames = {'All'};
structnames2 = {'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15',...
    'e16','e17','e18','e19','e20','e21','e22','e23','e24','e25','e26','e27','e28','e29','e30'};
if read_in
    [ta,years,ta_ensall,ensallyears,lon,lat] = Read_in_LENSfortrends(directory,files,'T');
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
dates3 = [2050,2080];
dateind3 = find(years2 >= dates3(1) & years2 <= dates3(2));

fsize = 20;
createfig('medium','on');
cbrew2 = cbrewer('qual','Paired',10);
LENSlonlines = permute(nanmean(cat(4,squeeze(cat(4,ta_alltogether(:,:,27,dateind(1)-1+9:12:dateind(length(dateind))))),squeeze(cat(4,ta_alltogether(:,:,27,dateind(1)-1+10:12:dateind(length(dateind))))),squeeze(cat(4,ta_alltogether(:,:,27,dateind(1)-1+11:12:dateind(length(dateind)))))),4),[2,3,1]);
test = reshape(LENSlonlines,288,size(LENSlonlines,2)*size(LENSlonlines,3));
plot(lon,test,'color',cbrew2(1,:))
hold on
amplitude = max(test)-min(test);
ampstd = std(amplitude);
ampmean = nanmean(amplitude);

LENSlonlines3 = permute(nanmean(cat(4,squeeze(cat(4,ta_alltogether(:,:,27,dateind3(1)-1+9:12:dateind3(length(dateind3))))),squeeze(cat(4,ta_alltogether(:,:,27,dateind3(1)-1+10:12:dateind3(length(dateind3))))),squeeze(cat(4,ta_alltogether(:,:,27,dateind3(1)-1+11:12:dateind3(length(dateind3)))))),4),[2,3,1]);
test3 = reshape(LENSlonlines3,288,size(LENSlonlines3,2)*size(LENSlonlines3,3));
hold on
plot(lon,test3,'color',cbrew2(3,:))

LENSlonlines2 = permute(nanmean(cat(4,squeeze(cat(4,ta_alltogether(:,:,27,dateind2(1)-1+9:12:dateind2(length(dateind2))))),squeeze(cat(4,ta_alltogether(:,:,27,dateind2(1)-1+10:12:dateind2(length(dateind2))))),squeeze(cat(4,ta_alltogether(:,:,27,dateind2(1)-1+11:12:dateind2(length(dateind2)))))),4),[2,3,1]);
test2 = reshape(LENSlonlines2,288,size(LENSlonlines2,2)*size(LENSlonlines2,3));
hold on
plot(lon,test2,'color',cbrew2(5,:))

ph(1) = plot(lon,nanmean(test,2),'LineWidth',5,'color',cbrew2(2,:));
ph(2) = plot(lon,nanmean(test2,2),'--','LineWidth',5,'color',cbrew2(6,:));
ph(3) = plot(lon,nanmean(test3,2),':','LineWidth',5,'color',cbrew2(4,:));
amplitude2 = max(test2)-min(test2);
ampstd2 = std(amplitude2);
ampmean2 = nanmean(amplitude2);
set(gca,'fontsize',fsize+2);
xlim([-5 365]);
xlabel('Longitude ({\circ}E)','fontsize',fsize+4);
ylabel('Temperature (K)','fontsize',fsize+4);
title('Austral spring temperatures at 65{\circ}S and 50 hPa','fontsize',fsize+6);
lh = legend(ph,'1995-2025 (high chlorine)','1955-1979 (low chlorine)','2050-2080 (high GHGs - low chlorine)');
set(lh,'fontsize',fsize+4,'box','off','location','south')

filename = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TalkFigures/LENSAmplitudes_65S';
export_fig(filename,'-png');

%% finding maximum and minimum for all members
lats = [-50,-55,-60,-65,-70,-75,-80,-85];
for i = 1:size(ta_alltogether)
    for k = 1:length(lats)
        for j = 1:12        
            [~,latind] = min(abs(lats(k)-lat));   
            [tempminvalue,tempminind] = min(squeeze(ta_alltogether(i,:,latind,j:12:end)));        
            [tempmaxvalue,tempmaxind] = max(squeeze(ta_alltogether(i,:,latind,j:12:end)));       
            LENS.minvalue_m(i).m(j,k,:) = tempminvalue;
            LENS.minind_m(i).m(j,k,:) = tempminind;
            LENS.maxvalue_m(i).m(j,k,:) = tempmaxvalue;
            LENS.maxind_m(i).m(j,k,:) = tempmaxind;
        end
        
        [tempminvalue,tempminind] = min(nanmean(cat(3,squeeze(ta_alltogether(i,:,latind,9:12:end)),...
            squeeze(ta_alltogether(i,:,latind,10:12:end)),...
            squeeze(ta_alltogether(i,:,latind,11:12:end))),3),[],1);        
        [tempmaxvalue,tempmaxind] = max(nanmean(cat(3,squeeze(ta_alltogether(i,:,latind,9:12:end)),...
            squeeze(ta_alltogether(i,:,latind,10:12:end)),...
            squeeze(ta_alltogether(i,:,latind,11:12:end))),3),[],1);        
        LENS.minvalue_m_spr(i).m(k,:) = tempminvalue;
        LENS.minind_m_spr(i).m(k,:) = tempminind;                                    
        LENS.maxvalue_m_spr(i).m(k,:) = tempmaxvalue;
        LENS.maxind_m_spr(i).m(k,:) = tempmaxind;
        
        
    end
    LENS.minlongitude_m(i).m = lon(LENS.minind_m(i).m);
    LENS.maxlongitude_m(i).m = lon(LENS.maxind_m(i).m);
    LENS.amplitude_m(i).m = LENS.maxvalue_m(i).m - LENS.minvalue_m(i).m;
    
    LENS.minlongitude_m_spr(i).m = lon(LENS.minind_m_spr(i).m);
    LENS.maxlongitude_m_spr(i).m = lon(LENS.maxind_m_spr(i).m);
    LENS.amplitude_m_spr(i).m = LENS.maxvalue_m_spr(i).m - LENS.minvalue_m_spr(i).m;
end

for i = 1:size(ta_alltogether,1)    
    for k = 1:length(lats)   
        for j = 1:12
            for l = 1:size(LENS.maxlongitude_m(i).m(j,k,:),3)
                if LENS.maxlongitude_m(i).m(j,k,l) > 300
                    LENS.maxlongitude_m(i).m(j,k,l) = LENS.maxlongitude_m(i).m(j,k,l)-360;
                end
                if LENS.minlongitude_m(i).m(j,k,l) < 100
                    LENS.minlongitude_m(i).m(j,k,l) = LENS.minlongitude_m(i).m(j,k,l)+360;
                end                                                
            end
        end
        for l = 1:size(LENS.maxlongitude_m_spr(i).m(k,:),2)
            if LENS.maxlongitude_m_spr(i).m(k,l) > 300
                LENS.maxlongitude_m_spr(i).m(k,l) = LENS.maxlongitude_m_spr(i).m(k,l)-360;
            end
            if LENS.minlongitude_m_spr(i).m(k,l) < 100
                LENS.minlongitude_m_spr(i).m(k,l) = LENS.minlongitude_m_spr(i).m(k,l)+360;                           
            end
        end
    end
end

ens = [1,31];
for i = 1:length(structnames)
    LENS.maxlongitude2.(structnames{i}) = nanmean(cat(4,LENS.maxlongitude_m(ens(i):ens(i+1)-1).m),4);
    LENS.minlongitude2.(structnames{i}) = nanmean(cat(4,LENS.minlongitude_m(ens(i):ens(i+1)-1).m),4);
    LENS.maxlongitude_spr2.(structnames{i}) = nanmean(cat(3,LENS.maxlongitude_m_spr(ens(i):ens(i+1)-1).m),3);
    LENS.minlongitude_spr2.(structnames{i}) = nanmean(cat(3,LENS.minlongitude_m_spr(ens(i):ens(i+1)-1).m),3);
    LENS.amplitude2.(structnames{i}) = nanmean(cat(4,LENS.amplitude_m(ens(i):ens(i+1)-1).m),4);
    LENS.amplitude_spr2.(structnames{i}) = nanmean(cat(3,LENS.amplitude_m_spr(ens(i):ens(i+1)-1).m),3);
end

LENS.longitude = lon;
LENS.latitude = lats;

%%
for i = 1:size(ta_ensalltogether,1)    
    for k = 1:length(lats)
        [~,latind] = min(abs(lats(k)-lat));            
        for j = 1:12                    
            [tempminvalue,tempminind] = min(squeeze(ta_ensalltogether(i,:,latind,j:12:end)));        
            [tempmaxvalue,tempmaxind] = max(squeeze(ta_ensalltogether(i,:,latind,j:12:end)));                   
            LENS.minvalue.(structnames{i})(j,k,:) = tempminvalue;
            LENS.minind.(structnames{i})(j,k,:) = tempminind;                                    
            LENS.maxvalue.(structnames{i})(j,k,:) = tempmaxvalue;
            LENS.maxind.(structnames{i})(j,k,:) = tempmaxind;            
        end
        
        [tempminvalue,tempminind] = min(nanmean(cat(3,squeeze(ta_ensalltogether(i,:,latind,9:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,10:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,11:12:end))),3),[],1);        
        [tempmaxvalue,tempmaxind] = max(nanmean(cat(3,squeeze(ta_ensalltogether(i,:,latind,9:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,10:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,11:12:end))),3),[],1);               
        LENS.minvalue_spr.(structnames{i})(k,:) = tempminvalue;
        LENS.minind_spr.(structnames{i})(k,:) = tempminind;                                    
        LENS.maxvalue_spr.(structnames{i})(k,:) = tempmaxvalue;
        LENS.maxind_spr.(structnames{i})(k,:) = tempmaxind;
    end    
    
    LENS.minlongitude.(structnames{i}) = lon(LENS.minind.(structnames{i}));
    LENS.maxlongitude.(structnames{i}) = lon(LENS.maxind.(structnames{i}));
    LENS.amplitude.(structnames{i}) = LENS.maxvalue.(structnames{i}) - LENS.minvalue.(structnames{i});
    
    LENS.minlongitude_spr.(structnames{i}) = lon(LENS.minind_spr.(structnames{i}));
    LENS.maxlongitude_spr.(structnames{i}) = lon(LENS.maxind_spr.(structnames{i}));
    LENS.amplitude_spr.(structnames{i}) = LENS.maxvalue_spr.(structnames{i}) - LENS.minvalue_spr.(structnames{i});
    
end
for i = 1:size(ta_ensalltogether,1)    
    for k = 1:length(lats)     
        for j = 1:12
            for l = 1:size(LENS.maxlongitude.(structnames{i})(j,k,:),3)
                if LENS.maxlongitude.(structnames{i})(j,k,l) > 300
                    LENS.maxlongitude.(structnames{i})(j,k,l) = LENS.maxlongitude.(structnames{i})(j,k,l)-360;
                end
                if LENS.minlongitude.(structnames{i})(j,k,l) < 100
                    LENS.minlongitude.(structnames{i})(j,k,l) = LENS.minlongitude.(structnames{i})(j,k,l)+360;
                end            
            end
        end

        for l = 1:size(LENS.maxlongitude_spr.(structnames{i})(k,:),2)
            if LENS.maxlongitude_spr.(structnames{i})(k,l) > 300
                LENS.maxlongitude_spr.(structnames{i})(k,l) = LENS.maxlongitude_spr.(structnames{i})(k,l)-360;
            end
            if LENS.minlongitude_spr.(structnames{i})(k,l) < 100
                LENS.minlongitude_spr.(structnames{i})(k,l) = LENS.minlongitude_spr.(structnames{i})(k,l)+360;                           
            end
        end        
    
    end
end
LENS.years = years;
LENS.ensyears = ensallyears;
LENS.ensnames = structnames;
LENS.names = structnames2;
save('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/LENS.mat','LENS');

clearvars dateind

%% plotting amplitude
createfig('medium','on')
contourf(amplitude.All,14)
%caxis([0,22]);
colorbar
title('LENS');


%% plotting
fsize = 18;
cbrew = cbrewer('qual','Set1',10);
createfig('large','on');
%subplot(2,1,1)
m_proj('Stereographic','lon',-180,'lat',-90,'rad',50)
m_coast('color','k','LineWidth',3);
m_grid('ytick',[-90 -80 -70 -60 -50 -40 -30 -20],'XaxisLocation','top','fontsize',fsize);

m_line(squeeze(minlongitude.(structnames{i})(10,:))',lats','color',cbrew(i,:),'LineWidth',3);    
m_line(squeeze(maxlongitude.(structnames{i})(10,:))',lats','color',cbrew(i,:),'LineWidth',3);
    
title([num2str(timeperiod(1)),'-',num2str(timeperiod(2))],'position',[0,.98,],'fontsize',20);