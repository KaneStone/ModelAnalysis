%CCMIWACCM_minmaxlongitudes
%% LENS minimum and maximum locations

% Read in Can2ESM Temperature data
clear all
preslev = '50hPa';
directory = ['/Volumes/MyBook/work/data/CESM-CCMI/T/',preslev,'/'];
files = dir([directory,'*.nc']);
read_in = 1;
plot_ens = 0;
plot_all = 0;
plot_ensall = 1;
line_plots = 0;
timeperiod = [1995,2014];
sims = {'REF-C2','RCP85','lowGHG','lowODS','REF-C1','REF-C1SD'};
structnames = {'REFC2','RCP85','lowGHG','lowODS','REFC1','REFC1SD'};
structnames2 = {'REFC2_1','REFC2_2','REFC2_3','RCP85_1','RCP85_2','RCP85_3','lowGHG_1','lowGHG_2','lowGHG_3','lowODS_1','lowODS_2','lowODS_3','REFC1_1',...
    'REFC1_2','REFC1_3','REFC1_4','REFC1_5','REFC1SD'};
if read_in
    [ta,years,ta_ensall,ensallyears,lon,lat] = Read_in_CCMIWACCM_fortrends(directory,files,'T');
end

%% restructurig temperature ensemble
for i = 1:size(ta_ensall,1)
    ta_ensalltogether(i,:,:,:) = permute(cat(3,ta_ensall(i,:).w),[3,1,2]);
end

for i = 1:size(ta,1)
    ta_alltogether(i,:,:,:) = permute(cat(3,ta(i,:).w),[3,1,2]);
end
dateind(1,:) = find(ensallyears(1).y >= timeperiod(1) & ensallyears(1).y <= timeperiod(2));         
%dateind2(1,:) = find(ensallyears(5).y >= timeperiod(1) & ensallyears(5).y <= timeperiod(2));         

% %% tsting
% 
% %for i = 1:size(
% amplitude = max(squeeze(ta_alltogether(12,:,17,dateind(1)+10-1:12:dateind(end)))) - min(squeeze(ta_alltogether(12,:,17,dateind(1)+10-1:12:dateind(end)))); 
% ampstd = std(amplitude);
% amplitude2 = max(squeeze(ta_alltogether(15,:,17,dateind2(1)+10-1:12:dateind2(end)))) - min(squeeze(ta_alltogether(15,:,17,dateind2(1)+10-1:12:dateind2(end)))); 
% ampstd2 = std(amplitude2);
% figure;
% plot(squeeze(ta_alltogether(2,:,17,dateind(1)+9-1:12:dateind(end))),'r')
% hold on
% plot(nanmean(squeeze(ta_alltogether(2,:,17,dateind(1)+9-1:12:dateind(end))),2),'k','LineWidth',3)
% 
% figure;
% plot(squeeze(ta_alltogether(11,:,17,dateind(1)+9-1:12:dateind(end))),'r')
% hold on
% plot(nanmean(squeeze(ta_alltogether(11,:,17,dateind(1)+9-1:12:dateind(end))),2),'k','LineWidth',3)
% 
% figure;
% plot(squeeze(ta_alltogether(15,:,17,dateind2(1)+10-1:12:dateind2(end))),'r')
% hold on
% plot(nanmean(squeeze(ta_alltogether(15,:,17,dateind2(1)+10-1:12:dateind2(end))),2),'k','LineWidth',3)


%% finding maximum and minimum for all members
lats = [-50,-55,-60,-65,-70,-75,-80,-85];
for i = 1:size(ta_alltogether)
    for k = 1:length(lats)
        for j = 1:12        
            [~,latind] = min(abs(lats(k)-lat));   
            [tempminvalue,tempminind] = min(squeeze(ta_alltogether(i,:,latind,j:12:end)));        
            [tempmaxvalue,tempmaxind] = max(squeeze(ta_alltogether(i,:,latind,j:12:end)));       
            WACCM_CCMI.minvalue_m(i).m(j,k,:) = tempminvalue;
            WACCM_CCMI.minind_m(i).m(j,k,:) = tempminind;
            WACCM_CCMI.maxvalue_m(i).m(j,k,:) = tempmaxvalue;
            WACCM_CCMI.maxind_m(i).m(j,k,:) = tempmaxind;
        end
        
        [tempminvalue,tempminind] = min(nanmean(cat(3,squeeze(ta_alltogether(i,:,latind,9:12:end)),...
            squeeze(ta_alltogether(i,:,latind,10:12:end)),...
            squeeze(ta_alltogether(i,:,latind,11:12:end))),3),[],1);        
        [tempmaxvalue,tempmaxind] = max(nanmean(cat(3,squeeze(ta_alltogether(i,:,latind,9:12:end)),...
            squeeze(ta_alltogether(i,:,latind,10:12:end)),...
            squeeze(ta_alltogether(i,:,latind,11:12:end))),3),[],1);        
        WACCM_CCMI.minvalue_m_spr(i).m(k,:) = tempminvalue;
        WACCM_CCMI.minind_m_spr(i).m(k,:) = tempminind;                                    
        WACCM_CCMI.maxvalue_m_spr(i).m(k,:) = tempmaxvalue;
        WACCM_CCMI.maxind_m_spr(i).m(k,:) = tempmaxind;
        
        
    end
    WACCM_CCMI.minlongitude_m(i).m = lon(WACCM_CCMI.minind_m(i).m);
    WACCM_CCMI.maxlongitude_m(i).m = lon(WACCM_CCMI.maxind_m(i).m);
    WACCM_CCMI.amplitude_m(i).m = WACCM_CCMI.maxvalue_m(i).m - WACCM_CCMI.minvalue_m(i).m;
    
    WACCM_CCMI.minlongitude_m_spr(i).m = lon(WACCM_CCMI.minind_m_spr(i).m);
    WACCM_CCMI.maxlongitude_m_spr(i).m = lon(WACCM_CCMI.maxind_m_spr(i).m);
    WACCM_CCMI.amplitude_m_spr(i).m = WACCM_CCMI.maxvalue_m_spr(i).m - WACCM_CCMI.minvalue_m_spr(i).m;
end

for i = 1:size(ta_alltogether,1)    
    for k = 1:length(lats)   
        for j = 1:12
            for l = 1:size(WACCM_CCMI.maxlongitude_m(i).m(j,k,:),3)
                if WACCM_CCMI.maxlongitude_m(i).m(j,k,l) > 300
                    WACCM_CCMI.maxlongitude_m(i).m(j,k,l) = WACCM_CCMI.maxlongitude_m(i).m(j,k,l)-360;
                end
                if WACCM_CCMI.minlongitude_m(i).m(j,k,l) < 100
                    WACCM_CCMI.minlongitude_m(i).m(j,k,l) = WACCM_CCMI.minlongitude_m(i).m(j,k,l)+360;
                end                                                
            end
        end
        for l = 1:size(WACCM_CCMI.maxlongitude_m_spr(i).m(k,:),2)
            if WACCM_CCMI.maxlongitude_m_spr(i).m(k,l) > 300
                WACCM_CCMI.maxlongitude_m_spr(i).m(k,l) = WACCM_CCMI.maxlongitude_m_spr(i).m(k,l)-360;
            end
            if WACCM_CCMI.minlongitude_m_spr(i).m(k,l) < 100
                WACCM_CCMI.minlongitude_m_spr(i).m(k,l) = WACCM_CCMI.minlongitude_m_spr(i).m(k,l)+360;                           
            end
        end
    end
    
    %averaging
    
    
end

ens = [1,4,7,10,13,18];
for i = 1:length(structnames)-1
    WACCM_CCMI.maxlongitude2.(structnames{i}) = nanmean(cat(4,WACCM_CCMI.maxlongitude_m(ens(i):ens(i+1)-1).m),4);
    WACCM_CCMI.minlongitude2.(structnames{i}) = nanmean(cat(4,WACCM_CCMI.minlongitude_m(ens(i):ens(i+1)-1).m),4);
    WACCM_CCMI.maxlongitude_spr2.(structnames{i}) = nanmean(cat(3,WACCM_CCMI.maxlongitude_m_spr(ens(i):ens(i+1)-1).m),3);
    WACCM_CCMI.minlongitude_spr2.(structnames{i}) = nanmean(cat(3,WACCM_CCMI.minlongitude_m_spr(ens(i):ens(i+1)-1).m),3);
    WACCM_CCMI.amplitude2.(structnames{i}) = nanmean(cat(4,WACCM_CCMI.amplitude_m(ens(i):ens(i+1)-1).m),4);
    WACCM_CCMI.amplitude_spr2.(structnames{i}) = nanmean(cat(3,WACCM_CCMI.amplitude_m_spr(ens(i):ens(i+1)-1).m),3);
end

WACCM_CCMI.longitude = lon;
WACCM_CCMI.latitude = lats;

%%
for i = 1:size(ta_ensalltogether,1)    
    for k = 1:length(lats)
        [~,latind] = min(abs(lats(k)-lat));            
        for j = 1:12                    
            [tempminvalue,tempminind] = min(squeeze(ta_ensalltogether(i,:,latind,j:12:end)));        
            [tempmaxvalue,tempmaxind] = max(squeeze(ta_ensalltogether(i,:,latind,j:12:end)));                   
            WACCM_CCMI.minvalue.(structnames{i})(j,k,:) = tempminvalue;
            WACCM_CCMI.minind.(structnames{i})(j,k,:) = tempminind;                                    
            WACCM_CCMI.maxvalue.(structnames{i})(j,k,:) = tempmaxvalue;
            WACCM_CCMI.maxind.(structnames{i})(j,k,:) = tempmaxind;            
        end
        
        [tempminvalue,tempminind] = min(nanmean(cat(3,squeeze(ta_ensalltogether(i,:,latind,9:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,10:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,11:12:end))),3),[],1);        
        [tempmaxvalue,tempmaxind] = max(nanmean(cat(3,squeeze(ta_ensalltogether(i,:,latind,9:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,10:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,11:12:end))),3),[],1);               
        WACCM_CCMI.minvalue_spr.(structnames{i})(k,:) = tempminvalue;
        WACCM_CCMI.minind_spr.(structnames{i})(k,:) = tempminind;                                    
        WACCM_CCMI.maxvalue_spr.(structnames{i})(k,:) = tempmaxvalue;
        WACCM_CCMI.maxind_spr.(structnames{i})(k,:) = tempmaxind;
    end    
    
    WACCM_CCMI.minlongitude.(structnames{i}) = lon(WACCM_CCMI.minind.(structnames{i}));
    WACCM_CCMI.maxlongitude.(structnames{i}) = lon(WACCM_CCMI.maxind.(structnames{i}));
    WACCM_CCMI.amplitude.(structnames{i}) = WACCM_CCMI.maxvalue.(structnames{i}) - WACCM_CCMI.minvalue.(structnames{i});
    
    WACCM_CCMI.minlongitude_spr.(structnames{i}) = lon(WACCM_CCMI.minind_spr.(structnames{i}));
    WACCM_CCMI.maxlongitude_spr.(structnames{i}) = lon(WACCM_CCMI.maxind_spr.(structnames{i}));
    WACCM_CCMI.amplitude_spr.(structnames{i}) = WACCM_CCMI.maxvalue_spr.(structnames{i}) - WACCM_CCMI.minvalue_spr.(structnames{i});
    
end
for i = 1:size(ta_ensalltogether,1)    
    for k = 1:length(lats)     
        for j = 1:12
            for l = 1:size(WACCM_CCMI.maxlongitude.(structnames{i})(j,k,:),3)
                if WACCM_CCMI.maxlongitude.(structnames{i})(j,k,l) > 300
                    WACCM_CCMI.maxlongitude.(structnames{i})(j,k,l) = WACCM_CCMI.maxlongitude.(structnames{i})(j,k,l)-360;
                end
                if WACCM_CCMI.minlongitude.(structnames{i})(j,k,l) < 100
                    WACCM_CCMI.minlongitude.(structnames{i})(j,k,l) = WACCM_CCMI.minlongitude.(structnames{i})(j,k,l)+360;
                end            
            end
        end

        for l = 1:size(WACCM_CCMI.maxlongitude_spr.(structnames{i})(k,:),2)
            if WACCM_CCMI.maxlongitude_spr.(structnames{i})(k,l) > 300
                WACCM_CCMI.maxlongitude_spr.(structnames{i})(k,l) = WACCM_CCMI.maxlongitude_spr.(structnames{i})(k,l)-360;
            end
            if WACCM_CCMI.minlongitude_spr.(structnames{i})(k,l) < 100
                WACCM_CCMI.minlongitude_spr.(structnames{i})(k,l) = WACCM_CCMI.minlongitude_spr.(structnames{i})(k,l)+360;                           
            end
        end        
    
    end
end
WACCM_CCMI.years = years;
WACCM_CCMI.ensyears = ensallyears;
WACCM_CCMI.ensnames = structnames;
WACCM_CCMI.names = structnames2;
save('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/WACCM_CCMI.mat','WACCM_CCMI');

%% plotting amplitude
createfig('medium','on')
contourf(amplitude.REFC2,14)
%caxis([0,22]);
colorbar
title('CESM-CCMI');

createfig('medium','on')
contourf(amplitude.REFC1SD,14)
%caxis([0,22]);
colorbar
title('CESM-CCMI - SD');

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
end
legend(mh,sims);

title([num2str(timeperiod(1)),'-',num2str(timeperiod(2))],'position',[0,.98,],'fontsize',20);