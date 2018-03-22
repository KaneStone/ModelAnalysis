%CCMIWACCM_minmaxlongitudes
%% LENS minimum and maximum locations

% Read in Can2ESM Temperature data
clear all
preslev = '50hPa';
directory = ['/Volumes/MyBook/work/data/CCMI/ACCESSforzonal/',preslev,'/T/'];
files = dir([directory,'*.nc']);
read_in = 1;
plot_ens = 0;
plot_all = 0;
plot_ensall = 1;
line_plots = 0;
timeperiod = [1990,2010];
sims = {'REF-C1','REF-C2','lowGHG','lowODS'};
structnames = {'REFC1','REFC2','lowGHG','lowODS'};
structnames2 = {'REFC1_1','REFC2_1','REFC2_2','lowGHG_1','lowGHG_2','lowODS_1','lowODS_2'};
if read_in
    [ta,years,ta_ensall,ensallyears,lon,lat] = Read_in_CCMIACCESS_fortrends(directory,files,'ta');
end

%% restructurig temperature ensemble
for i = 1:size(ta_ensall,1)
    ta_ensalltogether(i,:,:,:) = permute(cat(3,ta_ensall(i,:).w),[3,1,2]);
end

for i = 1:size(ta,1)
    ta_alltogether(i,:,:,:) = permute(cat(3,ta(i,:).w),[3,1,2]);
end

%% finding minimum and maximum longitudes
lats = [-50,-55,-60,-65,-70,-75,-80,-85];
for i = 1:size(ta_alltogether)
    for k = 1:length(lats)
        for j = 1:12        
            [~,latind] = min(abs(lats(k)-lat));   
            [tempminvalue,tempminind] = min(squeeze(ta_alltogether(i,:,latind,j:12:end)));        
            [tempmaxvalue,tempmaxind] = max(squeeze(ta_alltogether(i,:,latind,j:12:end)));       
            ACCESS.minvalue_m(i).m(j,k,:) = tempminvalue;
            ACCESS.minind_m(i).m(j,k,:) = tempminind;
            ACCESS.maxvalue_m(i).m(j,k,:) = tempmaxvalue;
            ACCESS.maxind_m(i).m(j,k,:) = tempmaxind;
        end
        
        [tempminvalue,tempminind] = min(nanmean(cat(3,squeeze(ta_alltogether(i,:,latind,9:12:end)),...
            squeeze(ta_alltogether(i,:,latind,10:12:end)),...
            squeeze(ta_alltogether(i,:,latind,11:12:end))),3),[],1);        
        [tempmaxvalue,tempmaxind] = max(nanmean(cat(3,squeeze(ta_alltogether(i,:,latind,9:12:end)),...
            squeeze(ta_alltogether(i,:,latind,10:12:end)),...
            squeeze(ta_alltogether(i,:,latind,11:12:end))),3),[],1);        
        ACCESS.minvalue_m_spr(i).m(k,:) = tempminvalue;
        ACCESS.minind_m_spr(i).m(k,:) = tempminind;                                    
        ACCESS.maxvalue_m_spr(i).m(k,:) = tempmaxvalue;
        ACCESS.maxind_m_spr(i).m(k,:) = tempmaxind;
        
        
    end
    ACCESS.minlongitude_m(i).m = lon(ACCESS.minind_m(i).m);
    ACCESS.maxlongitude_m(i).m = lon(ACCESS.maxind_m(i).m);
    ACCESS.amplitude_m(i).m = ACCESS.maxvalue_m(i).m - ACCESS.minvalue_m(i).m;
    
    ACCESS.minlongitude_m_spr(i).m = lon(ACCESS.minind_m_spr(i).m);
    ACCESS.maxlongitude_m_spr(i).m = lon(ACCESS.maxind_m_spr(i).m);
    ACCESS.amplitude_m_spr(i).m = ACCESS.maxvalue_m_spr(i).m - ACCESS.minvalue_m_spr(i).m;
end

for i = 1:size(ta_alltogether,1)    
    for k = 1:length(lats)   
        for j = 1:12
            for l = 1:size(ACCESS.maxlongitude_m(i).m(j,k,:),3)
                if ACCESS.maxlongitude_m(i).m(j,k,l) > 300
                    ACCESS.maxlongitude_m(i).m(j,k,l) = ACCESS.maxlongitude_m(i).m(j,k,l)-360;
                end
                if ACCESS.minlongitude_m(i).m(j,k,l) < 100
                    ACCESS.minlongitude_m(i).m(j,k,l) = ACCESS.minlongitude_m(i).m(j,k,l)+360;
                end                                                
            end
        end
        for l = 1:size(ACCESS.maxlongitude_m_spr(i).m(k,:),2)
            if ACCESS.maxlongitude_m_spr(i).m(k,l) > 300
                ACCESS.maxlongitude_m_spr(i).m(k,l) = ACCESS.maxlongitude_m_spr(i).m(k,l)-360;
            end
            if ACCESS.minlongitude_m_spr(i).m(k,l) < 100
                ACCESS.minlongitude_m_spr(i).m(k,l) = ACCESS.minlongitude_m_spr(i).m(k,l)+360;                           
            end
        end
    end
end

ens = [1,2,4,6,8];
for i = 1:length(structnames)
    ACCESS.maxlongitude2.(structnames{i}) = nanmean(cat(4,ACCESS.maxlongitude_m(ens(i):ens(i+1)-1).m),4);
    ACCESS.minlongitude2.(structnames{i}) = nanmean(cat(4,ACCESS.minlongitude_m(ens(i):ens(i+1)-1).m),4);
    ACCESS.maxlongitude_spr2.(structnames{i}) = nanmean(cat(3,ACCESS.maxlongitude_m_spr(ens(i):ens(i+1)-1).m),3);
    ACCESS.minlongitude_spr2.(structnames{i}) = nanmean(cat(3,ACCESS.minlongitude_m_spr(ens(i):ens(i+1)-1).m),3);
    ACCESS.amplitude2.(structnames{i}) = nanmean(cat(4,ACCESS.amplitude_m(ens(i):ens(i+1)-1).m),4);
    ACCESS.amplitude_spr2.(structnames{i}) = nanmean(cat(3,ACCESS.amplitude_m_spr(ens(i):ens(i+1)-1).m),3);
end

ACCESS.longitude = lon;
ACCESS.latitude = lats;

%%
for i = 1:size(ta_ensalltogether,1)    
    for k = 1:length(lats)
        [~,latind] = min(abs(lats(k)-lat));            
        for j = 1:12                    
            [tempminvalue,tempminind] = min(squeeze(ta_ensalltogether(i,:,latind,j:12:end)));        
            [tempmaxvalue,tempmaxind] = max(squeeze(ta_ensalltogether(i,:,latind,j:12:end)));                   
            ACCESS.minvalue.(structnames{i})(j,k,:) = tempminvalue;
            ACCESS.minind.(structnames{i})(j,k,:) = tempminind;                                    
            ACCESS.maxvalue.(structnames{i})(j,k,:) = tempmaxvalue;
            ACCESS.maxind.(structnames{i})(j,k,:) = tempmaxind;            
        end
        
        [tempminvalue,tempminind] = min(nanmean(cat(3,squeeze(ta_ensalltogether(i,:,latind,9:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,10:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,11:12:end))),3),[],1);        
        [tempmaxvalue,tempmaxind] = max(nanmean(cat(3,squeeze(ta_ensalltogether(i,:,latind,9:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,10:12:end)),...
            squeeze(ta_ensalltogether(i,:,latind,11:12:end))),3),[],1);               
        ACCESS.minvalue_spr.(structnames{i})(k,:) = tempminvalue;
        ACCESS.minind_spr.(structnames{i})(k,:) = tempminind;                                    
        ACCESS.maxvalue_spr.(structnames{i})(k,:) = tempmaxvalue;
        ACCESS.maxind_spr.(structnames{i})(k,:) = tempmaxind;
    end    
    
    ACCESS.minlongitude.(structnames{i}) = lon(ACCESS.minind.(structnames{i}));
    ACCESS.maxlongitude.(structnames{i}) = lon(ACCESS.maxind.(structnames{i}));
    ACCESS.amplitude.(structnames{i}) = ACCESS.maxvalue.(structnames{i}) - ACCESS.minvalue.(structnames{i});
    
    ACCESS.minlongitude_spr.(structnames{i}) = lon(ACCESS.minind_spr.(structnames{i}));
    ACCESS.maxlongitude_spr.(structnames{i}) = lon(ACCESS.maxind_spr.(structnames{i}));
    ACCESS.amplitude_spr.(structnames{i}) = ACCESS.maxvalue_spr.(structnames{i}) - ACCESS.minvalue_spr.(structnames{i});
    
end
for i = 1:size(ta_ensalltogether,1)    
    for k = 1:length(lats)     
        for j = 1:12
            for l = 1:size(ACCESS.maxlongitude.(structnames{i})(j,k,:),3)
                if ACCESS.maxlongitude.(structnames{i})(j,k,l) > 300
                    ACCESS.maxlongitude.(structnames{i})(j,k,l) = ACCESS.maxlongitude.(structnames{i})(j,k,l)-360;
                end
                if ACCESS.minlongitude.(structnames{i})(j,k,l) < 100
                    ACCESS.minlongitude.(structnames{i})(j,k,l) = ACCESS.minlongitude.(structnames{i})(j,k,l)+360;
                end            
            end
        end

        for l = 1:size(ACCESS.maxlongitude_spr.(structnames{i})(k,:),2)
            if ACCESS.maxlongitude_spr.(structnames{i})(k,l) > 300
                ACCESS.maxlongitude_spr.(structnames{i})(k,l) = ACCESS.maxlongitude_spr.(structnames{i})(k,l)-360;
            end
            if ACCESS.minlongitude_spr.(structnames{i})(k,l) < 100
                ACCESS.minlongitude_spr.(structnames{i})(k,l) = ACCESS.minlongitude_spr.(structnames{i})(k,l)+360;                           
            end
        end        
    
    end
end
ACCESS.years = years;
ACCESS.ensyears = ensallyears;
ACCESS.ensnames = structnames;
ACCESS.names = structnames2;
save('/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/ACCESS.mat','ACCESS');



%% plotting amplitude
createfig('medium','on')
contourf(amplitude.REFC1,14)
%caxis([0,22]);
colorbar

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