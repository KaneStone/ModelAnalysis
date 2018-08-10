% NO2 anomalies
clear all

%% import SWOOSH
Stimeperiod = [2000 2014];
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201611-latpress-10deg-L31.nc');
SWOOSHyears = repmat(1984:2015,[12,1]);
SWOOSHyears = [SWOOSHyears(:);ones(11,1)*2016];
SWOOSHextract = permute(SWOOSH.combinedo3q(:,:,SWOOSHyears >= Stimeperiod(1) & SWOOSHyears <= Stimeperiod(2)),[2,1,3]);

%% Read in WACCM NO2 data
timeperiod = [2000 2014];

commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[~, MolConc, ~, Pressure, ~, ~, Latitudes, Longitudes, datedata] = ReadWACCMvertical('NO2','monthly',commondir,0,1);

MolConc = rmfield(MolConc, 'MAM1990');
datedata = rmfield(datedata, 'MAM1990');
Pressure = rmfield(Pressure, 'MAM1990');
runs = fields(MolConc);
runs = {runs{2},runs{6}};

intpres = [SWOOSH.level;[.9,.8,.7,.6,.5,.4,.3,.2,.1]'];

for i = 1:2
    CCMIy = CCMI_years(datedata.(runs{i}).date,1);
    dateind = CCMIy >= timeperiod(1) & CCMIy <= timeperiod(2);
    dateindHF = CCMIy >= timeperiod(1)-1 & CCMIy <= timeperiod(2);
    for j = 1:size(MolConc.Chemonlynoleap,1)
        [SDWaccmData(i).NO2(:,j,:),~] = intRegPres(squeeze(MolConc.(runs{i})(j,:,dateind)),...
            squeeze(Pressure.(runs{i})(j,:,dateind))./100,intpres);
    end
    SDWaccmData(i).NO2nointerp = squeeze(MolConc.(runs{i})(:,:,dateind));
end

%% import Chem-only and MAM
SDtimeperiod = [2000 2014];
[SDWaccmDataO3,latitudes,longitudes,intpres] = ReadInSDWACCM(intpres,SDtimeperiod);

% read in level data

[~,levdata,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/O3/O3_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.1999.fsst.nl.004f.cam.h0zm_merged_c160407.nc');
WAClev = levdata.lev;

%% extract NO2 anomalies

lats = [-70 -60;60 70];
%lats = [80 90];

for j = 1:size(lats,1)
    latind = Latitudes >= lats(j,1) & Latitudes <= lats(j,2);
    SWOOSHlats = SWOOSH.lat >= lats(j,1) & SWOOSH.lat <= lats(j,2);
    for l = 1:length(SDWaccmData)    
        for i = 1:size(SDWaccmData(1).NO2,1)       
            SDweighted(j,l,i,:) = weightedaverage(squeeze(SDWaccmData(l).NO2(i,latind,:)),Latitudes(latind));                
            %for j = 1:size(SDWaccmData(1).O3,1)
                for k = 1:12
                    %SDanom(l,i,k,:) = squeeze(SDWaccmData(l).O3(j,i,k:12:end)) - nanmean(squeeze(SDWaccmData(l).O3(j,i,k:12:end)));
                    SDanom(j,l,i,k,:) = squeeze(SDweighted(j,l,i,k:12:end)) - nanmean(squeeze(SDweighted(j,l,i,96+k:12:144)));
                    for m = 1:size(SDWaccmData(1).NO2,2)       
                        SDanomalllats(l,i,m,k,:) = (squeeze(SDWaccmData(l).NO2(i,m,k:12:end)) - ...
                            nanmean(squeeze(SDWaccmData(l).NO2(i,m,96+k:12:144))));
                    end
                end
            %end
        end    
        for i = 1:size(SDWaccmData(1).NO2nointerp,2)       
            SDweighted_nointerp(j,l,i,:) = weightedaverage(squeeze(SDWaccmData(l).NO2nointerp(latind,i,:)),Latitudes(latind));        
            %for j = 1:size(SDWaccmData(1).O3,1)
                for k = 1:12
                    %SDanom(l,i,k,:) = squeeze(SDWaccmData(l).O3(j,i,k:12:end)) - nanmean(squeeze(SDWaccmData(l).O3(j,i,k:12:end)));
                    SDanom_nointerp(j,l,i,k,:) = squeeze(SDweighted_nointerp(j,l,i,k:12:end)) - nanmean(squeeze(SDweighted_nointerp(j,l,i,96+k:12:144)));
                    for m = 1:size(SDWaccmData(1).NO2nointerp,1)       
                        SDanomalllats_nointerp(l,i,m,k,:) = (squeeze(SDWaccmData(l).NO2nointerp(m,i,k:12:end)) - ...
                            nanmean(squeeze(SDWaccmData(l).NO2nointerp(m,i,96+k:12:144))));
                    end
                end
            %end
        end    
    end
end
    SDanom2 = SDanom(:,:,:,:)*1e9; %lats,run,pressure,time
    SDanom_nointerp2 = SDanom_nointerp(:,:,:,:)*1e9;

%% extract O3 anomalies
for m = 1:2
    latind = Latitudes >= lats(m,1) & Latitudes <= lats(m,2);
    SWOOSHlats = SWOOSH.lat >= lats(m,1) & SWOOSH.lat <= lats(m,2);
    for l = 1:length(SDWaccmDataO3)    
        for i = 1:size(SDWaccmDataO3(1).O3,1)       
            SDweighted_O3(m,l,i,:) = weightedaverage(squeeze(SDWaccmDataO3(l).O3(i,latind,:)),latitudes(latind));        
            %for j = 1:size(SDWaccmDataO3(1).O3,1)
                for k = 1:12
                    %SDanom(l,i,k,:) = squeeze(SDWaccmDataO3(l).O3(j,i,k:12:end)) - nanmean(squeeze(SDWaccmDataO3(l).O3(j,i,k:12:end)));
                    SDanom_O3(m,l,i,k,:) = squeeze(SDweighted_O3(m,l,i,k:12:end)) - nanmean(squeeze(SDweighted_O3(m,l,i,k:12:end)));
                    for j = 1:size(SDWaccmDataO3(1).NO2,2)       
                        SDanomalllats_O3(l,i,j,k,:) = (squeeze(SDWaccmDataO3(l).O3(i,j,k:12:end)) - ...
                            nanmean(squeeze(SDWaccmDataO3(l).O3(i,j,k:12:end))));
                    end
                end
            %end
        end    
        for i = 1:size(SDWaccmDataO3(1).O3nointerp,2)       
            SDweighted_nointerp_O3(m,l,i,:) = weightedaverage(squeeze(SDWaccmDataO3(l).O3nointerp(latind,i,:)),latitudes(latind));        
            %for j = 1:size(SDWaccmDataO3(1).O3,1)
                for k = 1:12
                    %SDanom(l,i,k,:) = squeeze(SDWaccmDataO3(l).O3(j,i,k:12:end)) - nanmean(squeeze(SDWaccmDataO3(l).O3(j,i,k:12:end)));
                    SDanom_nointerp_O3(m,l,i,k,:) = squeeze(SDweighted_nointerp(m,l,i,k:12:end)) - nanmean(squeeze(SDweighted_nointerp(m,l,i,k:12:end)));
                end
            %end
        end    
    end
    
    % SWOOSH
    for i = 1:size(SWOOSHextract,1)       
        SWOOSHweighted(m,i,:) = weightedaverage(squeeze(SWOOSHextract(i,SWOOSHlats,:)),SWOOSH.lat(SWOOSHlats)); 
        for k = 1:12
            SWOOSHanom_O3(m,i,k,:) = squeeze(SWOOSHweighted(m,i,k:12:end)) - nanmean(squeeze(SWOOSHweighted(m,i,k:12:end)));
            SWOOSHanom_O3_noweight(m,i,k,:) = squeeze(SWOOSHextract(i,SWOOSHlats,k:12:end)) - nanmean(squeeze(SWOOSHextract(i,SWOOSHlats,k:12:end)));
        end
    end
    
end
SWOOSHanom2_O3 = cat(2,SWOOSHanom_O3_noweight,zeros(size(SWOOSHanom_O3_noweight,1),9,size(SWOOSHanom_O3_noweight,3),size(SWOOSHanom_O3_noweight,4)));
SWOOSHanom2_O3 (SWOOSHanom2_O3 == 0) = NaN;
SWOOSHanom3_O3(:,1,:,:) = SWOOSHanom2_O3(:,:,:);
SDanom2_O3 = cat(2,SDanom_O3(:,:,:,:)*1e6,SWOOSHanom3_O3);
SDanom_nointerp2_O3 = SDanom_nointerp_O3(:,:,:,:)*1e6;


%%
lev = 3;
[~,levind] = min(abs(intpres-lev));
SDanomalllats_atlev = squeeze(SDanomalllats(:,levind,:,:))*1e9;
SDanomalllats_atlev_O3 = squeeze(SDanomalllats_O3(:,levind,:,:))*1e6;

%% contourplot 


save('/Volumes/MyBook/work/data/regressionOutput/NO2forregression2.mat','SDanomalllats','SDanomalllats_nointerp');

%% plot unusual regression results

for i = 1:size(SDanomalllats_O3,2)
    for j = 1:size(SDanomalllats_O3,3)
        for k = 1:size(SDanomalllats_O3,4)
            rNO2O3_chemonly(i,j,k) = corr(squeeze(SDanomalllats_O3(2,i,j,k,:)),squeeze(SDanomalllats(2,i,j,k,:)));
            rNO2O3_MAM(i,j,k) = corr(squeeze(SDanomalllats_O3(1,i,j,k,:)),squeeze(SDanomalllats(1,i,j,k,:)));
        end
    end
end

%% 
plotcorr = 0;
if plotcorr
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {'Chem-only NO2, O3 correlation','MAM NO2, O3 correlation'};

    for i = 1:12
        clearvars monthtoplot
        monthtoplot(1,:,:) = rNO2O3_chemonly(:,:,i);
        monthtoplot(2,:,:) = rNO2O3_MAM(:,:,i);

        monthtoplot = permute(monthtoplot,[1,3,2]);

        subplotmaps(monthtoplot,latitudes,log(intpres),{'div','RdBu'},1,[],22,titles,'Year','Pressure (hPa)','correlation','on',...
        [-1 1],22,-90:10:90,-90:10:90,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[-90 90],[log(.1) log(300)],1,'none',0,'');

        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/',sprintf('%02d',i),'_N2O_O3_Scorr'];

        export_fig(filename,'-png');

    end
    
    %%
    cbrewqual = cbrewer('qual','Set1',10);
    plot(squeeze(SDanomalllats_atlev_O3(2,80,1:12:end)),'color',cbrewqual(1,:))
    hold on
    plot(squeeze(SDanomalllats_atlev(2,80,1:12:end))./1e1,'--','color',cbrewqual(1,:))
    corr(squeeze(SDanomalllats_atlev_O3(2,80,1:12:end)),squeeze(SDanomalllats_atlev(2,80,1:12:end)))

    plot(squeeze(SDanomalllats_atlev_O3(2,88,1:12:end)),'color',cbrewqual(2,:))
    hold on
    plot(squeeze(SDanomalllats_atlev(2,88,1:12:end))./1e1,'--','color',cbrewqual(2,:))

    corr(squeeze(SDanomalllats_atlev_O3(2,95,1:12:end)),squeeze(SDanomalllats_atlev(2,95,1:12:end)))

    plot(squeeze(SDanomalllats_atlev_O3(2,95,1:12:end)),'color',cbrewqual(3,:))
    hold on
    plot(squeeze(SDanomalllats_atlev(2,95,1:12:end))./1e1,'--','color',cbrewqual(3,:))

end

      

%%
plotraw = 0;
if plotraw
    clearvars extract extract2

    extract(2,:,:) = squeeze(SDanom2(1,:,:));
    extract(1,:,:) = squeeze(SDanom2(2,:,:));
    extract(4,:,:) = squeeze(SDanom2_O3(1,:,:));
    extract(3,:,:) = squeeze(SDanom2_O3(2,:,:));
    extract = permute(extract,[1,3,2]);

    extract2(2,:,:) = squeeze(SDanomalllats_atlev(1,:,:));
    extract2(1,:,:) = squeeze(SDanomalllats_atlev(2,:,:));
    extract2(4,:,:) = squeeze(SDanomalllats_atlev_O3(1,:,:));
    extract2(3,:,:) = squeeze(SDanomalllats_atlev_O3(2,:,:));
    extract2 = permute(extract2,[1,3,2]);
    
    
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {['Chem-only NO2, ',num2str(lats(1)),' to ',num2str(lats(2)),'N'],['MAM NO2, ',num2str(lats(1)),' to ',num2str(lats(2)),'N'],...
        ['Chem-only O3, ',num2str(lats(1)),' to ',num2str(lats(2)),'N'],['MAM O3, ',num2str(lats(1)),' to ',num2str(lats(2)),'N']}; 

    subplotmaps(extract,1:size(extract,2),log(intpres),{'div','RdBu'},1,[],22,titles,'Year','Pressure (hPa)','ppmv anomaly','on',...
        [-1 1],22,1:24:size(extract,2),2000:2:2000+size(extract,2)/12,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[0 size(extract,2)+1],[log(.1) log(300)],1,'none',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','N2O_O3_SPEandchemonlyMAM_northpole'];

    export_fig(filename,'-png');
    %%
    titles = {['Chem-only NO2'],['MAM NO2'],...
        ['Chem-only O3'],['MAM O3']}; 

    subplotmaps(extract2,1:size(extract,2),Latitudes,{'div','RdBu'},1,[],22,titles,'Year','Latitude','ppmv anomaly','on',...
        [-1 1],22,1:24:size(extract,2),2000:2:2000+size(extract,2)/12,...
        -90:10:90,-90:10:90,{[num2str(lev),' hPa']} ,1,[0 size(extract,2)+1],[-90 90],0,'none',0,'');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','N2O_O3_SPEandchemonlyMAM_',num2str(lev),'hPa_latpres'];

    export_fig(filename,'-png');
end

%% plot NO2 in chem only and ozone overaly
cbrew2 = cbrewer('div','PRGn',12);
cbrew2 = [cbrew2(1:3,:);cbrew2(end-2:end,:)];
%extract3 = extract(1:2,:,:);
clearvars extract3 extract4
extract3(1:2,:,:) = squeeze(SDanom2(:,2,:,:));
extract3(3:4,:,:) = squeeze(SDanom2(:,1,:,:));
extract3(5:6,:,:) = squeeze(SDanom2(:,1,:,:));
extract3 = permute(extract3,[1,3,2]);

%extract4 = extract(3:4,:,:);

extract4(1:2,:,:) = squeeze(SDanom2_O3(:,2,:,:));
extract4(3:4,:,:) = squeeze(SDanom2_O3(:,1,:,:));
extract4(5:6,:,:) = squeeze(SDanom2_O3(:,3,:,:));
extract4 = permute(extract4,[1,3,2]);


prestick = [200,100,90:-10:10,9:-1:1,.9:-.1:.1];
presticklabel = {200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
    5,[],[],2,1,[],[],[],[],.5,[],[],.2,.1};
logprestick = log(prestick);

mtit = {['2000',char(8211),'2005 O3 and NO2 anomalies']};
titles = {['Chem-only, ',num2str(abs(lats(1,2))),char(8211),num2str(abs(lats(1,1))),'{\circ}','S'],...
    ['Chem-only, ',num2str(abs(lats(2,2))),char(8211),num2str(abs(lats(2,1))),'{\circ}','N'],...
    ['Chem-Dyn-Vol, ',num2str(abs(lats(1,1))),char(8211),num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['Chem-Dyn-Vol, ',num2str(abs(lats(1,1))),char(8211),num2str(abs(lats(1,2))),'{\circ}','N'],...
    ['SWOOSH, ',num2str(abs(lats(1,1))),char(8211),num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['SWOOSH, ',num2str(abs(lats(1,1))),char(8211),num2str(abs(lats(1,2))),'{\circ}','N']};

[fig,sh] = subplotmaps(extract4(:,1:60,:),1:size(extract3(1,1:60,:),2),log(intpres),{'div','RdBu'},1,[],20,titles,'Year','Pressure (hPa)','O3 ppmv anomaly','on',...
    [-1 1],22,1:12:size(extract3,2),2000:1:2000+size(extract3,2)/12,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[1 size(extract3(1,1:60),2)],[log(.1) log(300)],1,'none',0,'');
hold on

for i = 1:6
    shpos = get(sh(i),'position');
    if i > 2 && i <= 4 
        set(sh(i),'position',[shpos(1),shpos(2)+.02,shpos(3:4)]);
    elseif i > 4 
        set(sh(i),'position',[shpos(1),shpos(2)+.04,shpos(3:4)]);
    end
    axes2 = axes;
    axes2.Visible = 'off';
    ch2 = contour(1:size(extract3(1,1:60,:),2),log(intpres),squeeze(extract3(i,1:60,:))',[-3:1:3],'LineWidth',2);%[-.2:.2:2]
    set(gca,'ytick',fliplr(logprestick),'yticklabel',fliplr(presticklabel),'YDir','reverse');
    set(gca,'color','none');
    set(gca,'xtick',1:12:size(extract3,2)-1,'xticklabel',2000:1:2000+size(extract3,2)/12);
    set(axes2,'position',get(sh(i),'position'));
    colormap(axes2,cbrew2)
    caxis([-3 3])
    if i == 2
        ch = colorbar;
        set(ch,'orientation','horizontal','position',[.1 .062 .755 .02],'fontsize',22);
        set(ch,'Ticks',-3:1:3)  
        set(get(ch,'ylabel'),'string','NO2 ppbv anomaly','fontsize',24)
    end 
    axes2.XTick = [];
    axes2.YTick = [];
end

abc = get(gcf,'children');
set(abc(8),'position',[.87 .15 .022 .775])

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','N2O_O3_overlay_chemonlyMAM_',num2str(lats(1)),'-',num2str(lats(2)),'N'];
 
export_fig(filename,'-png');

