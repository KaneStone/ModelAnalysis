% Compare MLS and WACCM
clear all

pres = 100;
fsize = 18;
lwidth = 2;
numberdensity = 1;
sa = 1;
monthstoplot = [7,8,9,10,11,12,1,2,3,4,5,6];
daytitle = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',...
    'Oct','Nov','Dec'};

deviationlayer = 0;
vertical_deviation = 1;
temperature_deviations = 0;
line_plots = 0;
vertical_plots = 0;
comparetoWACCMSAD = 0;
include_SAD_SULFC = 1;
remove_tertiary = 1;
sulfdeviation = 0;
ext = 1;

addpath('/Volumes/ExternalOne/work/data/MLS/code/');
[MLS,presind,allyears,MLSPressure] = importMLSmonthly(pres,0);

%% importing OMI

dirOMI = dir(['/Users/kanestone/work/projects/OMI/output/','*.mat']);
for i = 1:length(dirOMI)
    OMIdata = load(['/Users/kanestone/work/projects/OMI/output/',dirOMI(i).name]);
    OMI_monthly(i,:,:,:) = OMIdata.data_monthly;
    if i == length(dirOMI)
        load(['/Users/kanestone/work/projects/OMI/output/',dirOMI(i).name]);        
    end
end

OMI_monthly_latmean = squeeze(nanmean(OMI_monthly,2));
OMI_monthly_latmean = permute(OMI_monthly_latmean,[1,3,2]);

Meanall = nanmean(OMI_monthly_latmean,1);
Meanall = repmat(Meanall,[size(OMI_monthly_latmean,2),1,1]);
if sa
    stdall = nanstd(OMI_monthly_latmean,0,1);    
    stdall = repmat(stdall,[size(OMI_monthly_latmean,2),1,1]);    
else
    stdall = 1;
end
OMIdeviationnorm = squeeze(OMI_monthly_latmean-Meanall)./stdall;

%OMIstd_OctNov = stdall(1,10:11,25);
%% importing WACCM TCO

TCOname = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};

TCOdirectory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/TCO/'];
TCOfiles = dir([TCOdirectory,'TOZ*']);

yearstoexclude = 2015;
deviation_year = 2015;
vert = 0;
for i = 1:length(TCOfiles)
    
    [~, WACCMTCO.(TCOname{i}), ~] = ...
        Read_in_netcdf([TCOdirectory,TCOfiles(i).name]);
end

MAM_TCO_latmean = squeeze(nanmean(WACCMTCO.MAM.toz));  
VCMAM_TCO_latmean = squeeze(nanmean(WACCMTCO.VCMAM.toz));  
for i = 1:12
    MAM_TCO_latmean2(i,1:17,:) = MAM_TCO_latmean(:,i:12:end-6)';
    VCMAM_TCO_latmean2(i,1:17,:) = VCMAM_TCO_latmean(:,i:12:end-6)';
end

MAM_TCO_final = MAM_TCO_latmean2(:,6:end,:);
VCMAM_TCO_final = VCMAM_TCO_latmean2(:,6:end,:);
MAM_TCO_final(1:9,1,:) = NaN;
VCMAM_TCO_final(1:9,1,:) = NaN;
% creating deviations

Meanall = nanmean(MAM_TCO_final,2);
Meanall = repmat(Meanall,[1, size(MAM_TCO_final,2),1]);
if sa
    stdall = nanstd(MAM_TCO_final,0,2);    
    stdall = repmat(stdall,[1,size(MAM_TCO_final,2),1]);    
else
    stdall = 1;
end
MAMTCO_deviationnorm = squeeze(MAM_TCO_final-Meanall)./stdall;
%MAMTCO_deviationnorm(1:9,1,:) = NaN;

%MAMstd_OctNov = stdall(10:11,1,14);

Meanall = nanmean(VCMAM_TCO_final,2);
Meanall = repmat(Meanall,[1, size(VCMAM_TCO_final,2),1]);
if sa
    stdall = nanstd(VCMAM_TCO_final,0,2);    
    stdall = repmat(stdall,[1,size(VCMAM_TCO_final,2),1]);    
else
    stdall = 1;
end
VCMAMTCO_deviationnorm = squeeze(VCMAM_TCO_final-Meanall)./stdall;
%VCMAMTCO_deviationnorm(1:9,1,:) = NaN;
%VCMAMstd_OctNov = stdall(10:11,1,14);

% %% plotting line plots at 70S and normalized anomaly at 70S
% figure
% subplot(2,1,1)
% plot(squeeze(nanmean(OMI_monthly_latmean(:,10,10:40),3)));
% hold on
% plot(squeeze(nanmean(MAM_TCO_final(10,:,6:22),3)));
% plot(squeeze(nanmean(VCMAM_TCO_final(10,:,6:22),3)));
% subplot(2,1,2)
% plot(squeeze(nanmean(OMIdeviationnorm(:,10,10:40),3)));
% hold on
% plot(squeeze(nanmean(MAMTCO_deviationnorm(10,:,6:22),3)));
% plot(squeeze(nanmean(VCMAMTCO_deviationnorm(10,:,6:22),3)));


%%
% cbrew_vert = cbrewer('div','RdBu',20);
% cbrew_vert1 = flipud(cbrew_vert);
% 
% vertintervals = -3:.3:3;
% figure;
% set(gcf,'color','white','position',[100 100 1000 1000]);
% for i = 1:3;
%     subplot(3,1,i)
%     if i == 1
%         contourf(1:12,latitudes(1:60),squeeze(OMIdeviationnorm(end,:,1:60))',vertintervals)    
%         title('OMI','fontsize',22)
%     elseif i == 2
%         contourf(1:12,Latitudes(1:32),squeeze(MAMTCO_deviationnorm(:,end,1:32))',vertintervals)
%         title('MAM','fontsize',22)
%     elseif i == 3
%         contourf(1:12,Latitudes(1:32),squeeze(VCMAMTCO_deviationnorm(:,end,1:32))',vertintervals)
%         title('VC-MAM','fontsize',22);
%          xlabel('month','fontsize',22);         
%     %colorbar     
%     end
%     colormap(cbrew_vert1)
%     caxis([-3 3])
%      ylabel('Latitude','fontsize',22);
% end


%% importing CALIPSO
CALIPSO = load('/Users/kanestone/work/projects/CALIPSO/MATLABoutput/CALIPSO_MLS_compare.mat');
CALIPSO.datafinal_latinterp_presinterp = cat(3,ones(7,55,4)*-9999,CALIPSO.datafinal_latinterp_presinterp,ones(7,55,3)*-9999);
CALIPSO.datafinal_latinterp_presinterp (CALIPSO.datafinal_latinterp_presinterp  == -9999) = NaN;
% year, month, lat_band
if numberdensity
    MLS_at_pressure = permute(squeeze(MLS.O3_monthmeanND(:,:,:,presind)),[2,1,3]);
else
    MLS_at_pressure = permute(squeeze(MLS.O3_monthmean(:,:,:,presind)),[2,1,3]);
end

%% import WACCM
% WACCM data is in the same grid as MLS and has taken into account MLS
% averaging kernels
load('/Users/kanestone/work/projects/WACCM/MatlabOutput/WACCM_MLS_compare.mat');
runs = fields(WACCM_MolC_convolve);
for i = 1:length(runs)
    if numberdensity
        WACCM_at_pressure.(runs{i}) = squeeze(WACCM_ND_convolve.(runs{i})(:,:,presind,:));
        WACCM_at_pressure.(runs{i}) = permute(WACCM_at_pressure.(runs{i}),[2,1,3]);
    else
        WACCM_at_pressure.(runs{i}) = squeeze(WACCM_MolC_convolve.(runs{i})(:,:,presind,:));
        WACCM_at_pressure.(runs{i}) = permute(WACCM_at_pressure.(runs{i}),[2,1,3]);
    end
        SAD_SULFC.(runs{i}) = permute(SAD_SULFClatMLSperiod.(runs{i}),[2,1,3,4]);
        STS.(runs{i}) = permute(STSlatMLSperiod.(runs{i}),[2,1,3,4]);
        HNO3GAS.(runs{i}) = permute(HNO3GASlatMLSperiod.(runs{i}),[2,1,3,4]);
        ICE.(runs{i}) = permute(SAD_ICElatMLSperiod.(runs{i}),[2,1,3,4]);
        Extinctionrearrange.(runs{i}) = permute(Extinctionrearrange.(runs{i}),[2,1,3,4]);
        SAD_SULFC_at_pressure.(runs{i}) = squeeze(SAD_SULFC.(runs{i})(:,:,presind,:));
        STS_at_pressure.(runs{i}) = squeeze(STS.(runs{i})(:,:,presind,:));
        HNO3GAS_at_pressure.(runs{i}) = squeeze(HNO3GAS.(runs{i})(:,:,presind,:));
        
        % calculating STS/GAS ratio
        STSGASRatio.(runs{i}) = STS.(runs{i})./HNO3GAS.(runs{i});
        STSGASRatio_at_pressure.(runs{i}) = squeeze(STSGASRatio.(runs{i})(:,:,presind,:));
        
        if remove_tertiary
            SAD_SULFC_at_pressure.(runs{i}) (STSGASRatio_at_pressure.(runs{i}) > .1) = NaN;            
            SAD_SULFC.(runs{i}) (STSGASRatio.(runs{i}) > .1) = NaN;
            Extinctionrearrange.(runs{i}) (STSGASRatio.(runs{i})(:,end-1:end,:,:) > .1) = NaN;
            Extinctionrearrange.(runs{i}) (ICE.(runs{i})(:,end-1:end,:,:) > 1e-7) = NaN;
            Extinctionrearrange.(runs{i})(:,:,1:10,:) = NaN;
        end
        %
end

%%
%import WACCM sulfc
commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[~, ~, ~, ~, ~, ~, Latitudes, ~] = ReadWACCMvertical('SAD_SULFC','monthly',commondir,0);

%% constructing deviations

%allyears = 2004:2015;

includeVCMAM = 1;
includeMLS = 1;

MLS.O3_monthmeanND = permute(MLS.O3_monthmeanND,[2,1,4,3]);
MLS.O3_monthmeanND(7:12,13,:,:) = NaN;

yearstoexclude = [2015];
deviation_year = 2015;
vert = 1;
[MLS_vertical_deviation_final, MLS_vertical_deviation_final_percent, MLSdeviationnorm,stdallMLS] = constructdeviations(MLS.O3_monthmeanND,allyears,yearstoexclude,deviation_year,vert,sa);     
%MLS_vertical_deviation_final = MLS_vertical_deviation_final(monthstoplot(1):monthstoplot(1)+length(monthstoplot)-1,:,:);
%MLS_vertical_deviation_final_percent = MLS_vertical_deviation_final_percent(monthstoplot(1):monthstoplot(1)+length(monthstoplot)-1,:,:);
for i = 1:length(runs)
    WACCM_ND_convolve.(runs{i}) = permute(WACCM_ND_convolve.(runs{i}),[2 1 3 4]);
    WACCM_Temp_convolve.(runs{i}) = permute(WACCM_Temp_convolve.(runs{i}),[2 1 3 4]);
    TemplatMLSperiod.(runs{i}) = permute(TemplatMLSperiod.(runs{i}),[2 1 3 4]);
    %SAD_SULFC.(runs{i}) = permute(SAD_SULFC.(runs{i}),[2 1 3 4]);
    [WACCM_vertical_deviation_final(i,:,:,:), WACCM_vertical_deviation_final_percent(i,:,:,:), WACCMdeviationnorm(i,:,:,:,:),stdallMAM(i,:,:,:,:)] = constructdeviations(...
        WACCM_ND_convolve.(runs{i}),allyears,yearstoexclude,deviation_year,vert,sa);     
    [WACCM_vertical_Temp_deviation_final(i,:,:,:), WACCM_vertical_Temp_deviation_final_percent(i,:,:,:), WACCM_Temp_deviationnorm(i,:,:,:,:)] = constructdeviations(...
        TemplatMLSperiod.(runs{i}),allyears,yearstoexclude,deviation_year,vert,sa);     
    if sulfdeviation
        [SAD_SULFC_vertical_deviation_final(i,:,:,:), SAD_SULFC_vertical_deviation_final_percent(i,:,:,:), SULFCdeviationnorm(i,:,:,:,:)] = constructdeviations(...
            SAD_SULFC.(runs{i}),allyears,yearstoexclude,deviation_year,vert,sa);                     
    else
        SAD_SULFC_vertical_deviation_final(i,:,:,:) = cat(1,SAD_SULFC.(runs{i})(:,end-1,:,:),SAD_SULFC.(runs{i})(:,end,:,:));
        Extinction_vertical_deviation_final(i,:,:,:) = cat(1,Extinctionrearrange.(runs{i})(:,1,:,:),Extinctionrearrange.(runs{i})(:,2,:,:))./50.*(550/532).^2.27;
    end
end

deviationindex = find(allyears == deviation_year);
WACCM_vertical_deviation_final2 = squeeze(WACCMdeviationnorm(:,:,deviationindex,:,:));
MLSdeviationnorm_final = squeeze(MLSdeviationnorm(:,deviationindex,:,:));

% MAM100hPastd = stdallMAM(2,10:11,1,13,14);
% VCMAM100hPastd = stdallMAM(3,10:11,1,13,14);
% MLSstd100hPa = stdallMLS(10:11,1,13,14);

%WACCM_vertical_Temp_deviation_final2 = WACCM_vertical_Temp_deviation_final;
WACCM_vertical_Temp_deviation_final2 = squeeze(WACCM_Temp_deviationnorm(:,:,deviationindex,:,:));
WACCM_vertical_deviation_final_percent2 = WACCM_vertical_deviation_final_percent;
Volcanoes_vertical = squeeze((WACCM_vertical_deviation_final_percent2(2,:,:,:) - WACCM_vertical_deviation_final_percent2(3,:,:,:)));

WACCM_vertical_deviation_final2(:,:,1:7,:) = NaN;
WACCM_vertical_Temp_deviation_final2(:,:,1:7,:) = NaN;
MLS_vertical_deviation_final(:,1:7,:) = NaN;

cbrewtemp = cbrewer('seq','Reds',20);

% for i = 1:20
% plot(MLSdeviationnorm_final(:,13,28-i),'color',cbrewtemp(i,:))
% hold on
% end
% 
% figure
% for i = 1:20
% plot(squeeze(WACCMdeviationnorm(2,:,12,13,28 - i)),'color',cbrewtemp(i,:))
% hold on
% end
% 
% figure
% for i = 1:20
% plot(squeeze(WACCMdeviationnorm(3,:,12,13,28 - i)),'color',cbrewtemp(i,:))
% hold on
% end

%%
% figure;
% 
% plot(squeeze(MLS.O3_monthmeanND(10,:,13,14)))
% hold on
% plot(squeeze(WACCM_ND_convolve.MAM(10,:,13,14)))
% plot(squeeze(WACCM_ND_convolve.VCMAM(10,:,13,14)))
% title('Absolute at 100 hPa and 65?S');
% xlabel('Year');
% ylabel('ND')
% legend('MLS','MAM','VC-MAM');
% 
% figure;
% plot(squeeze((MLSdeviationnorm(10,:,13,14))))
% hold on
% plot(squeeze(WACCMdeviationnorm(2,10,:,13,14)))
% plot(squeeze(WACCMdeviationnorm(3,10,:,13,14)))
% 
% title('Absolute anomalies at 100 hPa and 65?S');
% xlabel('Year');
% ylabel('ND')
% legend('MLS','MAM','VC-MAM');

%% percent differences
MAM_change_TCO = (MAM_TCO_final - VCMAM_TCO_final);
MAM_percentchange_TCO = (MAM_TCO_final - VCMAM_TCO_final)./VCMAM_TCO_final;
MAM_TCOdifference_toplot_percent = squeeze(MAM_percentchange_TCO(5:12,12,1:31))'*100;
MAM_TCOdifference_toplot = squeeze(MAM_change_TCO(5:12,12,1:31))';
MAM_change_atlayer =  (WACCM_ND_convolve.MAM - WACCM_ND_convolve.VCMAM);
MAM_percentchange_atlayer =  (WACCM_ND_convolve.MAM - WACCM_ND_convolve.VCMAM)./WACCM_ND_convolve.VCMAM*100;
MAM_layerdifference_toplot_percent = squeeze(MAM_percentchange_atlayer(5:12,12,13,1:31));
MAM_layerdifference_toplot = squeeze(MAM_change_atlayer(5:12,12,13,1:31))./1e11;

%% plotting differences
cmap = flipud(cbrewer('seq','YlOrBr',10));
createfig('medium','on');
sp(1) = subplot(2,2,1);
sp_pos(1,:) = get(sp(1),'position');
[~,h] = contourf(1:8,Latitudes(1:31),MAM_layerdifference_toplot',-10:1:0);
colormap(cmap);
ch = colorbar;
caxis([-10 0]);
set(ch,'Ticks',-10:1:0)                                     
set(get(ch,'ylabel'),'string','{\times} 10^1^1 molecules/cm^3','fontsize',20)
set(gca,'fontsize',16);
set(gca,'xtick',1:1:8,'xticklabel',{'M','J','J','A','S','O','N','D'});
ylabel('Latitude ({\circ}N)','fontsize',20);
title('MAM-VCMAM at 100 hPa','fontsize',22);

sp(2) = subplot(2,2,2)
sp_pos(2,:) = get(sp(2),'position');
[~,h] = contourf(1:8,Latitudes(1:31),MAM_TCOdifference_toplot,-20:2:0);
colormap(cmap);
ch = colorbar;
caxis([-20 0]);
set(get(ch,'ylabel'),'string','DU','fontsize',20)
set(ch,'Ticks',-20:2:0)
set(gca,'fontsize',16);
set(gca,'xtick',1:1:8,'xticklabel',{'M','J','J','A','S','O','N','D'});
title('TCO MAM-VCMAM','fontsize',22);

sp(1) = subplot(2,2,3);
sp_pos(1,:) = get(sp(1),'position');
[~,h] = contourf(1:8,Latitudes(1:31),MAM_layerdifference_toplot_percent',-80:8:0);
colormap(cmap);
ch = colorbar;
caxis([-80 0]);
set(get(ch,'ylabel'),'string','Percent','fontsize',20)
set(ch,'Ticks',-80:8:0)                                     
set(gca,'fontsize',16);
set(gca,'xtick',1:1:8,'xticklabel',{'M','J','J','A','S','O','N','D'});
xlabel('Month','fontsize',20);
ylabel('Latitude ({\circ}N)','fontsize',20);
title('Percent change at 100 hPa','fontsize',22);

sp(2) = subplot(2,2,4);
sp_pos(2,:) = get(sp(2),'position');
[~,h] = contourf(1:8,Latitudes(1:31),MAM_TCOdifference_toplot_percent,-10:1:0);
colormap(cmap);
ch = colorbar;
caxis([-10 0]);
set(get(ch,'ylabel'),'string','Percent','fontsize',20)
set(ch,'Ticks',-10:1:0)                                     
set(gca,'fontsize',16);
set(gca,'xtick',1:1:8,'xticklabel',{'M','J','J','A','S','O','N','D'});
xlabel('Month','fontsize',20);
title('TCO percent change','fontsize',22);

filename = '/Users/kanestone/Dropbox/Work_Share/MITwork/CalbucoPaper/NewFigures/MAMminusVCMAM/MAMminusVCMAM.pdf';
export_fig(filename,'-pdf');
%% layer deviations
if deviationlayer
    
    pressureplot = 100;  
    include_sulf = 0;
    layer_deviation(MLSdeviationnorm_final,WACCM_vertical_deviation_final2,...
        deviation_year,pressureplot,OMIdeviationnorm,MAMTCO_deviationnorm,VCMAMTCO_deviationnorm,...
        Extinction_vertical_deviation_final,...
        CALIPSO.datafinal_latinterp_presinterp,MLSPressure,Latitudes,latitudes,include_sulf,sa);
end
    
lats = [-90 -85; -85 -80; -80 -75; -75 -70; -70 -65; -65 -60; -60 -55; -55 -50;...
                -50 -45; -45 -40; -40 -35; -35 -30; -30 -25];

%% Differences between MLS and WACCM temperature
temp_diff = 0;
if temp_diff
    diffyear = 2005;
    diffind = find(allyears == diffyear);

Temp_diff = WACCM_Temp_convolve.MAM - MLS.temp_monthmean_vert;
% plotting
prestick = [1000,500,200,100,50,20,10,5,2,1,.5,.2,.1];
logprestick = log(prestick);
posind = [.05,.02,.01];    

cbrew_vert = cbrewer('div','RdBu',11);
cbrew_vert1 = flipud(cbrew_vert(:,:));
tempdiff_intervals = -11:2:11;

fig = figure;   
set(fig,'color','white','position',[100 100 500 1000])

monthplot = [9,10,11];

for i = 1:3;
    
    sp = subplot(3,1,i);  
    sp_pos(i,:) = get(sp,'position');
    MERRA = contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),squeeze(Temp_diff(diffind,monthplot(i),:,:)),tempdiff_intervals,...
        'LineStyle','--','LineColor',[.5 .5 .5]);                    
    colormap(cbrew_vert1);
    caxis([min(tempdiff_intervals) max(tempdiff_intervals)]);
    set(gca,'YDir','reverse');    
    set(gca,'ytick',fliplr(log(prestick)),'yticklabel',fliplr(prestick),'fontsize',fsize);
    ylim([log(1) log(300)]);
    set(gca,'color',[.8 .8 .8]);

    ylabel('Pressure (hPa)','fontsize',fsize+2);
    title(daytitle{monthplot(i)},'fontsize',fsize+4);

    set(sp,'position',[sp_pos(i,1),sp_pos(i,2)-i*posind(i),sp_pos(i,3)/1.15,sp_pos(i,4)]);
    sp_pos(i,:) = get(sp,'position');
    
    if i == 3
        xlabel(['Latitude (',char(176),'N)'],'fontsize',fsize+2);

        ch = colorbar;
        colormap(cbrew_vert1);           
        caxis([min(tempdiff_intervals) max(tempdiff_intervals)]);
        set(ch,'Ticks',tempdiff_intervals)
        set(get(ch,'ylabel'),'string','Degrees (K)','fontsize',fsize+2)
        cbaxloc = get(ch,'Position');
        subpos = get(sp,'OuterPosition');                    

        set(ch,'Position',[sp_pos(i,1)+sp_pos(i,3)+sp_pos(i,3)/30, cbaxloc(2),...
            cbaxloc(3), sp_pos(1,2)+sp_pos(1,4) - sp_pos(3,2)],...
            'Box','on','YAxisLocation','right','fontsize',fsize);

        set(ch,'YaxisLocation','right','fontsize',fsize);                
    end                   
end    
    mtit([' WACCM-MLS temperature difference ', num2str(diffyear)],'xoff',-.01,'yoff',.06,'fontsize',fsize);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
        num2str(diffyear),'_','WACCM_minus_MLS_temp_anomaly','_',daytitle{monthplot(1)},'-',daytitle{monthplot(end)},'.png'];
    export_fig(filename,'-png','-nofontswap');  
    
end

          
            
%% calculating vertical deviation


if vertical_deviation
    

%% Temperature deviations

if temperature_deviations
    
    prestick = [1000,500,200,100,50,20,10,5,2,1,.5,.2,.1];
    logprestick = log(prestick);
    posind = [.05,.02,.01];    
    
    cbrew_vert = cbrewer('div','RdBu',13);
    cbrew_vert1 = flipud(cbrew_vert);
    temp_intervals = -2.6:.4:2.6;
    
    fig = figure;   
    set(fig,'color','white','position',[100 100 500 1000])
 
    monthplot = [9,10,11];
    
    for i = 1:3;
        sp = subplot(3,1,i);  
        sp_pos(i,:) = get(sp,'position');
        MERRA = contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),squeeze(WACCM_vertical_Temp_deviation_final2(3,monthplot(i),:,:)),temp_intervals,...
            'LineStyle','--','LineColor',[.5 .5 .5]);                    
        colormap(cbrew_vert1);
        caxis([min(temp_intervals) max(temp_intervals)]);
        set(gca,'YDir','reverse');    
        set(gca,'ytick',fliplr(log(prestick)),'yticklabel',fliplr(prestick),'fontsize',fsize);
        ylim([log(10) log(300)]);
        set(gca,'color',[.8 .8 .8]);
        
        ylabel('Pressure (hPa)','fontsize',fsize+2);
        title(daytitle{monthplot(i)},'fontsize',fsize+4);
        
        set(sp,'position',[sp_pos(i,1),sp_pos(i,2)-i*posind(i),sp_pos(i,3)/1.15,sp_pos(i,4)]);
        sp_pos(i,:) = get(sp,'position');
        
        if i == 3
            xlabel(['Latitude (',char(176),'N)'],'fontsize',fsize+2);
            
            ch = colorbar;
            colormap(cbrew_vert1);           
            caxis([min(temp_intervals) max(temp_intervals)]);
            set(ch,'Ticks',temp_intervals)
            set(get(ch,'ylabel'),'string','Normalized anomaly','fontsize',fsize+2)
            cbaxloc = get(ch,'Position');
            subpos = get(sp,'OuterPosition');                    

            set(ch,'Position',[sp_pos(i,1)+sp_pos(i,3)+sp_pos(i,3)/30, cbaxloc(2),...
                cbaxloc(3), sp_pos(1,2)+sp_pos(1,4) - sp_pos(3,2)],...
                'Box','on','YAxisLocation','right','fontsize',fsize);

            set(ch,'YaxisLocation','right','fontsize',fsize);    
            
        end                   
        
    end
    mtit([num2str(deviation_year),' anomaly from 2004-2014 mean'],'xoff',-.01,'yoff',.06,'fontsize',fsize+4);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
        num2str(deviation_year),'_','MLS_WACCM_temp_anomaly','_',daytitle{monthplot(1)},'-',daytitle{monthplot(end)},'.pdf'];
    export_fig(filename,'-pdf','-nofontswap');   
end

%% plotting lines plots of normalised

montitle = {'January','February','March','April','May','June','July','August','September',...
       'October','November','December'};

latnorm = 17;
monthnorm = 9;
presnorm = 11;

WACMAMtemp = squeeze(nanmean(squeeze(WACCM_ND_convolve.MAM(:,1:12,presnorm,latnorm)),3));
WACVCMAMtemp = squeeze(nanmean(squeeze(WACCM_ND_convolve.VCMAM(:,1:12,presnorm,latnorm)),3));
MLStemp = squeeze(nanmean(squeeze(MLS.O3_monthmeanND(:,1:12,presnorm,latnorm)),3));

W = cosd(Latitudes(latnorm));
MAMweightedmean = zeros(size(WACMAMtemp,1),size(WACMAMtemp,2)-1);
VCMAMweightedmean = zeros(size(WACMAMtemp,1),size(WACMAMtemp,2)-1);
MLSweightedmean = zeros(size(WACMAMtemp,1),size(WACMAMtemp,2)-1);
for i = 1:size(WACMAMtemp,1)
    for j = 1:size(WACMAMtemp,2)
        tempMAM = squeeze(WACMAMtemp(i,j,:));
        tempVCMAM = squeeze(WACVCMAMtemp(i,j,:));
        tempMLS = squeeze(MLStemp(i,j,:));
        MAMweightedmean(i,j) = nansum(tempMAM(:).*W(:))./nansum(W(:));   
        VCMAMweightedmean(i,j) = nansum(tempVCMAM(:).*W(:))./nansum(W(:));   
        MLSweightedmean(i,j) = nansum(tempMLS(:).*W(:))./nansum(W(:));   
    end
end
MAMweightedmean (MAMweightedmean == 0) = NaN;
VCMAMweightedmean (VCMAMweightedmean == 0) = NaN;
MLSweightedmean (MLSweightedmean == 0) = NaN;

MeanallMAM = repmat(nanmean(MAMweightedmean,2),[1,12]);
MeanallVCMAM = repmat(nanmean(VCMAMweightedmean,2),[1,12]);
MeanallMLS = repmat(nanmean(MLSweightedmean,2),[1,12]);

StdallMAM = repmat(nanstd(MAMweightedmean,0,2),[1,12]);
StdallVCMAM = repmat(nanstd(VCMAMweightedmean,0,2),[1,12]);
StdallMLS = repmat(nanstd(MLSweightedmean,0,2),[1,12]);

MAMfinal = squeeze(MAMweightedmean-MeanallMAM)./StdallMAM;
VCMAMfinal = squeeze(VCMAMweightedmean-MeanallVCMAM)./StdallVCMAM;
MLSfinal = squeeze(MLSweightedmean-MeanallMLS)./StdallMLS;

%WACCMMAM_area = nanmean(squeeze(WACCMdeviationnorm(2,:,:,presnorm,latnorm)),4);
%WACCMMAM_area = nanmean(WACCMMAM_area,3);

%WACCMVCMAM_area = nanmean(squeeze(WACCMdeviationnorm(3,:,:,presnorm,latnorm)),4);
%WACCMVCMAM_area = nanmean(WACCMVCMAM_area,3);

%MLS_area = nanmean(squeeze(MLSdeviationnorm(:,:,presnorm,latnorm)),4);
%MLS_area = nanmean(MLS_area,3);

matstand = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

cbrew = cbrewer('qual','Paired',12);

fig = figure;
set(fig,'color','white','position',[100 100 700 1000]);

for i = 1:3
    
    subplot(3,1,i);
    plot(squeeze(MLSfinal(monthnorm,:)),'LineWidth',lwidth,'color',matstand(3,:),'LineStyle','-');    
    hold on
    plot(squeeze(MAMfinal(monthnorm,:)),'LineWidth',lwidth,'color',cbrew(2,:),'LineStyle','--');
    plot(squeeze(VCMAMfinal(monthnorm,:)),'LineWidth',lwidth,'color',matstand(2,:),'LineStyle','-.');
    
    setgca(fsize,'Normalized anomaly','Year',[0,0,1],[1,1,1],{montitle{monthnorm}},i)
    set(gca,'xtick',1:20,'xticklabel',{'2004',[],'2006',[],'2008',[],'2010',[],'2012',[],'2014',[],'2016',[]});
    ylim([-2.4 2.4]);
    xlim([0 14]);
    monthnorm = monthnorm+1;
    
end
mtit(['Normalized anomaly at 60{\circ}S and 150 hPa'],'xoff',-.01,'yoff',.04,'fontsize',fsize+2);
lh = legend('MLS','MAM','VC-MAM');
set(lh,'box','off','fontsize',fsize-4,'location','NorthEast');
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
        'Lineplots_normalized_60_150hPa','.pdf'];
export_fig(filename,'-pdf');
% figure;
% plot(squeeze(WACCM_ND_convolve.MAM(monthnorm,:,presnorm,latnorm)))
% hold on
% plot(squeeze(WACCM_ND_convolve.VCMAM(monthnorm,:,presnorm,latnorm)))
% plot(squeeze(MLS.O3_monthmeanND(monthnorm,:,presnorm,latnorm)))

 prestick = [1000:-100:300,200,150,100,90:-10:10,9:-1:3,2.5,1.5,1,.9:-.1:.1];
    presticklabel = {1000,[],[],[],[],500,[],[],200,150,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
        5,[],[],2.5,1.5,1,[],[],[],[],.5,[],[],.2,.1};
    logprestick = log(prestick);

%% plotting differnces between MAM and VC-MAM normalized anomalies
nadiff = squeeze(WACCM_vertical_deviation_final2(2,:,:,:)-WACCM_vertical_deviation_final2(3,:,:,:));
%nadiff = squeeze(WACCM_vertical_deviation_final2(2,:,:,:))-MLSdeviationnorm_final;
%nadiff = MLSdeviationnorm_final-squeeze(WACCM_vertical_deviation_final2(3,:,:,:));
nadifftit = {'September','October','November'};
cbrewfordiff = cbrewer('div','RdBu',27);
cbrewfordiff = flipud(cbrewfordiff(11:end,:));
diffintervals = -2.5:.2:.9;
figure;
for i = 1:3
    
    spdiff =  subplot(3,1,i);
    sp_pos_diff(i,:) = get(spdiff,'position');
    set(gcf,'color','white','position',[100 100 1200 1000]);
    contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),squeeze(nadiff(8+i,:,:)),...
        diffintervals);
    xlim([-90 -45])
    set(gca,'YDir','reverse');    
        set(gca,'ytick',fliplr(log(prestick)),'yticklabel',fliplr(presticklabel),'fontsize',fsize);
        ylim([log(10) log(300)]);
    title(nadifftit{i});
    ylabel('Pressure(hPa)');
    if i == 3;
        xlabel('Latitude');
    end
if i == 3
    ch = colorbar;
    colormap(cbrewfordiff);           
    caxis([min(diffintervals) max(diffintervals)]);
    set(ch,'YTick',diffintervals,'YTickLabel',diffintervals)
    %set(get(ch,'ylabel'),'string','\times 10^{11} molecules/cm^3','fontsize',fsize)
    set(get(ch,'ylabel'),'string','Normalized anomaly difference','fontsize',fsize)
    cbaxloc = get(ch,'Position');
    subpos = get(spdiff,'OuterPosition');                    

    set(ch,'Position',[sp_pos_diff(i,1)+sp_pos_diff(i,3)+sp_pos_diff(i,3)/60, cbaxloc(2),...
        cbaxloc(3), sp_pos_diff(1,2)+sp_pos_diff(1,4)-sp_pos_diff(3,2)],...
        'Box','on','YAxisLocation','right','fontsize',fsize);

    set(ch,'YaxisLocation','right','fontsize',fsize);                  
end

end

mtit('MAM - VC-MAM normalized anomaly','xoff',-.01,'yoff',.04,'fontsize',fsize+4);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
        'MAMminusVCMAM_nomalized_anomaly','.png'];
export_fig(filename,'-png');

%% plotting vertical deviations
    %monthplot = [6,7,8,9,10,11];
    %monthplot = [9,9,10,10,11,11];
    monthplot = [9,9,9,10,10,10,11,11,11];
    %monthplot = [10,10,10,11,11,11,12,12,12];
    monthplotCAL = [5,5,5,6,6,6,7,7,7];
    %daytitle = {'MAM-September','VCMAM-September','MAM-October','VCMAM-October','MAM-November','VCMAM-November'};
    %daytitle = {'MLS-September','MAM-September','MLS-October','MAM-October','MLS-November','MAM-November'};
    daytitle = {'MLS - September','MAM - September','VC-MAM - September','MLS - October',...
        'MAM - October','VC-MAM - October','MLS - November','MAM - November','VC-MAM - November'};
    montitle = {'January','February','March','April','May','June','July','August','September',...
       'October','November','December'};
    latindex = find(Latitudes >= -90 & Latitudes <= -30);
%     prestick = [1000:-100:300,250,150,100,90:-10:10,9:-1:3,2.5,1.5,1,.9:-.1:.1];
%     presticklabel = {1000,[],[],[],[],500,[],[],250,150,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
%         5,[],[],2.5,1.5,1,[],[],[],[],.5,[],[],.2,.1};           
    cbrew_vert = cbrewer('div','RdBu',10);
    cbrew_vert1 = flipud(cbrew_vert);
    fig = figure;
    set(fig,'color','white','position',[100 100 1500 1000])
    %vertintervals = -1:.2:1;
    vertintervals = -3:.6:3;
    
    if includeMLS && includeVCMAM
        nosublots = 9;
    else
        nosublots = 6;
    end
    for i = 1:nosublots;
    % plotting
        if includeMLS && includeVCMAM
            sp = subplot(3,3,i);
        else
            sp = subplot(3,2,i);
        end
        sp_pos(i,:) = get(sp,'position');
        sp_pos_outer(i,:) = get(sp,'OuterPosition');
        if includeVCMAM == 1 && includeMLS == 0
            if mod(i,2)
                contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),WACCM_vertical_deviation_final2(2,:,:,monthplot(i))',vertintervals,...
                    'LineStyle','--','LineColor',[.5 .5 .5]);    
            else contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),WACCM_vertical_deviation_final2(3,:,monthplot(i))',vertintervals,...%-2.3e12:2e11:9e11
                    'LineStyle','--','LineColor',[.5 .5 .5]);    
            end
        elseif includeMLS == 1 && includeVCMAM == 0
            if mod(i,2)
                contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),squeeze(MLS_vertical_deviation_final(monthplot(i),:,:)),vertintervals ,...
                    'LineStyle','--','LineColor',[.5 .5 .5]);    
            else contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),squeeze(WACCM_vertical_deviation_final2(2,monthplot(i),:,:)),vertintervals ,...
                'LineStyle','--','LineColor',[.5 .5 .5]);            
            end        
        elseif includeMLS && includeVCMAM
            if i == 1 || i == 4 || i == 7
                [a,c] = contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),squeeze(MLSdeviationnorm_final(monthplot(i),:,:)),vertintervals,...
                    'LineStyle','--','LineColor',[.5 .5 .5]);    
            elseif i == 2 || i == 5 || i == 8
                contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),squeeze(WACCM_vertical_deviation_final2(2,monthplot(i),:,:)),vertintervals,...
                'LineStyle','--','LineColor',[.5 .5 .5]);            
            elseif i == 3 || i == 6 || i == 9
                contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),squeeze(WACCM_vertical_deviation_final2(3,monthplot(i),:,:)),vertintervals,...
                'LineStyle','--','LineColor',[.5 .5 .5]);            
            end
        else
            contourf(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),squeeze(WACCM_vertical_deviation_final2(2,monthplot(i),:,:)),vertintervals,...
                'LineStyle','--','LineColor',[.5 .5 .5]);            
        end
        colormap(cbrew_vert1);
        caxis([min(vertintervals) max(vertintervals)]);
        set(gca,'YDir','reverse');    
        set(gca,'ytick',fliplr(log(prestick)),'yticklabel',fliplr(presticklabel),'fontsize',fsize);
        ylim([log(10) log(300)]);

        set(gca,'color',[.8 .8 .8]);

        if ~includeMLS && ~ includeVCMAM
            if i == 1 || i == 3 || i == 5
                ylabel(['Pressure (hPa)'],'fontsize',fsize+2);
            end
            if i == 5 || i == 6
                xlabel(['Latitude (',char(176),'N)'],'fontsize',fsize+2);
            end
            title(daytitle{i},'fontsize',fsize+4);

            if ~mod(i,2)
                set(sp,'position',[sp_pos(i,1)-.07 sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)]);
            else 
                set(sp,'position',[sp_pos(i,1)-.035 sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)]);
            end 
        else
            if i == 1 || i == 4 || i == 7
                ylabel(['Pressure (hPa)'],'fontsize',fsize+2);
            end
            if i == 7 || i == 8 || i == 9;
                xlabel(['Latitude (',char(176),'N)'],'fontsize',fsize+2);
            end
            title(daytitle{i},'fontsize',fsize+4);

            
            set(sp,'position',[sp_pos(i,1) sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)-.02]);
            
            sp_pos(i,:) = get(sp,'position');
            
            if i == 1 || i == 4 || i == 7
                set(sp,'position',[sp_pos(i,1)-.02 sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)]);
            elseif i == 2 || i == 5 || i == 8 
                set(sp,'position',[sp_pos(i,1)-.04 sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)]);
            elseif i == 3 || i == 6 || i == 9 
                set(sp,'position',[sp_pos(i,1)-.06 sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)]);
            end 
            
            sp_pos2(i,:) = get(sp,'Position');
            if i == 1 || i == 2 || i == 3        
                set(sp,'position',[sp_pos2(i,1) sp_pos2(i,2)-.010 sp_pos2(i,3) sp_pos2(i,4)]);            
            elseif i == 4 || i == 5 || i == 6
                set(sp,'position',[sp_pos2(i,1) sp_pos2(i,2)+.02 sp_pos2(i,3) sp_pos2(i,4)]);
            elseif i == 7 || i == 8 || i == 9
                set(sp,'position',[sp_pos2(i,1) sp_pos2(i,2)+.05 sp_pos2(i,3) sp_pos2(i,4)]);
            end            
        end
        set(gca,'LineWidth',1.5);
        

        outpos(i,:) = get(sp,'outerposition');

        sp_pos(i,:) = get(sp,'position'); 

        if i == 9
            ch = colorbar;
            colormap(cbrew_vert1);           
            caxis([min(vertintervals) max(vertintervals)]);
            set(ch,'YTick',vertintervals,'YTickLabel',vertintervals)
            %set(get(ch,'ylabel'),'string','\times 10^{11} molecules/cm^3','fontsize',fsize)
            set(get(ch,'ylabel'),'string','Normalized anomaly','fontsize',fsize)
            cbaxloc = get(ch,'Position');
            subpos = get(sp,'OuterPosition');                    

            set(ch,'Position',[sp_pos(i,1)+sp_pos(i,3)+sp_pos(i,3)/15, cbaxloc(2),...
                cbaxloc(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(9,2)],...
                'Box','on','YAxisLocation','right','fontsize',fsize);

            set(ch,'YaxisLocation','right','fontsize',fsize);                  
        end

        if include_SAD_SULFC 
%             if i == 6
%                 tic;
%                 toc;
%             end
            %cmapWACCM = repmat([166/255,86/255,40/255],8,1);
            cmapWACCM = cbrewer('seq','YlOrBr',12);
            cmapWACCM = cmapWACCM(2:10,:);
            if ext                
                sulfurtoplot = Extinction_vertical_deviation_final;
                %cmapWACCM = cbrewer('seq','YlOrBr',10);
                cmapWACCM = cbrewer('seq','YlOrBr',13);
                %cmapWACCM = cmapWACCM(2:13,:);
                cmapWACCM = cmapWACCM(3:12,:);
                %sulf_intervals = 0:.625e-6:5e-6;
                sulf_intervals = 0:.8e-8:8e-8;
                exttimes = 1e8;
            elseif ~ext && sulfdeviation
                sulfurtoplot = SAD_SULFC_vertical_deviation_final;
                cmapWACCM = cbrewer('seq','YlOrBr',10);
                cmapWACCM = cmapWACCM(2:10,:);
                sulf_intervals = -2.2e-8:2e-8:11.8e-8;
            elseif ~ext && ~sulfdeviation
                sulfurtoplot = SAD_SULFC_vertical_deviation_final;
                cmapWACCM = cbrewer('seq','YlOrBr',10);
                cmapWACCM = cmapWACCM(2:9,:);
                sulf_intervals = -2.2e-8:2e-8:13.8e-8;
                exttimes = 1e8;
            end
                        
            haxes1 = gca;
            haxes1_pos(i,:) = get(haxes1,'Position');
            haxes2 = axes('Position',haxes1_pos(i,:),...
                      'XAxisLocation','bottom',...
                      'YAxisLocation','left',...
                      'Color','none');
            hold on        
            if i == 2 || i == 5 || i == 8
                contour(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),...
                    squeeze(sulfurtoplot(2,monthplot(i),:,:)),sulf_intervals,'LineWidth',1.5);%,'ShowText','on');
            elseif i == 3 || i == 6 || i == 9  
                contour(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),...
                    squeeze(sulfurtoplot(3,monthplot(i),:,:)),sulf_intervals,'LineWidth',1.5);%,'ShowText','on');
            elseif i == 1 || i == 4 || i == 7  
                contour(Latitudes(1:size(MLS_at_pressure,3)),log(MLSPressure),...
                    squeeze(CALIPSO.datafinal_latinterp_presinterp(monthplotCAL(i),:,:)),sulf_intervals,'LineWidth',1.5);%,'ShowText','on');
            end
            set(gca,'YDir','reverse');    
            set(gca,'ytick',fliplr(log(prestick)),'yticklabel',fliplr(presticklabel),'fontsize',fsize);
            ylim([log(10) log(300)]);
            %xlim([-90 -27.5])
            if i == 9
                set(gca,'position',[haxes1_pos(3,1),haxes1_pos(8,2),haxes1_pos(3,3),haxes1_pos(8,4)]);
            end
            colormap(haxes2,cmapWACCM)
            %WACCMcount = WACCMcount+1;
            caxis([min(sulf_intervals) max(sulf_intervals)]);
            if i == 9;
                ch1 = colorbar('south');
                set(ch1,'ytick',sulf_intervals,'yticklabel',sulf_intervals*exttimes)
                cbaxloc1 = get(ch1,'Position');
                if ext                    
                    set(get(ch1,'ylabel'),'string','CALIPSO 532 nm backscatter and WACCM 532 nm backscatter estimation (\times 10^{-8} m^-^1)','fontsize',fsize);
                else
                    set(get(ch1,'ylabel'),'string','\times 10^{-8} WACCM absolute 2015 SAD SULFC (cm^2/cm^3)','fontsize',fsize);                    
                end
                set(ch1, 'YAxisLocation','bottom')

                %[subpos(1)+subpos(3)+(subpos(1)+subpos(3))/25, cbaxloc(2),...
                %    cbaxloc(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(9,2)]
                %set(ch1,'position',[.109 .06 .740 .02]);                
                set(ch1,'position',[sp_pos(1,1), .06, sp_pos(9,1)+sp_pos(1,3)-sp_pos(1,1), cbaxloc1(4)]);                
                set(ch1,'Ticks',sulf_intervals)             
            end
        end

    end
    if ext
        ext_append = 'ext1020nm_absolute';
    elseif ~ext && sulfdeviation
        ext_append = 'SAD_SULFC_deviation';
    elseif ~ext && ~sulfdeviation
        ext_append = 'SAD_SULFC_absolute';
    end    
        
    if includeMLS
    mtit([num2str(deviation_year),' anomaly from 2004-2014 mean'],'xoff',-.01,'yoff',.06,'fontsize',fsize+8);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
        num2str(deviation_year),'_MLS_WACCM_O3_anomaly_',ext_append,'_',montitle{monthplot(1)},'-',montitle{monthplot(end)},'.pdf'];
    filename2 = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
        num2str(deviation_year),'_MLS_WACCM_O3_anomaly_',ext_append,'_',montitle{monthplot(1)},'-',montitle{monthplot(end)},'.png'];
    else
    mtit('WACCM 2015 anomaly from 2004-2014 mean','xoff',-.01,'yoff',.04,'fontsize',fsize+8);
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
            'WACCM_O3_anomaly_',daytitle{monthplot(1)},'-',daytitle{monthplot(end)},'.png'];
    end


    export_fig(filename2,'-png','-nofontswap');   
    export_fig(filename,'-pdf','-nofontswap');   
end


%% line plots
if line_plots
    for i = 1:12
        createfig('medium','on');
        for j = 1:6 
            sp = subplot(3,2,j);
            sp_pos = get(sp,'position');
            mlsph = plot(allyears,MLS_at_pressure(i,:,j+1),'LineWidth',lwidth);
            hold on
            mamph = plot(allyears,WACCM_at_pressure.MAM(i,:,j+1),'LineWidth',lwidth);
            vcmamph = plot(allyears,WACCM_at_pressure.VCMAM(i,:,j+1),'LineWidth',lwidth);
            
            set(gca,'fontsize',fsize-2,'xtick',min(allyears):4:max(allyears));
            xlim([min(allyears)-1 max(allyears)+1]);
            
            %moving subplots
            if mod(j,2)
                set(sp,'position',[sp_pos(1)+.02,sp_pos(2:4)]);
            else
                set(sp,'position',[sp_pos(1)-.02,sp_pos(2:4)]);
            end
            
            if j == 1 || j == 3 || j == 5
                ylabel('vmr','fontsize',fsize+2);
            end
            if j == 5 || j == 6
                xlabel('Year','fontsize',fsize+2);
            end
            title([num2str(lats(j+1,1)),' to ',num2str(lats(j+1,2)),'N'],'fontsize',fsize+4);
            if j == 6
                lh = legend('MLS','MAM','VC-MAM');
                set(lh,'orientation','horizontal','position',[.5 .01 .05 .05],'fontsize',fsize+2,...
                    'box','off');
            end
        end
        mtit([daytitle{i},' at ',num2str(pres),'hPa'],'xoff',-.01,'yoff',.04,'fontsize',fsize+4);
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/Obs_comparison/',...
            'MLScomparison_withAK/',sprintf('%02d',i),'_','WACCM_MLS_comparison_at_',...
            num2str(pres),'hPa','.pdf'];
        export_fig(filename,'-pdf');
    end    
end

%% vertical plots
if vertical_plots
    % monthly mean vertical profiles through time
    meanall = 0;
    indyear = 2011;
    indyear_ind = find(allyears == indyear);
    if meanall
        MLS_vert_mean = squeeze(nanmean(MLS.O3_monthmeanND_vert));
        MAM_vert_mean = squeeze(nanmean(WACCM_ND_convolve.MAM));
        VCMAM_vert_mean = squeeze(nanmean(WACCM_ND_convolve.VCMAM));
    else
        MLS_vert_mean = squeeze(MLS.O3_monthmeanND_vert(indyear_ind ,:,:,:));
        MAM_vert_mean = squeeze(WACCM_ND_convolve.MAM(indyear_ind ,:,:,:));
        VCMAM_vert_mean = squeeze(WACCM_ND_convolve.VCMAM(indyear_ind ,:,:,:));
    end
    
    MAM_vert_mean(:,1:8,:) = NaN;
    VCMAM_vert_mean(:,1:8,:) = NaN;
    
    prestick = [1000,500,200,100,50,20,10,5,2,1,.5,.2,.1];
    logprestick = log(prestick);
    
     for i = 1:12
        fig = figure;
        set(fig,'color','white','position',[100 100 900 900]);
        for j = 1:6 
            sp = subplot(3,2,j);
            sp_pos = get(sp,'position');
            mlsph = plot(squeeze(MLS_vert_mean(i,:,j+1)),log(MLSPressure),'LineWidth',lwidth);
            hold on
            mamph = plot(squeeze(MAM_vert_mean(i,:,j+1)),log(MLSPressure),'LineWidth',lwidth);
            vcmmamph = plot(squeeze(VCMAM_vert_mean(i,:,j+1)),log(MLSPressure),'LineWidth',lwidth);
            set(gca,'Ydir','Reverse');
            set(gca,'ytick',fliplr(logprestick),'yticklabel',fliplr(prestick));
            set(gca,'fontsize',fsize-2);
            ylim([log(.9) log(1000)]);
            
            %moving subplots
            if mod(j,2)
                set(sp,'position',[sp_pos(1)+.02,sp_pos(2:4)]);
            else
                set(sp,'position',[sp_pos(1)-.02,sp_pos(2:4)]);
            end
            
            if j == 1 || j == 3 || j == 5
                ylabel('Pressure (hPa)','fontsize',fsize+2);
            end
            if j == 5 || j == 6
                xlabel('molecules/cm^3','fontsize',fsize+2);
            end
            title([num2str(lats(j+1,1)),' to ',num2str(lats(j+1,2)),' hPa'],'fontsize',fsize+4);
            if j == 6
                lh = legend('MLS','MAM','VC-MAM');
                set(lh,'orientation','horizontal','position',[.5 .01 .05 .05],'fontsize',fsize+2,...
                    'box','off');
            end
        end
        if meanall
            mtit([daytitle{i},' ',num2str(allyears(1)),'-',num2str(allyears(end)),' average'],'xoff',-.01,'yoff',.04,'fontsize',fsize+4);
            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/Obs_comparison/',...
                'MLScomparison_withAK/profiles/',sprintf('%02d',i),'_','WACCM_MLS_profiles_at_',...
                num2str(lats(j+1,1)),'to',num2str(lats(end,2)),'.pdf'];            
        else
            mtit([daytitle{i},' ',num2str(allyears(indyear_ind))],'xoff',-.01,'yoff',.04,'fontsize',fsize+4);
            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/Obs_comparison/',...
                'MLScomparison_withAK/profiles/',sprintf('%02d',i),'_','WACCM_MLS_profiles_at_',...
                num2str(lats(j+1,1)),'to',num2str(lats(end,2)),'_',num2str(allyears(indyear_ind)),'.pdf'];            
        end
        export_fig(filename,'-pdf');
        close all;
    end        
end
