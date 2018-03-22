% Read in and construct Chem-only CL/HCL
clear all
%% import SWOOSH
Stimeperiod = [2000 2014];
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201611-latpress-10deg-L31.nc');
SWOOSHyears = repmat(1984:2015,[12,1]);
SWOOSHyears = [SWOOSHyears(:);ones(11,1)*2016];
SWOOSHextract = permute(SWOOSH.combinedo3q(:,:,SWOOSHyears >= Stimeperiod(1) & SWOOSHyears <= Stimeperiod(2)),[2,1,3]);
%% import highCl
intpres = [SWOOSH.level;[.9,.8,.7,.6,.5,.4,.3,.2,.1]'];

%% Read in WACCM O3 data
timeperiod = [2000 2014];
commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[~, MolConc_HCL, ~, Pressure, ~, ~, Latitudes, Longitudes, datedata] = ReadWACCMvertical('HCL','monthly',commondir,0,1);
[~, MolConc_CLO, ~, ~,~,~,~,~,~] = ReadWACCMvertical('CLO','monthly',commondir,0,1);
[~, MolConc_CH4, ~, ~,~,~,~,~,~] = ReadWACCMvertical('CH4','monthly',commondir,0,1);
[~, MolConc_CLOBRO_loss, ~, ~,~,~,~,~,~] = ReadWACCMvertical('OddOx_CLOxBROx_Loss','monthly',commondir,0,0);
[~, MolConc_OddOx_Loss_Tot, ~, ~,~,~,~,~,~] = ReadWACCMvertical('OddOx_Loss_Tot','monthly',commondir,0,0);
MolConc_HCL = rmfield(MolConc_HCL, 'MAM1990');
MolConc_CLO = rmfield(MolConc_CLO, 'MAM1990');
MolConc_CH4 = rmfield(MolConc_CH4, 'MAM1990');
MolConc_CLOBRO_loss = rmfield(MolConc_CLOBRO_loss, 'MAM1990');
MolConc_OddOx_Loss_Tot = rmfield(MolConc_OddOx_Loss_Tot, 'MAM1990');

datedata = rmfield(datedata, 'MAM1990');
Pressure = rmfield(Pressure, 'MAM1990');
runs = fields(MolConc_HCL);
runs = {runs{2},runs{6}};

for i = 1:2
    CCMIy = CCMI_years(datedata.(runs{i}).date,1);
    dateind = CCMIy >= timeperiod(1) & CCMIy <= timeperiod(2);    
    for j = 1:size(MolConc_HCL.Chemonlynoleap,1)
        [SDWaccmData(i).HCL(:,j,:),~] = intRegPres(squeeze(MolConc_HCL.(runs{i})(j,:,dateind)),...
            squeeze(Pressure.(runs{i})(j,:,dateind))./100,intpres);
        [SDWaccmData(i).CLO(:,j,:),~] = intRegPres(squeeze(MolConc_CLO.(runs{i})(j,:,dateind)),...
            squeeze(Pressure.(runs{i})(j,:,dateind))./100,intpres);
        [SDWaccmData(i).CH4(:,j,:),~] = intRegPres(squeeze(MolConc_CH4.(runs{i})(j,:,dateind)),...
            squeeze(Pressure.(runs{i})(j,:,dateind))./100,intpres);        
        [SDWaccmData(i).ClBrOxloss(:,j,:),~] = intRegPres(squeeze(MolConc_CLOBRO_loss.(runs{i})(j,:,dateind)),...
            squeeze(Pressure.(runs{i})(j,:,dateind))./100,intpres);        
        [SDWaccmData(i).TotalOxLoss(:,j,:),~] = intRegPres(squeeze(MolConc_OddOx_Loss_Tot.(runs{i})(j,:,dateind)),...
            squeeze(Pressure.(runs{i})(j,:,dateind))./100,intpres);
        
    end
    SDWaccmData(i).HCLnointerp = squeeze(MolConc_HCL.(runs{i})(:,:,dateind));
    SDWaccmData(i).CLOnointerp = squeeze(MolConc_CLO.(runs{i})(:,:,dateind));
    SDWaccmData(i).CH4nointerp = squeeze(MolConc_CH4.(runs{i})(:,:,dateind));
    SDWaccmData(i).ClBrOxlossnointerp = squeeze(MolConc_CLOBRO_loss.(runs{i})(:,:,dateind));
    SDWaccmData(i).TotalOxLossnointerp = squeeze(MolConc_OddOx_Loss_Tot.(runs{i})(:,:,dateind));
    ratio(i,:,:,:) = SDWaccmData(i).CLO./SDWaccmData(i).HCL;
    ratio2(i,:,:,:) = SDWaccmData(i).CH4./SDWaccmData(i).HCL;
    CLyoxloss(i,:,:,:) = SDWaccmData(i).ClBrOxloss./SDWaccmData(i).TotalOxLoss;
end


%%
clearvars chemonlyratio MAMratio CH4_chemonly CH4_MAM ratiomean CH4_mean CH4ratio_chemonly CH4ratio_MAM ratiomean2
lats = [-80 -50; 50 80];

for i = 1:12
    for j = 1:2
        latind = Latitudes >= lats(j,1) & Latitudes <= lats(j,2);
        chemonlyratio(j,:,i,:) = squeeze(nanmean(ratio(2,:,latind,i:12:end),3));
        MAMratio(j,:,i,:) = squeeze(nanmean(ratio(1,:,latind,i:12:end),3));
        
        CH4_chemonly(j,:,i,:) = squeeze(nanmean(SDWaccmData(2).CH4(:,latind,i:12:end),2));
        CH4_MAM(j,:,i,:) = squeeze(nanmean(SDWaccmData(1).CH4(:,latind,i:12:end),2));
        
        CH4ratio_chemonly(j,:,i,:) = squeeze(nanmean(ratio2(2,:,latind,i:12:end),3));
        CH4ratio_MAM(j,:,i,:) = squeeze(nanmean(ratio2(1,:,latind,i:12:end),3));
        
        CLyoxloss_chemonly(j,:,i,:) = squeeze(nanmean(CLyoxloss(2,:,latind,i:12:end),3));
        CLyoxloss_MAM(j,:,i,:) = squeeze(nanmean(CLyoxloss(1,:,latind,i:12:end),3));
        
    end
end

ratiomean(1:2,:,:) = squeeze(nanmean(chemonlyratio,4));
ratiomean(3:4,:,:) = squeeze(nanmean(MAMratio,4));
ratiomean = permute(ratiomean,[1,3,2]);

ratiomean2(1:2,:,:) = squeeze(nanmean(CH4ratio_chemonly,4));
ratiomean2(3:4,:,:) = squeeze(nanmean(CH4ratio_MAM,4));
ratiomean2 = permute(ratiomean2,[1,3,2]);

CH4_mean(1:2,:,:) = squeeze(nanmean(CH4_chemonly,4));
CH4_mean(3:4,:,:) = squeeze(nanmean(CH4_MAM,4));
CH4_mean = permute(CH4_mean,[1,3,2]);

Ox_mean(1:2,:,:) = squeeze(nanmean(CLyoxloss_chemonly,4));
Ox_mean(3:4,:,:) = squeeze(nanmean(CLyoxloss_MAM,4));
Ox_mean = permute(Ox_mean,[1,3,2]);

%% Reading in 

SDtimeperiod = [2000 2014];
[SDWaccmDataOther] = ReadInSDWACCM(intpres,SDtimeperiod);

%% take regression of SD WACCM (chem-only and MAM)

[bchemonly,predictorschemonly,O3Anomalychemonly] = ...
    ozoneRegressionTrends(squeeze(ratio(2,:,:,:)),SDWaccmDataOther(2).U,SDWaccmDataOther(2).solar,SDWaccmDataOther(2).SPE,...
    SDWaccmDataOther(2).NINO34,SDWaccmDataOther(2).HFS,SDWaccmDataOther(2).HFN,SDWaccmDataOther(2).NO2,0,SDtimeperiod,1,0,0);%SDWaccmData(2).SPEintpres

%% Take regression of residuals
[btest,ballchemonly] = ozoneRegressionEnsAve(O3Anomalychemonly.percent_residuals_months,SDtimeperiod(2) - SDtimeperiod(1));
[btest2,ballchemonly2] = ozoneRegressionEnsAve(O3Anomalychemonly.percent,SDtimeperiod(2) - SDtimeperiod(1));

%%
clearvars watemp bmonth bmonth2
%lats = [70 90];
lats = [-80 -50; 50 80];

for i = 1:12
    for k = 1:2
        latind = Latitudes >= lats(k,1) & Latitudes <= lats(k,2);
        for j = 1:size(O3Anomalychemonly.percent_residuals_months,2)
            watemp(:,i,j,k) = weightedaverage(permute(squeeze(O3Anomalychemonly.percent_residuals_months(i:12:end,j,latind)),[2,1]),Latitudes(latind));
            watemp2(:,i,j,k) = weightedaverage(permute(squeeze(O3Anomalychemonly.percent(i:12:end,j,latind)),[2,1]),Latitudes(latind));
            bmonth(k,i,j,:) = regress(watemp(:,i,j,k),[ones(size(watemp,1),1),[1:size(watemp,1)]']);
            bmonth2(k,i,j,:) = regress(watemp2(:,i,j,k),[ones(size(watemp2,1),1),[1:size(watemp2,1)]']);
        end
    end    
end


%% plotting Ox loss
Ox_meantoplot(1,:,:) = Ox_mean(1,:,:);
Ox_meantoplot(2,:,:) = circshift(Ox_mean(2,:,:),[0,6,0]);

prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
    5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
logprestick = log(prestick);

mtit = {['2000',char(8211),'2014 average Ox fractional loss due to ClOx/BrOx']};
titles = {['Chem-only, ',num2str(abs(lats(1,1))),char(8211),num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['Chem-only, ',num2str(abs(lats(2,2))),char(8211),num2str(abs(lats(2,1))),'{\circ}','N'],...
    ['MAM, ',num2str(abs(lats(1,1))),char(8211),num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['MAM,  ',num2str(abs(lats(2,2))),char(8211),num2str(abs(lats(2,1))),'{\circ}','N']};

[fig,sh] = subplotmaps(Ox_meantoplot(1:2,:,:),1:12,log(intpres),{'seq','YlOrBr'},0,[],20,titles,'Month','Pressure (hPa)','Total O_x loss fraction','on',...
    [0 .6],13,1:12,{'J','F','M','A','M','J','J','A','S','O','N','D'},...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[1 12],[log(.1) log(30)],1,'--',0,'');

set(gca,'xticklabel',{'J','A','S','O','N','D','J','F','M','A','M','J'});

overlay = 1;

if overlay
    cbrew2 = cbrewer('seq','Greys',11);
    cbrew2 = cbrew2(6:10,:);
    h = get(gcf,'children');
    hold on

    CH4_meantoplot(1,:,:) = CH4_mean(1,:,:);
    CH4_meantoplot(2,:,:) = circshift(CH4_mean(2,:,:),[0,6,0]);

    for i = 1:2
        axes2(i).h = axes;
        axes2(i).h.Visible = 'off';
        [ch2,h2] = contour(1:12,log(intpres),squeeze(CH4_meantoplot(i,:,:))'*1e6,.2:.2:1,'LineWidth',3);%[-.2:.2:2]
        set(gca,'ytick',fliplr(logprestick),'yticklabel',fliplr(presticklabel),'YDir','reverse');
        set(gca,'color','none');
        set(gca,'xtick',1:12,'xticklabel',1:12);
        set(axes2(i).h,'position',get(sh(i),'position'));
        colormap(axes2(i).h,cbrew2)
        clabel(ch2,h2,'FontSize',14,'LabelSpacing',400);        
        axes2(i).h.XTick = [];
        axes2(i).h.YTick = [];
    end
    ch = colorbar;
    caxis([.2 1]);
    set(ch,'orientation','horizontal','position',[.1 .06 .78 .02],'fontsize',20,'Ticks',[.28:(.16):1],'TickLabels',[.2:.2:1]);
    %set(ch,'orientation','horizontal','position',[.1 .06 .78 .02],'fontsize',20,'Ticks',[.2:.2:1]);     
    set(get(ch,'ylabel'),'string','CH_4 (ppmv)','fontsize',20,'units','normalized','position',[.5,-1.3, .01]);
    
    

end


filename = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/CLOxBROx_OxLoss_2000-2014average.pdf';
 
%export_fig(filename,'-pdf');
print(filename,'-depsc'); 
 
 %% plotting Ox loss
Ox_meantoplot(1,:,:) = Ox_mean(3,:,:);
Ox_meantoplot(2,:,:) = circshift(Ox_mean(4,:,:),[0,6,0]);

prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
    5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
logprestick = log(prestick);

mtit = {['2000-2014 average Ox fractional loss due to ClOx/BrOx']};
titles = {['MAM, ',num2str(abs(lats(1,1))),'-',num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['MAM,  ',num2str(abs(lats(2,1))),'-',num2str(abs(lats(2,2))),'{\circ}','N'],...
    ['MAM, ',num2str(abs(lats(1,1))),'-',num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['MAM,  ',num2str(abs(lats(2,1))),'-',num2str(abs(lats(2,2))),'{\circ}','N']};

[fig,sh] = subplotmaps(Ox_meantoplot(1:2,:,:),1:12,log(intpres),{'seq','YlOrBr'},0,[],20,titles,'Month','Pressure (hPa)','Total Ox loss fraction','on',...
    [0 .6],13,1:12,1:12,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[1 12],[log(.1) log(30)],1,'--',0,'');

set(gca,'xticklabel',[7 8 9 10 11 12 1 2 3 4 5 6]);

overlay = 1;

if overlay
    cbrew2 = cbrewer('seq','Greys',10);
    cbrew2 = cbrew2(6:end,:);
    h = get(gcf,'children');
    hold on

    CH4_meantoplot(1,:,:) = CH4_mean(3,:,:);
    CH4_meantoplot(2,:,:) = circshift(CH4_mean(4,:,:),[0,6,0]);

    for i = 1:2
        axes2(i).h = axes;
        axes2(i).h.Visible = 'off';
        [ch2,h2] = contour(1:12,log(intpres),squeeze(CH4_meantoplot(i,:,:))'*1e6,0:.2:1,'LineWidth',3);%[-.2:.2:2]
        set(gca,'ytick',fliplr(logprestick),'yticklabel',fliplr(presticklabel),'YDir','reverse');
        set(gca,'color','none');
        set(gca,'xtick',1:12,'xticklabel',1:12);
        set(axes2(i).h,'position',get(sh(i),'position'));
        colormap(axes2(i).h,cbrew2)
        %clabel(ch2,h2);        
        
        axes2(i).h.XTick = [];
        axes2(i).h.YTick = [];
    end
    ch = colorbar;
    set(ch,'orientation','horizontal','position',[.1 .06 .78 .02],'fontsize',20,'Ticks',[0:.2:1]);
    caxis([0 1]);
    locate = get(ch,'title');
    pos = get(locate,'position');           
    set(get(ch,'ylabel'),'string','CH4 (ppmv)','fontsize',20)

end


filename = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/MAM_CLOxBROx_OxLoss_2000-2014average.pdf';
 
 export_fig(filename,'-pdf');
 
%% plotting

ratiomeantoplot(1,:,:) = ratiomean(1,:,:);
ratiomeantoplot(2,:,:) = circshift(ratiomean(2,:,:),[0,6,0]);

mtit = {['2000-2014 average ClO/HCl ratio']};
titles = {['Chem-only, ',num2str(abs(lats(1,1))),'-',num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['Chem-only,  ',num2str(abs(lats(2,1))),'-',num2str(abs(lats(2,2))),'{\circ}','N'],...
    ['MAM, ',num2str(abs(lats(1,1))),'-',num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['MAM,  ',num2str(abs(lats(2,1))),'-',num2str(abs(lats(2,2))),'{\circ}','N']};

[fig,sh] = subplotmaps(ratiomeantoplot(1:2,:,:),1:12,log(intpres),{'seq','YlOrBr'},0,[],22,titles,'Month','Pressure (hPa)','Cl/HCl','on',...
    [0 .3],16,1:12,1:12,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[1 12],[log(.1) log(30)],1,'--',0,'');

set(gca,'xticklabel',[7 8 9 10 11 12 1 2 3 4 5 6]);

filename = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/ClO_HClratio_2000-2014average.pdf';
 
 export_fig(filename,'-pdf');

%%

[fig,sh] = subplotmaps(bmonth(:,:,:,2),1:12,log(intpres),{'div','RdBu'},1,[],22,titles,'Month','Pressure (hPa)','Cl/HCl','on',...
    [-.5 .5],16,1:12,1:12,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[1 12],[log(.1) log(30)],1,'--',0,'');

[fig,sh] = subplotmaps(bmonth2(:,:,:,2),1:12,log(intpres),{'div','RdBu'},1,[],22,titles,'Month','Pressure (hPa)','Cl/HCl','on',...
    [-.5 .5],16,1:12,1:12,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[1 12],[log(.1) log(30)],1,'--',0,'');

%% overlay CH4
cbrew2 = cbrewer('seq','Blues',27);
cbrew2 = cbrew2(16:end,:);
h = get(gcf,'children');
hold on

CH4_meantoplot(1,:,:) = CH4_mean(1,:,:);
CH4_meantoplot(2,:,:) = circshift(CH4_mean(2,:,:),[0,6,0]);

for i = 1:2
    axes2(i).h = axes;
    axes2(i).h.Visible = 'off';
    ch2 = contour(1:12,log(intpres),squeeze(CH4_meantoplot(i,:,:))'*1e6,[0:.2:2],'LineWidth',2.5);%[-.2:.2:2]
    set(gca,'ytick',fliplr(logprestick),'yticklabel',fliplr(presticklabel),'YDir','reverse');
    set(gca,'color','none');
    set(gca,'xtick',1:12,'xticklabel',1:12);
    set(axes2(i).h,'position',get(sh(i),'position'));
    colormap(axes2(i).h,cbrew2)
    ch = colorbar;
    set(ch,'orientation','horizontal','position',[.1 .05 .78 .02],'fontsize',22);
    axes2(i).h.XTick = [];
    axes2(i).h.YTick = [];
end

% filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/',...
%    'Chem-only_and_MAM_ClHClratio'];
% 
% export_fig(filename,'-png');
% export_fig(filename,'-pdf');

%%

mtit = {['2000-2014 average ClO/HCl ratio']};
titles = {['Chem-only, ',num2str(abs(lats(1,1))),'-',num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['Chem-only,  ',num2str(abs(lats(2,1))),'-',num2str(abs(lats(2,2))),'{\circ}','N'],...
    ['MAM, ',num2str(abs(lats(1,1))),'-',num2str(abs(lats(1,2))),'{\circ}','S'],...
    ['MAM,  ',num2str(abs(lats(2,1))),'-',num2str(abs(lats(2,2))),'{\circ}','N']};

[fig,sh] = subplotmaps(ratiomean2,1:12,log(intpres),{'seq','YlOrBr'},0,[],22,titles,'Month','Pressure (hPa)','Cl/HCl','on',...
    [0 1e3],32,1:12,1:12,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[1 12],[log(.1) log(30)],1,'--',0,'');


