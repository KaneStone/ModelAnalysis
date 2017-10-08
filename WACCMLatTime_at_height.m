% WACCM maps

vertical_deviation = 1;
plot_deviation = 0;
includeVCMAM = 1;
lattime_at_height = 0;
comparetoWACCMSAD = 0;
lineplots = 0;

variable = 'O3';
pres = 150;
lats = [-90 -30];
runs = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
daytitle = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};
MLS_compare = 1;
%% Reading in WACCM data
[NumberDensity, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes, datedata] = ReadWACCMvertical(variable,'monthly');


%% MLS data
if MLS_compare == 1
    
    % a priori
    [MLSpressure,O3ap_monthmean,tempap_monthmean] = importMLSapriori;
    [~,MLSpresindex] = min(abs(MLSpressure - pres));
    % AK
    filename = '/Users/kanestone/work/projects/MLS/AK/MLS-Aura_L2AK-O3-LAT70N_v04-2x_0000d000.txt';
    filename_temp = '/Users/kanestone/work/projects/MLS/AK/MLS-Aura_L2AK-Temperature-LAT70N_v04-2x_0000d000.txt';
    startRow = 23;

    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%13s%13s%13s%13s%13s%s%[^\n\r]';
    
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

    fileID_temp = fopen(filename_temp,'r');
    dataArray_temp = textscan(fileID_temp, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    
    for i = 1:length(dataArray)-1
        AKdata(:,i) = str2double(dataArray{:,i});
        AKdata_temp(:,i) = str2double(dataArray_temp{:,i});
    end
    
    AKdata1 = reshape(AKdata',60,55);
    AKdatahold = AKdata1;
    AKdatahold (isnan(AKdata1)) = [];
    AKdata2 = reshape(AKdatahold,55,55);
    
    AKdata1_temp = reshape(AKdata_temp',60,55);
    AKdatahold_temp = AKdata1_temp;
    AKdatahold_temp (isnan(AKdata1_temp)) = [];
    AKdata2_temp = reshape(AKdatahold_temp,55,55);
    
    
end

%%
[~, SAD_SULFC_MolConc, SAD_SULFC_Temperature, SAD_SULFC_Pressure, SAD_SULFC_Altitude, SAD_SULFC_GMH, SAD_SULFC_Latitudes, SAD_SULFC_datedata] = ReadWACCMvertical('SAD_SULFC','monthly');

    for i = 1:length(runs)

        if MLS_compare
            %interp onto MLS grid
            for j = 1:size(MolConc.(runs{i}),1)
                for k = 1:size(MolConc.(runs{i}),3)
                    MolConcMLSgrid.(runs{i})(j,:,k) = interp1(log(Pressure.(runs{i})(j,:,k)./100),...
                        MolConc.(runs{i})(j,:,k),log(MLSpressure),'linear','extrap');
                    TemperatureMLSgrid.(runs{i})(j,:,k) = interp1(log(Pressure.(runs{i})(j,:,k)./100),...
                        Temperature.(runs{i}).T(j,:,k),log(MLSpressure),'linear','extrap');
                    NDMLSgrid.(runs{i})(j,:,k) = vmr2conc(MolConcMLSgrid.(runs{i})(j,:,k),...
                        TemperatureMLSgrid.(runs{i})(j,:,k),MLSpressure'*100,'O3','conc');                                        
                end
            end
            
            
            latforweight = [-90 -85; -85 -80; -80 -75; -75 -70; -70 -65; -65 -60; -60 -55; -55 -50;...
                -50 -45; -45 -40; -40 -35; -35 -30; -30 -25];
            for j = 1:size(latforweight,1)
                latsforweightindex = find(Latitudes >= latforweight(j,1) & Latitudes < latforweight(j,2));
                for k = 1:size(MolConcMLSgrid.(runs{i}),2)
                    MolClatweighted.(runs{i})(:,k,j) = weightedaverage(squeeze(MolConcMLSgrid.(runs{i})(latsforweightindex,k,:)),...
                        Latitudes(latsforweightindex));
                    Templatweighted.(runs{i})(:,k,j) = weightedaverage(squeeze(TemperatureMLSgrid.(runs{i})(latsforweightindex,k,:)),...
                        Latitudes(latsforweightindex));
                    NDlatweighted.(runs{i})(:,k,j) = weightedaverage(squeeze(NDMLSgrid.(runs{i})(latsforweightindex,k,:)),...
                        Latitudes(latsforweightindex));
                end                                                
            end                            
            
            %extract same time period as MLS
            %MLS years = 2004:2016
            MolClatweighted2.(runs{i}) = MolClatweighted.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
                datedata.(runs{i}).date <= 20170101,:,:);
            Templatweighted2.(runs{i}) = Templatweighted.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
                datedata.(runs{i}).date <= 20170101,:,:);
            NDlatweighted2.(runs{i}) = NDlatweighted.(runs{i})(datedata.(runs{i}).date >= 20040201 & ...
                datedata.(runs{i}).date <= 20170101,:,:);
            monthend = num2str(datedata.(runs{i}).date(end));
            yearend = num2str(datedata.(runs{i}).date(end));
            yearend = str2double(yearend(1:4));
            monthend = str2double(monthend(5:6))-1;
            yearappend = 2016-yearend;
            monthappend = (yearappend*12)+12-monthend;
            MolClatweighted2.(runs{i}) = cat(1,MolClatweighted2.(runs{i}),zeros(monthappend,...
                size(MolClatweighted2.(runs{i}),2),size(MolClatweighted2.(runs{i}),3)));
            Templatweighted2.(runs{i}) = cat(1,Templatweighted2.(runs{i}),zeros(monthappend,...
                size(MolClatweighted2.(runs{i}),2),size(MolClatweighted2.(runs{i}),3)));
            NDlatweighted2.(runs{i}) = cat(1,NDlatweighted2.(runs{i}),zeros(monthappend,...
                size(MolClatweighted2.(runs{i}),2),size(MolClatweighted2.(runs{i}),3)));
            MolClatweighted2.(runs{i}) (MolClatweighted2.(runs{i}) == 0) = NaN;
            Templatweighted2.(runs{i}) (Templatweighted2.(runs{i}) == 0) = NaN;
            NDlatweighted2.(runs{i}) (NDlatweighted2.(runs{i}) == 0) = NaN;
            for j = 1:12
                %for k = 1:17
                    MolClatweightedMLSperiod.(runs{i})(:,j,:,:) = MolClatweighted2.(runs{i})(j:12:end,:,:);
                    TemplatweightedMLSperiod.(runs{i})(:,j,:,:) = Templatweighted2.(runs{i})(j:12:end,:,:);
                    NDlatweightedMLSperiod.(runs{i})(:,j,:,:) = NDlatweighted2.(runs{i})(j:12:end,:,:);
                %end
            end
            % convolve with MLS averaging kernels
            
            for j = 1:size(MolClatweightedMLSperiod.(runs{i}),1)
                for k = 1:size(MolClatweightedMLSperiod.(runs{i}),2)
                    for l = 1:size(MolClatweightedMLSperiod.(runs{i}),4)
                        WACCM_MolC_convolve.(runs{i})(j,k,:,l) = squeeze(O3ap_monthmean(j,k,:,l))+...
                            AKdata2*(squeeze(MolClatweightedMLSperiod.(runs{i})(j,k,:,l))-squeeze(O3ap_monthmean(j,k,:,l)));
                        WACCM_Temp_convolve.(runs{i})(j,k,:,l) = squeeze(tempap_monthmean(j,k,:,l))+...
                            AKdata2_temp*(squeeze(TemplatweightedMLSperiod.(runs{i})(j,k,:,l))-squeeze(tempap_monthmean(j,k,:,l)));                        
                        WACCM_ND_convolve.(runs{i})(j,k,:,l) = vmr2conc(squeeze(WACCM_MolC_convolve.(runs{i})(j,k,:,l)),...
                            squeeze(WACCM_Temp_convolve.(runs{i})(j,k,:,l)),MLSpressure*100,'O3','conc');
                    end
                end
            end
            
            %comparing different grid interpolations
%             figure;
%             plot(squeeze(NDlatweighted.(runs{i})(26*12+10,:,2)),log(MLSpressure))            
%             hold on
%             plot(squeeze(NDlatweighted2.(runs{i})(22,:,2)),log(MLSpressure))
%             plot(squeeze(NDlatweightedMLSperiod.(runs{i})(2,10,:,2)),log(MLSpressure))
%             plot(squeeze(WACCM_ND_convolve.(runs{i})(2,10,:,2)),log(MLSpressure))
%             plot(squeeze(NDMLSgrid.(runs{i})(5,:,26*12+10)),log(MLSpressure))
%             plot(squeeze(NumberDensity.(runs{i})(5,:,26*12+10)),log(squeeze(Pressure.(runs{i})(6,:,27*12+10))./100))
%             legend('Lat weighted','Lat weighted2','Lat weighted MLS period','Convolved','MLS grid','WACCM grid')
%             figure;
%             plot(squeeze(MolClatweighted.(runs{i})(26*12+10,:,2)),log(MLSpressure))
%             hold on
%             plot(squeeze(MolClatweightedMLSperiod.(runs{i})(2,10,:,2)),log(MLSpressure))
%             plot(squeeze(WACCM_MolC_convolve.(runs{i})(2,10,:,2)),log(MLSpressure))
%             plot(squeeze(MolConcMLSgrid.(runs{i})(5,:,26*12+10)),log(MLSpressure))
%             plot(squeeze(MolConc.(runs{i})(5,:,26*12+10)),log(squeeze(Pressure.(runs{i})(5,:,27*12+10))./100))            
%             legend('Lat weighted','Lat weighted MLS period','Convolved','MLS grid','WACCM grid')
            
            %save('/Users/kanestone/work/projects/WACCM/MatlabOutput/WACCM_MLS_compare.mat','WACCM_MolC_convolve');
            %extracting height
            WACCM_MolC_convolve_atpres(i,:,:,:) = WACCM_MolC_convolve.(runs{i})(:,:,MLSpresindex,:);
            
        end
        [~,presind] = min(abs(Pressure.(runs{i})-pres*100),[],2);
        presind2 = floor(nanmean(presind(:)))+1;

        NDatheight.(runs{i}) = squeeze(NumberDensity.(runs{i})(:,presind2,:));
        Molcatheight.(runs{i}) = squeeze(MolConc.(runs{i})(:,presind2,:));
        for j = 1:size(NumberDensity.(runs{i}),1)
            for k = 1:size(NumberDensity.(runs{i}),3)
                NDatheightinterp.(runs{i})(j,k) = interp1(log(Pressure.(runs{i})(j,:,k)),...
                    NumberDensity.(runs{i})(j,:,k),log(150*100));
                Molcatheightinterp.(runs{i})(j,k) = interp1(log(Pressure.(runs{i})(j,:,k)),...
                    MolConc.(runs{i})(j,:,k),log(150*100));
            end
        end
        NDatheightatlat.(runs{i}) = NDatheightinterp.(runs{i})(Latitudes >= lats(1) & Latitudes <= lats(2),:);
        Molcatheightatlat.(runs{i}) = Molcatheightinterp.(runs{i})(Latitudes >= lats(1) & Latitudes <= lats(2),:);
        dateendstring = num2str(datedata.(runs{i}).date(end)); 
        dateendmonth = str2double(dateendstring(5:6))-1;
        NDatheightatlat.(runs{i}) = [NDatheightatlat.(runs{i}),zeros(size(NDatheightatlat.(runs{i}),1),12-dateendmonth)];
        Molcatheightatlat.(runs{i}) = [Molcatheightatlat.(runs{i}),zeros(size(Molcatheightatlat.(runs{i}),1),12-dateendmonth)];
        NDatheightatlat.(runs{i}) (NDatheightatlat.(runs{i}) == 0) = NaN;
        Molcatheightatlat.(runs{i}) (Molcatheightatlat.(runs{i}) == 0) = NaN;
    end

%% plotting

    [~,latind(1)] = min(abs(Latitudes - lats(1)));
    [~,latind(2)] = min(abs(Latitudes - lats(2)));
    latind(2) = latind(2) - 1;
    Volcanoes = NDatheightatlat.MAM - NDatheightatlat.VCMAM;
    Volcanoespercent = (NDatheightatlat.MAM - NDatheightatlat.VCMAM)./NDatheightatlat.VCMAM*100;

%     for_subplot = cat(3,NDatheightatlat.MAM(:,9:12:end),NDatheightatlat.VCMAM(:,9:12:end),...
%         NDatheightatlat.MAM(:,10:12:end),NDatheightatlat.VCMAM(:,10:12:end),...
%         NDatheightatlat.MAM(:,11:12:end),NDatheightatlat.VCMAM(:,11:12:end));
    
    for_subplot = cat(3,NDatheightatlat.MAM(:,6:12:end),NDatheightatlat.MAM(:,7:12:end),...
        NDatheightatlat.MAM(:,8:12:end),NDatheightatlat.MAM(:,9:12:end),...
        NDatheightatlat.MAM(:,10:12:end),NDatheightatlat.MAM(:,11:12:end));
    

    for_subplot_diff = cat(3,Volcanoes(:,9:12:end),Volcanoes(:,10:12:end),Volcanoes(:,11:12:end));
    for_subplot_diff_percent = cat(3,Volcanoespercent(:,9:12:end),Volcanoespercent(:,10:12:end),...
        Volcanoespercent(:,11:12:end));

    for_subplot = permute(for_subplot,[3,2,1]);
    for_diff = permute(for_subplot_diff,[3,2,1]);
    for_diff_percent = permute(for_subplot_diff_percent,[3,2,1]);

    %% Plotting
if lattime_at_height
    master_title = 'WACCM MAM ozone concentration at 140 hPa';
%     daytitle = {'September, MAM','September, VC-MAM','October, MAM',...
%         'October, VC-MAM','November, MAM','November, VC-MAM'};
    
    daytitle = {'June','July','August',...
        'September','October','November'};
    
    x_intervals = 1999:2015;

    subplotmaps(for_subplot,x_intervals,Latitudes(latind(1):latind(2)),{'div','RdBu'},1,[],14,daytitle,'Year',...
            'Latitude ({\circ}N)','molecules/cm^3','on',[6e11,3.4e12],16,...
            2000:2:2016,2000:2:2016,lats(1):5:lats(2),lats(1):5:lats(2),master_title,1,[2004 2016]);

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/atLevel/','TimeLat','_',num2str(lats(1)),'Sto',num2str(lats(2)),'N.png'];
    export_fig(filename,'-png','-nofontswap');   

    %% Plotting difference

    master_title_diff = ['MAM - VC-MAM ozone concentration at ',num2str(pres),' hPa'];
    daytitle_diff = {'September','October','November'};

    subplotmaps(for_diff,x_intervals,Latitudes(latind(1):latind(2)),{'seq','Blues'},0,[],14,daytitle_diff,'Year',...
            'Latitude ({\circ}N)','molecules/cm^3','on',[-7e11,0],15,...
            2000:2:2016,2000:2:2016,lats(1):5:lats(2),lats(1):5:lats(2),master_title_diff,1,[2004 2016]);    

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/atLevel/','TimeLat','_diff_',num2str(lats(1)),'Sto',num2str(lats(2)),'N.png'];
    export_fig(filename,'-png','-nofontswap');   

    %% plotting difference percent

    master_title_diff_percent = ['MAM - VC-MAM ozone percent at ',num2str(pres),' hPa'];
    daytitle_diff = {'September','October','November'};

    subplotmaps(for_diff_percent,x_intervals,Latitudes(latind(1):latind(2)),{'seq','Blues'},0,[],14,daytitle_diff,'Year',...
            'Latitude ({\circ}N)','Percent','on',[-56,0],15,...
            2000:2:2016,2000:2:2016,lats(1):5:lats(2),lats(1):5:lats(2),master_title_diff_percent,1,[2004 2016]);    

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/atLevel/','TimeLat','_diffpercent_',num2str(lats(1)),'Sto',num2str(lats(2)),'N.png'];
    export_fig(filename,'-png','-nofontswap'); 
end

%% Line plots
if lineplots
    load('/Users/kanestone/work/projects/MLS/code/MLS_O3_monthmeanND.mat');
    O3_monthmeanND = squeeze(O3_monthmeanND);
    O3_monthmean = squeeze(O3_monthmean);
    MLSlats = -87.5:5:-32.5;
   
    % take 5 degree latitude averages (weighted)
    latforweight = [-90 -85; -85 -80; -80 -75; -75 -70; -70 -65; -65 -60; -60 -55; -55 -50;...
       -50 -45; -45 -40; -40 -35; -35 -30];
    for i = 1:size(latforweight,1)
        latsforweightindex = find(Latitudes >= latforweight(i,1) & Latitudes < latforweight(i,2));
        layerMAMweighted(i,:) = weightedaverage(NDatheightatlat.MAM(latsforweightindex,:),Latitudes(latsforweightindex));
        layerVCMAMweighted(i,:) = weightedaverage(NDatheightatlat.VCMAM(latsforweightindex,:),Latitudes(latsforweightindex));
        MolclayerMAMweighted(i,:) = weightedaverage(Molcatheightatlat.MAM(latsforweightindex,:),Latitudes(latsforweightindex));
        MolclayerVCMAMweighted(i,:) = weightedaverage(Molcatheightatlat.VCMAM(latsforweightindex,:),Latitudes(latsforweightindex));
    end
    
   layerMAMweighted (layerMAMweighted == 0) = NaN;
   layerVCMAMweighted (layerVCMAMweighted == 0) = NaN;
   MolclayerMAMweighted (MolclayerMAMweighted == 0) = NaN;
   MolclayerVCMAMweighted (MolclayerVCMAMweighted == 0) = NaN;
   
    for i = 1:12
        layerMAM(:,i,:) = layerMAMweighted(:,i:12:end);
        layerVCMAM(:,i,:) = layerVCMAMweighted(:,i:12:end); 
        MolclayerMAM(:,i,:) = MolclayerMAMweighted(:,i:12:end);
        MolclayerVCMAM(:,i,:) = MolclayerVCMAMweighted(:,i:12:end); 
    end
    
    layerMAM = permute(layerMAM,[3,2,1]);
    layerVCMAM = permute(layerVCMAM,[3,2,1]);
    MolclayerMAM = permute(MolclayerMAM,[3,2,1]);
    MolclayerVCMAM = permute(MolclayerVCMAM,[3,2,1]);
      
    %% plotting
    
    plot_molc = 1;
    latstitle = {'-90{\circ}S to -85{\circ}S','-85{\circ}S to -80{\circ}S',...
        '-80{\circ}S to -75{\circ}S','-75{\circ}S to -70{\circ}S','-70{\circ}S to -65{\circ}S',...
        '-65{\circ}S to -60{\circ}S','-60{\circ}S to -55{\circ}S','-55{\circ}S to -50{\circ}S'};
    latstitle2 = {'90to85S','85to80S','80to75S','75to70S','70to65S','65to60S','60to55S','55to50S'};
    latstitle3 = {'87.5','82.5S','77.5S','72.5S','67.5S','62.5S','57.5S','52.5S'};
    
    monthsplot = 12;
    lattoplot = 4;    
    figure;
    
    fsize = 16;
    fig2 = gcf;   
    set(fig2,'color','white','position',[100 100 900 700],'Visible','on');
    if plot_molc
        ph = plot(2004:2016,squeeze(O3_monthmean(:,monthsplot,lattoplot))*1e6,'LineWidth',2);
        hold on
        phMAM = plot(1999:2016,squeeze(MolclayerMAM(:,monthsplot,lattoplot))*1e6,'LineWidth',2);
        phMAM = plot(1999:2016,squeeze(MolclayerVCMAM(:,monthsplot,lattoplot))*1e6,'LineWidth',2);
    else
        ph = plot(2004:2016,squeeze(O3_monthmeanND(:,monthsplot,lattoplot)),'LineWidth',2);
        hold on
        phMAM = plot(1999:2016,squeeze(layerMAM(:,monthsplot,lattoplot)),'LineWidth',2);
        phMAM = plot(1999:2016,squeeze(layerVCMAM(:,monthsplot,lattoplot)),'LineWidth',2);
    end
    title('November','fontsize',20)
    set(gca,'fontsize',18)
    xlabel('Year','fontsize',18)
    if plot_molc
        ylabel('ppmv','fontsize',18)        
    else
        ylabel('Number Density','fontsize',18)        
    end
    title([daytitle{monthsplot},'{ }',num2str(pres),'{ }','hPa','{ }',latstitle3{lattoplot}],'fontsize',20)
    legend('MLS','MAM','VCMAM')
    lh = legend('MLS','MAM','VCMAM');
    set(lh,'box','off','fontsize',18)
    
    if plot_molc
         filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/obscomp/MLScompareCalbuco/',...
            daytitle{monthsplot},'_WACCM_MLS_',num2str(pres),'hPa',latstitle2{lattoplot},'ppmv.pdf'];
    else
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/obscomp/MLScompareCalbuco/',...
            daytitle{monthsplot},'_WACCM_MLS_',num2str(pres),'hPa',latstitle2{lattoplot},'ND.pdf'];
    end
    export_fig(filename);
%     for i = 1:length(latstitle);
%         sp = subplot(length(latstitle)/2,2',i);
%         ph = plot(2004:2016,squeeze(O3_monthmeanND(:,monthsplot,i)),'LineWidth',2);
%         ax.ColorOrderIndex = 1;
%         hold on
%         phMAM = plot(1999:2016,squeeze(layerMAM(:,monthsplot,i)),'--','LineWidth',2);
%         sp_pos = get(sp,'position');
%         if mod(i,2)
%             ylabel('molec/cm^3','fontsize',fsize+2)
%             set(sp,'position',[sp_pos(1)+.025 sp_pos(2:4)]);
%         else
%             set(sp,'position',[sp_pos(1)-.025 sp_pos(2:4)]);
%         end
%         if i == size(O3_monthmeanND,3) || i == size(O3_monthmeanND,3)-1
%             xlabel('year','fontsize',fsize+2);
%         end
%         set(gca,'fontsize',fsize-2);
%         xlim([2003.5 2016.5]);
%         ylim([.9e12 2.6e12]);
%         title(latstitle{i},'fontsize',fsize+4);
%         if i == size(O3_monthmeanND,3)
%             lh = legend(daytitle,'fontsize',fsize-2);
%             set(lh,'Orientation','horizontal','position',[.5 .02 .01 .01],'box','off');            
%         end
%     end
%     mtit(['MLS Ozone at ',sprintf('%.2f',Pressure(presind)),' hPa'],'xoff',-.01,'yoff',.04,'fontsize',fsize+8);
% 
%     filename2 = ['/Users/kanestone/Dropbox/Work_Share/MITwork/MLS/TimeLatitude/',...
%         'MLS_TimeConc_at_',num2str(Pressure(presind)),'.pdf'];
%     export_fig(filename2,'-pdf');
    
    
end



%% deviation at height
if plot_deviation
    monthstoplot = [7,8,9,10,11,12,1,2,3,4,5,6];
    
    for i = 1:length(monthstoplot)
        daytitletemp = daytitle{monthstoplot(i)};
        daytitledev{i} = daytitletemp(1:3);
    end
    fsize = 18;
    if MLS_compare
            layerdeviationMAM(i,:,:) = WACCM_MolC_convolve_atpres(2,:,:,:);
            layerdeviationVCMAM(i,:,:) = WACCM_MolC_convolve_atpres(3,:,:,:); 
    else
        for i = 1:12
            layerdeviationMAM(i,:,:) = NDatheightatlat.MAM(:,i:12:end);
            layerdeviationVCMAM(i,:,:) = NDatheightatlat.VCMAM(:,i:12:end); 
        end
    end
                    
    without20062011 = 0;
    startrefdate = 2000;
    endrefdate = 2014;
    if without20062011
        Mean20002014 = squeeze(nanmean(cat(2,for_subplot(:,1:2,:),for_subplot(:,4:7,:),for_subplot(:,9:11,:)),2));
        Mean20002015 = squeeze(nanmean(cat(2,for_subplot(:,2:3,:),for_subplot(:,5:8,:),for_subplot(:,10:12,:)),2));
        deviation2015 = squeeze(for_subplot(:,12,:))-Mean20002014;
        deviation2016 = squeeze(for_subplot(:,13,:))-Mean20002015;
        deviation_final = [deviation2015(1:4,:); deviation2016(5:6,:)];
    else
        for i = 1:length(runs)
            Mean20002014MAM = squeeze(nanmean(layerdeviationMAM(:,:,6:16),3));
            Mean20002015MAM = squeeze(nanmean(layerdeviationMAM(:,:,6:17),3));                        
            deviation2015MAM = squeeze(layerdeviationMAM(:,:,17))-Mean20002014MAM;
            deviation2016MAM = squeeze(layerdeviationMAM(:,:,18))-Mean20002015MAM;
            deviation_finalMAM = [deviation2015MAM(7:12,:); deviation2016MAM(1:6,:)];
            
            Mean20002014VCMAM = squeeze(nanmean(layerdeviationVCMAM(:,:,6:16),3));
            Mean20002015VCMAM = squeeze(nanmean(layerdeviationVCMAM(:,:,6:17),3));                        
            deviation2015VCMAM = squeeze(layerdeviationVCMAM(:,:,17))-Mean20002014VCMAM;
            deviation2016VCMAM = squeeze(layerdeviationVCMAM(:,:,18))-Mean20002015VCMAM;
            deviation_finalVCMAM = [deviation2015VCMAM(7:12,:); deviation2016VCMAM(1:6,:)];
        end
    end


% plotting
createfig('medium');

if pres == 150
    minmax = [-7.5e11 4.5e11];
    contourinterval = 1e11;
    cbrew = cbrewer('div','RdBu',15);
    cbrew1 = flipud(cbrew(4:15,:));
elseif pres == 100
    minmax = [-15e11 5e11];
    contourinterval = 2e11;
    cbrew = cbrewer('div','RdBu',15);
    cbrew1 = flipud(cbrew(6:15,:));   
elseif pres == 70    
    minmax = [-2.1e12 7e11];
    contourinterval = 2e11;
    cbrew = cbrewer('div','RdBu',21);
    cbrew1 = flipud(cbrew(8:21,:));
end
[~,h] = contourf(1:size(deviation_finalMAM,1),Latitudes(1:size(deviation_finalMAM,2)),deviation_finalMAM',minmax(1):contourinterval:minmax(2));   

ch = colorbar;
cbaxloc = get(ch,'Position');
%subpos = get(sp,'OuterPosition');
set(gca,'color',[.8 .8 .8]);
set(get(ch,'ylabel'),'string','molecules/cm^3','fontsize',fsize+2)
colormap(cbrew1);  
caxis([minmax(1) minmax(2)]);
set(ch,'Ticks',minmax(1):contourinterval:minmax(2)) 

set(gca,'xtick',1:1:length(monthstoplot),'xticklabel',daytitledev,'fontsize',fsize);
set(gca,'ytick',-85:5:30,'yticklabel',-85:5:30,'fontsize',fsize);
ylabel(['Latitude (',char(176),'N)'],'fontsize',fsize+2);
xlabel('Month','fontsize',fsize+2);
xlim([1 length(monthstoplot)])
if without20062011
    title(['WACCM O3 2015-2016 deviation from 2004-(2014-2015) mean (excluding 2006, 2011) at ',sprintf('%.2f',pres),' hPa'],'fontsize',fsize);
else
    title(['WACCM O3 2015-2016 deviation from 2004-(2014-2015) mean at ',sprintf('%.2f',pres),' hPa'],'fontsize',fsize+2);
end

% if without20062011
%     filename3 = ['/Users/kanestone/Dropbox/Work_Share/MITwork/MLS/TimeLatitude/',...
%         'MLS_Conc_2015diff_at_',sprintf('%.2f',Pressure(presind)),'without20062011.png'];
% else
%     filename3 = ['/Users/kanestone/Dropbox/Work_Share/MITwork/MLS/TimeLatitude/',...
%         'MLS_O3_Conc_2015diff_at_',sprintf('%.2f',Pressure(presind)),'.png'];    
% end
% export_fig(filename3,'-png');
end

%% vertical deviation

unipres = [1000:-50:150 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.05:.1];

for i = 1:length(runs) 
    for l = 1:size(NumberDensity.(runs{i}),3)
        for k = 1:size(NumberDensity.(runs{i}),1)
            vertical_uniform_ozone.(runs{i})(k,:,l) = interp1(log(squeeze(Pressure.(runs{i})(k,:,l))/100),...
                squeeze(NumberDensity.(runs{i})(k,:,l)),log(unipres));
        end        
        %rawdata(i,:,:) = squeeze(var20002014.(names{i})(rawlatindex,:,rawmonth:12:end)); 
    end            
end

if vertical_deviation
    years = [2004, 2014];
    meanyears_mam = vertical_uniform_ozone.MAM(:,:,datedata.MAM.date > 20040301 & datedata.MAM.date <= 20150301);
    meanyears_vcmam = vertical_uniform_ozone.VCMAM(:,:,datedata.VCMAM.date > 20040301 & datedata.VCMAM.date <= 20150301);
    
    %meanyears_mam_pres = Pressure.MAM(:,:,datedata.MAM.date > 20040301 & datedata.MAM.date <= 20150301);
    %meanyears_vcmam_pres = Pressure.MAM(:,:,datedata.VCMAM.date > 20040301 & datedata.VCMAM.date <= 20150301);
    
    for i = 1:12
        clim_mam(:,:,i) = nanmean(meanyears_mam(:,:,i:12:end),3);
        clim_vcmam(:,:,i) = nanmean(meanyears_vcmam(:,:,i:12:end),3);
        
        %clim_mam_pres(:,:,i) = nanmean(meanyears_mam_pres(:,:,i:12:end),3);
        %clim_vcmam_pres(:,:,i) = nanmean(meanyears_vcmam_pres(:,:,i:12:end),3);
        
    end
    clim_mam = circshift(clim_mam,[0,0,2]);
    clim_vcmam = circshift(clim_vcmam,[0,0,2]);  
    
    latindex = find(Latitudes >= -90 & Latitudes <= -30);
    
    MAMdeviation =  vertical_uniform_ozone.MAM(latindex,:,datedata.MAM.date > 20150101 & datedata.MAM.date <= 20151201) - clim_mam(latindex,:,1:11);
    VCMAMdeviation =  vertical_uniform_ozone.VCMAM(latindex,:,datedata.VCMAM.date > 20150101 & datedata.VCMAM.date <= 20151201) - clim_vcmam(latindex,:,1:11);
    
    %MAMdeviation_pres =  Pressure.MAM(latindex,:,datedata.MAM.date > 20150101 & datedata.MAM.date <= 20151201) - clim_mam_pres(latindex,:,1:11);
    %VCMAMdeviation_pres =  Pressure.VCMAM(latindex,:,datedata.VCMAM.date > 20150101 & datedata.VCMAM.date <= 20151201) - clim_mam_pres(latindex,:,1:11);
               
end

if comparetoWACCMSAD
    load(['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
        'WACCM_SAD_SULFC_MAM_anomaly.mat']);            
end
    

%% plotting
createfig('large');
cbrew_vert = cbrewer('div','RdBu',25);
cbrew_vert1 = flipud(cbrew_vert(10:25,:));

fsize = 18;
prestick = [1000 500 200 100 50 20 10 5 2 1];
monthplot = [6,7,8,9,10,11];
%daytitle = {'MAM-September','VCMAM-September','MAM-October','VCMAM-October','MAM-November','VCMAM-November'};
daytitle = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};
latindex = find(Latitudes >= -90 & Latitudes <= -30);
for i = 1:6;

    sp = subplot(3,2,i);
    sp_pos(i,:) = get(sp,'position');
    sp_pos_outer(i,:) = get(sp,'OuterPosition');
    if includeVCMAM
        if mod(i,2)
            contourf(Latitudes(latindex),log(unipres),MAMdeviation(:,:,monthplot(i))',-2.3e12:2e11:9e11,...
                'LineStyle','--','LineColor',[.5 .5 .5]);    
        else contourf(Latitudes(latindex),log(unipres),VCMAMdeviation(:,:,monthplot(i))',-2.3e12:2e11:9e11,...
                'LineStyle','--','LineColor',[.5 .5 .5]);    
        end
    else
        contourf(Latitudes(latindex),log(unipres),MAMdeviation(:,:,monthplot(i))',-2.3e12:2e11:9e11,...
            'LineStyle','--','LineColor',[.5 .5 .5]);            
    end
    caxis([-2.3e12 9e11]);
    set(gca,'YDir','reverse');    
    set(gca,'ytick',fliplr(log(prestick)),'yticklabel',fliplr(prestick),'fontsize',fsize);
    ylim([log(10) log(1000)]);

    set(gca,'color',[.8 .8 .8]);
    
    if i == 1 || i == 3 || i == 5
        ylabel(['Pressure (hPa)'],'fontsize',fsize+2);
    end
    if i == 5 || i == 6
        xlabel(['Latitude (',char(176),'N)'],'fontsize',fsize+2);
    end
    title(daytitle{monthplot(i)},'fontsize',fsize+4);

    if ~mod(i,2)
        set(sp,'position',[sp_pos(i,1)-.07 sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)]);
    else 
        set(sp,'position',[sp_pos(i,1)-.035 sp_pos(i,2) sp_pos(i,3) sp_pos(i,4)]);
    end 

    sp_pos2(i,:) = get(sp,'Position');
    if i == 3 || i == 4        
        set(sp,'position',[sp_pos2(i,1) sp_pos2(i,2)+.025 sp_pos2(i,3) sp_pos2(i,4)]);
    elseif i == 5 || i == 6
        set(sp,'position',[sp_pos2(i,1) sp_pos2(i,2)+.05 sp_pos2(i,3) sp_pos2(i,4)]);
    end
        
    lastoutpos = get(sp,'outerposition');
  
    sp_pos(i,:) = get(sp,'position'); 
    
    if i == 6
        ch = colorbar;
        colormap(cbrew_vert1);           
        caxis([-2.3e12 9e11]);
        set(ch,'Ticks',-2.3e12:4e11:9e11)
        set(get(ch,'ylabel'),'string','molecules/cm^3','fontsize',fsize+2)
        cbaxloc = get(ch,'Position');
        subpos = get(sp,'OuterPosition');                    
            
        set(ch,'Position',[subpos(1)+subpos(3)+(subpos(1)+subpos(3))/25, cbaxloc(2),...
            cbaxloc(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(6,2)],...
            'Box','on','YAxisLocation','right','fontsize',fsize);
            
        set(ch,'YaxisLocation','right','fontsize',fsize);                  
    end
    
    if comparetoWACCMSAD  
        if i == 6
            tic;
            toc;
        end
        %cmapWACCM = repmat([166/255,86/255,40/255],8,1);
        cmapWACCM = cbrewer('seq','YlOrBr',10);
        cmapWACCM = cmapWACCM(2:8,:);
        
        haxes1 = gca;
        haxes1_pos = get(haxes1,'Position');
        haxes2 = axes('Position',haxes1_pos,...
                  'XAxisLocation','bottom',...
                  'YAxisLocation','left',...
                  'Color','none');
        hold on        
        contour(SULFCdeviations.Latitudes,log(SULFCdeviations.Pressure),...
            SULFCdeviations.MAMdeviation(:,:,monthplot(i))',-2.2e-8:2e-8:11.8e-8,'LineWidth',1.5);%,'ShowText','on');
        set(gca,'YDir','reverse');    
        set(gca,'ytick',fliplr(log(prestick)),'yticklabel',fliplr(prestick),'fontsize',fsize);
        ylim([log(10) log(1000)]);
        xlim([-90 -31.2632])
        colormap(haxes2,cmapWACCM)
        %WACCMcount = WACCMcount+1;
        caxis([-2.2e-8 11.8e-8]);
        if i == 3;
            ch1 = colorbar('south');
            set(get(ch1,'ylabel'),'string','WACCM MAM SAD SULFC (cm^2/cm^3)','fontsize',fsize)
            set(ch1,'position',[.095 .07 .740 .02]);
            set(ch1,'Ticks',-2.2e-8:2e-8:11.8e-8)             
        end
    end
    
end

mtit('WACCM 2015 anomaly (from 2004-2014 mean)','xoff',-.01,'yoff',.04,'fontsize',fsize+8);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
    'WACCM_O3_MAM-VCMAM_anomaly_',daytitle{monthplot(1)},'-',daytitle{monthplot(end)},'.png'];
export_fig(filename,'-png','-nofontswap');   
