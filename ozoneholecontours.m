% comparing calculated total column ozone to model derived total column ozone
clear all

fontsize = 26;
leveltoplot = 'plotTCO'; %plotTCO, plothPa
hPatoplot = 100;

matstand = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

if strcmp(leveltoplot,'plotTCO')
    contourlevel = 220;
    % Read in model derived total column ozone
    [~,WACCM_TCO,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/TCO/TOZ_f.e11.FWTREFC1SD.f19.f19.ccmi30.mam.004f_merged_c161106.nc');
    [~,WACCM_TCO_VC,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/TCO/TOZ_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.004f_merged_c161106.nc');

    % reading in OMI monthly average data

    dirOMI = dir(['/Users/kanestone/work/projects/OMI/output/','*.mat']);
    for i = 1:length(dirOMI)
        OMIdata = load(['/Users/kanestone/work/projects/OMI/output/',dirOMI(i).name]);
        OMI_monthly(i,:,:,:) = OMIdata.data_monthly;
        if i == length(dirOMI)
            load(['/Users/kanestone/work/projects/OMI/output/',dirOMI(i).name]);        
        end
    end
    lontoplot = [WACCM_TCO.lon;360];
    WACCM_TCO.toz = cat(1,WACCM_TCO.toz,WACCM_TCO.toz(1,:,:));
    WACCM_TCO_VC.toz = cat(1,WACCM_TCO_VC.toz,WACCM_TCO_VC.toz(1,:,:));
    OMI_clim = squeeze(nanmean(OMI_monthly));
elseif strcmp(leveltoplot,'plothPa')    
    % Reading in MLS
    pres = 100;
    addpath('/Users/kanestone/work/projects/MLS/code');
    [MLS,presind,allyears,MLSPressure] = importMLSmonthly(pres,0);
    % interpolating to hPa level
    % reading in WACCM
           
    MLSlon = [-175:10:185]';
    MLSlat = -90:5:-30;
    
    load('/Users/kanestone/work/projects/WACCM/output/latlon/WACCM_MolC_convolve.mat');
    [~,data,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/O3/O3_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.2015.mam.004f.cam.h0zm_new.nc');
    waclat1 = data.lat(1:34);
    
    waclon1 = 0:2.5:360;
    MLSclim = squeeze(nanmean(MLS.O3_monthmean_latlon));         
    MLSclim = cat(3,MLSclim,MLSclim(:,:,1,:));
    MLS.O3_monthmean_latlon = cat(4,MLS.O3_monthmean_latlon,MLS.O3_monthmean_latlon(:,:,:,1,:));
    WACCM_MolC_convolve.MAM = cat(1,WACCM_MolC_convolve.MAM,WACCM_MolC_convolve.MAM(1,:,:,:));
    WACCM_MolC_convolve.VCMAM = cat(1,WACCM_MolC_convolve.VCMAM,WACCM_MolC_convolve.VCMAM(1,:,:,:));
end


%% calculating
testing = 1;
if testing
    % interpolate to high resolution in latitude
    newlat = -90:.25:90;
    %finding month to plot
    datetoplot = 20151201;
    dateind = find(WACCM_TCO.date == datetoplot);        

    %% testing code to see if contouring is doing what it should.

    for j = 1:size(WACCM_TCO.toz,1)
        WACCM_toz_interp(j,:) = interp1(WACCM_TCO.lat,squeeze(WACCM_TCO.toz(j,:,dateind)),newlat);
        WACCM_toz_interp_VC(j,:) = interp1(WACCM_TCO_VC.lat,squeeze(WACCM_TCO_VC.toz(j,:,dateind)),newlat);
    end
    % finding area below 220 DU
    for i = 1:size(WACCM_TCO.toz,3)
        for j = 1:length(lontoplot)
            [xbelow(j,i).a, ybelow(j,i).a] = find(WACCM_TCO.toz(j,:,i) < 220);
            [xbelow_VC(j,i).a, ybelow_VC(j,i).a] = find(WACCM_TCO_VC.toz(j,:,i) < 220);
            if ~isempty(ybelow(j,i).a)
                maxlat(j,i) = max(ybelow(j,i).a);            
            else
                maxlat(j,i) = NaN;

            end
            if ~isempty(ybelow_VC(j,i).a)
                maxlat_VC(j,i) = max(ybelow_VC(j,i).a);
            else
                maxlat_VC(j,i) = NaN;
            end
        end
    end
    clearvars xbelow ybelow xbelow_VC ybelow_VC
        for j = 1:length(lontoplot)
    %     if i == 202
    %         a = 1;
    %     end
            [xbelow_interp(j).a, ybelow_interp(j).a] = find(WACCM_toz_interp(j,:) < 220);
            [xbelow_interp_VC(j).a, ybelow_interp_VC(j).a] = find(WACCM_toz_interp_VC(j,:) < 220);
            if ~isempty(ybelow_interp(j).a)
                maxlat_interp(j) = max(ybelow_interp(j).a);            
            else
                maxlat_interp(j) = NaN;

            end
            if ~isempty(ybelow_interp_VC(j).a)
                maxlat_interp_VC(j) = max(ybelow_interp_VC(j).a);
            else
                maxlat_interp_VC(j) = NaN;
            end
        end
        clearvars xbelow_interp ybelow_interp xbelow_interp_VC ybelow_interp_VC
    %end

    
%%
    figure;
    contour(double(WACCM_TCO.lat),double(lontoplot),double(WACCM_TCO.toz(:,:,dateind)),[220 220],'color','b','LineWidth',2);
    hold on
    contour(double(WACCM_TCO.lat),double(lontoplot),double(WACCM_TCO_VC.toz(:,:,dateind)),[220 220],'color','r','LineWidth',2);
    contour(double(WACCM_TCO.lat),double(lontoplot),double(WACCM_TCO.toz(:,:,dateind)),[219 219],'color','b','LineWidth',2);
    contour(double(WACCM_TCO.lat),double(lontoplot),double(WACCM_TCO_VC.toz(:,:,dateind)),[219 219],'color','r','LineWidth',2);
    plot(double(WACCM_TCO.lat(maxlat(:,dateind))),double(lontoplot),'--b','LineWidth',2)
    plot(double(WACCM_TCO.lat(maxlat_VC(:,dateind))),double(lontoplot),'--r','LineWidth',2)
    plot(newlat(maxlat_interp),double(lontoplot),'color','k','LineWidth',2);
    plot(newlat(maxlat_interp_VC),double(lontoplot),'color',[.7 .7 .7],'LineWidth',2);
    xlim([-80 -60])
end


%% TCO
if strcmp(leveltoplot,'plotTCO');
    cbrew = cbrewer('qual','Set1',8);
    year = 2015;
    month = 10;
    omimonth = month-1;
    figure;
    fig  = gcf;
    set(fig,'color','white','position',[100 100 1800 600],'Visible','on');

    for i = 1:3
        datetoplot = str2double([num2str(year),sprintf('%02d',month),'01']);
        dateind = find(WACCM_TCO.date == datetoplot);
        sp(i) = subplot(1,3,i);
        sp_pos(i,:) = get(sp(i),'position');
        if i == 1

        elseif i == 3

        end
        %axesm('MapProjection','stereo','MapLatLimit',[-90 -50])       
        axesm('MapProjection','eqaazim','MapLatLimit',[-90 -50])
        %h = geoshow('landareas.shp', 'FaceColor',[178,223,138]./255);%[0.6 .9 0.6]);
        h = geoshow('landareas.shp', 'FaceColor',[.80 .80 .80]);%[0.6 .9 0.6]);
        hold on
        [C0,OMI_c] = contourm(latitudes,longitudes,OMI_clim(:,:,omimonth)',[contourlevel contourlevel],'color','k','LineWidth',2);
        [C,OMI] = contourm(latitudes,longitudes,data_monthly(:,:,omimonth)',[contourlevel contourlevel],'color',matstand(3,:),'LineWidth',2);
        [C1,MAM] = contourm(double(WACCM_TCO.lat),double(lontoplot),double(WACCM_TCO.toz(:,:,dateind))',[contourlevel contourlevel],'--','color',cbrew(2,:),'LineWidth',2);
        [C2,VCMAM] = contourm(double(WACCM_TCO.lat),double(lontoplot),double(WACCM_TCO_VC.toz(:,:,dateind))',[contourlevel contourlevel],'-.','color',cbrew(1,:),'LineWidth',2);

        lonwidthMAM = diff(C1(1,:));
        lonwidthVCMAM = diff(C2(1,:));
        
        % integrating area
        earthellipsoid = referenceSphere('earth','km');
        for k = 1:length(lonwidthMAM(2:end))
            
            aquadtech(k) = areaquad(-90,0,C1(2,k+1),abs(lonwidthMAM(k+1)),earthellipsoid);%*4*pi*6371000^2;
            
            lat = C1(2,k+1);
            lat2 = -90-lat;
            latstep = abs(lat2/100);
            latint = -90+latstep:latstep:lat;
            for j = 1:length(latint)
                intstep(j) = cosd(latint(j));               
            end
            cosdint = cosd(-90)-cosd(lat);
            MAMintegral(k) = 2*pi*abs((6371000^2).*sum(intstep)./(360./lonwidthMAM(k+1)));
            MAMintegral2(k) = 2*pi*abs((6371000^2).*cosdint./(360./lonwidthMAM(k+1)));
        end
        MAMcompleteintegral = sum(MAMintegral);
        MAMcompleteintegral2 = sum(MAMintegral2);
        aquadtechcomplete = sum(aquadtech)+aquadtech(k);
        for k = 1:length(lonwidthVCMAM(2:end))
            
            aquadtechVCMAM(k) = areaquad(-90,0,C2(2,k+1),abs(lonwidthVCMAM(k+1)),earthellipsoid);%*4*pi*6371000^2;
            
            lat = C2(2,k+1);
            lat2 = -90-lat;
            latstep = abs(lat2/100);
            latint = -90+latstep:latstep:lat;
            for j = 1:length(latint)
                intstep(j) = cosd(latint(j));               
            end
            cosdint = cosd(-90)-cosd(lat);
            
            VCMAMintegral(k) = 2*pi*abs((6371000^2).*sum(intstep)*.01./(360./lonwidthVCMAM(k+1)));
            VCMAMintegral2(k) = 2*pi*abs((6371000^2).*cosdint./(360./lonwidthVCMAM(k+1)));
        end
        VCMAMcompleteintegral = sum(VCMAMintegral);
        VCMAMcompleteintegral2 = sum(VCMAMintegral2);
        aquadtechVCMAMcomplete = sum(aquadtechVCMAM)+(aquadtechVCMAM(k)); 
        
        
        areadiff(i) = MAMcompleteintegral - VCMAMcompleteintegral;
        areadiff2(i) = MAMcompleteintegral2 - VCMAMcompleteintegral2;
        quadtechdiff(i) = aquadtechcomplete - aquadtechVCMAMcomplete;
        
        waclat = [-90,-70.7,-69,-68.6,-64.2]; % South Pole, Neumayer, Syowa, Davis, Marambio
        waclon = [158,351.7,39.6,78,303.4];


        %sites = [-70.7,351.7;-69.0,39.54;-68.6,77.97;-64.2,303.4];
        sites = [-70.7,351.7;-69.0,39.54;-68.6,77.97];

        h(1) = plotm(sites(1,1),sites(1,2),'Marker','pentagram','linestyle','none',...
            'MarkerEdgeColor',cbrew(3,:),'MarkerFaceColor',cbrew(3,:),'MarkerSize',18.5);
            set(gca,'fontsize',fontsize);

        h(2) = plotm(sites(2,1),sites(2,2),'Marker','s','linestyle','none',...
            'MarkerEdgeColor',cbrew(5,:),'MarkerFaceColor',cbrew(5,:),'MarkerSize',15);
            set(gca,'fontsize',fontsize);

        h(3) = plotm(sites(3,1),sites(3,2),'Marker','v','linestyle','none',...
            'MarkerEdgeColor',cbrew(4,:),'MarkerFaceColor',cbrew(4,:),'MarkerSize',15);
            set(gca,'fontsize',fontsize);

%         h(4) = plotm(sites(4,1),sites(4,2),'Marker','d','linestyle','none',...
%             'MarkerEdgeColor',cbrew(8,:),'MarkerFaceColor',cbrew(8,:),'MarkerSize',13.5);
%             set(gca,'fontsize',fontsize);

        setm(gca,'frame','on')
        set(gca','Visible','off')
        setm(gca,'parallellabel','on')
        %setm(gca,'meridianlabel','on','mlabellocation',360,'fontsize',fontsize-2)
        gridm('on');
        gridm('mlinelocation',90,'plinelocation',10,'plabellocation',[-90 -80 -70 -60],'fontsize',fontsize-10,'GLineWidth',2,'gcolor',[.7 .7 .7]);
        if i == 3 
            lh = legend([OMI_c,OMI,MAM,VCMAM,h(1),h(2),h(3)],'MLS / OMI 2004-2014 mean','MLS / OMI','MAM','VC-MAM','Neumayer','Syowa','Davis');
            set(lh,'orientation','horizontal','position',[.49 .18 .05 .05],'box','off');
        end

        if month ~= 12
            month = month+1;
            omimonth = omimonth+1;
        else month = 1;
            year = year+1;
            omimonth = omimonth+1;
        end        

    end

    set(sp(1),'position',[sp_pos(1,1)+.1, sp_pos(1,2:4)])
    set(sp(3),'position',[sp_pos(3,1)-.1, sp_pos(3,2:4)])
    mtit('September','xoff',-.295,'yoff',-.2,'fontsize',fontsize+2)
    mtit('October','xoff',.005,'yoff',-.2,'fontsize',fontsize+2)
    mtit('November','xoff',+.305,'yoff',-.2,'fontsize',fontsize+2)
    mtit('2015 monthly mean ozone hole contours','xoff',0,'yoff',-.1,'fontsize',fontsize+6);

    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
    'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

    text(.225, .37,'TCO < 200 DU','HorizontalAlignment'...
        ,'left','color','k','VerticalAlignment','top','rotation',90,'fontsize',fontsize);

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/areabelow220/','OctNov',num2str(year),'.png'];
    %export_fig(filename,'-png')
else
    
    cbrew = cbrewer('qual','Set1',8);
    year = 2015;
    month = 9;   
    figure;
    fig  = gcf;
    set(fig,'color','white','position',[100 100 1800 600],'Visible','on');
    mlspres = 2;    
    if mlspres == 1
        wacpresind = 11;
        contourlevel = .275e-6;        
    elseif mlspres == 2
        wacpresind = 13;
        contourlevel = .275e-6;
    elseif mlspres == 3
        wacpresind = 17;
        contourlevel = 1.5e-6;
    end
    for i = 1:3
        datetoplot = str2double([num2str(year),sprintf('%02d',month),'01']);
        dateind = 12*11+month;
        sp(i) = subplot(1,3,i);
        sp_pos(i,:) = get(sp(i),'position');
        if i == 1

        elseif i == 3

        end
        axesm('MapProjection','stereo','MapLatLimit',[-90 -50])
        %h = geoshow('landareas.shp', 'FaceColor',[178,223,138]./255);%[0.6 .9 0.6]);
        h = geoshow('landareas.shp', 'FaceColor',[.80 .80 .80]);%[0.6 .9 0.6]);
        hold on
        [C0,MLS_c] = contourm(MLSlat(1:9),MLSlon,double(squeeze(MLSclim(month,1:9,:,mlspres))),[contourlevel contourlevel],'color','k','LineWidth',2);
        [C,MLS_s] = contourm(MLSlat(1:9),MLSlon,double(squeeze(MLS.O3_monthmean_latlon(12,month,1:9,:,mlspres))),[contourlevel contourlevel]','color',matstand(3,:),'LineWidth',2);
        [C1,MAM] = contourm(waclat1(1:22),waclon1,double(squeeze(WACCM_MolC_convolve.MAM(:,1:22,wacpresind,dateind)))',[contourlevel contourlevel],'--','color',cbrew(2,:),'LineWidth',2);
        [C2,VCMAM] = contourm(waclat1(1:22),waclon1,double(squeeze(WACCM_MolC_convolve.VCMAM(:,1:22,wacpresind,dateind)))',[contourlevel contourlevel],'-.','color',cbrew(1,:),'LineWidth',2);

        waclat = [-90,-70.7,-69,-68.6,-64.2]; % South Pole, Neumayer, Syowa, Davis, Marambio
        waclon = [158,351.7,39.6,78,303.4];


        %sites = [-70.7,351.7;-69.0,39.54;-68.6,77.97;-64.2,303.4];
        sites = [-70.7,351.7;-69.0,39.54;-68.6,77.97];

        h(1) = plotm(sites(1,1),sites(1,2),'Marker','pentagram','linestyle','none',...
            'MarkerEdgeColor',cbrew(3,:),'MarkerFaceColor',cbrew(3,:),'MarkerSize',18.5);
            set(gca,'fontsize',fontsize);

        h(2) = plotm(sites(2,1),sites(2,2),'Marker','s','linestyle','none',...
            'MarkerEdgeColor',cbrew(5,:),'MarkerFaceColor',cbrew(5,:),'MarkerSize',15);
            set(gca,'fontsize',fontsize);

        h(3) = plotm(sites(3,1),sites(3,2),'Marker','v','linestyle','none',...
            'MarkerEdgeColor',cbrew(4,:),'MarkerFaceColor',cbrew(4,:),'MarkerSize',15);
            set(gca,'fontsize',fontsize);

%         h(4) = plotm(sites(4,1),sites(4,2),'Marker','d','linestyle','none',...
%             'MarkerEdgeColor',cbrew(8,:),'MarkerFaceColor',cbrew(8,:),'MarkerSize',13.5);
%             set(gca,'fontsize',fontsize);

        setm(gca,'frame','on')
        set(gca','Visible','off')
        setm(gca,'parallellabel','on')
        %setm(gca,'meridianlabel','on','mlabellocation',360,'fontsize',fontsize-2)
        gridm('on');
        gridm('mlinelocation',90,'plinelocation',10,'plabellocation',[-90 -80 -70 -60],'fontsize',fontsize-10,'GLineWidth',2,'gcolor',[.7 .7 .7]);
        if i == 3
            lh = legend([MLS_c,MLS_s,MAM,VCMAM,h(1),h(2),h(3)],'MLS / OMI 2004-2014 mean','MLS / OMI','MAM','VC-MAM','Neumayer','Syowa','Davis');
            set(lh,'orientation','horizontal','position',[.49 .18 .05 .05],'box','off','fontsize',fontsize-4);
        end
        %GLineWidth
        
            month = month+1;
            %omimonth = omimonth+1;        

    end

    if mlspres == 1
        prestit = '146.78 hPa';
    elseif mlspres == 2
        prestit = '100 hPa';
    elseif mlspres == 3
        prestit = '46.42 hPa';
    end
    
    set(sp(1),'position',[sp_pos(1,1)+.1, sp_pos(1,2:4)]);
    set(sp(3),'position',[sp_pos(3,1)-.1, sp_pos(3,2:4)]);
    mtit('September','xoff',-.295,'yoff',-.2,'fontsize',fontsize+2);
    mtit('October','xoff',.005,'yoff',-.2,'fontsize',fontsize+2);
    mtit('November','xoff',+.305,'yoff',-.2,'fontsize',fontsize+2);
    mtit('2015 monthly mean ozone hole contours','xoff',0,'yoff',-.1,'fontsize',fontsize+6);

    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
    'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

    text(.225, .28,[sprintf('%.1f', MLSpressure(wacpresind)),' hPa < ',num2str(contourlevel*1e6),' ppmv'],'HorizontalAlignment'...
        ,'left','color','k','VerticalAlignment','top','rotation',90,'fontsize',fontsize);
   
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/areabelow220/','SepOctNov_',prestit,'_',num2str(year),'.png'];
    %export_fig(filename,'-pdf')
end