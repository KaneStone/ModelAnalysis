% compare ozonesonde vertical profiles to WACCM

% read in ozonesonde raw data
% remove datapoints with the same pressure values
% least square fit data to
clear all
% Read in WACCM
variable = 'O3';
commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
%commondir = '/Volumes/My Book for Mac/work/data/WACCM/';

[NumberDensity, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes,...
    Longitudes,datedata] = ReadWACCMvertical(variable,'monthly',commondir,0);

waclat2 = [-90,-70.7,-69,-68.6]; % South Pole, Neumayer, Syowa, Davis, Marambio
waclon2 = [158,351.7,39.6,78,303.4];
waclat = [-90,-70];

years = [1992,1993,1994,1996:2000,2012:2014,2015];

latmean = 0;
latlon = 1;
month = 10;
monthtitle = {'January','February','March','April','May','June','July','August',...
    'September','October','November','December'};

    runs = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
    month = 10;
    for i = 1:length(runs)
        for j = 1:length(waclat)
            [~,latindex(j)] = min(abs(Latitudes-waclat(j)));
            latextract(j).(runs{i}) = NumberDensity.(runs{i})(latindex(j),:,10:12:end);        
            latextractpres(j).(runs{i}) = Pressure.(runs{i})(latindex(j),:,10:12:end);        
        end    
    end

    for j = 1:2
        CCMItoplot(1,j,:) = squeeze(nanmean(latextract(j).CCMI(1,:,14:15),3));
        CCMItoplot(2,j,:) = squeeze(nanmean(latextract(j).CCMI(1,:,18:22),3));
        pCCMItoplot(1,j,:) = squeeze(nanmean(latextractpres(j).CCMI(1,:,14:15),3));    
        pCCMItoplot(2,j,:) = squeeze(nanmean(latextractpres(j).CCMI(1,:,18:22),3));
        MAMtoplot(1,j,:) = squeeze(nanmean(latextract(j).MAM(1,:,14:16),3));
        MAMtoplot(2,j,:) = squeeze(nanmean(latextract(j).MAM(1,:,17),3));
        pMAMtoplot(1,j,:) = squeeze(nanmean(latextractpres(j).MAM(1,:,14:16),3));
        pMAMtoplot(2,j,:) = squeeze(nanmean(latextractpres(j).MAM(1,:,17),3));
        %VCMAMtoplot(1,j,:) = squeeze(nanmean(latextract(j).VCMAM(1,:,14:16),3));
        %VCMAMtoplot(2,j,:) = squeeze(nanmean(latextract(j).VCMAM(1,:,15),3));
    end


if latlon
    %% read in daily data

    [MAMinfo, MAMdata, MAMatt] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/newfromDoug/latlon/O3_f.e11.FWTREFC1SD.f19.f19.ccmi30.mam.004f.cam.h5_merged_latlon.nc');
    [VCMAMinfo, VCMAMdata, VCMAMatt] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/newfromDoug/latlon/O3_f.e11.FWTREFC1SD.f19.f19.ccmi30.vc.mam.004f.cam.h5_merged_latlon.nc');
    [MAMTinfo, MAMTdata, MAMTatt] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/newfromDoug/latlon/T_f.e11.FWTREFC1SD.f19.f19.ccmi30.mam.004f.cam.h5_merged_latlon.nc');
    szMAM = size(MAMdata.O3);
    PressureMAM = permute(repmat(MAMdata.ap*100000,[1,szMAM([1,2,4])]),[2,3,1,4]) + permute(repmat(MAMdata.b,[1,szMAM([1,2,4])]),[2,3,1,4]).*permute(repmat(MAMdata.PS,[1,1,1,88]),[1,2,4,3]);      

    MAMlatlonconc = vmr2conc(MAMdata.O3,MAMTdata.T,PressureMAM,'03','conc');
    VCMAMlatlonconc = vmr2conc(VCMAMdata.O3,MAMTdata.T,PressureMAM,'03','conc');

    [MAM1990info, MAM1990data, MAM1990att] = Read_in_netcdf('/Volumes/MyBook/work/data/WACCM/O3/T_O3_f.e11.FWTREFC1SD.f19.f19.ccmi30.1990-2000.mam.004f.cam.h0.nc');
    sz1990MAM = size(MAM1990data.O3);
    Pressure1990MAM = permute(repmat(MAM1990data.hyam*100000,[1,sz1990MAM([1,2,4])]),[2,3,1,4]) + permute(repmat(MAM1990data.hybm,[1,sz1990MAM([1,2,4])]),[2,3,1,4]).*permute(repmat(MAM1990data.PS,[1,1,1,88]),[1,2,4,3]);      
    
    MAM1990latlonconc = vmr2conc(MAM1990data.O3,MAM1990data.T,Pressure1990MAM,'03','conc');
    
    MAMall = cat(4,MAM1990latlonconc(:,:,:,1:108),MAMlatlonconc);
    presall = cat(4,Pressure1990MAM(:,:,:,1:108),PressureMAM);
    datesall = [MAM1990data.date(1:108);MAMdata.date];
    for j = 1:length(waclat2)
            [~,latindex2(j)] = min(abs(MAMdata.lat-waclat2(j)));
            [~,lonindex2(j)] = min(abs(MAMdata.lon-waclon2(j)));
            
            latlonextract(1,j,:,:,:) = squeeze(nanmean(MAMall(lonindex2(j),latindex2(j),:,[34,46]),4));        
            latlonextract(2,j,:,:,:) = squeeze(nanmean(MAMall(lonindex2(j),latindex2(j),:,[82,94,106,118,130]),4));        
            latlonextract(3,j,:,:,:) = squeeze(nanmean(MAMall(lonindex2(j),latindex2(j),:,[274,286,298]),4));        
            latlonextract(4,j,:,:,:) = squeeze(nanmean(MAMall(lonindex2(j),latindex2(j),:,310),4));        
            
            latlonextractpres(1,j,:,:,:) = squeeze(nanmean(presall(lonindex2(j),latindex2(j),:,[34,46]),4));        
            latlonextractpres(2,j,:,:,:) = squeeze(nanmean(presall(lonindex2(j),latindex2(j),:,[82,94,106,118,130]),4));                                            
            latlonextractpres(3,j,:,:,:) = squeeze(nanmean(presall(lonindex2(j),latindex2(j),:,[274,286,298]),4));                    
            latlonextractpres(4,j,:,:,:) = squeeze(nanmean(presall(lonindex2(j),latindex2(j),:,310),4));        
    end    
    
end
%%
if latlon
    
titles = {'South Pole (90.0{\circ}S)','Neumayer (70.7{\circ}S)','Syowa (69.0{\circ}S)','Davis (68.6{\circ}S)'};   
fsize = 18;
standtickpres = [1000, 700, 500,300, 200, 150, 100, 70, 50, 30, 20, 10];

fig = figure;
set(fig,'color','white','position',[100 100 900 800]);
linestyles = {':','--','-.','-'};
markers = {'s','d','^','o'};
matstand = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

cbrew = cbrewer('qual','Paired',12);
ozcolors = [matstand(4,:);matstand(3,:);matstand(2,:);cbrew(4,:);];

for i = 1:4
    sh(i) = subplot(2,2,i);
    sh_pos(i,:) = get(sh(i),'position');
    hold on        
%     if i == 1
%             plot(squeeze(CCMItoplot(1,1,:))./1e12,log(squeeze(pCCMItoplot(1,1,:))./100),linestyles{1},...
%                 'LineWidth',3,'color',ozcolors(1,:));
%             plot(squeeze(CCMItoplot(2,1,:))./1e12,log(squeeze(pCCMItoplot(2,1,:))./100),linestyles{2},...
%                 'LineWidth',3,'color',ozcolors(2,:));
%     end
    for j = 1:4                
        plot(squeeze(latlonextract(j,i,:))./1e12,log(squeeze(latlonextractpres(j,i,:))./100),linestyles{j},...
            'LineWidth',3,'color',ozcolors(j,:));
        %plot(ozonetoplot(i,j).a./1e12,log(prestomatch(i,j).a),'Marker',markers{j},...
        %    'LineWidth',1,'color',ozcolors(j,:),'MarkerFaceColor',ozcolors(j,:),...
        %    'MarkerEdgeColor',ozcolors(j,:)./1.3,'MarkerSize',6)
%         hold on
%         plot(squeeze(latlonextractvc(j,i,:))./1e12,log(squeeze(latlonextractpres(j,i,:))./100),linestyles{j},...
%             'LineWidth',3,'color',ozcolors(j,:));
        set(gca,'Ydir','reverse','ytick',fliplr(log(standtickpres)),'yticklabel',fliplr(standtickpres))
        ylim([log(9) log(1000)]);
        
        xlim([-.1 4])                
    end
    set(gca,'fontsize',fsize-2,'box','on')
    title(titles{i},'fontsize',fsize+4)
    if i == 1 || i == 2
        set(sh(i),'position',[sh_pos(i,1),sh_pos(i,2)-.04,sh_pos(i,3:4)]);
    elseif i == 3 || i == 4
        set(sh(i),'position',[sh_pos(i,1),sh_pos(i,2)+.01,sh_pos(i,3:4)]);
    end
    
    sh_pos(i,:) = get(sh(i),'position');
    if i == 1 || i == 3
        set(sh(i),'position',[sh_pos(i,1)+.015,sh_pos(i,2),sh_pos(i,3:4)]);
    elseif i == 2 || i == 4
        set(sh(i),'position',[sh_pos(i,1)-.015,sh_pos(i,2),sh_pos(i,3:4)]);
    end
    
    if i == 1 || i == 3
        ylabel('Pressure (hPa)','fontsize',fsize+2)
    end
    if i == 3 || i == 4
        xlabel('Molecules/cm^3 ({\times}10^{12})','fontsize',fsize+2)
    end
    if i == 1
        lh = legend('1992-1993','1996-2000','2012-2014','2015');
        set(lh,'location','SouthEast','fontsize',fsize-2,'box','off');
    end
end
end

tb = annotation('textbox',[.40 .88 .1 .1],'String',[char(monthtitle{month}),' model'],'LineStyle','none','fontsize',fsize+10);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/Ozonesondes/','Modelcomparison_',char(monthtitle{month}),'_','Verticalprofiles','.pdf'];
export_fig(filename,'-pdf')

%% plotting
if latmean

titles = {'90.0 {\circ}S','69 {\circ}S'};   
fsize = 18;
standtickpres = [1000, 700, 500,300, 200, 150, 100, 70, 50, 30, 20, 10];

fig = figure;
set(fig,'color','white','position',[100 100 900 800]);
linestyles = {':','--','-.','-'};
markers = {'s','d','^','o'};
matstand = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

cbrew = cbrewer('qual','Paired',12);
ozcolors = [matstand(4,:);matstand(3,:);matstand(2,:);cbrew(4,:);];

sh(1) = subplot(2,1,1);
sh_pos(1,:) = get(sh(1),'position');

plot(squeeze(CCMItoplot(1,1,:))./1e12,log(squeeze(pCCMItoplot(1,1,:))./100),linestyles{1},...
    'LineWidth',3,'color',ozcolors(1,:));
hold on
plot(squeeze(CCMItoplot(2,1,:))./1e12,log(squeeze(pCCMItoplot(2,1,:))./100),linestyles{2},...
    'LineWidth',3,'color',ozcolors(2,:));
plot(squeeze(MAMtoplot(1,1,:))./1e12,log(squeeze(pMAMtoplot(1,1,:))./100),linestyles{3},...
    'LineWidth',3,'color',ozcolors(3,:));
plot(squeeze(MAMtoplot(2,1,:))./1e12,log(squeeze(pMAMtoplot(2,1,:))./100),linestyles{4},...
    'LineWidth',3,'color',ozcolors(4,:));
%plot(ozonetoplot(i,j).a./1e12,log(prestomatch(i,j).a),'Marker',markers{j},...
%    'LineWidth',1,'color',ozcolors(j,:),'MarkerFaceColor',ozcolors(j,:),...
%    'MarkerEdgeColor',ozcolors(j,:)./1.3,'MarkerSize',6)
hold on
set(gca,'Ydir','reverse','ytick',fliplr(log(standtickpres)),'yticklabel',fliplr(standtickpres))
ylim([log(9) log(1000)]);
%if month == 10
    xlim([-.1 4])
    
if i == 1
    lh = legend('1992-1993','1996-2000','2012-2014','2015');
    set(lh,'location','SouthEast','fontsize',fsize-2,'box','off');
end
    
sh(2) = subplot(2,1,2);    

sh_pos(2,:) = get(sh(2),'position');

plot(squeeze(CCMItoplot(1,2,:))./1e12,log(squeeze(pCCMItoplot(1,2,:))./100),linestyles{1},...
    'LineWidth',3,'color',ozcolors(1,:));
hold on
plot(squeeze(CCMItoplot(2,2,:))./1e12,log(squeeze(pCCMItoplot(2,2,:))./100),linestyles{2},...
    'LineWidth',3,'color',ozcolors(2,:));
plot(squeeze(MAMtoplot(1,2,:))./1e12,log(squeeze(pMAMtoplot(1,2,:))./100),linestyles{3},...
    'LineWidth',3,'color',ozcolors(3,:));
plot(squeeze(MAMtoplot(2,2,:))./1e12,log(squeeze(pMAMtoplot(2,2,:))./100),linestyles{4},...
    'LineWidth',3,'color',ozcolors(4,:));
%plot(ozonetoplot(i,j).a./1e12,log(prestomatch(i,j).a),'Marker',markers{j},...
%    'LineWidth',1,'color',ozcolors(j,:),'MarkerFaceColor',ozcolors(j,:),...
%    'MarkerEdgeColor',ozcolors(j,:)./1.3,'MarkerSize',6)
hold on
set(gca,'Ydir','reverse','ytick',fliplr(log(standtickpres)),'yticklabel',fliplr(standtickpres))
ylim([log(9) log(1000)]);
%if month == 10
    xlim([-.1 4])
    
    lh = legend('1992-1993','1996-2000','2012-2014','2015');
%         set(lh,'location','SouthEast','fontsize',fsize-2,'box','off');
end
    
%elseif month == 11            
%    xlim([-.1 4.5])
%end

% for i = 1:4
%     sh(i) = subplot(2,2,i);
%     sh_pos(i,:) = get(sh(i),'position');
%     for j = 1:4        
%         plot(ozonetoplot(i,j).a./1e12,log(prestomatch(i,j).a),linestyles{j},...
%             'LineWidth',3,'color',ozcolors(j,:));
%         %plot(ozonetoplot(i,j).a./1e12,log(prestomatch(i,j).a),'Marker',markers{j},...
%         %    'LineWidth',1,'color',ozcolors(j,:),'MarkerFaceColor',ozcolors(j,:),...
%         %    'MarkerEdgeColor',ozcolors(j,:)./1.3,'MarkerSize',6)
%         hold on
%         set(gca,'Ydir','reverse','ytick',fliplr(log(standtickpres)),'yticklabel',fliplr(standtickpres))
%         ylim([log(9) log(1000)]);
%         if month == 10
%             xlim([-.1 4])
%         elseif month == 11            
%             xlim([-.1 4.5])
%         end
%     end
%     set(gca,'fontsize',fsize-2)
%     title(titles{i},'fontsize',fsize+4)
%     if i == 1 || i == 2
%         set(sh(i),'position',[sh_pos(i,1),sh_pos(i,2)-.04,sh_pos(i,3:4)]);
%     elseif i == 3 || i == 4
%         set(sh(i),'position',[sh_pos(i,1),sh_pos(i,2)+.01,sh_pos(i,3:4)]);
%     end
%     
%     sh_pos(i,:) = get(sh(i),'position');
%     if i == 1 || i == 3
%         set(sh(i),'position',[sh_pos(i,1)+.015,sh_pos(i,2),sh_pos(i,3:4)]);
%     elseif i == 2 || i == 4
%         set(sh(i),'position',[sh_pos(i,1)-.015,sh_pos(i,2),sh_pos(i,3:4)]);
%     end
%     
%     if i == 1 || i == 3
%         ylabel('Pressure (hPa)','fontsize',fsize+2)
%     end
%     if i == 3 || i == 4
%         xlabel('Molecules/cm^3 ({\times}10^{12})','fontsize',fsize+2)
%     end
%     if i == 1
%         lh = legend('1992-1993','1996-2000','2012-2014','2015');
%         set(lh,'location','SouthEast','fontsize',fsize-2,'box','off');
%     end
% end
% tb = annotation('textbox',[.33 .88 .1 .1],'String',[char(monthtitle{month}),' ozonesondes'],'LineStyle','none','fontsize',fsize+10);
% 
% filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/Ozonesondes/',char(monthtitle{month}),'_','Verticalprofiles','.pdf'];
% export_fig(filename,'-pdf')