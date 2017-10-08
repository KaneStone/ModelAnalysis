% SAD_SULFC time_series
% Read in and plot SAD_SULFC only


[~, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes, datedata] = ReadWACCMvertical('SAD_SULFC','monthly');

runs = fields(MolConc);

%% vertical deviation

unipres = [1000:-50:150 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.05:.1];

for i = 1:length(runs) 
    for l = 1:size(MolConc.(runs{i}),3)
        for k = 1:size(MolConc.(runs{i}),1)
            vertical_uniform.(runs{i})(k,:,l) = interp1(log(squeeze(Pressure.(runs{i})(k,:,l))/100),...
                squeeze(MolConc.(runs{i})(k,:,l)),log(unipres));
        end        
        %rawdata(i,:,:) = squeeze(var20002014.(names{i})(rawlatindex,:,rawmonth:12:end)); 
    end            
end

years = [2004, 2014];
meanyears_mam = vertical_uniform.MAM(:,:,datedata.MAM.date > 20040301 & datedata.MAM.date <= 20150301);
meanyears_vcmam = vertical_uniform.VCMAM(:,:,datedata.VCMAM.date > 20040301 & datedata.VCMAM.date <= 20150301);

for i = 1:12
        clim_mam(:,:,i) = nanmean(meanyears_mam(:,:,i:12:end),3);
        clim_vcmam(:,:,i) = nanmean(meanyears_vcmam(:,:,i:12:end),3);                        
end

clim_mam = circshift(clim_mam,[0,0,2]);
clim_vcmam = circshift(clim_vcmam,[0,0,2]);  

latindex = find(Latitudes >= -90 & Latitudes <= -30);

MAMdeviation =  vertical_uniform.MAM(latindex,:,datedata.MAM.date > 20150101 & datedata.MAM.date <= 20151201) - clim_mam(latindex,:,1:11);
VCMAMdeviation =  vertical_uniform.VCMAM(latindex,:,datedata.VCMAM.date > 20150101 & datedata.VCMAM.date <= 20151201) - clim_vcmam(latindex,:,1:11);

%% plotting
createfig('medium');
set(gcf,'position',[10 10 900 650]);
cbrew_vert = cbrewer('seq','YlOrBr',14);
%cbrew_vert1 = flipud(cbrew_vert(10:25,:));

fsize = 18;
prestick = [1000 500 200 100 50 20 10 5 2 1];
%monthplot = [9,9,10,10,11,11];
monthplot = [6,6,7,7,8,8];
%daytitle = {'MAM-September','VCMAM-September','MAM-October','VCMAM-October','MAM-November','VCMAM-November'};
daytitle = {'MAM-June','VCMAM-June','MAM-July','VCMAM-July','MAM-August','VCMAM-August'};

plot_deviation = 0;

if plot_deviation
    toplotMAM = MAMdeviation;
    toplotVCMAM = VCMAMdeviation;
    clim = -2.2e-8:1e-8:11.8e-8;
else
    toplotMAM = vertical_uniform.MAM(latindex,:,datedata.MAM.date > 20150101 & datedata.MAM.date <= 20151201);
    toplotVCMAM = vertical_uniform.VCMAM(latindex,:,datedata.VCMAM.date > 20150101 & datedata.VCMAM.date <= 20151201);
    clim = -2.2e-8:2e-8:25.8e-8;
end

for i = 1:6;

    sp = subplot(3,2,i);
    sp_pos(i,:) = get(sp,'position');
    sp_pos_outer(i,:) = get(sp,'OuterPosition');

    if mod(i,2)
        contourf(Latitudes(latindex),log(unipres),toplotMAM(:,:,monthplot(i))',clim);    
    else contourf(Latitudes(latindex),log(unipres),toplotVCMAM(:,:,monthplot(i))',clim);    
    end
    caxis([min(clim) max(clim)]);
    set(gca,'YDir','reverse');    
    set(gca,'ytick',fliplr(log(prestick)),'yticklabel',fliplr(prestick),'fontsize',fsize);
    ylim([log(1) log(1000)]);

    set(gca,'color',[.8 .8 .8]);
    
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

    lastoutpos = get(sp,'outerposition');
  
    sp_pos(i,:) = get(sp,'position'); 
    
    if i == 6
        ch = colorbar;
        colormap(cbrew_vert);           
        caxis([min(clim) max(clim)]);
        set(ch,'Ticks',clim)
        set(get(ch,'ylabel'),'string','cm^2/cm^3','fontsize',fsize+2)
        cbaxloc = get(ch,'Position');
        subpos = get(sp,'OuterPosition');                    
            
        set(ch,'Position',[subpos(1)+subpos(3)+(subpos(1)+subpos(3))/25, cbaxloc(2),...
            cbaxloc(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(6,2)],...
            'Box','on','YAxisLocation','right','fontsize',fsize);
            
        set(ch,'YaxisLocation','right','fontsize',fsize);                  
    end
end

if plot_deviation
    mtit('WACCM 2015 anomaly (from 2004-2014 mean)','xoff',-.01,'yoff',.04,'fontsize',fsize+8);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
    'WACCM_SAD_SULFC_MAM-VCMAM_anomaly-SepNov','.png'];
else
    mtit('WACCM 2015 SAD SULFC','xoff',-.01,'yoff',.04,'fontsize',fsize+8);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
    'WACCM_SAD_SULFC_2015_JulAug','.png'];
end

export_fig(filename,'-png','-nofontswap','-r250');   

SULFCdeviations = struct('MAMdeviation',MAMdeviation,'VCMAMdeviation',VCMAMdeviation,...
    'Latitudes',Latitudes(latindex),'Pressure',unipres);
%save(['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
%    'WACCM_SAD_SULFC_MAM_anomaly.mat'],'SULFCdeviations');

%% plotting MAM-VCMAM
sulflog = 0;

MAMmVCMAM = (toplotMAM-toplotVCMAM)./toplotVCMAM*100;
MAMmVCMAM (MAMmVCMAM < 0) = NaN;
if sulflog
    MAMmVCMAM = log(MAMmVCMAM);
end

MAMmVCMAM (MAMmVCMAM > 240) = 240;

createfig('medium');
set(gcf,'position',[10 10 900 650]);
monthplot = [1,2,3,4,5,6];
%monthplot = [6,7,8,9,10,11];
daytitle = {'January','February','March','April','May','June','July','August','September','October','November','December'};
%clim = 0:75:900;
clim = [0:20:240];
cbrew_vert2 = cbrewer('seq','YlOrBr',length(clim)-1);
for i = 1:6;

    sp = subplot(3,2,i);
    sp_pos(i,:) = get(sp,'position');
    sp_pos_outer(i,:) = get(sp,'OuterPosition');

    contourf(Latitudes(latindex),log(unipres),MAMmVCMAM(:,:,monthplot(i))',clim);    
    caxis([min(clim) max(clim)]);
    set(gca,'YDir','reverse');    
    set(gca,'ytick',fliplr(log(prestick)),'yticklabel',fliplr(prestick),'fontsize',fsize);
    ylim([log(1) log(1000)]);

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

    lastoutpos = get(sp,'outerposition');
  
    sp_pos(i,:) = get(sp,'position'); 
    
    if i == 6
        ch = colorbar;
        colormap(cbrew_vert2);           
        caxis([min(clim) max(clim)]);
        set(ch,'Ticks',clim)
        set(get(ch,'ylabel'),'string','Percent','fontsize',fsize+2)
        cbaxloc = get(ch,'Position');
        subpos = get(sp,'OuterPosition');                    
            
        set(ch,'Position',[subpos(1)+subpos(3)+(subpos(1)+subpos(3))/25, cbaxloc(2),...
            cbaxloc(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(6,2)],...
            'Box','on','YAxisLocation','right','fontsize',fsize);
            
        set(ch,'YaxisLocation','right','fontsize',fsize);                  
    end
end


mtit('WACCM MAM-VCMAM 2015 SAD SULFC','xoff',-.01,'yoff',.04,'fontsize',fsize+8);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/2015anomaly/',...
    'WACCM_SAD_SULFC_2015_MAM-VCMAM_JanJun','.png'];


export_fig(filename,'-png','-nofontswap','-r250');   

