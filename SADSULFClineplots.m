% plot SADSULFC 1992 1993 and 2015 profiles


directory = '/Users/kanestone/work/projects/WACCM/netcdffiles/';

[~, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes, Longitudes,datedata] = ReadWACCMvertical('SAD_SULFC','monthly',directory,0);

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

% %%
% 
% [~,data1990,~] = Read_in_netcdf([directory,'SAD_SULFC_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.1990-2000.mam.004f.cam.h0zm_c160314.nc']);
% [~,data2015,~] = Read_in_netcdf([directory,'SAD_SULFC_species_f.e11.FWTREFC1SD.f19.f19.ccmi30.2015.mam.004f.cam.h0zm_new.nc']);
% 
%%

lats = [-50,-60,-70,-80,-90];
pressures = [100,150,70,50,10];
% 
for i = 1:length(lats)
    [~,latind(i)] = min(abs(Latitudes-lats(i)));
end

for i = 1:length(pressures)
    [~,presind(i)] = min(abs(unipres-pressures(i)));
end

%% plotting
fsize = 18;
for i = 1:length(lats)
    for j = 1:length(pressures)
        createfig('small','off');
        hold on
        plot(squeeze(vertical_uniform.MAM1990(latind(i),presind(j),13:28)),'LineWidth',3)
        plot(squeeze(vertical_uniform.MAM1990(latind(i),presind(j),25:40)),':','LineWidth',3)
        plot(squeeze(vertical_uniform.MAM1990(latind(i),presind(j),37:52)),'LineWidth',3)
        plot(squeeze(vertical_uniform.MAM(latind(i),presind(j),193:208)),'-.','LineWidth',3)
        set(gca,'fontsize',fsize);
        xlabel('Month','fontsize',fsize+2);
        set(gca,'xtick',1:1:100,'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});        
        ylabel('cm^2/cm^3','fontsize',fsize+2);
        title(['Chemical sulfur surface area density at ',num2str(abs(lats(i))),char(176),'S and ',num2str(pressures(j)),'hPa'],'fontsize',fsize+2);
        lh = legend('Pinatubo - 1991','Pinatubo - 1992','Pinatubo - 1993','Calbuco - 2015','location','NorthWest');
        set(lh,'fontsize',fsize+2,'box','off');
        set(gca,'box','on');
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/CalbucoPaper/NewFigures/CalbucoPinatuboNewFigures/',num2str(abs(lats(i))),'S','and',num2str(pressures(j))];
        export_fig(filename,'-pdf');
    end
end
