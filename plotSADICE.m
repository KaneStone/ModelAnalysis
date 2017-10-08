% plot SAD ICE

commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[~,WACCMdata,~] = Read_in_netcdf('/Users/kanestone/work/projects/WACCM/netcdffiles/Extinction/Extinction_WACCM_MAM_201501-201609_zm.nc');
%[NumberDensity, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes, Longitudes,datedata] = ReadWACCMvertical(variable,'monthly',commondir,0);
[~, Extinction, Temperature, Pressure, Altitude, GMH, Latitudes, Longitudes,datedata] = ReadWACCMvertical('Extinction','monthly',commondir,0);

load('/Users/kanestone/work/projects/WACCM/MatlabOutput/WACCM_MLS_compare.mat');
runs = fields(WACCM_MolC_convolve);

MLSpressure = [1000;825.40417;681.29205;562.34131;464.15887;383.11868;316.22775;261.01572;...
    215.44347;177.82794;146.77992;121.15276;100;82.540421;68.129204;56.234131;46.415890;...
    38.311867;31.622776;26.101572;21.544348;17.782795;14.677993;12.115276;10;8.25404170...
    ;6.8129206;5.6234131;4.6415887;3.8311868;3.1622777;2.6101573;2.1544347;1.7782794;...
    1.4677993;1.2115277;1;0.68129206;0.46415889;0.31622776;0.21544346;0.14677992;...
    0.10000000;0.046415888;0.021544347;0.0099999998;0.0046415888;0.00215443480...
    ;0.0010000000;0.00046415889;0.00021544346;9.9999997e-05;4.6415887e-05;2.1544347e-05;9.9999997e-06];

for i = 1:length(runs)
    ICE.(runs{i}) = permute(SAD_ICElatMLSperiod.(runs{i}),[2,1,3,4]);        
end

figure
%contourf(Latitudes(1:34),log(MLSpressure),squeeze(ICE.MAM(10,end-1,:,:)))
contourf(squeeze(ICE.MAM(10,end-1,:,:)),10)
% figure
% contourf(squeeze(ICE.MAM(11,end-1,:,:)))
% figure
% contourf(squeeze(ICE.MAM(12,end-1,:,:)))
% figure
% contourf(squeeze(ICE.MAM(9,end-1,:,:)))

% cmapWACCM = cbrewer('seq','YlOrBr',17);
% cmapWACCM = cmapWACCM(2:17,:);
% sulf_intervals = 0:.25e-6:4e-6;
% sulf_intervals2 = 0:.25e-6:4e-6;
% %sulf_intervals = 0:.12e-6:.96e-6;
% exttimes = 1e6;
% exttimes2 = 1e3;
% fsize = 16;
% figure;
% set(gcf,'color','white','position',[100 100 1000 600]);
% count = 1;
% count1 = 1;
% titles = {'CALIPSO October','WACCM October','CALIPSO November','WACCM November'};
% for i = 1:4
%     sp = subplot(2,2,i);
%     
%     sp_pos(i,:) = get(sp,'position');
%     if i == 1 || i == 3
%         [~,spCAL(count1)] = contourf(Latitudes(1:34),log(MLSpressure),squeeze(datafinal_latinterp_presinterp(count1+4,:,:)),sulf_intervals);
%         set(gca,'color',[.8 .8 .8]);
%         ylim([log(10) log(300)])
%         set(gca,'ydir','reverse','ytick',fliplr(log([1000 500 250 150 50 20 10])),'yticklabel',fliplr([1000 500 250 150 50 20 10]));
%     end
% end