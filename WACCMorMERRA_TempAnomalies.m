% Temperature normalized anomalies

load('/Users/kanestone/work/projects/WACCM/MatlabOutput/WACCM_MLS_compare.mat');
runs = fields(WACCM_MolC_convolve);


allyears = 2004:2015;
includeVCMAM = 1;
includeMLS = 1;

yearstoexclude = [2015];
deviation_year = 2015;
vert = 1;

for i = 1:length(runs)
    WACCM_Temp_convolve.(runs{i}) = permute(WACCM_Temp_convolve.(runs{i}),[2 1 3 4]);
    [WACCM_vertical_deviation_final(i,:,:,:), WACCM_vertical_deviation_final_percent(i,:,:,:), WACCMdeviationnorm(i,:,:,:,:)] = constructdeviations(...
        WACCM_Temp_convolve.(runs{i}),allyears,yearstoexclude,deviation_year,vert);     
end


%% plotting
prestick = [1000,500,200,100,50,20,10,5,2,1,.5,.2,.1];
logprestick = log(prestick);
posind = [.05,.02,.01];    

cbrew_vert = cbrewer('div','RdBu',9);
cbrew_vert1 = flipud(cbrew_vert);
temp_intervals = -2.7:.6:2.7;

fig = figure;   
set(fig,'color','white','position',[100 100 500 1000])

monthplot = [9,10,11];    

for i = 1:3
    sp = subplot(1,3,i);
    h = contourf(squeeze(WACCMdeviationnorm(2,i+8,13,:,:)));
end
