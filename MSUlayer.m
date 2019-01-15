% Read in and create MSU layer temperature
clear variables
run = 'highCl';
var = 'T';
%% Read in MSU weighting function
TMTlanddir = '/Volumes/ExternalOne/work/data/MSU/weightingfunctions/std_atmosphere_wt_function_chan_tmt_land.txt';
TMToceandir = '/Volumes/ExternalOne/work/data/MSU/weightingfunctions/std_atmosphere_wt_function_chan_tmt_ocean.txt';

[TMT.land,TMT.landsurface,TMT.headers] = ReadinMSUweightingfunctions(TMTlanddir);
[TMT.ocean,TMT.oceansurface,~] = ReadinMSUweightingfunctions(TMToceandir);
%% Read in landfrac
[~,landfrac,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/predruns/landfrac/LANDFRAC_b.e11.BWTREFC2.f19_g16.ccmi32.001.cam.h0.2097-08.nc');

%% Read in standardized highCl data individually
directory = ['/Volumes/ExternalOne/work/data/predruns/T/',run,'/raw/'];
TSdirectory = ['/Volumes/ExternalOne/work/data/predruns/TS/',run,'/'];

files = dir([directory,'*.nc']);
TSfiles = dir([TSdirectory,'*.nc']);

for i = 1:length(files)
    %read in temperature data
    [~,Tdata,~] = Read_in_netcdf([directory,files(i).name]);
    
    %read in surface pressure data
    
    pressure = permute(repmat(Tdata.hyai.*Tdata.P0,[1,size(Tdata.(var),1),size(Tdata.(var),2),...
        size(Tdata.(var),4)]) + repmat(Tdata.hybi,[1,size(Tdata.(var),1),size(Tdata.(var),2),...
        size(Tdata.(var),4)]) .* permute(repmat(Tdata.PS,[1,1,1,length(Tdata.ilev)]),[4,1,2,3]),[2,3,1,4])./100;    
    %read in surface temperature
    [~,TSdata,~] = Read_in_netcdf([TSdirectory,TSfiles(i).name]);
    
    
    
%% calculate weight over each layer    
%Wmp = zeros(size(Tdata.T));
Tweight = zeros(size(Tdata.T,1),size(Tdata.T,2),size(Tdata.T,4));
% w2_ind = zeros(size(Tdata.T));
% w1_ind = zeros(size(Tdata.T));

for j = 1:size(Tdata.T,1) %lons        
    tic;
    for k = 1:size(Tdata.T,2) %lats
        island = 0;
        island (landfrac.LANDFRAC(j,k) > .5) = 1;     
        if island
            w = TMT.land(:,6);
            wsurface = TMT.landsurface;
            wp = TMT.land(:,4)./100;
        else
            w = TMT.ocean(:,6);
            wsurface = TMT.oceansurface;
            wp = TMT.ocean(:,4)./100;
        end        
        for m = 1:size(Tdata.T,4)
            Wmp = zeros(size(pressure,3)-1,1);
            for l = 1:size(pressure,3)-1 %levels - 1                                                        
                [~,w2_ind] = min(abs(wp-pressure(j,k,l,m)));
                [~,w1_ind] = min(abs(wp-pressure(j,k,l+1,m)));
                if w2_ind == w1_ind
                    %Wmp(j,k,l,m) = 0;
                    Wmp(l) = 0;
                elseif w1_ind == 1
                    %Wmp(j,k,l,m) = w(1);
                    Wmp(l) = w(1);
                else  
                    %Wmp(j,k,l,m) = 1./(wp(w1_ind)-wp(w2_ind)+1).*sum(w(w1_ind:w2_ind).*(wp(w1_ind-1:w2_ind-1)-wp(w1_ind:w2_ind)));
                    Wmp(l) = 1./(wp(w1_ind)-wp(w2_ind)+1).*sum(w(w1_ind:w2_ind).*(wp(w1_ind-1:w2_ind-1)-wp(w1_ind:w2_ind)));
                end            
            end
            %Tweight(j,k,m) = squeeze(sum(Wmp(j,k,:,m).*Tdata.T(j,k,:,m),3))./sum(Wmp(j,k,:,m))+wsurface.*TSdata.TS(j,k,m); 
            Tweight(i,j,k,m) = sum(Wmp.*squeeze(Tdata.T(j,k,:,m)))./sum(Wmp)+wsurface.*TSdata.TS(j,k,m); 
        end        
    end
    toc;
end

end
save('/Volumes/ExternalOne/work/data/predruns/output/MSU/TMT.mat','Tweight');

    %%              
            
    
    
