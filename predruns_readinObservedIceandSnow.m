function [data,obvar] = predruns_readinObservedIceandSnow(inputs,modeldim,ERAtemp)

if strcmp(inputs.obstouse,'seaiceindex')
    directory = '/Volumes/ExternalOne/work/data/NSIDC_seaice/seaiceindex/';
    files = dir([directory,'*.csv']);
    for i = 1:length(files)
        datatemp = importdata([directory,files(i).name]);
        data(:,:,i) = datatemp.data;
    end
    data (data == -9999) = NaN;
    obvar = [];
elseif strcmp(inputs.obstouse,'Goddard')
    directory = '/Volumes/ExternalOne/work/data/NSIDC_seaice/';
    files = dir([directory,'*.nc']);
    obvar = 'siconc';
    if ~exist('/Volumes/ExternalOne/work/data/NSIDC_seaice/output/GoddardMergedStandard_bootstrapped.mat')    
        for i = 1:length(files)
            [goddard(i).(obvar)] = ncread([directory,files(i).name],'goddard_bt_seaice_conc_monthly'); %goddard_merged_seaice_conc_monthly
            [goddard(i).time] = ncread([directory,files(i).name],'time');
            [goddard(i).latitude] = ncread([directory,files(i).name],'latitude');
            [goddard(i).longitude] = ncread([directory,files(i).name],'longitude');
            goddard(i).(obvar) (goddard(i).(obvar) < 0) = NaN;        
        end    
        
        %goddard = predruns_readinNDISC(datatype);
        data.years = repmat(1979:2017,[12,1]);
        data.timeperiod = 1979:2017;
        standardlatsbins = 45:1.5:90;
        standardlonsbins = 0:1.5:360;
        standardlats = 45.75:1.5:90;
        standardlons = .75:1.5:360;
        for k = 1:length(goddard)
            testdata = goddard(k).siconc(:);
            testlats = goddard(k).latitude(:);            
            testlons = goddard(k).longitude(:);
            testlonsind = (testlons < 0);
            testlons(testlonsind) = testlons(testlonsind) + 360;            
            for i = 1:length(standardlatsbins)-1
                for j = 1:length(standardlonsbins)-1
                    data.(obvar)(i,j,k) = nanmean(testdata(testlats >= standardlatsbins(i) & testlats < standardlatsbins(i+1) & testlons >= standardlonsbins(j) & testlons < standardlonsbins(j+1)));
                end
            end
        end
        data.latitude = standardlats';
        data.longitude = standardlons';
        data.(obvar) = permute(data.(obvar),[2,1,3]);
        %adding in missing years
        data.(obvar) = cat(3,data.(obvar)(:,:,1:107),ones(size(data.(obvar)(:,:,1:2)))*-9999,data.(obvar)(:,:,108:end));
        data.(obvar) (data.(obvar) == -9999) = NaN;
        data.(obvar) = cat(2,ones(size(data.(obvar),1),91,size(data.(obvar),3))*-9999,data.obvar);
        save('/Volumes/ExternalOne/work/data/NSIDC_seaice/output/GoddardMergedStandard_bootstrapped.mat','data');
        
    else
        load('/Volumes/ExternalOne/work/data/NSIDC_seaice/output/GoddardMergedStandard_bootstrapped.mat');        
    end
    
    data.latitude = -90:1.5:90;
    
elseif strcmp(inputs.obstouse,'ERA')
    data.timeperiod = 1979:2018;
    if strcmp(inputs.var,'ICEFRAC')
        directory = '/Volumes/ExternalOne/work/data/ERA-Interim/SEAICE/';   
        obvar = 'siconc';
    elseif strcmp(inputs.var,'SNOWHICE')
        directory = '/Volumes/ExternalOne/work/data/ERA-Interim/SNOWD/';
        obvar = 'sd';
    end
    files = dir([directory,'*.nc']);
    for i = 1:length(files)
        [data(i).(obvar)] = ncread([directory,files(i).name],obvar);
        [data(i).time] = ncread([directory,files(i).name],'time');
        [data(i).latitude] = ncread([directory,files(i).name],'latitude');
        [data(i).longitude] = ncread([directory,files(i).name],'longitude');
        data(i).(obvar) (data(i).(obvar) < 0) = NaN;
        data(i).(obvar) (data(i).(obvar) > 1) = NaN;
    end 
elseif strcmp(inputs.obstouse,'MERRA')
    data.timeperiod = 1980:2018;
    if strcmp(inputs.var,'ICEFRAC')
        directory = '/Volumes/ExternalOne/work/data/MERRA/ICEFRAC/';   
        obvar = 'FRSEAICE';
    elseif strcmp(inputs.var,'SNOWHICE')
        directory = '/Volumes/ExternalOne/work/data/MERRA/snowmass/';
        obvar = 'SNOMAS';
    elseif strcmp(inputs.var,'PRECSL')
        directory = '/Volumes/ExternalOne/work/data/MERRA/snowfall/';
        obvar = 'PRECSN';
    end
    files = dir([directory,'*.nc']);
    for i = 1:length(files)
        [data(i).(obvar)] = ncread([directory,files(i).name],obvar);
        [data(i).time] = ncread([directory,files(i).name],'time');
        [data(i).latitude] = ncread([directory,files(i).name],'lat');
        [data(i).longitude] = ncread([directory,files(i).name],'lon');
        data(i).(obvar) (data(i).(obvar) < 0) = NaN;
        if strcmp(inputs.var,'SNOWHICE')
            data(i).(obvar) = data(i).(obvar)./1000; % convert from kgm-2 to water equivalent snow depth
            %data(i).(obvar) = sqrt(1000/data(i).(obvar)); % convert from kgm-2 to water equivalent snow depth
        elseif strcmp(inputs.var,'PRECSL')
            data(i).(obvar) = data(i).(obvar)./1000; % convert from kgm-2 to water equivalent snow depth
            %data(i).(obvar) = sqrt(1000/data(i).(obvar)); % convert from kgm-2 to water equivalent snow depth
        end
    end    
else    
    
    directory = '/Volumes/ExternalOne/work/data/NSIDC_seaice/G10010_SIBT1850_v1.1/';
    files = dir([directory,'*.nc']);
    for i = 1:length(files)
        [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    end    
    obvar = [];
end

% put into yearly array
for i = 1:12
    if inputs.detrend
        for j = 1:length(data.latitude)
            %data.montharrange(i,:,j,:) = detrend(data.(obvar)(:,j,i:12:end))'+nanmean(squeeze(data.(obvar)(:,j,i:12:end)),2);    
            data.montharrange(i,:,j,:) = squeeze(detrendNaN3(data.(obvar)(:,j,i:12:end)))+nanmean(squeeze(data.(obvar)(:,j,i:12:end)),2);    
        end        
    else
        data.montharrange(i,:,:,:) = permute(data.(obvar)(:,:,i:12:end),[3,1,2]);    
    end
end
if inputs.detrend
    data.montharrange = permute(data.montharrange,[1,4,2,3]);
end

% Read in Bodeker total column ozone
%% READ in Bodeker 
tcolat = inputs.lats;
tozmonth = inputs.tozmonth;
data.toz.timeperiod = [1980,2016];
[~,BSdata,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/BodekerScientific/TCO/Bodeker_TCO_monavg.nc');

BSyears = 1980:2016;
BSyear_vector = repmat(BSyears,12,1);
BSdata.years = [BSyear_vector(:);ones(length(BSdata.time)-length(BSyear_vector(:)),1)*max(BSyear_vector(:))+1];

ozonedateindex(1) = find(BSdata.years == data.toz.timeperiod(1),1,'first');
ozonedateindex(2) = find(BSdata.years == data.toz.timeperiod(2),1,'last');

latindex_bs = find(BSdata.lat >= tcolat(1) & BSdata.lat <= tcolat(2));

data.toz.zm = detrend(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+...
    tozmonth-1:12:ozonedateindex(2)),1)),BSdata.lat(latindex_bs))) + ...
    nanmean(weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+...
    tozmonth-1:12:ozonedateindex(2)),1)),BSdata.lat(latindex_bs)));
data.toz.zm_nodetrend = weightedaverage(squeeze(nanmean(BSdata.tco(:,latindex_bs,ozonedateindex(1)+...
    tozmonth-1:12:ozonedateindex(2)),1)),BSdata.lat(latindex_bs));   

lowperc = prctile(data.toz.zm,inputs.percentile);
highperc = prctile(data.toz.zm,100 - inputs.percentile);

data.toz.lowind = find(data.toz.zm <= lowperc);
data.toz.highind = find(data.toz.zm >= highperc);


%% extract surface data time period

data.montharrange = data.montharrange(:,data.timeperiod >= data.toz.timeperiod(1) & data.timeperiod <= data.toz.timeperiod(2),:,:);

%% calculate ENSO
if inputs.removeENSO && ~exist(['/Volumes/ExternalOne/work/data/predruns/output/NINO34/',inputs.var,'_removeENSO_',inputs.obstouse,'.mat'])
    
    datamonthall = permute(data.montharrange,[3,4,1,2]);
    datamonthall = datamonthall(:,:,:);

    varmonth2 = inputs.varmonth;
    
    %ERAdateindex(1) = find(ERAdata.years == obstimeperiod(1),1,'first');
    %ERAdateindex(2) = find(ERAdata.years == obstimeperiod(2),1,'last');

    latlimits = [-5 5];
    lonlimits = [190 240];

    latindex = find(data.latitude >= latlimits(1) & data.latitude <= latlimits(2));
    lonindex = find(data.longitude >= lonlimits(1) & data.longitude <= lonlimits(2));

    for j = 1:12
        NINO_mn(:,j,:) = squeeze(nanmean(ERAtemp(lonindex,latindex,j:12:end),1));
        NINO_mn2(j,:) = squeeze(nanmean(NINO_mn(:,j,:),1));
        NINOmonth(j,:) = (NINO_mn2(j,:) - nanmean(NINO_mn2(j,:),2))./std(NINO_mn2(j,:),1,2);
    end   
    laglength = 3;
    NINO34all = NINOmonth(:);
    ext = 0;
    for k = 1:length(varmonth2)
        k
        for l = 1:size(datamonthall,1) % longitudes
            for m = 1:size(datamonthall,2) % latitudes
                for lag = 1:laglength % lag
               
                    regressors = [ones(length(squeeze(NINO34all(varmonth2(k)-lag+1:12:end-ext))),1),...
                        squeeze(NINO34all(varmonth2(k)-lag+1:12:end-ext))];



                    [b(k,l,m,lag,:)] = regress(squeeze(datamonthall(l,m,varmonth2(k):12:end-ext)),...
                        regressors); 

                                 
                    % finding largest lag correlation
                end
                
                [~,llc(k,l,m)] = max(abs(squeeze(b(k,l,m,:,2))));

      
                blag(k,l,m,:) = b(k,l,m,llc(k,l,m),:);
          
          
                dataVarMonth(:,k,l,m) = squeeze(datamonthall(l,m,varmonth2(k):12:end-ext)) - ...
                    squeeze(b(k,l,m,llc(k,l,m),2))*...
                    squeeze(NINO34all(varmonth2(k)-llc(k,l,m)+1:12:end-ext));                        



            end
        end
    end       
        
    dataVarMonth = permute(dataVarMonth,[2,1,3,4]);    
    ENSO = NINO34all;
    montharrange_re = cat(1,dataVarMonth(1:2,:,:,:),dataVarMonth);
    data.montharrange_re = cat(1,dataVarMonth(1:2,:,:,:),dataVarMonth);
    data.ENSO = NINO34all;
    save(['/Volumes/ExternalOne/work/data/predruns/output/NINO34/',inputs.var,'_removeENSO_',inputs.obstouse,'.mat'],'montharrange_re','ENSO');
elseif inputs.removeENSO && exist(['/Volumes/ExternalOne/work/data/predruns/output/NINO34/',inputs.var,'_removeENSO_',inputs.obstouse,'.mat'])
    load(['/Volumes/ExternalOne/work/data/predruns/output/NINO34/',inputs.var,'_removeENSO_',inputs.obstouse,'.mat']);
    data.montharrange_re = montharrange_re;
    data.ENSO = ENSO;
end

%% convert to model dims
for i = 1:size(data.montharrange,1)
    for j = 1:size(data.montharrange,2)
        for k = 1:size(data.montharrange,3)
            templatinterp(i,j,k,:) = interp1(data.latitude,squeeze(data.montharrange(i,j,k,:)),modeldim.lat);
        end
        for k = 1:size(templatinterp,4)
            data.montharrange_modeldim(i,j,:,k) = interp1(data.longitude,squeeze(templatinterp(i,j,:,k)),modeldim.lon);
        end
    end
end


end
