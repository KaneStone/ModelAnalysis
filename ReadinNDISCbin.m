function [output] = ReadinNDISCbin

% Read in NSIDC daily data

if ~exist('/Volumes/ExternalOne/work/data/NSIDC_seaice/output/Smithrecreate.mat','file')
    years = 1979:2017;

    parentdirectory = '/Volumes/ExternalOne/work/data/NSIDC_seaice/daily/';

    datatemp = zeros(304,448,366);
    datatemp (datatemp == 0) = NaN;

    tic;
    for i = 1:length(years)
        toc;
        count = 1;
        directory = [parentdirectory,num2str(years(i)),'/'];
        files = dir([directory,'*.bin']);
        data(i).d = datatemp;
        for j = 1:length(files)
            % calculating day of year
            mon = files(j).name(8:9);
            da = files(j).name(10:11);
            year = files(j).name(4:7);
            dateve = [year,'-',mon,'-',da];
            DV = datevec(dateve);
            DV2 = DV;
            DV2(2:3) = 0;
            doy = datenum(DV) - datenum(DV2);                

            fileID = fopen([directory,files(j).name]);
            data(i).d(:,:,doy) = fread(fileID,[304,448],'int16',0,'l');
            fclose(fileID);

            count = count + 1;
        end
        %fill in missing day data
        datafilled(i).d = datatemp;
        for k = 1:size(data(i).d,1)
            datafilled(i).d(k,:,:) = fillmissing(squeeze(data(i).d(k,:,:))','linear')';
        end
    end

    % filling in missing data from Dec 3rd 87 to 12th Jan 88 from 1979-1992
    % daily average.
    dailyaverage_19791997 = nanmean(cat(4,datafilled(1:14).d),4);
    datafilled(9).d(:,:,337:end) = dailyaverage_19791997(:,:,337:end);
    datafilled(10).d(:,:,1:12) = dailyaverage_19791997(:,:,1:12);
        
    %%
    seaicedata2 = predruns_readinNDISC('Goddard');
    
    %%
    standardlatsbins = 45:1.5:90;
    standardlonsbins = 0:1.5:360;
    standardlats = 45.75:1.5:90;
    standardlons = .75:1.5:360;    
        
    tempstandardgrid = zeros(30,240,366);    
    for i = 1:length(datafilled)
        i
        tic;
        standardgrid(i).d = tempstandardgrid;
        %for k = 1:size(datafilled(i).d,3)
            datafilled(i).d (datafilled(i).d == 1200) = NaN;
            datafilled(i).d (datafilled(i).d == 1100) = NaN;
            testdata = permute(datafilled(i).d,[3,1,2]);
            testdata = testdata(:,:);
            testlats = seaicedata2(1).latitude(:);            
            testlons = seaicedata2(1).longitude(:);
            testlonsind = (testlons < 0);
            testlons(testlonsind) = testlons(testlonsind) + 360;

            for l = 1:length(standardlatsbins)-1
                for j = 1:length(standardlonsbins)-1
                    standardgrid(i).d(l,j,:) = nanmean(testdata(:,testlats >= standardlatsbins(l) & testlats < standardlatsbins(l+1) & testlons >= standardlonsbins(j) & testlons < standardlonsbins(j+1)),2);
                end
            end

            
        %end
        %standardgrid(i).d = uint16(standardgrid(i).d); 
        toc;
    end
    standardlats = standardlats';
    standardlons = standardlons';
    
    
    
    
    %%
    save('/Volumes/ExternalOne/work/data/NSIDC_seaice/output/Smithrecreate.mat','datafilled','standardgrid','years','standardlats','standardlons','-v7.3')
       
    %%
    
    output.datafilled = datafilled;
    output.standardgrid = standardgrid;
    output.years = years;
    output.latitude = standardlats;
    output.longitude = standardlons';
else
    output = load('/Volumes/ExternalOne/work/data/NSIDC_seaice/output/Smithrecreate.mat');
end

end