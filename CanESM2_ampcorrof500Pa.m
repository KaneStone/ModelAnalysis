
clear all
directory = '/Volumes/MyBook/work/data/CanESM2/400hPa/090S/strat/';
files = dir([directory,'*.nc']);
strat = 1;
%%
for i = 1:length(files)
    tic;
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    toc;
    tempspring(i).t = squeeze(nanmean(cat(5,data(i).ta(:,:,:,9:12:end),data(i).ta(:,:,:,10:12:end),data(i).ta(:,:,:,11:12:end)),5));
end

% GHG,highcl,lowcl
tempspringcat(3).m = cat(3,tempspring(1:50).t);
tempspringcat(1).m = cat(3,tempspring(51:100).t);
tempspringcat(2).m = cat(3,tempspring(101:150).t);

%% reading in amplitudes
directory2 = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/';
files2 = dir([directory2,'*.mat']);
Canyears = 1950:2100;
for i = 1:length(files2)
    nameofdata{i} = files2(i).name(1:end-4);
    if i < length(files2)
        temp = load([directory2,files2(i).name]);
        data2.(nameofdata{i}) = temp.(nameofdata{i});
    else
        WACCM_predruns = load([directory2,files2(i).name]);
    end
    nameofdatafortitles{i} = nameofdata{i}; 
    nameofdatafortitles{i} (nameofdatafortitles{i} == '_') = ' ';
end

%%
years = [1995,2024];
years2 = [1955,1979];
years3 = [2050,2079];
yearsind(1).m = find(Canyears >= years(1) & Canyears <= years(2));
yearsind(2).m = find(Canyears >= years2(1) & Canyears <= years2(2));
yearsind(3).m = find(Canyears >= years3(1) & Canyears <= years3(2));
if strat
    count = 102:2:200;
else
    count = 51:100;
end
for i = 1:length(yearsind)
    for j  = 1:50
        ampexttemp(j).m = data2.Can2ESM.amplitude_m_spr(count(j)).m(:,yearsind(i).m);       
        %ampext_strat(i,j).m = data2.Can2ESM.amplitude_m_spr(count2(i)).m(:,yearsind(i,:));       
    end
    ampextcat(i).m = cat(2,ampexttemp(:).m);
    clearvars ampexttemp
end

%ampextcat_strat = cat(2,ampext_strat(:).m);

%% taking correlations
latind = 4;
for k = 1:length(yearsind)
    for i = 1:size(tempspringcat(k).m)
        for j = 1:size(tempspringcat(k).m,2)
            [r(k,i,j),p(k,i,j)] = corr(squeeze(tempspringcat(k).m(i,j,:)),ampextcat(k).m(latind,:)');
            %[r(2,i,j),p(2,i,j)] = corr(squeeze(tempspringcatstrat(i,j,:)),ampextcat_strat(latind,:)');
        end
    end
end

%% plotting

%r = permute(r,[1,3,2]);
%rCanESM2 = reshape(rCanESM2,[1,size(rCanESM2)]);
cbrew = cbrewer('div','RdBu',16);         

contourtitle = {'400 hPa temperature correlations with 50 hPa 65S temperature amplitudes'};       

%contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','2050-2080 (high GHG - low chlorine)'};
contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','2050-2080 (low chlorine)'};

subplotmaps(r,data(1).lon,data(1).lat,{'div','RdBu'},1,[],18,contourtitle2,'Longitude','Pressure (hPa)','Correlation','on',...
    [-1,1],18,0:20:360,0:20:360,-90:10:90,-90:10:90,contourtitle,1,[0 360],[-90 0],1,'-',0,'none');

 filename = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/FinalPosterFigures/strat_Can400hPacorrelations';
 export_fig(filename,'-pdf');