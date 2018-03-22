clear all

directory = '/Volumes/MyBook/work/data/LENS/500hPa/highcl/';
files = dir([directory,'*.nc']);

%%

for i = 1:length(files)
    tic;
    [~,data(i),~] = Read_in_netcdf([directory,files(i).name]);
    toc;
    tempspring(i).t = squeeze(nanmean(cat(5,data(i).T(:,:,:,9:12:end),data(i).T(:,:,:,10:12:end),data(i).T(:,:,:,11:12:end)),5));
end

tempspringcat = cat(3,tempspring(:).t);

%% reading in amplitudes

directory2 = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/';
files2 = dir([directory2,'*.mat']);
LENSyears = 1920:2100;
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
yearsind = find(LENSyears >= years(1) & LENSyears <= years(2));
for i  = 1:30
    ampext(i).m = data2.LENS.amplitude_m_spr(i).m(:,yearsind);
        
end

ampextcat = cat(2,ampext(:).m);

%% taking correlations
latind = 4;
for i = 1:size(tempspringcat,1)
    for j = 1:size(tempspringcat,2)
        r(1,i,j) = corr(squeeze(tempspringcat(i,j,:)),ampextcat(latind,:)');
    end
end
  

%% plotting

%r = permute(r,[1,3,2]);
%rCanESM2 = reshape(rCanESM2,[1,size(rCanESM2)]);
cbrew = cbrewer('div','RdBu',16);         

contourtitle = {'test'};       

contourtitle2 = {'1995-2024 (high chlorine)'};

subplotmaps(r,data(1).lon,data(1).lat,{'div','RdBu'},1,[],18,contourtitle2,'Longitude','Pressure (hPa)','Correlation','on',...
    [-1,1],18,0:20:360,0:20:360,-90:10:90,-90:10:90,contourtitle,1,[0 360],[-90 -30],1,'-',0,'none');