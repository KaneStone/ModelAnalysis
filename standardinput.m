clear all
directory = '/Volumes/MyBook/work/data/predruns/T/lowGHG/4060S/';
files = dir([directory,'*.nc']);
var = 'T';
tic;

%%

for i = 1:length(files)
    if i == 1
        [standardData(i).(var),regpres,latitudes,longitudes,years] = standardizedata([directory,files(i).name],var);
    else
        [standardData(i).(var),~,~,~,~] = standardizedata([directory,files(i).name],var);
    end
        
end

%%

for i = 1:length(files)
    for j = 1:size(standardData(i).(var),1)
        for k = 1:size(standardData(i).(var),3)
            standardwa(i).(var)(j,k,:) = weightedaverage(squeeze(standardData(i).(var)(j,:,k,:)),latitudes);
        end
    end
end

%%

outfile1 = ['/Volumes/MyBook/work/data/predruns/output/zonalanalysis/','lowGHG_standard.mat'];
outfile2 = ['/Volumes/MyBook/work/data/predruns/output/zonalanalysis/','lowGHG_standard_4060Swa.mat'];

save(outfile1,'standardData','latitudes','longitudes','years','regpres');
save(outfile2,'standardwa','longitudes','years','regpres');

toc;