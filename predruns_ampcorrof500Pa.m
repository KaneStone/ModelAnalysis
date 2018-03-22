directory = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/data/zonalAnalysis/';
files = dir([directory,'*.mat']);
% user inputs

%% reading in data
for i = 1:length(files)
    nameofdata{i} = files(i).name(1:end-4);
    if i < length(files)
        temp = load([directory,files(i).name]);
        data.(nameofdata{i}) = temp.(nameofdata{i});
    else
        WACCM_predruns = load([directory,files(i).name]);
    end
    %nameofdatafortitles{i} = nameofdata{i}; 
    %nameofdatafortitles{i} (nameofdatafortitles{i} == '_') = ' ';
end

%%
% ClLevel = 'highCl';
% timeperiodhigh = [1995,2016];%[1955,1975]
% 
% vardirectory = ['/Volumes/MyBook/work/data/predruns/','TS','/',ClLevel,'/'];
% varfiles = dir([vardirectory,'*.nc']);
% 
% [highcl,~,~,~]= predruns_ReadInlayer(vardirectory,varfiles,'TS',timeperiodhigh,[-90,-75],0);

%% Read in 500 hPa temperature

temp500.lowcl = load('/Volumes/MyBook/work/data/predruns/output/T_500hPa_T.LowCl..mat');
temp500.highcl = load('/Volumes/MyBook/work/data/predruns/output/T_500hPa_THighCl..mat');
%%
fieldnames = fields(temp500);

for i = 1:length(fieldnames)
    for j = 1:12
        temp500month.(fieldnames{i})(:,:,:,:,j) = temp500.(fieldnames{i}).data_at_level(:,:,:,j:12:end);
    end
    temp500spring.(fieldnames{i}) = reshape(permute(nanmean(temp500month.(fieldnames{i})(:,:,:,:,9:11),5),[2,3,4,1]),...
        [144,96,size(temp500month.(fieldnames{i}),4)*size(temp500month.(fieldnames{i}),1)]);
    amplitudes.(fieldnames{i}) = WACCM_predruns.var_maxvalue.(fieldnames{i}) - WACCM_predruns.var_minvalue.(fieldnames{i});
    ampspring.(fieldnames{i}) = squeeze(nanmean(amplitudes.(fieldnames{i})(:,9:11,:),2));
end

%% find correlations
for k = 1:length(fieldnames)
    for i = 1:size(temp500spring.highcl,1)
        for j = 1:size(temp500spring.highcl,2)
            [r(k,i,j),p(k,i,j)] = corr(squeeze(temp500spring.(fieldnames{k})(i,j,:)),ampspring.(fieldnames{k})(6,:)');
        end
    end
end

%% plotting

lats = [-90;-88.1052631578947;-86.2105263157895;-84.3157894736842;-82.4210526315790;-80.5263157894737;...
    -78.6315789473684;-76.7368421052632;-74.8421052631579;-72.9473684210526;-71.0526315789474;...
    -69.1578947368421;-67.2631578947369;-65.3684210526316;-63.4736842105263;-61.5789473684211;...
    -59.6842105263158;-57.7894736842105;-55.8947368421053;-54;-52.1052631578947;-50.2105263157895;...
    -48.3157894736842;-46.4210526315790;-44.5263157894737;-42.6315789473684;-40.7368421052632;...
    -38.8421052631579;-36.9473684210526;-35.0526315789474;-33.1578947368421;-31.2631578947368;...
    -29.3684210526316;-27.4736842105263;-25.5789473684211;-23.6842105263158;-21.7894736842105;...
    -19.8947368421053;-18;-16.1052631578947;-14.2105263157895;-12.3157894736842;-10.4210526315789;...
    -8.52631578947369;-6.63157894736843;-4.73684210526317;-2.84210526315790;-0.947368421052630;...
    0.947368421052630;2.84210526315789;4.73684210526315;6.63157894736841;8.52631578947369;...
    10.4210526315789;12.3157894736842;14.2105263157895;16.1052631578947;18;19.8947368421053;...
    21.7894736842105;23.6842105263158;25.5789473684210;27.4736842105263;29.3684210526316;...
    31.2631578947368;33.1578947368421;35.0526315789474;36.9473684210526;38.8421052631579;...
    40.7368421052632;42.6315789473684;44.5263157894737;46.4210526315789;48.3157894736842;...
    50.2105263157895;52.1052631578947;54;55.8947368421053;57.7894736842105;59.6842105263158;...
    61.5789473684210;63.4736842105263;65.3684210526316;67.2631578947368;69.1578947368421;...
    71.0526315789474;72.9473684210526;74.8421052631579;76.7368421052632;78.6315789473684;...
    80.5263157894737;82.4210526315789;84.3157894736842;86.2105263157895;88.1052631578947;90];

%r = permute(r,[1,3,2]);
%rCanESM2 = reshape(rCanESM2,[1,size(rCanESM2)]);
cbrew = cbrewer('div','RdBu',16);         

contourtitle = {'test'};       

contourtitle2 = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)'};

subplotmaps(r,data.WACCM_CCMI.longitude,lats,{'div','RdBu'},1,[],18,contourtitle2,'Longitude','Pressure (hPa)','Correlation','on',...
    [-1,1],18,0:20:360,0:20:360,-90:10:90,-90:10:90,contourtitle,1,[0 360],[-90 -30],1,'-',0,'none');

% filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/',...
%     'CanESM2_Tempdiff_',num2str(abs(lats(latind))),'S_',num2str(perc),'_perc_spring'];
% export_fig(filename,'-pdf');
