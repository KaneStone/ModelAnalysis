% plot toz and temperature longitude lines
clear all
%% Read in temperature and toz
tozvar = 'toz';
lats = [-90 -60];
detrend_ozone = 0;
sig = .05;

%% plotting parameters
fsize = 18;

tic;

%% Read in TOZ highcl
ClLevel = 'highCl';
tozdates = [1995,2024];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfiles = dir([directory,'*.nc']);
[toz_data.highcl,toz_years.highcl,toz_varweighted.highcl,toz_composite.highcl,toz_dataMonthArrange.highcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfiles,tozvar,tozdates,lats,detrend_ozone);

%% Read in TOZ lowGHG
ClLevel = 'lowGHG';
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfilesGHG = dir([directory,'*.nc']);
[toz_data.lowGHG,toz_years.lowGHG,toz_varweighted.lowGHG,toz_composite.lowGHG,toz_dataMonthArrange.lowGHG] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfilesGHG,tozvar,tozdates,lats,detrend_ozone);

latitudes = toz_data.lowGHG.lat;
longitudes = toz_data.lowGHG.lon;

%% Read in TOZ lowcl
ClLevel = 'lowCl';
tozpastdates = [1955,1979];
directory = ['/Volumes/ExternalOne/work/data/predruns/',tozvar,'/',ClLevel,'/'];
tozfilespast = dir([directory,'*.nc']);
[toz_data.lowcl,toz_years.lowcl,toz_varweighted.lowcl,toz_composite.lowcl,toz_dataMonthArrange.lowcl] = ...
    predruns_ReadInlayer_areaaverage(directory,tozfilespast,tozvar,tozpastdates,lats,detrend_ozone);

latitudes = toz_data.lowcl.lat;
longitudes = toz_data.lowcl.lon;

%% Read in highcl Temperature at 50hPa
var = 'T';
ClLevel = 'highCl';
timeperiodhigh = [1995,2024];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/50hPa/'];
varfilespast = dir([vardirectory,'*.nc']);
[data.highcl,years.highcl,composite.highcl,dataMonthArrange.highcl]...
    = predruns_ReadInlayer(vardirectory,varfilespast,var,timeperiodhigh,lats,1);

%% Read in lowGHG Temperature at 50hPa
var = 'T';
ClLevel = 'lowGHG';
timeperiodhigh = [1995,2024];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/50hPa/'];
varfilesGHG = dir([vardirectory,'*.nc']);
[data.lowGHG,years.lowGHG,composite.lowGHG,dataMonthArrange.lowGHG]...
    = predruns_ReadInlayer(vardirectory,varfilesGHG,var,timeperiodhigh,lats,1);

%% Read in lowcl Temperature at 50hPa
var = 'T';
ClLevel = 'lowCl';
timeperiodhigh = [1955,1979];%[1955,1975]

vardirectory = ['/Volumes/ExternalOne/work/data/predruns/',var,'/',ClLevel,'/50hPa/'];
varfilespast = dir([vardirectory,'*.nc']);
[data.lowcl,years.lowcl,composite.lowcl,dataMonthArrange.lowcl]...
    = predruns_ReadInlayer(vardirectory,varfilespast,var,timeperiodhigh,lats,1);

fieldnames = fields(data);

%% Read in zonal average 4060S data

wa4060S.highcl = load('/Volumes/ExternalOne/work/data/predruns/output/zonalanalysis/highCl_standard_4060Swa.mat');
wa4060S.lowcl = load('/Volumes/ExternalOne/work/data/predruns/output/zonalanalysis/lowCl_standard_4060Swa.mat');
wa4060S.lowGHG = load('/Volumes/ExternalOne/work/data/predruns/output/zonalanalysis/lowGHG_standard_4060Swa.mat');
regpres = wa4060S.lowcl.regpres;

%% create composite for correlations
month = 10;
for i = 1:length(fieldnames)
    % month arrange first
    for j = 1:length(wa4060S.(fieldnames{i}).standardwa)
        for k = 1:12
            montharrange.(fieldnames{i})(j).T(:,:,k,:) = wa4060S.(fieldnames{i}).standardwa(j).T(:,:,k:12:end);
        end
    end
    wa4060Scomposite.(fieldnames{i}) = cat(4,montharrange.(fieldnames{i})(:).T);    
end

%% extracting data to plot
lats = [-50,-55,-60,-65,-70,-75,-80,-85];


for i = 1:length(fieldnames)
    month = 1;
for m = 1:12

    % extracting latitude and month from toz and temperature data
    for l = 1:length(lats)
        [~,latind] = min(abs(lats(l)-latitudes));
        for j = 1:length(toz_data.(fieldnames{i}))
            toz_extract.(fieldnames{i})(j).d(l,:,:) = squeeze(toz_data.(fieldnames{i})(j).toz(:,latind,month:12:end));
            var_extract.(fieldnames{i})(j).d(l,:,:) = squeeze(data.(fieldnames{i})(j).(var)(:,latind,1,month:12:end));
            toz_extracteddy.(fieldnames{i})(j).d(l,:,:) = squeeze(toz_data.(fieldnames{i})(j).toz(:,latind,month:12:end)) ...
                - repmat(nanmean(squeeze(toz_data.(fieldnames{i})(j).toz(:,latind,month:12:end))),[144,1]);
            var_extracteddy.(fieldnames{i})(j).d(l,:,:) = squeeze(data.(fieldnames{i})(j).(var)(:,latind,1,month:12:end)) ...
                - repmat(nanmean(squeeze(data.(fieldnames{i})(j).(var)(:,latind,1,month:12:end))),[144,1]);    
            
            var_extract_2.(fieldnames{i})(j).d(l,:,m,:) = squeeze(data.(fieldnames{i})(j).(var)(:,latind,1,month:12:end));
            
        end
    end
        
        % creating composites
        toz_extract_all.(fieldnames{i}) = cat(3,toz_extract.(fieldnames{i})(:).d);        
        toz_extract_allmean(i,:,:) = nanmean(toz_extract_all.(fieldnames{i}),3);
        
        var_extract_all.(fieldnames{i}) = cat(3,var_extract.(fieldnames{i})(:).d);
        var_extract_allmean(i,:,:) = nanmean(var_extract_all.(fieldnames{i}),3);        
    
        for j = 1:length(var_extract.(fieldnames{i}))
            [var_minvalue_m.(fieldnames{i})(j,:,m,:),~] = min(var_extract.(fieldnames{i})(j).d,[],2); 
            [var_maxvalue_m.(fieldnames{i})(j,:,m,:),~] = max(var_extract.(fieldnames{i})(j).d,[],2); 
        end        
        
        % find min longitude values for composite
        [minvalue.(fieldnames{i})(:,m,:),minlongind] = min(toz_extract_all.(fieldnames{i}),[],2); 
        minlong.(fieldnames{i})(:,m,:) = longitudes(squeeze(minlongind));  
        minlong_std.(fieldnames{i})(m,:) = std(squeeze(longitudes(minlongind)),0,2);
        
        [var_minvalue.(fieldnames{i})(:,m,:),var_minlongind] = min(var_extract_all.(fieldnames{i}),[],2);    
        var_minlong.(fieldnames{i})(:,m,:) = longitudes(squeeze(var_minlongind));
        var_minlong_std.(fieldnames{i})(m,:) = std(squeeze(longitudes(var_minlongind)),0,2);

        for j = 1:length(minlong.(fieldnames{i}))
            for l = 1:length(lats)
                if minlong.(fieldnames{i})(l,m,j) < 150 
                    minlong.(fieldnames{i})(l,m,j) = minlong.(fieldnames{i})(l,m,j)+360;
                end      
                if var_minlong.(fieldnames{i})(l,m,j) < 150 
                    var_minlong.(fieldnames{i})(l,m,j) = var_minlong.(fieldnames{i})(l,m,j)+360;
                end        
            end
        end

        % find max longitude values for composite
        
        [maxvalue.(fieldnames{i})(:,m,:),maxlongind] = max(toz_extract_all.(fieldnames{i}),[],2); 
        maxlong.(fieldnames{i})(:,m,:) = longitudes(squeeze(maxlongind));  
        maxlong_std.(fieldnames{i})(m,:) = std(squeeze(longitudes(maxlongind)),0,2);
        
        [var_maxvalue.(fieldnames{i})(:,m,:),var_maxlongind] = max(var_extract_all.(fieldnames{i}),[],2);    
        var_maxlong.(fieldnames{i})(:,m,:) = longitudes(squeeze(var_maxlongind));
        var_maxlong_std.(fieldnames{i})(m,:) = std(squeeze(longitudes(var_maxlongind)),0,2);
        
        % take mean of individual minimums and maximums
        minlong_mean.(fieldnames{i}) = nanmean(minlong.(fieldnames{i}),3);
        var_minlong_mean.(fieldnames{i}) = nanmean(var_minlong.(fieldnames{i}),3);
        maxlong_mean.(fieldnames{i}) = nanmean(maxlong.(fieldnames{i}),3);
        var_maxlong_mean.(fieldnames{i}) = nanmean(var_maxlong.(fieldnames{i}),3);
        
        % taking min and max of mean of composite        
        [minvalue2.(fieldnames{i})(m,:),minlong_mean2_ind] = min(toz_extract_allmean(i,:,:),[],3);
        minlong_mean2.(fieldnames{i})(m,:) = longitudes(minlong_mean2_ind);
        [var_minvalue2.(fieldnames{i})(m,:),var_minlong_mean2_ind] = min(var_extract_allmean(i,:,:),[],3);
        var_minlong_mean2.(fieldnames{i})(m,:) = longitudes(var_minlong_mean2_ind);
                
        [maxvalue2.(fieldnames{i})(m,:),maxlong_mean2_ind] = max(toz_extract_allmean(i,:,:),[],3);
        maxlong_mean2.(fieldnames{i})(m,:) = longitudes(maxlong_mean2_ind);
        [var_maxvalue2.(fieldnames{i})(m,:),var_maxlong_mean2_ind] = max(var_extract_allmean(i,:,:),[],3);
        var_maxlong_mean2.(fieldnames{i})(m,:) = longitudes(var_maxlong_mean2_ind);

        month = month + 1;
end

end

%% save output so I can compare to others
textfilename = '/Volumes/ExternalOne/work/data/matlabOutput/zonalAnalysis/WACCM_predruns.mat';
save(textfilename,'minvalue','maxvalue','minlong','maxlong','minlong_std','maxlong_std',...
    'minlong_mean','maxlong_mean','minvalue2','maxvalue2','minlong_mean2','maxlong_mean2',...
    'var_minvalue','var_maxvalue','var_minlong','var_maxlong','var_minlong_std','var_maxlong_std',...
    'var_minlong_mean','var_maxlong_mean','var_minvalue2','var_maxvalue2','var_minlong_mean2','var_maxlong_mean2');


%% plot amplitudes
latindamp = 4;
cbrew2 = cbrewer('qual','Paired',10);
highclcat = cat(4,var_extract_2.highcl(:).d);
lowclcat = cat(4,var_extract_2.lowcl(:).d);
lowGHGcat = cat(4,var_extract_2.lowGHG(:).d);

highclcat_spring = squeeze(nanmean(cat(5,highclcat(latindamp,:,9,:),highclcat(latindamp,:,10,:),highclcat(4,:,11,:)),5));
lowclcat_spring = squeeze(nanmean(cat(5,lowclcat(latindamp,:,9,:),lowclcat(latindamp,:,10,:),lowclcat(4,:,11,:)),5));
lowGHGcat_spring = squeeze(nanmean(cat(5,lowGHGcat(latindamp,:,9,:),lowGHGcat(latindamp,:,10,:),lowGHGcat(4,:,11,:)),5));

%% plotting
fsize = 20;
createfig('medium','on');

plot(longitudes,highclcat_spring,'color',cbrew2(1,:))
hold on
plot(longitudes,lowclcat_spring,'color',cbrew2(3,:))
plot(longitudes,lowGHGcat_spring,'color',cbrew2(5,:))
ph(1) = plot(longitudes,nanmean(highclcat_spring,2),'color',cbrew2(2,:),'LineWidth',5);
ph(3) = plot(longitudes,nanmean(lowGHGcat_spring,2),'--','color',cbrew2(4,:),'LineWidth',5);
ph(2) = plot(longitudes,nanmean(lowclcat_spring,2),':','color',cbrew2(6,:),'LineWidth',5);
set(gca,'fontsize',fsize+2);
xlim([-5 365]);
xlabel('Longitude ({\circ}E)','fontsize',fsize+4);
ylabel('Temperature (K)','fontsize',fsize+4);
title(['Austral spring temperatures at ',num2str(abs(lats(latindamp))),'{\circ}S and 50 hPa'],'fontsize',fsize+6);
lh = legend(ph,'1995-2025 (high chlorine)','1955-1979 (low chlorine)','1995-2025 (high chlorine - low GHG)');
set(lh,'fontsize',fsize+4,'box','off','location','south')

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TalkFigures/predrunsAmplitudes_',num2str(abs(lats(latindamp))),'S'];
export_fig(filename,'-pdf');
%%
amplowGHG = max(nanmean(lowGHGcat_spring,2)) - min(nanmean(lowGHGcat_spring,2));
amplowcl = max(nanmean(lowclcat_spring,2)) - min(nanmean(lowclcat_spring,2));
amphighcl = max(nanmean(highclcat_spring,2)) - min(nanmean(highclcat_spring,2));

stdlowGHG = std(max(lowGHGcat_spring,[],1) - min(lowGHGcat_spring,[],1),0);
stdhighcl = std(max(highclcat_spring,[],1) - min(highclcat_spring,[],1),0);
stdlowcl = std(max(lowclcat_spring,[],1) - min(lowclcat_spring,[],1),0);

%% plot min and max locations and standard deviations on a map
cbrew = cbrewer('qual','Set1',10);
createfig('large','on');
m_proj('Stereographic','lon',-180,'lat',-90,'rad',45)
m_coast('color','k','LineWidth',3);
m_grid('ytick',[-90 -80 -70 -60 -50 -40 -30 -20],'XaxisLocation','top','fontsize',fsize);
for i = 1:length(fieldnames)
    m_line(squeeze(minlong_mean2.(fieldnames{i})(10,:))',lats','color',cbrew(i,:),'LineWidth',3);    
    m_line(squeeze(maxlong_mean2.(fieldnames{i})(10,:))',lats','color',cbrew(i,:),'LineWidth',3);
    
    m_line(squeeze(var_minlong_mean2.(fieldnames{i})(10,:))',lats','color',cbrew(i,:),'LineWidth',3,'LineStyle','--');    
    m_line(squeeze(var_maxlong_mean2.(fieldnames{i})(10,:))',lats','color',cbrew(i,:),'LineWidth',3,'LineStyle','--');
    
    %m_line(squeeze(minlong_mean.(fieldnames{i})(:,10)),lats','color',cbrew(i,:),'LineWidth',3);    
    %m_line(squeeze(maxlong_mean.(fieldnames{i})(:,10)),lats','color',cbrew(i,:),'LineWidth',3);
    
end

%% plotting minimum locations


for i = 1:length(fieldnames)    
    plot(squeeze(minlong.(fieldnames{i})(latindamp,10,:)),'LineWidth',3,'color',cbrew(i,:))    
end

%% correlate with 4060S zonal T composite;
month = 10;
% for i = 1:length(fieldnames)
%     amplitude.(fieldnames{i}) = maxvalue.(fieldnames{i}) - minvalue.(fieldnames{i});
%     for j = 1:size(wa4060Scomposite.(fieldnames{i}),2)
%         r(i,:,j) = corr(squeeze(wa4060Scomposite.(fieldnames{i})(:,j,month,:))',squeeze(minlong.(fieldnames{i})(7,10,:)));
%         rmax(i,:,j) = corr(squeeze(wa4060Scomposite.(fieldnames{i})(:,j,month,:))',squeeze(maxlong.(fieldnames{i})(7,10,:)));
%         ramp(i,:,j) = corr(squeeze(wa4060Scomposite.(fieldnames{i})(:,j,month,:))',squeeze(amplitude.(fieldnames{i})(latindamp,10,:)));
%         %r_eddy(i,:,j) = corr(squeeze(wa4060Scomposite.(fieldnames{i})(:,j,month,:))',minlong_eddy.(fieldnames{i}));
%     end
% end

for i = 1:length(fieldnames)
    var_amplitude.(fieldnames{i}) = var_maxvalue.(fieldnames{i}) - var_minvalue.(fieldnames{i});
    for j = 1:size(wa4060Scomposite.(fieldnames{i}),2)
        [rvar(i,:,j),pvar(i,:,j)] = corr(squeeze(nanmean(wa4060Scomposite.(fieldnames{i})(:,j,9:11,:),3))',squeeze(nanmean(var_minlong.(fieldnames{i})(latindamp,9:11,:),2)));
        [rmaxvar(i,:,j),pmaxvar(i,:,j)] = corr(squeeze(nanmean(wa4060Scomposite.(fieldnames{i})(:,j,9:11,:),3))',squeeze(nanmean(var_maxlong.(fieldnames{i})(latindamp,9:11,:),2)));
        [rampvar(i,:,j),pampvar(i,:,j)] = corr(squeeze(nanmean(wa4060Scomposite.(fieldnames{i})(:,j,9:11,:),3))',squeeze(nanmean(var_amplitude.(fieldnames{i})(latindamp,9:11,:),2)));
        %rvar_eddy(i,:,j) = corr(squeeze(wa4060Scomposite.(fieldnames{i})(:,j,month,:))',var_minlong_eddy.(fieldnames{i}));
    end
end

%% correlate with 4060S zonal T individuals;

for i = 1:length(fieldnames)            
    for k = 1:size(var_maxvalue_m.(fieldnames{i}),1)
        %var_amplitudetemp = var_maxvalue.(fieldnames{i})(:,:,count:count+countadd) - var_minvalue.(fieldnames{i})(:,:,count:count+countadd);    
        var_amplitudetemp = squeeze(var_maxvalue_m.(fieldnames{i})(k,:,:,:) - var_minvalue_m.(fieldnames{i})(k,:,:,:));
        for j = 1:size(wa4060Scomposite.(fieldnames{i}),2)        
            %rvar(i,:,j) = corr(squeeze(wa4060Scomposite.(fieldnames{i})(:,j,month,:))',squeeze(var_minlong.(fieldnames{i})(latindamp,10,:)));
            %rmaxvar(i,:,j) = corr(squeeze(wa4060Scomposite.(fieldnames{i})(:,j,month,:))',squeeze(var_maxlong.(fieldnames{i})(latindamp,10,:)));
            temp = squeeze(nanmean(squeeze(montharrange.(fieldnames{i})(k).T(:,j,9:11,:)),2));
            tempamp = nanmean(squeeze(var_amplitudetemp(latindamp,9:11,:)),1);
            if strcmp(fieldnames{i},'highcl')
                rampvar_ind_highcl(k,:,j) = corr(temp',tempamp');
            elseif strcmp(fieldnames{i},'lowcl')    
                rampvar_ind_lowcl(k,:,j) = corr(temp',tempamp');
            elseif strcmp(fieldnames{i},'lowGHG')    
                rampvar_ind_lowGHG(k,:,j) = corr(temp',tempamp');
            end
            %rvar_eddy(i,:,j) = corr(squeeze(wa4060Scomposite.(fieldnames{i})(:,j,month,:))',var_minlong_eddy.(fieldnames{i}));            
        end
    end
end

%% quick contour plots
quick_contour_plots = 1;
if quick_contour_plots
    sig = .05;
    cbrew = cbrewer('div','RdBu',16);         

    prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];

    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
        5,[],[],2,1};


    logprestick = log(prestick);
    % subplotmaps(btoplot,1:12,log(regPres.highcl(1,:)),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
    %     contourlims,22,1:12,{'M','J','J','A','S','O','N','D','J','F','M','A'},...
    %     fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(maxheight) log(700)],1,'-',0);

    contourtitle = {'1995-2024 (high chlorine)','1955-1979 (low chlorine)','1995-2024 (high chlorine - low GHG)'};       

%     subplotmaps(r,longitudes,log(regpres),{'div','RdBu'},1,[],12,contourtitle,'Longitude','','Correlation','on',...
%         [-.75,.75],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),'',1,[0 360],[log(1) log(1000)],1,'none',0,'none');
% 
%     subplotmaps(rvar,longitudes,log(regpres),{'div','RdBu'},1,[],12,contourtitle,'Longitude','','Correlation','on',...
%         [-.75,.75],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),'',1,[0 360],[log(1) log(1000)],1,'none',0,'none');
%     
%     subplotmaps(rmax,longitudes,log(regpres),{'div','RdBu'},1,[],12,contourtitle,'Longitude','','Correlation','on',...
%         [-.75,.75],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),'',1,[0 360],[log(1) log(1000)],1,'none',0,'none');
% 
%     subplotmaps(rmaxvar,longitudes,log(regpres),{'div','RdBu'},1,[],12,contourtitle,'Longitude','','Correlation','on',...
%         [-.75,.75],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),'',1,[0 360],[log(1) log(1000)],1,'none',0,'none');
    
%     subplotmaps(ramp,longitudes,log(regpres),{'div','RdBu'},1,[],12,contourtitle,'Longitude','','Correlation','on',...
%         [-.75,.75],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),'',1,[0 360],[log(1) log(1000)],1,'none',0,'none');
    %rampvar2 = cat(1,rampvar(1,:,:),rampvar(3,:,:),rampvar(2,:,:));
    
    %pampvar (pampvar <= sig) = 0;
    %pampvar (pampvar > sig) = .1;
    
    subplotmaps(rampvar,longitudes,log(regpres),{'div','RdBu'},1,[],18,contourtitle,'Longitude','Pressure (hPa)','Correlation','on',...
        [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),['Correlation between ',num2str(abs(lats(latindamp))),'{\circ}S amplitude and 40-60{\circ}S temperature'],1,[0 360],[log(1) log(1000)],1,'-',0,'none');
    
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/WACCM_predrunsamp4060Scorr_',num2str(abs(lats(latindamp))),'S'];
    export_fig(filename,'-pdf');
    
    
 
    rminvar = cat(1,rvar(1,:,:),rvar(3,:,:),rvar(2,:,:));
%     
%     subplotmaps(rminvar,longitudes,log(regpres),{'div','RdBu'},1,[],18,contourtitle,'Longitude','Pressure','Correlation','on',...
%         [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),['Correlation between ',num2str(abs(lats(latindamp))),'{\circ}S min long and 40-60{\circ}S temperature'],1,[0 360],[log(1) log(1000)],1,'-',0,'none');
%     
%     filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/WACCM_predrunsminlon4060Scorr_',num2str(abs(lats(latindamp))),'S'];
%     export_fig(filename,'-pdf');
%     
%     
%     
%     rmaxvar = cat(1,rmaxvar(1,:,:),rmaxvar(3,:,:),rmaxvar(2,:,:));
%     
%     subplotmaps(r,longitudes,log(regpres),{'div','RdBu'},1,[],18,contourtitle,'Longitude','Pressure','Correlation','on',...
%         [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),['Correlation between ',num2str(abs(lats(latindamp))),'{\circ}S max long and 40-60{\circ}S temperature'],1,[0 360],[log(1) log(1000)],1,'-',0,'none');
%     
%     filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/WACCM_predrunsmaxlon4060Scorr_',num2str(abs(lats(latindamp))),'S'];
%     export_fig(filename,'-pdf');
    
    %%
    for i = 1:10
        contourtitle = {'No .1','No .2','No .3','No .4','No. 5','No. 6','No. 7','No. 8','No. 9','No. 10'};       
        subplotmaps(rampvar_ind_lowGHG(i,:,:),longitudes,log(regpres),{'div','RdBu'},1,[],18,contourtitle,'Longitude','Pressure','Correlation','on',...
        [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),['1995-2024 Correlation between ',num2str(abs(lats(latindamp))),'{\circ}S min long and 40-60{\circ}S temperature'],1,[0 360],[log(1) log(1000)],1,'-',0,'none');
    
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/',contourtitle{i},'_lowGHG_','WACCM_predrunsminlon4060Scorr_',num2str(abs(lats(latindamp))),'S'];
        export_fig(filename,'-pdf');
    
    end
    
    
end



%% plotting all data point
plot_ozone = 0;
if plot_ozone
    cbrew = cbrewer('qual','Set1',10);
    createfig('medium','on')

    subplot(2,1,1);    
    for j = 1:size(toz_extract_allmean,1)
        plot(longitudes,toz_extract_alleddymean(j,:)','LineWidth',3,'color',cbrew(j,:));
        hold on       
    end

    subplot(2,1,2);    
    for j = 1:size(toz_extract_allmean,1)
        plot(longitudes,std(toz_extract_all.(fieldnames{j}),0,2),'LineWidth',3,'color',cbrew(j,:));           
        hold on       
    end    


    %% plotting individuals
    createfig('medium','on')
    for i = 1:length(fieldnames)
    %     plot(minlong_ind.(fieldnames{i})','color',cbrew(i,:),'LineWidth',1);
    %     hold on
        ph(i) = plot(nanmean(minlong_ind.(fieldnames{i}),1),'color',cbrew(i,:),'LineWidth',3);
        hold on
        legend(fieldnames);
    end
end

%% plotting all data points temperature
plot_temperature = 1;
if plot_temperature
    cbrew = cbrewer('qual','Set1',10);
    createfig('medium','on')
    linesty = {'-','--','-.'};
    subplot(2,1,1);    
    for j = 1:size(var_extract_allmean,1)
        ph(j,:) = plot(longitudes,var_extract_alleddymean(j,:)','LineWidth',3,'color',cbrew(j,:),...
            'LineStyle',linesty{j});
        hold on          
    end
    lh = legend(ph,fieldnames,'box','off','fontsize',fsize+2);
    set(gca,'fontsize',fsize)
    xlim([-5 365]);
    xlabel('Longitude ({\circ}E)','fontsize',fsize+2)
    ylabel('K','fontsize',fsize+2)
    title(['Ensemble average temperatures at ',num2str(lat),'{\circ}S'],'fontsize',fsize+4) 
    
    subplot(2,1,2);    
    for j = 1:size(var_extract_allmean,1)
        plot(longitudes,std(var_extract_all.(fieldnames{j}),0,2),'LineWidth',3,'color',cbrew(j,:),...
            'LineStyle',linesty{j});           
        hold on       
    end    
    set(gca,'fontsize',fsize)
    xlim([-5 365]);
    xlabel('Longitude ({\circ}E)','fontsize',fsize+2)
    ylabel('K','fontsize',fsize+2)
    title(['Ensemble average temperature standard deviations at ',num2str(lat),'{\circ}S'],'fontsize',fsize+4) 

%     %% plotting individuals
%     createfig('medium','on')
%     for i = 1:length(fieldnames)
%     %     plot(var_minlong_ind.(fieldnames{i})','color',cbrew(i,:),'LineWidth',1);
%     %     hold on
%         ph(i) = plot(nanmean(var_minlong_ind.(fieldnames{i}),1),'color',cbrew(i,:),'LineWidth',3);
%         hold on
%         legend(fieldnames);
%     end
end

