% plot Surface temperature trends
clear all

pastfutureens = 0;
variable = 'TS';
pastfuturemembers = 1;
line_plots = 0;
trendyears = [1979,1999;2000,2014];
pastdates = [num2str(trendyears(1,1)),'-',num2str(trendyears(1,2))];
futuredates = [num2str(trendyears(2,1)),'-',num2str(trendyears(2,2))];
contourlimspast = [-2 2];
contourlimsfuture = [-2 2];

if strcmp(variable,'T') || strcmp(variable,'TS')
    ctitle = 'K/decade';
elseif strcmp(variable,'Z3')
    ctitle = 'm/decade';
end

directory = ['/Volumes/ExternalOne/work/data/CESM-CCMI/',variable,'/'];
datedirectory = ['/Volumes/ExternalOne/work/data/CESM-CCMI/date/'];
datefiles = dir([datedirectory,'*.nc*']);
files = dir([directory,'/','*nc*']);
filenames = {files.name}';
%% Reading in data
for i = 1:length(files);    
    [~, data(i),~] = Read_in_netcdf([directory,files(i).name]);
    [~, date(i),~] = Read_in_netcdf([datedirectory,datefiles(i).name]);    
    data(i).(variable) = squeeze(data(i).(variable));

    years = num2str(date(i).date);
    for j = 1:length(years)
        yearsfinal(i).y(j,:) = str2double(years(j,1:4));
    end
    yearsfinal(i).y = circshift(yearsfinal(i).y,1);
    yearsfinal(i).y(1) = yearsfinal(i).y(2); 
    TSonly(i).TS = squeeze(nanmean(data(i).(variable),1));
end

TSonly(length(TSonly)+1).TS = nanmean(cat(3,TSonly(1).TS,TSonly(2).TS,TSonly(3).TS),3);
yearsfinal(length(yearsfinal)+1).y = yearsfinal(1).y;

TSonly(length(TSonly)+1).TS = nanmean(cat(3,TSonly(10).TS,TSonly(11).TS,TSonly(12).TS,TSonly(13).TS,...
    TSonly(14).TS),3);
yearsfinal(length(yearsfinal)+1).y = yearsfinal(1).y;
yearsfinal(length(yearsfinal)+1).y = yearsfinal(10).y;

%% Read in ERA
includeERA = 0;
if trendyears(1,1) > 1979 || trendyears(2,2) < 2016
    includeERA = 1;
    scale = 1;
    eravar = 't2m';
    ERAdata = ncread(['/Volumes/ExternalOne/work/data/ERA-Interim/',variable,'/',variable,'_ERA-Interim.nc'],eravar)./scale;
    ERAlatitude = ncread('/Volumes/ExternalOne/work/data/ERA-Interim/Z3/Z3_ERA-Interim.nc','latitude');    
    yearsERA = repmat([1979:2016],[12,1]);
    yERA.y = yearsERA(:); 
    ERAInt.wa = squeeze(nanmean(ERAdata));
    ERAInt_interp.wa = interp1(ERAlatitude,ERAInt.wa,data(1).lat);
end


%%
x = zeros(1,length(TSonly));
for i = 1:length(TSonly)
    if ~isempty(find(yearsfinal(i).y == trendyears(1,1), 1)) && ~isempty(find(yearsfinal(i).y == trendyears(2,2), 1))
        x(i) = 1;    
    end
end

TSonly(~x) = [];
yearsfinal(~x) = [];
files(~x) = [];

%% take linear trends over apropriate times
[b,bint] = regress_heighttime(TSonly,yearsfinal,[trendyears(1,1),trendyears(1,2)],1);
[b2,bint2] = regress_heighttime(TSonly,yearsfinal,[trendyears(2,1),trendyears(2,2)],1);

if includeERA
    [bERA,bintERA] = regress_heighttime(ERAInt_interp,yERA,[trendyears(1,1),trendyears(1,2)],1);
    [b2ERA,bint2ERA] = regress_heighttime(ERAInt_interp,yERA,[trendyears(2,1),trendyears(2,2)],1);
end

%% taking ensemble averages

REFC2_exist = regexp(filenames,'REFC2');
REFC1_exist = regexp(filenames,'REFC1');
REFC1SD_exist = regexp(filenames,'REFC1SD');
GHG_exist = regexp(filenames,'ODS');
ODS_exist = regexp(filenames,'GHG');
REFC2_exist_ind = zeros(1,length(filenames));
REFC1_exist_ind = zeros(1,length(filenames));
REFC1SD_exist_ind = zeros(1,length(filenames));
GHG_exist_ind = zeros(1,length(filenames));
ODS_exist_ind = zeros(1,length(filenames));
for i = 1:length(filenames)
    REFC2temp = REFC2_exist{i};
    REFC1temp = REFC1_exist{i};
    REFC1SDtemp = REFC1SD_exist{i};
    GHGtemp = GHG_exist{i};
    ODStemp = ODS_exist{i};
    if ~isempty(REFC2temp)
        REFC2_exist_ind(i) = REFC2temp/REFC2temp;        
    end
    if ~isempty(REFC1temp)
        REFC1_exist_ind(i) = REFC1temp/REFC1temp;
    end
    if ~isempty(REFC1SDtemp)
        REFC1SD_exist_ind(i) = REFC1SDtemp/REFC1SDtemp;
    end
    if ~isempty(GHGtemp)
        GHG_exist_ind(i) = GHGtemp/GHGtemp;
    end
    if ~isempty(ODStemp)
        ODS_exist_ind(i) = ODStemp/ODStemp;
    end
end

%% now I can plot

monthnames = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};
    
%% plot members of future simulations (past)
btoplot = cat(1,b(find(REFC2_exist_ind),:,:,2),b(find(GHG_exist_ind),:,:,2),...
    b(find(ODS_exist_ind),:,:,2));   
titles = {'All forcings no.1' ,'All forcings no.2' ,'All forcings no.3' ,...
    'GHGs only no.1' ,'GHGs only no.2' ,'GHGs only no.3' ,...
        'ODSs only no.1' ,'ODSs only no.2' ,'ODSs only no.3' };    

plotmtitle = [variable,' trends (',pastdates,')'];

btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

btoplot (btoplot < contourlimspast(1)) = contourlimspast(1)+1;
btoplot (btoplot > contourlimspast(2)) = contourlimspast(2)-1;

subplotmaps(btoplot,1:12,data(1).lat,{'div','RdBu'},1,[],16,titles,'Month','Latitude ({\circ}N)',ctitle,'on',...
    contourlimspast,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    -90:20:90,-90:20:90,plotmtitle ,1,[1 12],[-90 90],0);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/',variable,'_individualfuture_',...
    pastdates,'.pdf'];
export_fig(filename,'-pdf');        

%% future members of future simulations (future)
btoplot = cat(1,b2(find(REFC2_exist_ind),:,:,2),b2(find(GHG_exist_ind),:,:,2),...
    b2(find(ODS_exist_ind),:,:,2));   
titles = {'All forcings no.1' ,'All forcings no.2' ,'All forcings no.3' ,...
    'GHGs only no.1' ,'GHGs only no.2' ,'GHGs only no.3' ,...
        'ODSs only no.1' ,'ODSs only no.2' ,'ODSs only no.3' };    

plotmtitle = [variable,' trends (',futuredates,')'];

btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

btoplot (btoplot < contourlimsfuture(1)) = contourlimsfuture(1)+1;
btoplot (btoplot > contourlimsfuture(2)) = contourlimsfuture(2)-1;

subplotmaps(btoplot,1:12,data(1).lat,{'div','RdBu'},1,[],16,titles,'Month','Latitude ({\circ}S)',ctitle,'on',...
    contourlimsfuture,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    -90:20:90,-90:20:90,plotmtitle ,1,[1 12],[-90 90],0);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/',variable,'_individualfuture_',...
    futuredates,'.pdf'];
export_fig(filename,'-pdf');

%% plot REF-C1 with specified dynamics (past)
btoplot = cat(1,b(find(REFC1_exist_ind),:,:,2));
titles = {'REF-C1 no.1' ,'REF-C1 no.2' ,'REF-C1 no.3' ,...
    'REF-C1 no.4' ,'REF-C1 no.5' ,'REF-C1SD'};        

plotmtitle = [variable,' trends (',pastdates,')'];

btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

btoplot (btoplot < contourlimspast(1)) = contourlimspast(1)+1;
btoplot (btoplot > contourlimspast(2)) = contourlimspast(2)-1;

subplotmaps(btoplot,1:12,data(1).lat,{'div','RdBu'},1,[],16,titles,'Month','Latitude ({\circ}S)',ctitle,'on',...
    contourlimspast,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    -90:20:90,-90:20:90,plotmtitle ,1,[1 12],[-90 90],0);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/REF-C1',variable,'_individualfuture_',...
    pastdates,'_','.pdf'];
export_fig(filename,'-pdf');

%% plot REF-C1 with specified dynamics (future)
btoplot = cat(1,b2(find(REFC1_exist_ind),:,:,2));   
titles = {'REF-C1 no.1' ,'REF-C1 no.2' ,'REF-C1 no.3' ,...
    'REF-C1 no.4' ,'REF-C1 no.5' ,'REF-C1SD'};        

plotmtitle = [variable,' trends (',futuredates,')'];

btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

btoplot (btoplot < contourlimsfuture(1)) = contourlimsfuture(1)+1;
btoplot (btoplot > contourlimsfuture(2)) = contourlimsfuture(2)-1;

subplotmaps(btoplot,1:12,data(1).lat,{'div','RdBu'},1,[],16,titles,'Month','Latitude ({\circ}S)',ctitle,'on',...
    contourlimsfuture,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    -90:20:90,-90:20:90,plotmtitle ,1,[1 12],[-90 90],0);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/REF-C1',variable,'_individualfuture_',...
    futuredates,'_','.pdf'];
export_fig(filename,'-pdf');

%% plot REF-C1SD with specified dynamics and ERA-Interim(past)
btoplot = cat(1,b(find(REFC1SD_exist_ind),:,:,2),bERA(:,:,:,2));
titles = {'REF-C1SD','ERA-Interim'};        

plotmtitle = [variable,' trends (',pastdates,')'];

btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

btoplot (btoplot < contourlimspast(1)) = contourlimspast(1)+1;
btoplot (btoplot > contourlimspast(2)) = contourlimspast(2)-1;

subplotmaps(btoplot,1:12,data(1).lat,{'div','RdBu'},1,[],16,titles,'Month','Latitude ({\circ}S)',ctitle,'on',...
    contourlimspast,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    -90:20:90,-90:20:90,plotmtitle ,1,[1 12],[-90 90],0);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/Re-analysis_',variable,'_individualfuture_',...
    pastdates,'.pdf'];
export_fig(filename,'-pdf');

%% plot REF-C1SD with specified dynamics and ERA-Interim(future)
btoplot = cat(1,b2(find(REFC1SD_exist_ind),:,:,2),b2ERA(:,:,:,2));   
titles = {'REF-C1SD','ERA-Interim'};        

plotmtitle = [variable,' trends (',futuredates,')'];

btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

btoplot (btoplot < contourlimsfuture(1)) = contourlimsfuture(1)+1;
btoplot (btoplot > contourlimsfuture(2)) = contourlimsfuture(2)-1;

subplotmaps(btoplot,1:12,data(1).lat,{'div','RdBu'},1,[],16,titles,'Month','Latitude ({\circ}S)',ctitle,'on',...
    contourlimsfuture,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    -90:20:90,-90:20:90,plotmtitle ,1,[1 12],[-90 90],0);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/Re-analysis_',variable,'_individualfuture_',...
    futuredates,'.pdf'];
export_fig(filename,'-pdf');

%% MERRA vs free running
btoplot = cat(1,b(find(REFC1SD_exist_ind),:,:,2),b(end-1:end,:,:,2));
titles = {'MERRA','Free running','Free running (fixed SSTs)'};        

plotmtitle = [variable,' trends (',pastdates,')'];

btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

btoplot (btoplot < contourlimspast(1)) = contourlimspast(1)+1;
btoplot (btoplot > contourlimspast(2)) = contourlimspast(2)-1;

subplotmaps(btoplot,1:12,data(1).lat,{'div','RdBu'},1,[],16,titles,'Month','Latitude ({\circ}S)',ctitle,'on',...
    contourlimsfuture,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    -90:5:90,-90:5:90,plotmtitle ,1,[1 12],[-90 -50],0);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/MERRAvFR_',variable,'_individualfuture_',...
    pastdates,'.pdf'];
export_fig(filename,'-pdf');

%% plot REF-C1SD with specified dynamics and ERA-Interim(future)
btoplot = cat(1,b2(find(REFC1SD_exist_ind),:,:,2),b2(end-1:end,:,:,2));
titles = {'MERRA','Free running','Free running (prescribed SSTs)'};        

plotmtitle = [variable,' trends (',futuredates,')'];

btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

btoplot (btoplot < contourlimsfuture(1)) = contourlimsfuture(1)+1;
btoplot (btoplot > contourlimsfuture(2)) = contourlimsfuture(2)-1;

subplotmaps(btoplot,1:12,data(1).lat,{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
    contourlimsfuture,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
    -90:5:90,-90:5:90,plotmtitle ,1,[1 12],[-90 -50],0);
filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/MERRAvFR_',variable,'_individualfuture_',...
    futuredates,'.pdf'];
export_fig(filename,'-pdf');

%%
if line_plots
for month = [8,10,12]
    for pres = [2,100]    
        sim = 'REF-C2';
        %REF-C2
        Z3lineplots(month,pres,REFC2,b.REFC2,b2.REFC2,b.REFC2_ensmean,b2.REFC2_ensmean,yearsREFC2,ERApressure,trendyears,sim)

        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/Z3_',sim,'_',...
            num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
        export_fig(filename,'-pdf');

        sim = 'SEN-C2-fGHG';
        %REF-C2
        Z3lineplots(month,pres,GHG,b.GHG,b2.GHG,b.GHG_ensmean,b2.GHG_ensmean,yearsGHG,ERApressure,trendyears,sim)

        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/Z3_',sim,'_',...
            num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
        export_fig(filename,'-pdf');

        sim = 'SEN-C2-fODS';
        %REF-C2
        Z3lineplots(month,pres,ODS,b.ODS,b2.ODS,b.ODS_ensmean,b2.ODS_ensmean,yearsODS,ERApressure,trendyears,sim)

        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/Z3_',sim,'_',...
            num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
        export_fig(filename,'-pdf');

        if ~future
            % line plots REF-C1
            sim = 'REF-C1';
            %REF-C2
            Z3lineplots(month,pres,REFC1,b.REFC1,b2.REFC1,b.REFC1_ensmean,b2.REFC1_ensmean,yearsREFC1,ERApressure,trendyears,sim)

            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/Z3_',sim,'_',...
                num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
            export_fig(filename,'-pdf');
        end
    end
end
end
close all
