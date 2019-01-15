% plot WACCM CCMI Z3 line and contour plots
clear all

variable = 'Z3';
trendyears = [1970,1999;2000,2025];
pastdates = [num2str(trendyears(1,1)),'-',num2str(trendyears(1,2))];
futuredates = [num2str(trendyears(2,1)),'-',num2str(trendyears(2,2))];
area = '7590S';

%contourlimspast = [-9 9];
%contourlimsfuture = [-4.5 4.5];

if strcmp(variable,'T')
    ctitle = 'K/decade';
    linetitle = 'Temperature (K)';
    contourlimspast = [-9 9];
    contourlimsfuture = [-4.5 4.5];
elseif strcmp(variable,'Z3')
    ctitle = 'm/decade';
    linetitle = 'Height (m)';
    contourlimspast = [-600 600];
    contourlimsfuture = [-600 600];
end

% different plot groups
MERRAvEnsemble = 0;
MERRAvERA = 0;
past_members = 0;
future_members = 1;
line_plots = 1;

%% loading in data

load(['/Volumes/ExternalOne/work/data/CESM-CCMI/',variable,'/output/',variable,'_',area,'.mat']);
load(['/Volumes/ExternalOne/work/data/CESM-CCMI/',variable,'/output/',variable,'_ERA-Interim_',area,'.mat']);

SDsave = regp(end);
SDsaveyears = yearsfinal(end);
x = zeros(1,length(regp));
for i = 1:length(regp)
    if ~isempty(find(yearsfinal(i).y == trendyears(1,1), 1)) && ~isempty(find(yearsfinal(i).y == trendyears(2,2), 1))
        x(i) = 1;    
    end
end

regp(~x) = [];
yearsfinal(~x) = [];
filenames(~x) = [];

if trendyears(1,1) < 1979 || trendyears(2,2) > 2016
    includeERA = 0;
else
    includeERA = 1;
    ERAInt.wa = ERAInterim;
    yERA.y = yearsERA;
end
filenamestemp = filenames;
if length(filenames) > 0
    regp(length(regp)+1).wa = nanmean(cat(3,regp(1).wa,regp(2).wa,regp(3).wa),3);    
    regp(length(regp)+1).wa = nanmean(cat(3,regp(4).wa,regp(5).wa,regp(6).wa),3);    
    regp(length(regp)+1).wa = nanmean(cat(3,regp(7).wa,regp(8).wa,regp(9).wa),3);
    yearsfinal(length(yearsfinal)+1).y = yearsfinal(1).y;
    yearsfinal(length(yearsfinal)+1).y = yearsfinal(4).y;
    yearsfinal(length(yearsfinal)+1).y = yearsfinal(7).y;
    filenames{length(filenames)+1} = 'r2ens';
    filenames{length(filenames)+1} = 's2odsens';
    filenames{length(filenames)+1} = 's2ghgens';
end
if length(filenamestemp) > 9
    regp(length(regp)+1).wa = nanmean(cat(3,regp(10).wa,regp(11).wa,regp(12).wa,...
        regp(13).wa,regp(14).wa),3);
    yearsfinal(length(yearsfinal)+1).y = yearsfinal(10).y;
    filenames{length(filenames)+1} = 'r1ens';
end
%% take linear trends over apropriate times
[b,bint] = regress_heighttime(regp,yearsfinal,[trendyears(1,1),trendyears(1,2)],1);
[b2,bint2] = regress_heighttime(regp,yearsfinal,[trendyears(2,1),trendyears(2,2)],1);

if includeERA
    [bERA,bintERA] = regress_heighttime(ERAInt,yERA,[trendyears(1,1),trendyears(1,2)],1);
    [b2ERA,bint2ERA] = regress_heighttime(ERAInt,yERA,[trendyears(2,1),trendyears(2,2)],1);
end

%% taking ensemble averages

REFC2_exist = regexp(filenames,'REFC2');
REFC1_exist = regexp(filenames,'REFC1');
REFC1SD_exist = regexp(filenames,'REFC1SD');
GHG_exist = regexp(filenames,'ODS');
ODS_exist = regexp(filenames,'GHG');

refc2ens_exist = regexp(filenames,'r2ens');
refc1ens_exist = regexp(filenames,'r1ens');
ghgens_exist = regexp(filenames,'s2ghgens');
odsens_exist = regexp(filenames,'s2odsens');

REFC2_exist_ind = zeros(1,length(filenames));
REFC1_exist_ind = zeros(1,length(filenames));
REFC1SD_exist_ind = zeros(1,length(filenames));
GHG_exist_ind = zeros(1,length(filenames));
ODS_exist_ind = zeros(1,length(filenames));

refc2ens_exist_ind = zeros(1,length(filenames));
refc1ens_exist_ind = zeros(1,length(filenames));
ghgens_exist_ind = zeros(1,length(filenames));
odsens_exist_ind = zeros(1,length(filenames));

for i = 1:length(filenames)
    REFC2temp = REFC2_exist{i};
    REFC1temp = REFC1_exist{i};
    REFC1SDtemp = REFC1SD_exist{i};
    GHGtemp = GHG_exist{i};
    ODStemp = ODS_exist{i};
    
    r2enstemp = refc2ens_exist{i};
    r1enstemp = refc1ens_exist{i};   
    ghgenstemp = ghgens_exist{i};
    odsenstemp = odsens_exist{i};
    
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
    
    if ~isempty(r2enstemp)
        refc2ens_exist_ind(i) = r2enstemp/r2enstemp;        
    end
    if ~isempty(r1enstemp)
        refc1ens_exist_ind(i) = r1enstemp/r1enstemp;
    end
    if ~isempty(ghgenstemp)
        ghgens_exist_ind(i) = ghgenstemp/ghgenstemp;
    end
    if ~isempty(odsenstemp)
        odsens_exist_ind(i) = odsenstemp/odsenstemp;
    end    
end

REFC2_exist_ind = find(REFC2_exist_ind);
REFC1_exist_ind = find(REFC1_exist_ind);
REFC1SD_exist_ind = find(REFC1SD_exist_ind);
GHG_exist_ind = find(GHG_exist_ind);
ODS_exist_ind = find(ODS_exist_ind);
refc2ens_exist_ind = find(refc2ens_exist_ind);
refc1ens_exist_ind = find(refc1ens_exist_ind);
ghgens_exist_ind = find(ghgens_exist_ind);
odsens_exist_ind = find(odsens_exist_ind);

%% now I can plot

pressure = flip(double(ERApressure));
prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];
    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],5,[],[],2,1};
    logprestick = log(prestick);

monthnames = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};
    
%% plot members of future simulations (past)
if future_members
    btoplot = cat(1,b((REFC2_exist_ind),:,:,2),b((GHG_exist_ind),:,:,2),...
        b((ODS_exist_ind),:,:,2));   
    titles = {'All forcings no.1' ,'All forcings no.2' ,'All forcings no.3' ,...
        'GHGs only no.1' ,'GHGs only no.2' ,'GHGs only no.3' ,...
            'ODSs only no.1' ,'ODSs only no.2' ,'ODSs only no.3' };    

    plotmtitle = [area,' ',variable,' trends (',pastdates,')'];

    btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

    btoplot (btoplot < contourlimspast(1)) = contourlimspast(1)+1;
    btoplot (btoplot > contourlimspast(2)) = contourlimspast(2)-1;

    subplotmaps(btoplot,1:12,log(pressure),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
        contourlimspast,14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(700)],1,'-',0);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/',variable,'_individualfuture_',...
        pastdates,'_',area,'.png'];    
    export_fig(filename,'-png');   

    %% future members of future simulations (future)
    btoplot = cat(1,b2((REFC2_exist_ind),:,:,2),b2((GHG_exist_ind),:,:,2),...
        b2((ODS_exist_ind),:,:,2));   

    plotmtitle = [area,' ',variable,' trends (',futuredates,')'];

    btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

    btoplot (btoplot < contourlimsfuture(1)) = contourlimsfuture(1)+1;
    btoplot (btoplot > contourlimsfuture(2)) = contourlimsfuture(2)-1;

    subplotmaps(btoplot,1:12,log(pressure),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
        contourlimsfuture,14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(700)],1,'-',0);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/',variable,'_individualfuture_',...
        futuredates,'_',area,'.png'];
    export_fig(filename,'-png');
end

%% plot REF-C1 with specified dynamics (past)
if past_members
    btoplot = cat(1,b((REFC1_exist_ind),:,:,2));
    titles = {'REF-C1 no.1' ,'REF-C1 no.2' ,'REF-C1 no.3' ,...
        'REF-C1 no.4' ,'REF-C1 no.5' ,'REF-C1SD (MERRA)'};        

    plotmtitle = [area,' ',variable,' trends (',pastdates,')'];

    btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

    btoplot (btoplot < contourlimspast(1)) = contourlimspast(1)+1;
    btoplot (btoplot > contourlimspast(2)) = contourlimspast(2)-1;

    subplotmaps(btoplot,1:12,log(pressure),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
        contourlimspast,14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(700)],1,'-',0);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/REF-C1',variable,'_individualfuture_',...
        pastdates,'_',area,'.png'];
    export_fig(filename,'-png');

    %% plot REF-C1 with specified dynamics (future)
    btoplot = cat(1,b2((REFC1_exist_ind),:,:,2));       

    plotmtitle = [area,' ',variable,' trends (',futuredates,')'];

    btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

    btoplot (btoplot < contourlimsfuture(1)) = contourlimsfuture(1)+1;
    btoplot (btoplot > contourlimsfuture(2)) = contourlimsfuture(2)-1;

    subplotmaps(btoplot,1:12,log(pressure),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
        contourlimsfuture,14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(700)],1,'-',0);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/REF-C1',variable,'_individualfuture_',...
        futuredates,'_',area,'.png'];
    export_fig(filename,'-png');
end
%% plot REF-C1SD with specified dynamics and ERA-Interim(past)
if MERRAvERA
    btoplot = cat(1,b((REFC1SD_exist_ind),:,:,2),bERA(:,:,:,2));
    titles = {'REF-C1SD (MERRA)','ERA-Interim'};        

    plotmtitle = [area,' ',variable,' trends (',pastdates,')'];

    btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

    btoplot (btoplot < contourlimspast(1)) = contourlimspast(1)+1;
    btoplot (btoplot > contourlimspast(2)) = contourlimspast(2)-1;

    subplotmaps(btoplot,1:12,log(pressure),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
        contourlimspast,14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(700)],1,'-',0);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/Re-analysis_',variable,'_individualfuture_',...
        pastdates,'_',area,'.png'];
    export_fig(filename,'-png');

    %% plot REF-C1SD with specified dynamics and ERA-Interim(future)
    btoplot = cat(1,b2((REFC1SD_exist_ind),:,:,2),b2ERA(:,:,:,2));      

    plotmtitle = [area,' ',variable,' trends (',futuredates,')'];

    btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

    btoplot (btoplot < contourlimsfuture(1)) = contourlimsfuture(1)+1;
    btoplot (btoplot > contourlimsfuture(2)) = contourlimsfuture(2)-1;

    subplotmaps(btoplot,1:12,log(pressure),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
        contourlimsfuture,14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(700)],1,'-',0);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/Re-analysis_',variable,'_individualfuture_',...
        futuredates,'_',area,'.png'];
    export_fig(filename,'-png');
end

%% MERRA vs free running (past)
if MERRAvEnsemble
    btoplot = cat(1,b((REFC1SD_exist_ind),:,:,2),b([refc2ens_exist_ind,refc1ens_exist_ind],:,:,2));
    titles = {'MERRA','Free running (REF-C2)','Free running (REF-C1 - prescribed SSTs)'};        

    plotmtitle = [area,' ',variable,' trends (',pastdates,')'];

    btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

    btoplot (btoplot < contourlimspast(1)) = contourlimspast(1)+1;
    btoplot (btoplot > contourlimspast(2)) = contourlimspast(2)-1;

    subplotmaps(btoplot,1:12,log(pressure),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
        contourlimspast,14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(700)],1,'-',0);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/MERRAvFR_',variable,'_individualfuture_',...
        pastdates,'_',area,'.png'];
    export_fig(filename,'-png');

    %% MERRA vs free running (future)
    btoplot = cat(1,b2((REFC1SD_exist_ind),:,:,2),b2([refc2ens_exist_ind,refc1ens_exist_ind],:,:,2));
    titles = {'MERRA','Free running (REF-C2)','Free running (REF-C1 - prescribed SSTs)'};        

    plotmtitle = [area,' ',variable,' trends (',futuredates,')'];

    btoplot = circshift(flip(btoplot,3),[0,7,0])*10;

    btoplot (btoplot < contourlimsfuture(1)) = contourlimsfuture(1)+1;
    btoplot (btoplot > contourlimsfuture(2)) = contourlimsfuture(2)-1;

    subplotmaps(btoplot,1:12,log(pressure),{'div','RdBu'},1,[],16,titles,'Month','Pressure (hPa)',ctitle,'on',...
        contourlimsfuture,14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(30) log(700)],1,'-',0);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/MERRAvFR_',variable,'_individualfuture_',...
        futuredates,'_',area,'.png'];
    export_fig(filename,'-png');
end

%%

if line_plots
    for month = [8,10,12]
        for pres = [2,100]    
            sim = 'REF-C2';
            %REF-C2
            Z3lineplots(month,pres,regp([REFC2_exist_ind,refc2ens_exist_ind]),...
                b([REFC2_exist_ind,refc2ens_exist_ind],:,:,:),...
                b2([REFC2_exist_ind,refc2ens_exist_ind],:,:,:),...
                yearsfinal([REFC2_exist_ind,refc2ens_exist_ind]),...
                ERApressure,trendyears,SDsave,SDsaveyears,sim,linetitle,area,[1949 2101]);

            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/',variable,'_',sim,'_',...
                num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
            export_fig(filename,'-pdf');

            if ~isempty(regp([REFC1_exist_ind,refc1ens_exist_ind]))
                sim = 'REF-C1';
                %REF-C2
                Z3lineplots(month,pres,regp([REFC1_exist_ind,refc1ens_exist_ind]),...
                    b([REFC1_exist_ind,refc1ens_exist_ind],:,:,:),...
                    b2([REFC1_exist_ind,refc1ens_exist_ind],:,:,:),...
                    yearsfinal([REFC1_exist_ind,refc1ens_exist_ind]),...
                    ERApressure,trendyears,SDsave,SDsaveyears,sim,linetitle,area,[1954 2016]);

                filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/',variable,'_',sim,'_',...
                    num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
                export_fig(filename,'-pdf');
            end
            
            sim = 'GHGs-only';
                %REF-C2
            Z3lineplots(month,pres,regp([GHG_exist_ind,ghgens_exist_ind]),...
                b([GHG_exist_ind,ghgens_exist_ind],:,:,:),...
                b2([GHG_exist_ind,ghgens_exist_ind],:,:,:),...
                yearsfinal([GHG_exist_ind,ghgens_exist_ind]),...
                ERApressure,trendyears,SDsave,SDsaveyears,sim,linetitle,area,[1959 2101]);

            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/',variable,'_',sim,'_',...
                num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
            export_fig(filename,'-pdf');
            
            sim = 'ODSs-only';
                %REF-C2
            Z3lineplots(month,pres,regp([ODS_exist_ind,odsens_exist_ind]),...
                b([ODS_exist_ind,odsens_exist_ind],:,:,:),...
                b2([ODS_exist_ind,odsens_exist_ind],:,:,:),...
                yearsfinal([ODS_exist_ind,odsens_exist_ind]),...
                ERApressure,trendyears,SDsave,SDsaveyears,sim,linetitle,area,[1959 2101]);

            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/',variable,'_',sim,'_',...
                num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
            export_fig(filename,'-pdf');
        end
    end
end
close all
