% plot WACCM CCMI T line and contour plots
clear all

pastfutureens = 1;
REFC1only = 0;
pastfuturemembers = 0;
line_plots = 0;
trendyears = [1960,1998;1999,2099];
contourlims = [-10 10];
if trendyears(2,2) > 2014
    future = 1;
else 
    future = 0;
end

if trendyears(1,1) < 1979
    past = 1;
else 
    past = 0;
end

%% loading in data

load('/Volumes/My Book for Mac/work/data/CESM-CCMI/T/output/REFC1_6090S_ERApressure.mat');
load('/Volumes/My Book for Mac/work/data/CESM-CCMI/T/output/REFC2_6090S_ERApressure.mat');
load('/Volumes/My Book for Mac/work/data/CESM-CCMI/T/output/SENC2fGHG_6090S_ERApressure.mat');
load('/Volumes/My Book for Mac/work/data/CESM-CCMI/T/output/SENC2fODS_6090S_ERApressure.mat');
load('/Volumes/My Book for Mac/work/data/CESM-CCMI/T/output/ERA-Interim_6090S_ERApressure.mat');

ERAInt.wa = ERAInterim;

%% take linear trends over apropriate times
if ~future
    [b.REFC1,bint.REFC1] = regress_heighttime(REFC1,yearsREFC1,[trendyears(1,1),trendyears(1,2)],1);
end
if ~future && ~past
    [b.ERA,bint.ERA] = regress_heighttime(ERAInt,yearsERA,[trendyears(1,1),trendyears(1,2)],1);
end
[b.REFC2,bint.REFC2] = regress_heighttime(REFC2,yearsREFC2,[trendyears(1,1),trendyears(1,2)],1);
[b.GHG,bint.GHG] = regress_heighttime(GHG,yearsGHG,[trendyears(1,1),trendyears(1,2)],1);
[b.ODS,bint.ODS] = regress_heighttime(ODS,yearsODS,[trendyears(1,1),trendyears(1,2)],1);

if ~future
    [b2.REFC1,bint2.REFC1] = regress_heighttime(REFC1,yearsREFC1,[trendyears(2,1),trendyears(2,2)],1);
    
end
if ~future && ~past
    [b2.ERA,bint2.ERA] = regress_heighttime(ERAInt,yearsERA,[trendyears(2,1),trendyears(2,2)],1);
end

[b2.REFC2,bint2.REFC2] = regress_heighttime(REFC2,yearsREFC2,[trendyears(2,1),trendyears(2,2)],1);
[b2.GHG,bint2.GHG] = regress_heighttime(GHG,yearsGHG,[trendyears(2,1),trendyears(2,2)],1);
[b2.ODS,bint2.ODS] = regress_heighttime(ODS,yearsODS,[trendyears(2,1),trendyears(2,2)],1);

%% taking ensemble means
if ~future
    b.REFC1_ensmean = squeeze(nanmean(b.REFC1,1));
end
b.REFC2_ensmean = squeeze(nanmean(b.REFC2,1));
b.GHG_ensmean = squeeze(nanmean(b.GHG,1));
b.ODS_ensmean = squeeze(nanmean(b.ODS,1));

if ~future
    b2.REFC1_ensmean = squeeze(nanmean(b2.REFC1,1));
end
b2.REFC2_ensmean = squeeze(nanmean(b2.REFC2,1));
b2.GHG_ensmean = squeeze(nanmean(b2.GHG,1));
b2.ODS_ensmean = squeeze(nanmean(b2.ODS,1));

%% now I can plot

pressure = flip(double(ERApressure));
prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];
    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],5,[],[],2,1};
    logprestick = log(prestick);

monthnames = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};
    
%% past and future ensembles
if pastfutureens
   if ~future && past
        benstoplot = permute(cat(3,b.REFC1_ensmean(:,:,2),...
            b2.REFC1_ensmean(:,:,2),b.REFC2_ensmean(:,:,2),b2.REFC2_ensmean(:,:,2),...
            b.GHG_ensmean(:,:,2),b2.GHG_ensmean(:,:,2),b.ODS_ensmean(:,:,2),b2.ODS_ensmean(:,:,2)),[3,1,2]);
        titles = {'REF-C1 1979-1998','REF-C1 1999-2014',...
            'REF-C2 1979-1998','REF-C2 1999-2014','SEN-C2-fGHG 1978-1998','SEN-C2-fGHG 1999-2014',...
            'SEN-C2-fODS 1979-1998','SEN-C2-fODS 1999-2014'};    
    elseif future && past
        benstoplot = permute(cat(3,b.REFC2_ensmean(:,:,2),b2.REFC2_ensmean(:,:,2),...
            b.GHG_ensmean(:,:,2),b2.GHG_ensmean(:,:,2),b.ODS_ensmean(:,:,2),b2.ODS_ensmean(:,:,2)),[3,1,2]);
        titles = {'REF-C2 1979-1998','REF-C2 1999-2099','SEN-C2-fGHG 1978-1998','SEN-C2-fGHG 1999-2099',...
            'SEN-C2-fODS 1979-1998','SEN-C2-fODS 1999-2099'};    
   else
        contourlims = [-10 10];
        benstoplot = permute(cat(3,squeeze(b.ERA(:,:,:,2)),squeeze(b2.ERA(:,:,:,2)),b.REFC1_ensmean(:,:,2),...
            b2.REFC1_ensmean(:,:,2),b.REFC2_ensmean(:,:,2),b2.REFC2_ensmean(:,:,2),...
            b.GHG_ensmean(:,:,2),b2.GHG_ensmean(:,:,2),b.ODS_ensmean(:,:,2),b2.ODS_ensmean(:,:,2)),[3,1,2]);
        titles = {'ERA-Interim 1979-1998','ERA-Interim 1999-2014','REF-C1 1979-1998','REF-C1 1999-2014',...
            'REF-C2 1979-1998','REF-C2 1999-2014','SEN-C2-fGHG 1978-1998','SEN-C2-fGHG 1999-2014',...
            'SEN-C2-fODS 1979-1998','SEN-C2-fODS 1999-2014'};    
    end
            
    plotmtitle = 'Ensemble average Temperature trends';

    benstoplot = circshift(flip(benstoplot,3),[0,7,0])*10;

    benstoplot (benstoplot < contourlims(1)) = contourlims(1)+1;
    benstoplot (benstoplot > contourlims(2)) = contourlims(2)-1;
    
    subplotmaps(benstoplot,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','T/decade','on',...
        contourlims,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/Temperature_ensemble_',...
        num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
    export_fig(filename,'-pdf');
end

%% past and future member
if pastfuturemembers
   contourlims = [-10 10];
    benstoplot_members = permute(cat(3,squeeze(b.REFC2(1,:,:,2)),squeeze(b.REFC2(2,:,:,2)),squeeze(b.REFC2(3,:,:,2)),...
        squeeze(b.GHG(1,:,:,2)),squeeze(b.GHG(2,:,:,2)),squeeze(b.GHG(3,:,:,2)),...
        squeeze(b.ODS(1,:,:,2)),squeeze(b.ODS(2,:,:,2)),squeeze(b.ODS(3,:,:,2))),[3,1,2]);
    titles = {'REF-C2 no.1','REF-C2 no.2','REF-C2 no.3','SEN-C2-fGHG no.1','SEN-C2-fGHG no.2',...
        'SEN-C2-fGHG no.3','SEN-C2-fODS no.1','SEN-C2-fODS no.2','SEN-C2-fODS no.3'};       
            
    plotmtitle = ['Ensemble member Temperature trends ',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2))];

    benstoplot_members = circshift(flip(benstoplot_members,3),[0,7,0])*10;
    
    benstoplot_members (benstoplot_members < contourlims(1)) = contourlims(1)+1;
    benstoplot_members (benstoplot_members > contourlims(2)) = contourlims(2)-1;
    
    subplotmaps(benstoplot_members,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','K/decade','on',...
        contourlims,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/Temperature_members_',...
        num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'hPa.pdf'];
    export_fig(filename,'-pdf');
    
    % FUTURE
    contourlims = [-5 5];
    benstoplot_members = permute(cat(3,squeeze(b2.REFC2(1,:,:,2)),squeeze(b2.REFC2(2,:,:,2)),squeeze(b2.REFC2(3,:,:,2)),...
        squeeze(b2.GHG(1,:,:,2)),squeeze(b2.GHG(2,:,:,2)),squeeze(b2.GHG(3,:,:,2)),...
        squeeze(b2.ODS(1,:,:,2)),squeeze(b2.ODS(2,:,:,2)),squeeze(b2.ODS(3,:,:,2))),[3,1,2]);
    titles = {'REF-C2 no.1','REF-C2 no.2','REF-C2 no.3','SEN-C2-fGHG no.1','SEN-C2-fGHG no.2',...
        'SEN-C2-fGHG no.3','SEN-C2-fODS no.1','SEN-C2-fODS no.2','SEN-C2-fODS no.3'};       
            
    plotmtitle = ['Ensemble member Temperature trends ',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2))];

    benstoplot_members = circshift(flip(benstoplot_members,3),[0,7,0])*10;
    
    benstoplot_members (benstoplot_members < contourlims(1)) = contourlims(1)+1;
    benstoplot_members (benstoplot_members > contourlims(2)) = contourlims(2)-1;
    
    subplotmaps(benstoplot_members,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','T/decade','on',...
        contourlims,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/Temperature_members_',...
        num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
    export_fig(filename,'-pdf');
    
end


%% only REF-C1 ensembles (with ERA-Interim)

if REFC1only
    benstoplot_REFC1_past = permute(cat(3,squeeze(b.ERA(:,:,:,2)),squeeze(b.REFC1(1,:,:,2)),squeeze(b.REFC1(2,:,:,2)),...
        squeeze(b.REFC1(3,:,:,2)),squeeze(b.REFC1(4,:,:,2)),squeeze(b.REFC1(5,:,:,2))),[3,1,2]);
    benstoplot_REFC1_past = circshift(flip(benstoplot_REFC1_past,3),[0,7,0])*10;

    benstoplot_REFC1_future = permute(cat(3,squeeze(b2.ERA(:,:,:,2)),squeeze(b2.REFC1(1,:,:,2)),squeeze(b2.REFC1(2,:,:,2)),...
        squeeze(b2.REFC1(3,:,:,2)),squeeze(b2.REFC1(4,:,:,2)),squeeze(b2.REFC1(5,:,:,2))),[3,1,2]);
    benstoplot_REFC1_future = circshift(flip(benstoplot_REFC1_future,3),[0,7,0])*10;

    titles = {'ERA-Interim','REF-C1 - ens 1','REF-C1 - ens 2','REF-C1 - ens 3','REF-C1 - ens 4',...
        'REF-C1 - ens 5'};
    plotmtitle = '1979-1998 Temperature trends';

    benstoplot_REFC1_past (benstoplot_REFC1_past < contourlims(1)) = contourlims(1)+1;
    benstoplot_REFC1_past (benstoplot_REFC1_past > contourlims(2)) = contourlims(2)-1;
    
    subplotmaps(benstoplot_REFC1_past,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','T/decade','on',...
        contourlims,22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);

    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/Temperature_REFC1_1979-1998.pdf'],'-pdf');

    plotmtitle = '1999-2014 Temperature trends';

    benstoplot_REFC1_future (benstoplot_REFC1_future < contourlims(1)) = contourlims(1)+1;
    benstoplot_REFC1_future (benstoplot_REFC1_future > contourlims(2)) = contourlims(2)-1;
    
    subplotmaps(benstoplot_REFC1_future,1:12,log(pressure),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','T/decade','on',...
        [-10 10],22,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(1) log(1000)],1);

    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/trends/Temperaturecl_REFC1_1999-2014.pdf'],'-pdf');
end

%% line plots
if line_plots
for month = [8,10,12]
    for pres = [2,100]    
        sim = 'REF-C2';
        %REF-C2
        Tlineplots(month,pres,REFC2,b.REFC2,b2.REFC2,b.REFC2_ensmean,b2.REFC2_ensmean,yearsREFC2,ERApressure,trendyears,sim)

        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/T_',sim,'_',...
            num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
        export_fig(filename,'-pdf');

        sim = 'SEN-C2-fGHG';
        %REF-C2
        Tlineplots(month,pres,GHG,b.GHG,b2.GHG,b.GHG_ensmean,b2.GHG_ensmean,yearsGHG,ERApressure,trendyears,sim)

        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/T_',sim,'_',...
            num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
        export_fig(filename,'-pdf');

        sim = 'SEN-C2-fODS';
        %REF-C2
        Tlineplots(month,pres,ODS,b.ODS,b2.ODS,b.ODS_ensmean,b2.ODS_ensmean,yearsODS,ERApressure,trendyears,sim)

        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/T_',sim,'_',...
            num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
        export_fig(filename,'-pdf');

        if ~future
            % line plots REF-C1
            sim = 'REF-C1';
            %REF-C2
            Tlineplots(month,pres,REFC1,b.REFC1,b2.REFC1,b.REFC1_ensmean,b2.REFC1_ensmean,yearsREFC1,ERApressure,trendyears,sim)

            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/lineplots/T_',sim,'_',...
                num2str(month),'_',num2str(pres),'_',num2str(trendyears(1,1)),'-',num2str(trendyears(1,2)),'and',num2str(trendyears(2,1)),'-',num2str(trendyears(2,2)),'hPa.pdf'];
            export_fig(filename,'-pdf');
        end
    end
end
end
close all
