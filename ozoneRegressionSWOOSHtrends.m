% SWOOSH trends
clear all
%% import SWOOSH
Stimeperiod = [1999 2014];

%[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedanomfillo3q_swoosh-v02.6-198401-201712-latpress-2.5deg-L31.nc');
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedeqfillo3q_swoosh-v02.6-198401-201712-latpress-2.5deg-L31.nc');
sfields = fieldnames(SWOOSH);
SWOOSH.(sfields{1}) = cat(3,SWOOSH.(sfields{1}),zeros(size(SWOOSH.(sfields{1}),1),size(SWOOSH.(sfields{1}),2),1));
SWOOSH.(sfields{1}) (SWOOSH.(sfields{1}) == 0) = NaN;
SWOOSHyears = repmat(1984:2017,[12,1]);
%SWOOSHyears = [SWOOSHyears(:);ones(12,1)*2016];
SWOOSHextract = permute(SWOOSH.(sfields{1})(:,:,SWOOSHyears >= Stimeperiod(1) & SWOOSHyears <= Stimeperiod(2)),[2,1,3]);

% removing missing lats
SWOOSHextract = SWOOSHextract(:,4:end-3,:);
SWOOSH.lat = SWOOSH.lat(4:end-3);

%read in SWOOSH regression functions
[SWOOSHregfun] = ozoneRegression_SWOOSHregfun(Stimeperiod);

%% Regress SWOOSH
linear = 1;
[bSWOOSH,predictorsSWOOSH,O3AnomalySWOOSH,uncertainty] = ...
    ozoneRegressionTrends(SWOOSHextract,SWOOSHregfun.singaporedata,...
    SWOOSHregfun.solardata,0,SWOOSHregfun.MEIdata,0,0,0,SWOOSHregfun.aerosols,Stimeperiod,0,linear,0);%SDWaccmData(2).SPEintpres,, SDWaccmData(1).NO2*1e9

%%
if linear
    ball = bSWOOSH.percent(:,:,1:2);
else    
    [~,ball] = ozoneRegressionEnsAve(O3AnomalySWOOSH.percent_residuals,Stimeperiod(2) - Stimeperiod(1));
end

[~,bmonthall,pval,umonthraw] = ozoneRegressionEnsAve(O3AnomalySWOOSH.percent_residuals_months,Stimeperiod(2) - Stimeperiod(1));
%[~,ball,pval] = ozoneRegressionEnsAve(O3AnomalySWOOSH.percent_residuals,Stimeperiod(2) - Stimeperiod(1));
%[~,bmonthallARgone] = ozoneRegressionEnsAve(O3AnomalySWOOSH.percent_residuals_months_ARgone,Stimeperiod(2) - Stimeperiod(1));
[~,braw,pvalraw,uraw] = ozoneRegressionEnsAve(O3AnomalySWOOSH.percent,Stimeperiod(2) - Stimeperiod(1));

%%

clearvars watemp bmonth
lats = [-80 -50];

latind = SWOOSH.lat >= lats(1) & SWOOSH.lat <= lats(2);

for i = 1:12
    for j = 1:size(O3AnomalySWOOSH.percent_residuals_months,2)
        watemp(:,i,j) = weightedaverage(permute(squeeze(O3AnomalySWOOSH.percent_residuals_months(i:12:end,j,latind)),[2,1]),SWOOSH.lat(latind));
        watemp2(:,i,j) = weightedaverage(permute(squeeze(O3AnomalySWOOSH.percent_residuals(i:12:end,j,latind)),[2,1]),SWOOSH.lat(latind));
        bmonth(i,j,:) = regress(watemp(:,i,j),[ones(size(watemp,1),1),[1:size(watemp,1)]']);
        bmonth2(i,j,:) = regress(watemp2(:,i,j),[ones(size(watemp2,1),1),[1:size(watemp2,1)]']);
    end
    
end

%% plotting

 %% testing with contour plot
climits = [-8 8];
clearvars btoplot pval2plot
btoplot(1,:,:) = ball(:,:,2)*120;    
%btoplot(2,:,:) = bmonthall(:,:,2)*120;
btoplot(2,:,:) = braw(:,:,2)*120;    
%btoplot(4,:,:) = bmonthallARgone(:,:,2)*120;
btoplot = permute(btoplot,[1,3,2]);

%pval2plot(1,:,:) = pval;
pval2plot(1,:,:) = uncertainty.sig;
%pval2plot(2,:,:) = umonthraw.sig;
pval2plot(2,:,:) = uraw.sig;
pval2plot = permute(pval2plot,[1,3,2]);

%btoplot (btoplot < climits(1)) = climits(1);
%btoplot(:,[1,2,end-1,end],:) = NaN;

prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
    5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
logprestick = log(prestick);

titles = {['SWOOSH ', num2str(Stimeperiod(1)), char(8211), num2str(Stimeperiod(2)), ' yearly ozone trends (multiple linear regression)'],...
    ['SWOOSH ', num2str(Stimeperiod(1)), char(8211), num2str(Stimeperiod(2)), ' yearly ozone trends (simple linear trend)'],'3','4'};%'SWOOSH 1984-1999 yearly ozone trends (full regression)',    

%bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
% bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
%bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
%bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
subplotmaps(btoplot,SWOOSH.lat,log(SWOOSH.level),{'div','RdBu'},2,pval2plot,20,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
    [-8 8],24,-90:10:90,-90:10:90,...
    fliplr(logprestick),fliplr(presticklabel),{''} ,1,[-90 90],[log(1) log(200)],1,'none',0,'');
%set(gcf,'position',[100 100 900 700])

abc = get(gcf,'children');
set(abc(1), 'FontName', 'Times');
set(abc(2), 'FontName', 'Times');
set(abc(3), 'FontName', 'Times');
set(abc(4), 'FontName', 'Times');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','SWOOSHmonthregression_eq_',num2str(Stimeperiod(1)),'-',num2str(Stimeperiod(2))];
print(filename,'-depsc');
%export_fig(filename,'-pdf');



    %%
    clearvars bmonthtoplot
    bmonthtoplot(1,:,:) = permute(bmonth(:,:,2),[3,1,2])*10;
    %bmonthtoplot(2,:,:) = permute(bmonth2(:,:,2),[3,1,2])*10;
    prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {300,200,100,[],[],[],[],50,[],30,20,10,[],[],[],[],...
        5,[],3,2,1,[],[],[],[],.5,[],.3,.2,.1};
    logprestick = log(prestick);

    titles = {['Monthly ozone trends between ', num2str(lats(1)),' and ',num2str(lats(2)),'{\circ}','N']};

    %bforplot = permute(reshape(squeeze(nanmean(b(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    % bforplot = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    %bforplot(1,:,:) = permute(reshape(squeeze(byearly(:,:,10,8)),[1,31,96]),[1,3,2])*12*10;
    %bforplot(2,:,:) = permute(reshape(squeeze(nanmean(blinearQBOrm(:,:,:,10,2),1)),[1,31,96]),[1,3,2])*1e6*10;
    subplotmaps(bmonthtoplot,1:12,log(SWOOSH.level),{'div','RdBu'},1,[],22,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
        [-4 4],22,1:12,1:12,...
        fliplr(logprestick),fliplr(presticklabel),{''} ,1,[1 12],[log(.1) log(30)],1,'-',0,'');

%     filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/',...
%         'SWOOSHmonthregression_',num2str(SDtimeperiod(1)),'-',num2str(SDtimeperiod(2)),num2str(lats(1)),'and',num2str(lats(2)),'N'];
% 
%     export_fig(filename,'-png');