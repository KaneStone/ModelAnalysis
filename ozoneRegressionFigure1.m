% produce Paper figure 1

clear all
%% Read in All HighCl data

SWOOSH.level = [316.22775;261.01572;215.44347;177.82794;146.77992;121.15276;100;82.540421;...
    68.129204;56.234131;46.415890;38.311867;31.622776;26.101572;21.544348;17.782795;14.677993;...
    12.115276;10;8.2540417;6.8129206;5.6234131;4.6415887;3.8311868;3.1622777;2.6101573;2.1544347;...
    1.7782794;1.4677993;1.2115277;1];
%% import highCl
highCllevel = [SWOOSH.level;[.9,.8,.7,.6,.5,.4,.3,.2,.1]'];
highClTimperiod = [1998 2016];
[highClData,WAClat] = ReadInHighClRegression(highCllevel,highClTimperiod,'highCl');

%% import Chem-only and MAM
SDtimeperiod = [2000 2014];
[SDWaccmData] = ReadInSDWACCM(highCllevel,SDtimeperiod);

%% take regression of SD WACCM (chem-only and MAM)

[bchemonly,predictorschemonly,O3Anomalychemonly,uncertaintychemonly] = ...
    ozoneRegressionTrends(SDWaccmData(2).O3,SDWaccmData(2).U,SDWaccmData(2).solar,SDWaccmData(2).SPE,...
    SDWaccmData(2).NINO34,SDWaccmData(2).HFS,SDWaccmData(2).HFN,SDWaccmData(2).NO2,0,SDtimeperiod,1,0,0);%SDWaccmData(2).SPEintpres 

%% Take regression of residuals
[btest,ballchemonly,ballchem_pvalue,uncall_chemonly] = ozoneRegressionEnsAve(O3Anomalychemonly.percent_residuals_months,SDtimeperiod(2) - SDtimeperiod(1));
[btest2,ballchemonly2] = ozoneRegressionEnsAve(O3Anomalychemonly.percent_residuals,SDtimeperiod(2) - SDtimeperiod(1));

sig = .01;
ballchem_pvalue_toplot = ballchem_pvalue;
ballchem_pvalue_toplot (ballchem_pvalue_toplot <= sig) = 0;
ballchem_pvalue_toplot (ballchem_pvalue_toplot > sig) = 1;


%% take regression of highcl
clearvars O3Anomaly predictors b
for i = 1:length(highClData)
    tic;
   [b(i),predictors(i,:,:),O3Anomaly(i),uncertainty(i,:,:)] = ...
       ozoneRegressionTrends(highClData(i).O3,highClData(i).U,0,0,highClData(i).NINO34,...
       highClData(i).HFsouth,highClData(i).HFnorth,0,0,highClTimperiod,0,0,0); % 
   toc;
end


O3Anomaly(11).percent_residuals_months = nanmean(cat(4,O3Anomaly(1:9).percent_residuals_months),4);

% O3Anomaly(11).percent_residuals = nanmean(cat(4,O3Anomaly(1:9).percent_residuals),4);

%% second time period

%% import highCl second time period
highClTimperiod2 = [1998 2024];
[highClData2,WAClat2] = ReadInHighClRegression(highCllevel,highClTimperiod2,'highCl');


%% take regression of highcl second time period
clearvars O3Anomaly2 predictors2 b2
for i = 1:length(highClData2)
    tic;
   [b2(i),predictors2(i,:,:),O3Anomaly2(i),uncertainty2(i,:,:)] = ...
       ozoneRegressionTrends(highClData2(i).O3,highClData2(i).U,0,0,highClData2(i).NINO34,...
       highClData2(i).HFsouth,highClData2(i).HFnorth,0,0,highClTimperiod2,0,0,0); % 
   toc;
end

O3Anomaly2(11).percent_residuals_months = nanmean(cat(4,O3Anomaly2(1:9).percent_residuals_months),4);
O3Anomaly2(11).percent_residuals= nanmean(cat(4,O3Anomaly2(1:9).percent_residuals),4);

%%

[~,ball,~,uncall] = ozoneRegressionEnsAve(O3Anomaly(11).percent_residuals_months,highClTimperiod(2) - highClTimperiod(1));
[~,ball2,~,uncall2] = ozoneRegressionEnsAve(O3Anomaly2(11).percent_residuals_months,highClTimperiod2(2) - highClTimperiod2(1));

% sig = .01;
% ball_pvalue_toplot = ball_pvalue;
% ball_pvalue_toplot (ball_pvalue_toplot <= sig) = 0;
% ball_pvalue_toplot (ball_pvalue_toplot > sig) = 1;

%% individual global trends
for i = 1:9
    [bhclind(i).b,ballhclind(i).b,pvalue(i).b,uncallhighcl(i)] = ozoneRegressionEnsAve(O3Anomaly(i).percent_residuals,highClTimperiod(2) - highClTimperiod(1));
end

for i = 1:9
    [bhclind2(i).b,ballhclind2(i).b,pvalue(i).b,uncallhighcl2(i)] = ozoneRegressionEnsAve(O3Anomaly2(i).percent_residuals,highClTimperiod2(2) - highClTimperiod2(1));
end

%% plot select
ball1cat = cat(4,ballhclind(:).b);
ball2cat = cat(4,ballhclind2(:).b);
bindall = permute(squeeze(ball1cat(:,:,2,:)*120),[3,2,1]);
bindall2 = permute(squeeze(ball2cat(:,:,2,:)*120),[3,2,1]);

bstd = std(bindall,1,1);
bstd2 = std(bindall2,1,1);

bcat = cat(4,ballhclind([2,3,5,9]).b);%,ball(:,:,2),ball2(:,:,2)):%,bstd,bstd2);
ucat = squeeze(cat(3,uncallhighcl([2,3,5,9]).sig));%,ball(:,:,2),ball2(:,:,2)):%,bstd,bstd2);
bindtoplottemp = permute(squeeze(bcat(:,:,2,:)*120),[3,2,1]);
uindtoplottemp = permute(squeeze(ucat),[3,2,1]);

ball1temp = permute(ball(:,:,2)*120,[3,2,1]);
ball2temp = permute(ball2(:,:,2)*120,[3,2,1]);

uall1temp(1,:,:) = uncall.sig;
uall2temp(1,:,:) = uncall2.sig;

bindtoplot = cat(1,bindtoplottemp,ball1temp,ball2temp,bstd,bstd2);
uindtoplot = cat(1,uindtoplottemp,ones(size(bstd))*-1,ones(size(bstd))*-1,ones(size(bstd))*-1,ones(size(bstd))*-1);
%uindtoplot = cat(1,uindtoplottemp,permute(uall1temp,[1,3,2]),permute(uall2temp,[1,3,2]),ones(size(bstd))*-1,ones(size(bstd))*-1);

%%

prestick = [300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
    presticklabel = {[],[],100,[],[],[],[],[],[],30,[],10,[],[],[],[],...
        [],[],3,[],1,[],[],[],[],[],[],.3,[],.1};
    logprestick = log(prestick);

titles = {['No. 2, ', num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))],...
    ['No. 3, ', num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))],...
    ['No. 5, ',num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))],...
    ['No. 9, ',num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))],...
    ['Ens-ave, ',num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))],...
    ['Ens-ave, ',num2str(highClTimperiod2(1)),char(8211),num2str(highClTimperiod2(2))],...,
    ['Ens-ave std of linear trends, ',num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))],...
    ['Ens-ave std of linear trends, ',num2str(highClTimperiod2(1)),char(8211),num2str(highClTimperiod2(2))]};

mtit = {['Individual ensemble member global ozone linear trends over ', num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))]};

fig = subplotmaps(bindtoplot(1:8,:,:),WAClat,log(highCllevel),{'div','RdBu'},1,uindtoplot(1:8,:,:),16,titles,'Latitude','Pressure (hPa)','Percent/decade','on',...
    [-8 8],22,-90:30:90,-90:30:90,...
    fliplr(logprestick),fliplr(presticklabel),mtit ,1,[-90 90],[log(.1) log(200)],1,'-',0,'');

% abc = get(gcf,'children');
% for i = 1:10
%     set(abc(i), 'FontName', 'Times');
% end

fig.Renderer='Painters';
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/Paper/Draft/Figures/','Figure1_withsig_onlysingle_lines'];

export_fig(filename,'-pdf');
print(filename,'-depsc');
%saveas(fig,filename,'epsc')
%%
subplotmaps(bstd,WAClat,log(highCllevel),{'div','RdBu'},1,[],20,{['Ens-ave standard deviation in linear trends ',num2str(highClTimperiod(1)),char(8211),num2str(highClTimperiod(2))]},'Latitude','Pressure (hPa)','Percent/decade','on',...
    [-4 4],22,-90:30:90,-90:30:90,...
    fliplr(logprestick),fliplr(presticklabel),{''},1,[-90 90],[log(1) log(200)],1,'-',0,'');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/','Ensave_indmem_std',num2str(highClTimperiod(1)),'-',num2str(highClTimperiod(2))];

set(gcf,'position',[100 100 900 700])

export_fig(filename,'-pdf');