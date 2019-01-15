% Compare WACCM to SWOOSH gas phase trends

clear all

fsize = 18;
percent = 0;
% import MLS data
lats = [60 70];
compareyears = [2000,2015];
highclcompareyears = [1998,2024];

for p = 1:1
clearvars -except lats fsize plot_MLScomp plot_WACCMonly percent compareyears plotSWOOSHwaccm highclcompareyears

%% import SWOOSH
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201611-latpress-10deg-L31.nc');
SWOOSHlatindex = find(SWOOSH.lat > lats(1) & SWOOSH.lat < lats(2));
SWOOSH.combinedo3q = SWOOSH.combinedo3q(:,:,1:end-11);
SWOOSHyears = 1984:2015;
SWOOSHtimeindex = find(SWOOSHyears >= compareyears(1) & SWOOSHyears <= compareyears(2));

%% import WACCM data
commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[NumberDensity, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes, Longitudes,datedata] = ReadWACCMvertical('O3','monthly',commondir,0);
latindex = find(Latitudes >= lats(1) & Latitudes <= lats(2));

MolConc = rmfield(MolConc, 'MAM1990');
datedata = rmfield(datedata, 'MAM1990');
Pressure = rmfield(Pressure, 'MAM1990');
runs = fields(MolConc);

%% Read in highcl runs
directory = '/Volumes/ExternalOne/work/data/predruns/O3/highCl/zonalmean/';
files = dir([directory,'*.nc']);
for i = 1:length(files)
    [~,highcldata(i),~] = Read_in_netcdf([directory,files(i).name]);
    highcldataPressure(i).p =  permute(repmat(highcldata(i).hyam*100000,[1,size(highcldata(i).PS)]),[2,1,3]) + ...
        permute(repmat(highcldata(i).hybm,[1,size(highcldata(i).PS)]),[2,1,3]) .* ...
        double(permute(repmat(highcldata(i).PS,[1,1,length(highcldata(i).lev)]),[1,3,2]));       
    if i == 1
        highclyears = 1995:2024;
    end
    for j = 1:size(highcldata(i).O3,1)
        [higcl_dataRegPres(i).h(:,j,:),~] = intRegPres(squeeze(highcldata(i).O3(j,:,:)),...
            squeeze(highcldataPressure(i).p(j,:,:))./100,SWOOSH.level);
    end
    for j = 1:12
        higcl_ra(i).h(j,:,:,:) = higcl_dataRegPres(i).h(:,:,j:12:end);
    end
        
    hcldateindex = find(highclyears >= highclcompareyears(1) & highclyears <= highclcompareyears(2));    
    higcl_ra_te(i).h = permute(circshift(squeeze(nanmean(higcl_ra(i).h(:,:,latindex,hcldateindex),3)),[7,0,0]),[3,1,2]);
end
higcl_ra_te(i+1).h = nanmean(cat(4,higcl_ra_te(:).h),4);

%% Read in Uwind
directory = '/Volumes/ExternalOne/work/data/predruns/U/highCl/';
filesU = dir([directory,'*.nc']);
for i = 1:length(filesU)
    [~,highclU(i),~] = Read_in_netcdf([directory,filesU(i).name]);
    highclUzonalmean(i).U = squeeze(nanmean(highclU(i).U(:,:,[2,7],:),1));
    for j = 1:12
        highclU_ra(i).U(j,:,:) = highclUzonalmean(i).U(:,j:12:end);
    end
    highclU_ra_te(i).U = permute(circshift(highclU_ra(i).U(:,:,hcldateindex),[7,0,0]),[3,1,2]);
end

highclU_ra_te(i+1).U = nanmean(cat(4,highclU_ra_te(:).U),4);

%% testing something els

for i = 1:9
    if mod(i,2)
        higcl_ra_te_shift(:,:,:,i) = higcl_ra_te(i).h(1:end-1,:,:);    
    else
        higcl_ra_te_shift(:,:,:,i) = higcl_ra_te(i).h(2:end,:,:);
    end
end

higcl_ra_te_shift = nanmean(higcl_ra_te_shift,4);
plot(higcl_ra_te_shift(:,7,27));

%% remove QBO
for i = 1:length(higcl_ra_te)
     [higcl_ra_te_QBOremoved(i).h,higcl_ra_te_QBOremoved_nl(i).h,b] = removeQBO(higcl_ra_te(i).h,highclU_ra_te(i).U);
end

%% testing
abc = nanmean(cat(4,higcl_ra_te_QBOremoved(1:9).h),4);
figure
plot(higcl_ra_te(10).h(:,9,27),'k')
hold on
plot(higcl_ra_te_QBOremoved_nl(10).h(:,9,27),'r')
plot(higcl_ra_te_QBOremoved(10).h(:,9,27),'b')
plot(higcl_ra_te_shift(:,9,27),'g');
%plot(abc(:,7,27),'g');
%plot(abc2(:,7,27),'c');
%%

for j = 1:length(runs)
    
    CCMIy.(runs{j}) = CCMI_years(datedata.(runs{j}).date,1);
    
    for i = 1:size(MolConc.(runs{j}),1)
        [dataRegPres.(runs{j})(:,i,:),~] = intRegPres(squeeze(MolConc.(runs{j})(i,:,:)),...
            squeeze(Pressure.(runs{j})(i,:,:))./100,SWOOSH.level);
    end
    if j == 1 || j > 3
        dataRegPres.(runs{j}) = cat(3,dataRegPres.(runs{j}),...
            zeros(size(dataRegPres.(runs{j}),1),size(dataRegPres.(runs{j}),2),1));
        %MolConc.(runs{j}) = cat(3,MolConc.(runs{j}),zeros(96,88,1));
    %    Pressure.(runs{j}) = cat(3,Pressure.(runs{j}),zeros(96,88,1));
    else
        %MolConc.(runs{j}) = cat(3,MolConc.(runs{j}),zeros(96,88,6));
        dataRegPres.(runs{j}) = cat(3,dataRegPres.(runs{j}),...
            zeros(size(dataRegPres.(runs{j}),1),size(dataRegPres.(runs{j}),2),6));
    %    Pressure.(runs{j}) = cat(3,Pressure.(runs{j}),zeros(96,88,6));
    end    
    
    MolConcra.(runs{j}) = zeros(12,size(dataRegPres.(runs{j}),1),size(dataRegPres.(runs{j}),2),ceil(size(dataRegPres.(runs{j}),3)/12));
    for i = 1:12    
        MolConcra.(runs{j})(i,:,:,:) = dataRegPres.(runs{j})(:,:,i:12:end);
    end
    MolConcra.(runs{j}) (MolConcra.(runs{j}) == 0) = NaN;
    MolConcra.(runs{j}) = permute(MolConcra.(runs{j}),[4,1,2,3]);
    MolConcralatmean.(runs{j}) = circshift(squeeze(nanmean(MolConcra.(runs{j})(:,:,:,latindex),4)),[0,7,0]);
    WAdateindex = find(CCMIy.(runs{j})(1:12:end) >= compareyears(1) & CCMIy.(runs{j})(1:12:end) <= compareyears(2));    
    MolConcralatmean.(runs{j}) = MolConcralatmean.(runs{j})(WAdateindex,:,:,:);    
end

%% interpolate data onto regular pressure

%%
for i = 1:12
    SWOOSHmontharrange(i,:,:,:) = SWOOSH.combinedo3q(:,:,i:12:end); 
end
SWOOSHlatmean = permute(circshift(squeeze(SWOOSHmontharrange(:,SWOOSHlatindex,:,SWOOSHtimeindex)),[7,0,0]),[3,1,2]);

%% regression to produce linear trends

for j = 1:size(MolConcralatmean.MAM,2)
    for k = 1:size(MolConcralatmean.MAM,3)      
        for l = 1:length(runs)
            [bWAC.(runs{l})(j,k,:),bintWAC.(runs{l})(j,k,:,:),~,~,bWACstats.(runs{l})(j,k).s] = regress(squeeze(MolConcralatmean.(runs{l})(:,j,k)),[ones(1,size(MolConcralatmean.(runs{l}),1));1:size(MolConcralatmean.(runs{l}),1)]');
        end
    end
end

for j = 1:size(higcl_ra_te(1).h,2)
    for k = 1:size(higcl_ra_te(1).h,3)      
        for l = 1:length(higcl_ra_te)
            [bhighcl(l).h(j,k,:),binthighcl(l).h(j,k,:,:),~,~,bhighclstats(i).h(j,k).s] = regress(squeeze(higcl_ra_te(l).h(:,j,k)),[ones(1,size(higcl_ra_te(l).h,1));1:size(higcl_ra_te(l).h,1)]');
            [bhighcl_QBOrm(l).h(j,k,:),binthighcl_QBOrm(l).h(j,k,:,:),~,~,bhighclstats_QBOrm(i).h(j,k).s] = regress(squeeze(higcl_ra_te_QBOremoved_nl(l).h(:,j,k)),[ones(1,size(higcl_ra_te_QBOremoved_nl(l).h,1));1:size(higcl_ra_te_QBOremoved_nl(l).h,1)]');
        end
    end
end

for i = 1:length(runs)
    bWACpercent.(runs{i})(:,:,:) = bWAC.(runs{i})(:,:,2)./squeeze(nanmean(MolConcralatmean.(runs{i}),1))*100;
end
% Volcanoes MAM - VC-MAM
b_Volcanoes = bWAC.MAM - bWAC.VCMAM;
b_Volcanoespercent = bWACpercent.MAM - bWACpercent.VCMAM;
% Dynamics VCMAM - ChemOnly
b_dynamics = bWAC.VCMAM - bWAC.Chemonly;
b_dynamicspercent = bWACpercent.VCMAM - bWACpercent.Chemonly;

for j = 1:size(SWOOSHlatmean,2)        
    for k = 1:size(SWOOSHlatmean,3)
         [bSWOOSH(j,k,:),bSWOOSH_int(j,k,:,:),~,~,SWOOSHstats(j,k).s] = regress(squeeze(SWOOSHlatmean(:,j,k)),[ones(1,size(SWOOSHlatmean,1));1:size(SWOOSHlatmean,1)]');
    end
end    
bpercentSWOOSH(:,:,:) = bSWOOSH(:,:,2)./squeeze(nanmean(SWOOSHlatmean(:,:,:),1))*100;

% %% if plotting swoosh WACCM comparison
% 
% if plotSWOOSHwaccm
%     if percent
%         bswooshforplot = reshape(bpercentSWOOSH,[1,size(bpercentSWOOSH)]);
%     else
%         bswooshforplot = reshape(bSWOOSH(:,:,2),[1,size(bSWOOSH(:,:,2))])*10;
%     end
% end
% 
% prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
%     presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
%         5,[],[],2,1,[],[],[],[],.5,[],[],.2,.1};
%         logprestick = log(prestick);
% plotmtitle = ['Ozone trends over 2005-2015 and ',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'{\circ}S'];
% subplotmaps(bswooshforplot,1:12,log(SWOOSH.level),{'div','RdBu'},1,[],18,{'SWOOSH'},'Month','Pressure (hPa)','%/decade','on',...
%         [-.6 .6],14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
%         fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(.1) log(300)],1,'none',0,'');    
%     
% % plot WACCM
%% PLOTTING HIGHCL
prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
    5,[],[],2,1,[],[],[],[],.5,[],[],.2,.1};
logprestick = log(prestick);

titles = {'No .1','No. 2','No. 3','No. 4','No. 5','No. 6','No. 7','No. 8','No. 9','Ensave'};

bhighclforplot = cat(4,bhighcl_QBOrm(:).h);
bhighclforplot = permute(squeeze(bhighclforplot(:,:,2,:)),[3,1,2])*10*1e6;

subplotmaps(bhighclforplot,1:12,log(SWOOSH.level),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','ppmv/decade','on',...
        [-.4 .4],26,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),'' ,1,[1 12],[log(.1) log(300)],1,'none',0,'');
    
export_fig(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/solarcycleandtrends/swooshWACCMcompare/percenthighcl_rmQBO_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S','.pdf'],'-pdf');
    
%% plotting MLS WACCM comparison 2005 - 2015
if percent
    bforplot = permute(cat(3,bpercentSWOOSH(:,:),bWACpercent.MAM(:,:), bWACpercent.VCMAM(:,:), bWACpercent.Chemonly(:,:),...
        b_Volcanoespercent(:,:), b_dynamicspercent(:,:)),[3,1,2])*10;
    %bforplot = permute(cat(3,bpercent(:,:),bpercentSWOOSH(:,:),bWACpercent.MAM(:,:), bWACpercent.VCMAM(:,:), bWACpercent.Chemonly(:,:),...
    %    b_dynamicspercent(:,:)),[3,1,2])*10;
else
    bforplot = permute(cat(3,bSWOOSH(:,:,2)./1e6,bWAC.MAM(:,:,2), bWAC.VCMAM(:,:,2), bWAC.Chemonly(:,:,2),...
        b_Volcanoes(:,:,2), b_dynamics(:,:,2)),[3,1,2])*1e6*10;
end

prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1,.9:-.1:.1];
presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
    5,[],[],2,1,[],[],[],[],.5,[],[],.2,.1};
logprestick = log(prestick);

titles = {'SWOOSH','MAM','VC-MAM','Chem-only','Volcanoes (MAM - VC-MAM)','Dynamics (VC-MAM - Chem-only)'};
if lats(1) < 0 && lats(2) <= 0
    plotmtitle = ['Ozone trends over ',num2str(compareyears(1)),'-',num2str(compareyears(2)),' and ',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'{\circ}S'];
elseif lats(1) >= 0 && lats(2) > 0
    plotmtitle = ['Ozone trends over ',num2str(compareyears(1)),'-',num2str(compareyears(2)),' and ',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'{\circ}N'];
else
    plotmtitle = ['Ozone trends over ',num2str(compareyears(1)),'-',num2str(compareyears(2)),' and ',num2str(abs(lats(1))),'{\circ}S','-',num2str(abs(lats(1))),'{\circ}N'];
end
if percent
    subplotmaps(bforplot,1:12,log(SWOOSH.level),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','%/decade','on',...
        [-20 20],12,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(.1) log(300)],1,'none',0,'');
else
    subplotmaps(bforplot,1:12,log(SWOOSH.level),{'div','RdBu'},1,[],18,titles,'Month','Pressure (hPa)','ppmv/decade','on',...
        [-.6 .6],14,1:12,{'J','J','A','S','O','N','D','J','F','M','A','M'},...
        fliplr(logprestick),fliplr(presticklabel),plotmtitle ,1,[1 12],[log(.1) log(300)],1,'none',0,'');
end
if percent    
    if lats(1) < 0 && lats(2) <= 0
        export_fig(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/solarcycleandtrends/swooshWACCMcompare/percentSWOOSH_WACCM_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S','.pdf'],'-pdf');
    elseif lats(1) >= 0 && lats(2) > 0
        export_fig(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/solarcycleandtrends/swooshWACCMcompare/percentSWOOSH_WACCM_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'N','.pdf'],'-pdf');
    else
        export_fig(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/solarcycleandtrends/swooshWACCMcompare/percentSWOOSH_WACCM_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'SN','.pdf'],'-pdf');
    end
else
    if lats(1) < 0 && lats(2) <= 0
        export_fig(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/solarcycleandtrends/swooshWACCMcompare/SWOOSH_WACCM_',num2str(compareyears(1)),'-',num2str(compareyears(2)),'_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S','.pdf'],'-pdf');
    elseif lats(1) >= 0 && lats(2) > 0
        export_fig(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/solarcycleandtrends/swooshWACCMcompare/SWOOSH_WACCM_',num2str(compareyears(1)),'-',num2str(compareyears(2)),'_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'N','.pdf'],'-pdf');
    else
        export_fig(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/solarcycleandtrends/swooshWACCMcompare/SWOOSH_WACCM_',num2str(compareyears(1)),'-',num2str(compareyears(2)),'_',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'SN','.pdf'],'-pdf');
    end
end
%%
%lats = lats+10;
close all
end

%% now find locations and plot line plots

lev = 2;
monind = 9;
monind2 = monind - 5;
if monind > 12
    monind = monind -12;
end
[~,levind] = min(abs(SWOOSH.level - lev));

for i = 1:length(runs)
    WACCMlines(:,i) = squeeze(MolConcralatmean.(runs{i})(:,monind2,levind))*1e6;
end

SWOOSHlines = squeeze(SWOOSHlatmean(:,monind2,levind));

createfig('medium','on');
ph = plot(compareyears(1):compareyears(2),WACCMlines(:,[2,3,4]),'LineWidth',2);
hold on
sph = plot(compareyears(1):compareyears(2),SWOOSHlines,'LineWidth',2);
highcllines = cat(4,higcl_ra_te(:).h);
highcllines = squeeze(highcllines(:,monind2,levind,:))*1e6;
highcllines2 = cat(4,higcl_ra_te_QBOremoved_nl(:).h);
highcllines2 = squeeze(highcllines2(:,monind2,levind,:))*1e6;
%plot(compareyears(1):compareyears(2),highcllines,'r');
hph = plot(highclcompareyears(1):highclcompareyears(2),nanmean(highcllines(:,1:9),2),'k','LineWidth',3);
hpl = plot(highclcompareyears(1)+1:highclcompareyears(2),nanmean(highcllines2(:,1:9),2),'color',[.7 .7 .7],'LineWidth',3);
xlim([1997,2025])
lh = legend([ph;sph;hph],'MAM','VC-MAM','Chem-only','SWOOSH','HighCl ensemble');
set(lh,'box','off','fontsize',18,'location','SouthEast');
xlabel('Year','fontsize',20);
ylabel('ppmv','fontsize',20);
title(['Ozone at ',num2str(lev),' hPa - ',monthnames(monind,0,0),' - ',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'N'],'fontsize',20);
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/solarcycleandtrends/swooshWACCMcompare/',...
    num2str(lev),'hPa-',monthnames(monind,0,0),'-',num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'N'];
export_fig(filename,'-pdf');
