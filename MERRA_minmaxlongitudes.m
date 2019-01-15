%MERRA_minmaxlongitudes

% MERRA-interim minimum and maximum
clear all 
[~,MERRAdata,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/CESM-CCMI/T/50hPa/3090S_50hPa_T_f.e11.FWTREFC1SD.f19.f19.ccmi30.001.nc');
years(1).y = repmat(1979:2014,12,1);
%years(1).y = [years(1).y(:);2017;2017];

[~,MERRAdata4060S,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/CESM-CCMI/T/4060S/wgt_4060S_T_f.e11.FWTREFC1SD.f19.f19.ccmi30.001.nc');
plotlines = 1;
plotcorr = 1;
%% finding minimum and maximum longitudes
timeperiod = [1995,2014];
lat = -50;
[~,latind] = min(abs(MERRAdata.lat - lat));
yeartoremove = 2002;

dateind = find(years(1).y >= timeperiod(1) & years(1).y <= timeperiod(2));         
%dateind = find(years(1).y >= timeperiod(1) & years(1).y < yeartoremove | years(1).y > yeartoremove & years(1).y <= timeperiod(2));         
%is2002 = find(years(1).y == dateind);

    MERRAtempMERRAture = squeeze(MERRAdata.T);
    % testing
%     spring = squeeze(nanmean(cat(4,MERRAtempMERRAture(:,24,dateind(1)+9-1:12:dateind(end)),...
%         MERRAtempMERRAture(:,24,dateind(1)+10-1:12:dateind(end)),...
%         MERRAtempMERRAture(:,24,dateind(1)+11-1:12:dateind(end))),4));
    
    spring = squeeze(nanmean(cat(4,MERRAtempMERRAture(:,latind,dateind(9:12:end)),...
        MERRAtempMERRAture(:,latind,dateind(10:12:end)),...
        MERRAtempMERRAture(:,latind,dateind(11:12:end))),4));

    amplitude = max(spring) - min(spring); 
    ampmean = nanmean(amplitude);
    ampstd = std(amplitude);
if plotlines
    fsize = 20;
    createfig('medium','on');
    cbrew2 = cbrewer('qual','Paired',10);
    plot(MERRAdata.lon,spring,'color',cbrew2(1,:),'LineWidth',2)
    hold on
    %plot(MERRAdata.longitude,spring(:,24),'color',cbrew2(6,:),'LineWidth',2)
    plot(MERRAdata.lon,nanmean(spring,2),'color',cbrew2(2,:),'LineWidth',5);


    set(gca,'fontsize',fsize+2);
    xlim([-5 365]);
    xlabel('Longitude ({\circ}E)','fontsize',fsize+4);
    ylabel('TempMERRAture (K)','fontsize',fsize+4);
    title(['Austral spring temperatures at ',num2str(abs(lat)),'{\circ}S'],'fontsize',fsize+6);
    %lh = legend(ph,'1995-2025 (high chlorine)','1955-1979 (low chlorine)','1995-2025 (strat ozone only)');
    %set(lh,'fontsize',fsize+4,'box','off','location','south')

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TalkFigures/MERRAAmplitudes',num2str(abs(lat)),'S_',num2str(timeperiod(1)),'-',num2str(timeperiod(2))];
    export_fig(filename,'-pdf');
end

%%
if plotcorr
    MERRAspring = nanmean(cat(4,MERRAdata4060S.T(:,:,dateind(9:12:length(dateind))),MERRAdata4060S.T(:,:,dateind(10:12:length(dateind))),MERRAdata4060S.T(:,:,dateind(11:12:length(dateind)))),4);
    for i = 1:size(MERRAspring,1)
        tic;
        for j = 1:size(MERRAspring,2)
           r(1,i,j) = corr(squeeze(MERRAspring(i,j,:)),amplitude');
        end
    end
end

%%
    cbrew = cbrewer('div','RdBu',16);         

    prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];

    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
        5,[],[],2,1};


    logprestick = log(prestick);

    %contourtitle = {'Correlation between 65{\circ}S Amplitude and 40-60{\circ}S tempMERRAture'};       
    contourtitle2 = {[num2str(timeperiod(1)),'-',num2str(timeperiod(2)),' (high chlorine)']};


    subplotmaps(r,MERRAdata.lon,log(double(MERRAdata4060S.lev)),{'div','RdBu'},1,[],18,contourtitle2,'Longitude {\circ}E','Pressure (hPa)','Correlation','on',...
        [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),'',1,[0 360],[log(1) log(1000)],1,'-',0,'none');

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/MERRAamp4060Scorr_spr',num2str(abs(lat)),'S_',num2str(timeperiod(1)),'-',num2str(timeperiod(2))];
    export_fig(filename,'-pdf');

%%
