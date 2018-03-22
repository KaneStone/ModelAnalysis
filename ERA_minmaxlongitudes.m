% ERA-interim minimum and maximum
clear all 
ERAdata = ReadinERA('/Volumes/MyBook/work/data/ERA-Interim/T/3090S_50hPa_T_ERA.nc');
years(1).y = repmat(1979:2016,12,1);
years(1).y = [years(1).y(:);2017;2017];

ERAdata4060S = ReadinERA('/Volumes/MyBook/work/data/ERA-Interim/T/wgt4060S_T_ERA-Interim.nc');
plotlines = 1;
plotcorr = 1;
%% finding minimum and maximum longitudes
timeperiod = [2003,2016];
lats = [-50,-55,-60,-65,-70,-75,-80,-85];
yeartoremove = 2002;
dateind = find(years(1).y >= timeperiod(1) & years(1).y <= timeperiod(2));         
%dateind = find(years(1).y >= timeperiod(1) & years(1).y < yeartoremove | years(1).y > yeartoremove & years(1).y <= timeperiod(2));         
%is2002 = find(years(1).y == dateind);

    ERAtemperature = squeeze(ERAdata.t);
    % testing
%     spring = squeeze(nanmean(cat(4,ERAtemperature(:,24,dateind(1)+9-1:12:dateind(end)),...
%         ERAtemperature(:,24,dateind(1)+10-1:12:dateind(end)),...
%         ERAtemperature(:,24,dateind(1)+11-1:12:dateind(end))),4));
    
    spring = squeeze(nanmean(cat(4,ERAtemperature(:,24,dateind(9:12:end)),...
        ERAtemperature(:,24,dateind(10:12:end)),...
        ERAtemperature(:,24,dateind(11:12:end))),4));

    amplitude = max(spring) - min(spring); 
    ampmean = nanmean(amplitude);
    ampstd = std(amplitude);
if plotlines
    fsize = 20;
    createfig('medium','on');
    cbrew2 = cbrewer('qual','Paired',10);
    plot(ERAdata.longitude,spring,'color',cbrew2(1,:),'LineWidth',2)
    hold on
    %plot(ERAdata.longitude,spring(:,24),'color',cbrew2(6,:),'LineWidth',2)
    plot(ERAdata.longitude,nanmean(spring,2),'color',cbrew2(2,:),'LineWidth',5);


    set(gca,'fontsize',fsize+2);
    xlim([-5 365]);
    xlabel('Longitude ({\circ}E)','fontsize',fsize+4);
    ylabel('Temperature (K)','fontsize',fsize+4);
    title('Austral spring temperatures at 65{\circ}S','fontsize',fsize+6);
    %lh = legend(ph,'1995-2025 (high chlorine)','1955-1979 (low chlorine)','1995-2025 (strat ozone only)');
    %set(lh,'fontsize',fsize+4,'box','off','location','south')

    filename = '/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TalkFigures/ERAAmplitudes_65S';
    export_fig(filename,'-png');
end

%%
if plotcorr
    ERAspring = nanmean(cat(4,ERAdata4060S.t(:,:,dateind(1)+9-1:12:dateind(length(dateind))),ERAdata4060S.t(:,:,dateind(1)+10-1:12:dateind(length(dateind))),ERAdata4060S.t(:,:,dateind(1)+11-1:12:dateind(length(dateind)))),4);
    for i = 1:size(ERAspring,1)
        tic;
        for j = 1:size(ERAspring,2)
           r(1,i,j) = corr(squeeze(ERAspring(i,j,:)),amplitude');
        end
    end
end

%%
    cbrew = cbrewer('div','RdBu',16);         

    prestick = [1000:-100:300,200,100,90:-10:10,9:-1:1];

    presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
        5,[],[],2,1};


    logprestick = log(prestick);

    %contourtitle = {'Correlation between 65{\circ}S Amplitude and 40-60{\circ}S temperature'};       
    contourtitle2 = {'1995-2016 (high chlorine)'};


    subplotmaps(r,ERAdata.longitude,log(double(ERAdata4060S.level)),{'div','RdBu'},1,[],18,contourtitle2,'Longitude','','Correlation','on',...
        [-1,1],18,0:20:360,0:20:360,fliplr(logprestick),fliplr(presticklabel),'',1,[0 360],[log(1) log(1000)],1,'-',0,'none');

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/ZonalAssym/TempforTalk/ERAamp4060Scorr_spr',];
    export_fig(filename,'-pdf');

%%

for k = 1:length(lats)
    for j = 1:12
        [~,latind] = min(abs(lats(k)-ERAdata.latitude));            
        [tempminvalue,tempminind] = min(squeeze(ERAtemperature(:,latind,j:12:end-2)));        
        [tempmaxvalue,tempmaxind] = max(squeeze(ERAtemperature (:,latind,j:12:end-2)));                
        ERA.minvalue(j,k,:) = tempminvalue;
        ERA.minind(j,k,:) = tempminind;              
        ERA.maxvalue(j,k,:) = tempmaxvalue;
        ERA.maxind(j,k,:) = tempmaxind;

    end
    [tempminvalue,tempminind] = min(nanmean(cat(3,squeeze(ERAtemperature(:,latind,9:12:end)),...
        squeeze(ERAtemperature(:,latind,10:12:end)),...
        squeeze(ERAtemperature(:,latind,11:12:end))),3),[],1);        
    [tempmaxvalue,tempmaxind] = max(nanmean(cat(3,squeeze(ERAtemperature(:,latind,9:12:end)),...
        squeeze(ERAtemperature(:,latind,10:12:end)),...
        squeeze(ERAtemperature(:,latind,11:12:end))),3),[],1);               
    ERA.minvalue_spr(k,:) = tempminvalue;
    ERA.minind_spr(k,:) = tempminind;                                    
    ERA.maxvalue_spr(k,:) = tempmaxvalue;
    ERA.maxind_spr(k,:) = tempmaxind;
end    

ERA.minlongitude = ERAdata.longitude(ERA.minind);
ERA.maxlongitude = ERAdata.longitude(ERA.maxind);
ERA.amplitude = ERA.maxvalue - ERA.minvalue;

ERA.minlongitude_spr = ERAdata.longitude(ERA.minind_spr);
ERA.maxlongitude_spr = ERAdata.longitude(ERA.maxind_spr);
ERA.amplitude_spr = ERA.maxvalue_spr - ERA.minvalue_spr;

for k = 1:length(lats)     
    for j = 1:12
        for l = 1:size(ERA.maxlongitude(j,k,:),3)
            if ERA.maxlongitude(j,k,l) > 300
                ERA.maxlongitude(j,k,l) = ERA.maxlongitude(j,k,l)-360;
            end
            if ERA.minlongitude(j,k,l) < 100
                ERA.minlongitude(j,k,l) = ERA.minlongitude(j,k,l)+360;
            end            
        end
    end

    for l = 1:size(ERA.maxlongitude_spr(k,:),2)
        if ERA.maxlongitude_spr(k,l) > 300
            ERA.maxlongitude_spr(k,l) = ERA.maxlongitude_spr(k,l)-360;
        end
        if ERA.minlongitude_spr(k,l) < 100
            ERA.minlongitude_spr(k,l) = ERA.minlongitude_spr(k,l)+360;                           
        end
    end        

end
ERA.years = years;
ERA.latitude = lats;
ERA.longitude = ERAdata.longitude;
save('/Volumes/MyBook/work/data/matlabOutput/zonalAnalysis/ERA.mat','ERA');

clearvars dateind

%% plotting amplitude
createfig('medium','on')
contourf(amplitude,14)
%caxis([0,22]);
colorbar
title('ERA-Interim');

%% plotting
fsize = 18;
cbrew = cbrewer('qual','Set1',10);
createfig('large','on');
%subplot(2,1,1)
m_proj('Stereographic','lon',-180,'lat',-90,'rad',50)
m_coast('color','k','LineWidth',3);
m_grid('ytick',[-90 -80 -70 -60 -50 -40 -30 -20],'XaxisLocation','top','fontsize',fsize);

m_line(squeeze(minlongitude(10,:))',lats','color',cbrew(1,:),'LineWidth',3);    
m_line(squeeze(maxlongitude(10,:))',lats','color',cbrew(1,:),'LineWidth',3);
    
title([num2str(timeperiod(1)),'-',num2str(timeperiod(2))],'position',[0,.98,],'fontsize',20);
