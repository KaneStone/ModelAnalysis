function [] = predruns_cryregression(differences,Seasplot,...
    modlat,modlon,obslat,obslon,inputs,mons,surfacedata,observations,modeltoz,obstoz)
% regression predictions of last 4 and 6 years

% need raw data
% calculate sea ice extent for different areas
% take out last 4 years
% construct regression model using remaining

%% calculate sea ice extent for all years before detrending.

areacell = gridcellarea(modlat,modlon);
areacellobs = gridcellarea(obslat,obslon);

if strcmp(inputs.var,'ICEFRAC')
    modeldata = surfacedata.highCl.dataMonthArrange;
    modeldata (modeldata < inputs.fraction) = 0;
    modeldata (modeldata >= inputs.fraction) = 1;

    obsdata = observations.montharrange;
    obsdata (obsdata < inputs.fraction) = 0;
    obsdata (obsdata >= inputs.fraction) = 1;
else
    modeldata = surfacedata.highCl.dataMonthArrange;
    obsdata = observations.montharrange;
end

seasons = [[4,5];[8,9];[11,12]];
%% extract last 4 years
yearstoextractmodel = 1:26;
yearstoextractobs = 1:31;

for i = 1:size(Seasplot.lon,1)
    if strcmp(inputs.var,'ICEFRAC')
        latind = modlat >= Seasplot.lat(i,1) & modlat <= Seasplot.lat(i,2);
        latindobs = obslat >= Seasplot.lat(i,1) & obslat <= Seasplot.lat(i,2);
        if i == 3
            lonind = modlon >= Seasplot.lon(i,1) | modlon < Seasplot.lon(i,2)-360;
            lonindobs = obslon >= Seasplot.lon(i,1) | obslon < Seasplot.lon(i,2)-360;
        else
            lonind = modlon >= Seasplot.lon(i,1) & modlon < Seasplot .lon(i,2);
            lonindobs = obslon >= Seasplot.lon(i,1) & obslon < Seasplot.lon(i,2);
        end
        alldataextent(i).sea = permute(areacell(latind,lonind).*permute(modeldata(:,:,yearstoextractmodel,latind,lonind),[4,5,1,2,3]),[3,4,5,1,2]);
        restdataextent(i).sea = permute(areacell(latind,lonind).*permute(modeldata(:,:,yearstoextractmodel(end)+1:end,latind,lonind),[4,5,1,2,3]),[3,4,5,1,2]);
        seaiceextenttemp(i).sea = nansum(alldataextent(i).sea(:,:,:,:),4);
        restseaiceextenttemp(i).sea = nansum(restdataextent(i).sea(:,:,:,:),4);

        obsdataextent(i).sea = permute(areacellobs(latindobs,lonindobs).*permute(obsdata(:,yearstoextractobs,lonindobs,latindobs),[4,3,1,2]),[3,4,1,2]);
        restobsdataextent(i).sea = permute(areacellobs(latindobs,lonindobs).*permute(obsdata(:,yearstoextractobs(end)+1:end,lonindobs,latindobs),[4,3,1,2]),[3,4,1,2]);
        obsseaiceextenttemp(i).sea = nansum(obsdataextent(i).sea(:,:,:),3);
        restobsseaiceextenttemp(i).sea = nansum(restobsdataextent(i).sea(:,:,:),3);
    else
        latind = modlat >= Seasplot.lat(i,1) & modlat <= Seasplot.lat(i,2);
        latindobs = obslat >= Seasplot.lat(i,1) & obslat <= Seasplot.lat(i,2);

        lonind = modlon >= Seasplot.lon(i,1) & modlon < Seasplot .lon(i,2);
        lonindobs = obslon >= Seasplot.lon(i,1) & obslon < Seasplot.lon(i,2);

        alldataextent(i).sea = modeldata(:,:,:,latind,lonind);
        
        obsdataextent(i).sea = obsdata(:,:,lonindobs,latindobs);        

        snowextenttemp(i).sea = nanmean(alldataextent(i).sea(:,:,:,:),4);
        obssnowextenttemp(i).sea = nanmean(obsdataextent(i).sea(:,:,:),3);
        
        seaiceextenttemp(i).sea = snowextenttemp(i).sea(:,:,yearstoextractmodel);
        restseaiceextenttemp(i).sea = snowextenttemp(i).sea(:,:,yearstoextractmodel(end)+1:end);
        
        obsseaiceextenttemp(i).sea = obssnowextenttemp(i).sea(:,yearstoextractobs);
        restobsseaiceextenttemp(i).sea = obssnowextenttemp(i).sea(:,yearstoextractobs(end)+1:end);
        
    end
    
    %% remove enso now observations
    laglength = 3;
    if inputs.removeENSO
        for m = 1:12
            for lag = 1:laglength                                                            
                ENSOyearind = observations.ENSO;             
                if m == 1 || m == 2
                    m2 = 3;
                else
                    m2 = m;
                end
                ensopredictors = [squeeze(ENSOyearind(m2-lag+1:12:end-6*12,:)),ones(length(squeeze(ENSOyearind(m2-lag+1:12:end-6*12,:))),1)];
                benso(lag,:,:) = ensopredictors\squeeze(obsseaiceextenttemp(i).sea(m,:))';

                %                     ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonth(m)-lag+1:12:end,:))',ones(size(alldata,3)-ext,1)];
                %                     benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedataind(k,:,m,j,:));                        
                %find max lag
            end
            [~,bensomax_ind] = max(abs(benso),[],1);                
            %[~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   

            bensomax = benso(bensomax_ind(1),1);
            %bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
            obsseaiceextent(i).sea(m,:) = squeeze(obsseaiceextenttemp(i).sea(m,:))' - bensomax.*squeeze(ENSOyearind(m2-bensomax_ind(1)+1:12:end-6*12,:));   
            restobsseaiceextent(i).sea(m,:) = squeeze(restobsseaiceextenttemp(i).sea(m,:))' - bensomax.*squeeze(ENSOyearind(31*12+m2-bensomax_ind(1)+1:12:end,:));   
            
        end
    else
        obsseaiceextent(i).sea =  obsseaiceextenttemp(i).sea;
    end
    
    % remove model
    if inputs.removeENSO
        for k = 1:size(surfacedata.highCl.ENSO,1)
            for m = 1:12
                for lag = 1:laglength                                                            
                    ENSOyearind = surfacedata.highCl.ENSO(k,:);             
                    if m == 1 || m == 2
                        m2 = 3;
                    else
                        m2 = m;
                    end
                    ensopredictors = [squeeze(ENSOyearind(m2-lag+1:12:end-4*12))',ones(length(squeeze(ENSOyearind(m2-lag+1:12:end-4*12))),1)];
                    benso(lag,:,:) = ensopredictors\squeeze(seaiceextenttemp(i).sea(k,m,:));

                    %                     ensopredictors_alldata = [squeeze(ENSOall(k,inputs.varmonth(m)-lag+1:12:end,:))',ones(size(alldata,3)-ext,1)];
                    %                     benso_alldata(lag,:,:) = ensopredictors_alldata\squeeze(surfacedataind(k,:,m,j,:));                        
                    %find max lag
                end
                [~,bensomax_ind] = max(abs(benso),[],1);                
                %[~,bensomax_ind_alldata] = max(abs(benso_alldata),[],1);   

                bensomax = benso(bensomax_ind(1),1);
                %bensomax_alldata(li) = benso_alldata(bensomax_ind_alldata(1,1,li),1,li);
                seaiceextent(i).sea(k,m,:) = squeeze(seaiceextenttemp(i).sea(k,m,:))' - bensomax.*squeeze(ENSOyearind(m2-bensomax_ind(1)+1:12:end-4*12));   
                restseaiceextent(i).sea(k,m,:) = squeeze(restseaiceextenttemp(i).sea(k,m,:))' - bensomax.*squeeze(ENSOyearind(26*12+m2-bensomax_ind(1)+1:12:end));   

            end
        end
    else
        seaiceextent(i).sea =  seaiceextenttemp(i).sea;
    end
    
    % now remove linear
    linearpredictors = [[1:26]',ones(26,1)];
    for k = 1:size(surfacedata.highCl.ENSO,1)
        for m  = 1:12
            blinear(i,k,m,:) = linearpredictors\squeeze(seaiceextent(i).sea(k,m,:));            
            seaiceextent(i).seadetrend(k,m,:) = squeeze(seaiceextent(i).sea(k,m,:)) - blinear(i,k,m,1).*linearpredictors(:,1);
            restseaiceextent(i).seadetrend(k,m,:) = squeeze(restseaiceextent(i).sea(k,m,:)) - blinear(i,k,m,1).*[27:30]';
        end
    end
    
    obslinearpredictors = [[1:31]',ones(31,1)];

    for m  = 1:12
        obsblinear(i,m,:) = obslinearpredictors\squeeze(obsseaiceextent(i).sea(m,:))';
        obsseaiceextent(i).seadetrend(m,:) = squeeze(obsseaiceextent(i).sea(m,:)) - obsblinear(i,m,1).*obslinearpredictors(:,1)';
        restobsseaiceextent(i).seadetrend(m,:) = squeeze(restobsseaiceextent(i).sea(m,:)) - obsblinear(i,m,1).*[32:37];        
    end
    
    
    % take mean of seasons
    for l = 1:3
        obsseaiceextent(i).seadetrendseasons(l,:) = nanmean(obsseaiceextent(i).seadetrend(seasons(l,:),:),1);
        restobsseaiceextent(i).seadetrendseasons(l,:) = nanmean(restobsseaiceextent(i).seadetrend(seasons(l,:),:),1);
        for k = 1:9
            seaiceextent(i).seadetrendseasons(k,l,:) = nanmean(seaiceextent(i).seadetrend(k,seasons(l,:),:),2);
            restseaiceextent(i).seadetrendseasons(k,l,:) = nanmean(restseaiceextent(i).seadetrend(k,seasons(l,:),:),2);
        end    
    end
    % now construct model using ozone
    
    % model
    
    
    
    for k = 1:size(surfacedata.highCl.ENSO,1)
        tozb(i,k,:) = linearpredictors\squeeze(modeltoz(k,inputs.tozmonth,1:26));
        tozdetrend = squeeze(modeltoz(k,inputs.tozmonth,1:26)) - tozb(i,k,1).*linearpredictors(:,1) - tozb(i,k,2);
        resttozdetrend = squeeze(modeltoz(k,inputs.tozmonth,27:30)) - tozb(i,k,1).*[27:30]' - tozb(i,k,2);
        ozonepredictors = [squeeze(tozdetrend),ones(26,1)];
        
        for m  = 1:12           
            btoz(i,k,m,:) = ozonepredictors\squeeze(seaiceextent(i).seadetrend(k,m,:));
            modelfit(i,k,m,:) = squeeze(tozdetrend).*btoz(i,k,m,1) + btoz(i,k,m,2) - median(tozdetrend.*btoz(i,k,m,1) + btoz(i,k,m,2));% - btoz(i,k,m,2);
            modelprediction(i,k,m,:) = resttozdetrend.*btoz(i,k,m,1) + btoz(i,k,m,2) - median(tozdetrend.*btoz(i,k,m,1) + btoz(i,k,m,2));% - btoz(i,k,m,2);
            surface(i,k,m,:) = squeeze(seaiceextent(i).seadetrend(k,m,:)) - median(squeeze(seaiceextent(i).seadetrend(k,m,:)));
            surfaceleft(i,k,m,:) = squeeze(restseaiceextent(i).seadetrend(k,m,:)) - median(squeeze(seaiceextent(i).seadetrend(k,m,:)));
            modelcorr(i,k,m) = corr(squeeze(seaiceextent(i).seadetrend(k,m,:)),tozdetrend);
        end
        for m  = 1:3            
            btoz_seasons(i,k,m,:) = ozonepredictors\squeeze(seaiceextent(i).seadetrendseasons(k,m,:));
            modelfit_seasons(i,k,m,:) = squeeze(tozdetrend).*btoz_seasons(i,k,m,1) + btoz_seasons(i,k,m,2) - median(tozdetrend.*btoz_seasons(i,k,m,1) + btoz_seasons(i,k,m,2));% - btoz(i,k,m,2);
            modelprediction_seasons(i,k,m,:) = resttozdetrend.*btoz_seasons(i,k,m,1) + btoz_seasons(i,k,m,2) - median(tozdetrend.*btoz_seasons(i,k,m,1) + btoz_seasons(i,k,m,2));% - btoz(i,k,m,2);
            surface_seasons(i,k,m,:) = squeeze(seaiceextent(i).seadetrendseasons(k,m,:)) - median(squeeze(seaiceextent(i).seadetrendseasons(k,m,:)));
            surfaceleft_seasons(i,k,m,:) = squeeze(restseaiceextent(i).seadetrendseasons(k,m,:)) - median(squeeze(seaiceextent(i).seadetrendseasons(k,m,:)));
            modelcorr_seasons(i,k,m,:) = corr(squeeze(seaiceextent(i).seadetrendseasons(k,m,:)),tozdetrend);
        end
    end
    
    obstozb(i,:) = obslinearpredictors\obstoz(1:31)';
    obstozdetrend = obstoz(1:31)' - obstozb(i,1).*obslinearpredictors(1,:) - obstozb(i,2);
    restobstozdetrend = obstoz(32:end)' - obstozb(i,1).*[32:37]' - obstozb(i,2);
    obsozonepredictors = [squeeze(obstozdetrend(1:31))',ones(31,1)];
    obsozonerest = restobstozdetrend;
    % observations
    for m  = 1:12
       
        btozobs(i,m,:) = obsozonepredictors\squeeze(obsseaiceextent(i).seadetrend(m,:))';
        obsmodelfit(i,m,:) = btozobs(i,m,1).*obsozonepredictors(:,1) + btozobs(i,m,2) - median(squeeze(btozobs(i,m,1).*obsozonepredictors(:,1) + btozobs(i,m,2)));
        obsmodelprediction(i,m,:) = btozobs(i,m,1).*obsozonerest + btozobs(i,m,2) - median(squeeze(btozobs(i,m,1).*obsozonepredictors(:,1) + btozobs(i,m,2)));
        obssurface(i,m,:) = squeeze(obsseaiceextent(i).seadetrend(m,:))' - median(squeeze(obsseaiceextent(i).seadetrend(m,:))');
        obssurfaceleft(i,m,:) = squeeze(restobsseaiceextent(i).seadetrend(m,:))' - median(squeeze(obsseaiceextent(i).seadetrend(m,:))');
        obscorr(i,m) = corr(obsseaiceextent(i).seadetrend(m,:)',squeeze(obstozdetrend(1:31))');
    end    
    for m  = 1:3        
        btozobs_seasons(i,m,:) = obsozonepredictors\squeeze(obsseaiceextent(i).seadetrendseasons(m,:))';
        obsmodelfit_seasons(i,m,:) = btozobs_seasons(i,m,1).*obsozonepredictors(:,1) + btozobs(i,m,2) - median(squeeze(btozobs(i,m,1).*obsozonepredictors(:,1) + btozobs(i,m,2)));
        obsmodelprediction_seasons(i,m,:) = btozobs_seasons(i,m,1).*obsozonerest + btozobs(i,m,2) - median(squeeze(btozobs(i,m,1).*obsozonepredictors(:,1) + btozobs(i,m,2)));
        obssurface_seasons(i,m,:) = squeeze(obsseaiceextent(i).seadetrendseasons(m,:))' - median(squeeze(obsseaiceextent(i).seadetrendseasons(m,:))');
        obssurfaceleft_seasons(i,m,:) = squeeze(restobsseaiceextent(i).seadetrendseasons(m,:))' - median(squeeze(obsseaiceextent(i).seadetrendseasons(m,:))');
        obscorr_seasons(i,m) = corr(obsseaiceextent(i).seadetrendseasons(m,:)',squeeze(obstozdetrend(1:31))');
    end    
    
end

%% plotting
cbrew = cbrewer('qual','Set1',10);
%modelnums = [NaN,NaN,3,3,7,7,8,8];
%modelnums = [NaN,NaN,1,1,2,2,3,3];
%modelnums = [NaN,NaN,4,4,5,5,6,6];
modelnums = [NaN,NaN,7,7,8,8,9,9];
createfig('largeportrait','on');
obsloctoplot = [1,2];
modloctoplot = [NaN,NaN,1,2,1,2,1,2];
obsmontoplot = [4,4];
%obsmontoplot = [1,2];
modmontoplot = [NaN,NaN,4,4,4,4,4,4];
%modmontoplot = [NaN,NaN,1,2,1,2,1,2];
if strcmp(inputs.var,'ICEFRAC')
    titles = {'Obs, OBCLEB, May','Obs, GN, July',['No. ',sprintf('%02d',modelnums(3)),', OBCLEB, May'],['No. ',sprintf('%02d',modelnums(4)),', GN, September'],...
        ['No. ',sprintf('%02d',modelnums(5)),', OBCLEB, May'],['No. ',sprintf('%02d',modelnums(6)),', GN, September'],...
        ['No. ',sprintf('%02d',modelnums(7)),', OBCLEB, May'],['No. ',sprintf('%02d',modelnums(8)),', GN, September']};
else
    titles = {['Obs, Northern Siberia, ',monthnames(obsmontoplot(1),0,'long')],['Obs, Southern Siberia, ', monthnames(obsmontoplot(2),0,'long')],['No. ',sprintf('%02d',modelnums(3)),', Northern Siberia, ', monthnames(modmontoplot(3),0,'long')],['No. ',sprintf('%02d',modelnums(4)),', Southern Siberia, ', monthnames(modmontoplot(4),0,'long')],...
        ['No. ',sprintf('%02d',modelnums(5)),', Northern Siberia, ', monthnames(modmontoplot(5),0,'long')],['No. ',sprintf('%02d',modelnums(6)),', Southern Siberia, ', monthnames(modmontoplot(6),0,'long')],...
        ['No. ',sprintf('%02d',modelnums(7)),', Northern Siberia, ', monthnames(modmontoplot(7),0,'long')],['No. ',sprintf('%02d',modelnums(8)),', Southern Siberia, ',monthnames(modmontoplot(8),0,'long')]};
end
fsize = 18;

seasons = 0;
if seasons 
    obstoplot = obssurface_seasons;
    obstoplotleft = obssurfaceleft_seasons;
    obsmodetoplot = obsmodelfit_seasons;
    obsmodtoplotleft = obsmodelprediction_seasons;
    modtoplot = surface_seasons;
    modtoplotleft = surfaceleft_seasons;
    modmodtoplot = modelfit_seasons;
    modmodtoplotleft = modelprediction_seasons;
    modc = modelcorr_seasons;
    obsc = obscorr_seasons;
else
    obstoplot = obssurface;
    obstoplotleft = obssurfaceleft;
    obsmodetoplot = obsmodelfit;
    obsmodtoplotleft = obsmodelprediction;
    modtoplot = surface;
    modtoplotleft = surfaceleft;
    modmodtoplot = modelfit;
    modmodtoplotleft = modelprediction;
    modc = modelcorr;
    obsc = obscorr;
end

for i  = 1:8
    sp(i) = subplot(4,2,i);
    box on
    if i < 3
        phobs(i) = plot(1:37,[squeeze(obstoplot(obsloctoplot(i),obsmontoplot(i),:));squeeze(obstoplotleft(obsloctoplot(i),obsmontoplot(i),:))],...
            'color',cbrew(1,:),'LineWidth',2,'Marker','o','MarkerFaceColor',cbrew(1,:),'MarkerEdgeColor',cbrew(1,:)*.8);
        hold on
        phobsmod(i) = plot(1:37,[squeeze(obsmodetoplot(obsloctoplot(i),obsmontoplot(i),:));squeeze(obsmodtoplotleft(obsloctoplot(i),obsmontoplot(i),:))],...
            'color',cbrew(2,:),'LineWidth',2,'Marker','o','MarkerFaceColor',cbrew(2,:),'MarkerEdgeColor',cbrew(2,:)*.8);
        plot([-1 39],[0,0],'--k');
        plot([31 31],[-10e5 10e5],'--k');
        if i == 1
            if strcmp(inputs.var,'ICEFRAC')
                ylim([-7e-9,7e9]);
            else
                ylim([-.1,.1]);
            end
                
        elseif i == 2
            if strcmp(inputs.var,'ICEFRAC')
                ylim([-5e5,5e5]);
            else
                ylim([-.08,.08]);
            end
        end
        xlim([-1 39]);
        th(i) = title(titles{i},'fontsize',fsize+2);
        thpos(i,:) = get(th(i),'position');
        %set(th(i),'position',[thpos(i,1),thpos(i,2)+thpos(i,2)./5,thpos(i,3)]);
        
        set(gca,'xtick',1:10:37,'xticklabel',1980:10:2030,'fontsize',fsize-2);
        if strcmp(inputs.var,'ICEFRAC')
            ylabel('Sea ice extent (km^2)','fontsize',fsize);
        else
            ylabel('Snow depth (m)','fontsize',fsize);
        end
        xlabel('Year','fontsize',fsize);
        annotation('textbox',get(sp(i),'position'),'String',['r = ', sprintf('%.2f',obsc(obsloctoplot(i),obsmontoplot(i)))],'FitBoxToText','on','LineStyle','none','Fontsize',fsize);
    else
        phobs(i) = plot(1:30,[squeeze(modtoplot(modloctoplot(i),modelnums(i),modmontoplot(i),:));squeeze(modtoplotleft(modloctoplot(i),modelnums(i),modmontoplot(i),:))],...
            'color',cbrew(1,:),'LineWidth',2,'Marker','o','MarkerFaceColor',cbrew(1,:),'MarkerEdgeColor',cbrew(1,:)*.8);
        hold on
        phobsmod(i) = plot(1:30,[squeeze(modmodtoplot(modloctoplot(i),modelnums(i),modmontoplot(i),:));squeeze(modmodtoplotleft(modloctoplot(i),modelnums(i),modmontoplot(i),:))],...
            'color',cbrew(2,:),'LineWidth',2,'Marker','o','MarkerFaceColor',cbrew(2,:),'MarkerEdgeColor',cbrew(2,:)*.8);
        plot([-1 39],[0,0],'--k');
        plot([26 26],[-10e6 10e6],'--k');
        if mod(i,2)
            if strcmp(inputs.var,'ICEFRAC')
                ylim([-12e5,12e5]);
            else
                ylim([-.1,.1]);
            end
        else
            if strcmp(inputs.var,'ICEFRAC')
                ylim([-8e5,8e5]);
            else
                ylim([-.05,.05]);
            end
        end
        xlim([-1 32]);
        th(i) = title(titles{i},'fontsize',fsize+2);
        thpos(i,:) = get(th(i),'position');
        %set(th(i),'position',[thpos(i,1),thpos(i,2)+thpos(i,2)./5,thpos(i,3)]);
        
        set(gca,'xtick',1:5:30,'xticklabel',1995:5:2030,'fontsize',fsize-2);
        
        if strcmp(inputs.var,'ICEFRAC')
            ylabel('Sea ice extent (km^2)','fontsize',fsize);
        else
            ylabel('Snow depth (m)','fontsize',fsize);
        end
        xlabel('Year','fontsize',fsize);        
        annotation('textbox',get(sp(i),'position'),'String',['r = ', sprintf('%.2f',modc(modloctoplot(i),modelnums(i),modmontoplot(i)))],'FitBoxToText','on','LineStyle','none','Fontsize',fsize);
        
    end
    
end
filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/';
filename = [filedir,'Nos.',num2str(modelnums(3)),'_',num2str(modelnums(5)),'_',num2str(modelnums(7)),'_NEW_',inputs.var,'_regressionprediction_from_',...
    monthnames(inputs.tozmonth,1,'long'),'_',inputs.obstouse,'_Arcticozoneextremes_over_',...
    num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
export_fig(filename,'-pdf');
    
    % now extract regression model period
    
end