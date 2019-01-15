function [] = predruns_enscorr(dataMonthArrange,toz_dataMonthArrange,varmonth,tozmonth,var,...
    lats,lons,thigh,lats2,removeENSO,ClLevel)

% This function plots all individual and ensemble mean correlations

fields = fieldnames(dataMonthArrange);

%[~,ENSAVE,~] = Read_in_netcdf('/Volumes/ExternalOne/work/data/predruns/TS/lowcl/ensave/TS_b.e11.BWTREFC2.f19_g16.ccmi34.LowCl.ENSAVE.cam.h0.nc');

%% calculate ENSO

latlimits = [-5 5];
lonlimits = [190 240];

latindex = lats >= latlimits(1) & lats <= latlimits(2);
lonindex = lons >= lonlimits(1) & lons <= lonlimits(2);

for i = 1:length(fields)

    
    for j = 1:12
        NINO_mn(:,j,:,:) = squeeze(nanmean(dataMonthArrange.(fields{i})(:,j,:,latindex,lonindex),5));
        NINO_mn2(:,j,:) = squeeze(nanmean(NINO_mn(:,j,:,:),4));
        NINOmonth.(fields{i})(:,j,:) = (squeeze(NINO_mn2(:,j,:)) - nanmean(squeeze(NINO_mn2(:,j,:)),2))./std(squeeze(NINO_mn2(:,j,:)),1,2);
    end   

    NINO34all.(fields{i}) = NINOmonth.(fields{i})(:,:);    
             
    save('/Volumes/ExternalOne/work/data/predruns/output/NINO34/highcl_1995_2014','NINO34all','NINOmonth');
end

%Calculate correlations between ENSO and TCO



%% take ensemble averages
varmonth2 = varmonth;
isless = varmonth2 <= 2;
varmonth2 = (varmonth2 + isless*12);

laglength = 3;
for i = 1:length(fields)
    ensave.(fields{i}).(var) = squeeze(nanmean(dataMonthArrange.(fields{i}),1));
    ensave.(fields{i}).toz = squeeze(nanmean(toz_dataMonthArrange.(fields{i}),1));
                                    
    %varmonths
    datamonthall.(fields{i}) = permute(dataMonthArrange.(fields{i}),[1,4,5,2,3]);
    datamonthall.(fields{i}) = datamonthall.(fields{i})(:,:,:,:);
           
    % Here I am regressing with the of the seasonal average
    for j = 1:length(varmonth)
        together.(fields{i})(:,:,:,:,j) = squeeze(datamonthall.(fields{i})(:,:,:,varmonth2(j):12:end-varmonth(1)+varmonth(end)));
    end
    temp.(fields{i}) = nanmean(together.(fields{i}),5);
        
    if removeENSO        
        % finding maximum regression lag time
        for j = 1:size(datamonthall.(fields{i}),1) % members
            for k = 1:length(varmonth2)
                for l = 1:size(datamonthall.(fields{i}),2) % latitudes
                    for m = 1:size(datamonthall.(fields{i}),3) % longitudes
                        for lag = 1:laglength % lag
                            regressors = [ones(length(squeeze(NINO34all.(fields{i})(j,varmonth2(k)-lag+1:12:end-lag))),1),...
                                squeeze(NINO34all.(fields{i})(j,varmonth2(k)-lag+1:12:end-lag))'];
                            [b(j,k,l,m,lag,:)] = regress(squeeze(datamonthall.(fields{i})(j,l,m,varmonth2(k):12:end-lag)),...
                                regressors);                        
                            % finding largest lag correlation
                        end
                        [~,llc.(fields{i})(j,k,l,m)] = max(abs(squeeze(b(j,k,l,m,:,2))));
                        blag(j,k,l,m,:) = b(j,k,l,m,llc.(fields{i})(j,k,l,m),:);
                        dataVarMonth.(fields{i})(j,:,k,l,m) = squeeze(datamonthall.(fields{i})(j,l,m,varmonth2(k):12:end-laglength)) - ...
                            squeeze(b(j,k,l,m,llc.(fields{i})(j,k,l,m),2))*...
                            squeeze(NINO34all.(fields{i})(j,varmonth2(k)-llc.(fields{i})(j,k,l,m)+1:12:end-laglength))';                        
                    end
                end
            end        
        end
        dataVarMonthAve.(fields{i}) = squeeze(nanmean(dataVarMonth.(fields{i}),3));
        dataVarMonthAve_ensmean.(fields{i}) = squeeze(nanmean(dataVarMonthAve.(fields{i}),1));        
        
%         for j = 1:size(temp.(fields{i}),1)
%             for k = 1:size(temp.(fields{i}),3)
%                 for l = 1:size(temp.(fields{i}),4)
%                     [b(j,k,l,:)] = regress(squeeze(temp.(fields{i})(j,:,k,l))',[ones(length(NINO34.(fields{i}))-1,1),NINO34.(fields{i})(1:end-1,j)]);                        
%                     dataVarMonthAve.(fields{i})(j,:,k,l) = squeeze(temp.(fields{i})(j,:,k,l))' - squeeze(b(j,k,l,2))*NINO34.(fields{i})(1:end-1,j);            
%                 end
%             end
%         end        
%         dataVarMonthAve_ensmean.(fields{i}) = squeeze(nanmean(dataVarMonthAve.(fields{i}),1));        
    else
        dataVarMonthAve_ensmean.(fields{i}) = squeeze(nanmean(temp.(fields{i}),1));        
        dataVarMonthAve_ensmean.(fields{i}) = permute(dataVarMonthAve_ensmean.(fields{i}),[3,1,2]);
        dataVarMonthAve.(fields{i}) = permute(temp.(fields{i}),[1,4,2,3]);
    end        

    %composites
    for k = 1:size(dataVarMonthAve.(fields{i}),1)
        dataVarMonthAve2.(fields{i})(k).t = squeeze(dataVarMonthAve.(fields{i})(k,:,:,:));
        toz_dataVarMonthAve2.(fields{i})(k).t = squeeze(toz_dataMonthArrange.(fields{i})(k,:,:));
    end
    composite.(fields{i}).(var) = cat(1,dataVarMonthAve2.(fields{i}).t);
    composite.(fields{i}).toz = cat(2,toz_dataVarMonthAve2.(fields{i}).t)';
%     szvar = size(dataVarMonthAve.(fields{i}));
%     sztoz = size(toz_dataMonthArrange.(fields{i}));
%     composite.(fields{i}).(var) = reshape(dataVarMonthAve.(fields{i}),[szvar(2)*szvar(1),szvar([3,4])]); % I think this is right
%     composite.(fields{i}).toz = reshape(permute(toz_dataMonthArrange.(fields{i}),[1,3,2]),...
%         [sztoz(3)*sztoz(1),sztoz(2)]); % I think this is right    
end

%% taking varmonth averages
for j = 1:length(fields)    
    for i = 1:length(varmonth)
        if varmonth(i) <= 2
            ensave.(fields{j}).([(var),'_montemp'])(i,:,:,:) = ensave.(fields{j}).(var)(varmonth(i),2:end,:,:);                    
        else
            ensave.(fields{j}).([(var),'_montemp'])(i,:,:,:) = ensave.(fields{j}).(var)(varmonth(i),1:end-1,:,:);                    
        end     
        if i == length(varmonth)
            ensave.(fields{j}).([(var),'_monave']) = squeeze(nanmean(ensave.(fields{j}).([(var),'_montemp']),1));                        
        end
    end
    ensave.(fields{j}).(['toz','_monave']) = ensave.(fields{j}).toz(tozmonth,:);
    composite.(fields{j}).(['toz','_monave']) = composite.(fields{j}).toz(:,tozmonth);
end

%% taking correlations
for i = 1:length(fields)
    for j = 1:size(ensave.(fields{i}).([(var),'_monave']),2)
        for k = 1:size(ensave.(fields{i}).([(var),'_monave']),3)
                             
            [r.ensave(i,k,j),p.ensave(i,k,j)] = corr(dataVarMonthAve_ensmean.(fields{i})(:,j,k),ensave.(fields{i}).(['toz','_monave'])');
            [r.composite(i,k,j),p.composite(i,k,j)] = corr(composite.(fields{i}).(var)(:,j,k),composite.(fields{i}).toz(:,tozmonth));
        end
    end
end

maxminloclats = [-50 0];

maxminlatind = lats >= maxminloclats(1) & lats <= maxminloclats(2);

maxminlats = lats(maxminlatind);

%% taking correlations of individual members
rollwindow = 10;

if removeENSO
    %ENSOext = 'rmENSO';
    ENSOext = 'rmENSO_lagged';
else
    ENSOext = [];
end
calcroll = 1;
if calcroll
for i = 1:length(fields)
    for l = 1:size(dataVarMonthAve.(fields{i}),1)
        for j = 1:size(dataVarMonthAve.(fields{i}),3)
            for k = 1:size(dataVarMonthAve.(fields{i}),4)                             
                [r.ind.(fields{i})(l,k,j),p.ind.(fields{i})(l,k,j)] = corr(squeeze(dataVarMonthAve.(fields{i})(l,:,j,k))',squeeze(toz_dataMonthArrange.(fields{i})(l,tozmonth,:)));                
                [bpred.ind.(fields{i})(l,k,j,:),bpred_int.ind.(fields{i})(l,k,j,:,:)] = ...
                    regress(squeeze(dataVarMonthAve.(fields{i})(l,:,j,k))',...
                    [ones(length(squeeze(toz_dataMonthArrange.(fields{i})(l,tozmonth,:))),1),squeeze(toz_dataMonthArrange.(fields{i})(l,tozmonth,:))]);                
            end
        end
        [rollcorr.ind.(fields{i})(l,:,:,:),rolltrend.ind.(fields{i})(l,:,:,:)] = rollingCorrelations(permute(squeeze(dataVarMonthAve.(fields{i})(l,:,:,:)),[3,2,1]),...
            squeeze(toz_dataMonthArrange.(fields{i})(l,tozmonth,:)),rollwindow);
    end
         
    p.ind.(fields{i}) (p.ind.(fields{i}) <= .05) = 0;
    p.ind.(fields{i}) (p.ind.(fields{i}) > .05) = 1;
end
    if exist(['/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefs',monthnames(tozmonth,0,0),'toz','_TS_detrend',...
        num2str(thigh(1)),'-',num2str(thigh(2)),'_',num2str(abs(lats2(1))),'-',...
        num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,1),'_',ENSOext,'_','roll',num2str(rollwindow),'_',ClLevel],'file')
    else
        save(['/Volumes/ExternalOne/work/data/predruns/output/regression/regcoefs',monthnames(tozmonth,0,0),'toz','_TS_detrend',...
            num2str(thigh(1)),'-',num2str(thigh(2)),'_',num2str(abs(lats2(1))),'-',...
            num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,1),'_',ENSOext,'_','roll',num2str(rollwindow),'_',ClLevel],'r','p','bpred','bpred_int','lats','lons','rollcorr','rolltrend');
    end
    %calculate standard deviation of ensemble members
    rstd(i,:,:) = std(r.ind.(fields{i}),0,1);
end


%% taking differences of upper and lower percentiles in each member


%% plotting
if lats2(2) < 0
    hemext = 'S';
    xlims = [-90 0];
    clims = [-.75 .75];
    clims2 = [-.5 .5];
else
    hemext = 'N';
    xlims = [0 90];    
    clims = [-.75 .75];
    clims2 = [-.5 .5];
end

cbrew = cbrewer('div','RdBu',16);    
  

%% plotting ensembles
plotens = 1;
if plotens

    p.ensave (p.ensave <= .05) = 0;
    p.ensave (p.ensave > .05) = 1;
    p.composite (p.composite <= .05) = 0;
    p.composite (p.composite > .05) = 1;

    % contourtitle = {[monthnames(varmonth,1,shortnames),'{ }',var,'{ }',...
    %     monthnames(tozmonth,0,0),'{ }','NINO34',' correlations ',...
    %     num2str(abs(lats(1))),'-',num2str(abs(lats(2))),'S, ',num2str(timeperiodhigh(1)),'-',num2str(timeperiodhigh(2))]};       

%     titles = {[num2str(tlow(1)),'-',num2str(tlow(2)),' (low chlorine)'],...
%         [num2str(thigh(1)),'-',num2str(thigh(2)),' (high chlorine)'],...
%         [num2str(thigh(1)),'-',num2str(thigh(2)),' (high chlorine - low GHG)']};
    if strcmp(ClLevel,'highCl')        
        %titles = {[num2str(thigh(1)),'-',num2str(thigh(2)),' (high chlorine)']};
        titles = {[num2str(thigh(1)),'-',num2str(thigh(2)),' (high chlorine)'],'absolute ensave - std'};
    elseif strcmp(ClLevel,'lowCl')
         titles = {[num2str(thigh(1)),'-',num2str(thigh(2)),' (low chlorine)']};
    elseif strcmp(ClLevel,'lowGHG')
         titles = {[num2str(thigh(1)),'-',num2str(thigh(2)),' (low GHGs - high chlorine)']};
    end
    
    mtitle = '';
    %mtitle = ['Esemble mean correlations of ',num2str(abs(lats2(1))),'-',num2str(abs(lats2(2))),hemext,' toz and ',var];

%     subplotmaps(r.ensave(1,:,:),lons,lats,{'div','RdBu'},1,p.ensave(1,:,:),16,titles,'Longitude','latitude','Correlation','on',...
%         clims,18,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],mtitle,1,[0 360],xlims,0,'none',1,'Miller Cylindrical');
% 
%     filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/','correlations/maps/Ensmean2_',monthnames(tozmonth,0,0),'toz','_',var,'_detrend',...
%             num2str(thigh(1)),'-',num2str(thigh(2)),'_',num2str(abs(lats2(1))),'-',...
%             num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,1),'_',ENSOext,'_',ClLevel];
% 
%     export_fig(filename,'-png');
    %export_fig(filename,'-pdf');

    %mtitle = ['Composite mean correlations of ',num2str(abs(lats2(1))),'-',num2str(abs(lats2(2))),hemext,' toz and ',var];
    mtitle = '';
    
    comptoplot = cat(1,r.composite,abs(r.composite)-rstd);

    subplotmaps(comptoplot,lons,lats,{'div','RdBu'},1,cat(1,p.composite(1,:,:),ones(size(p.composite(1,:,:)))),16,titles,'Longitude','Latitude','Correlation','on',...
        clims2,22,[-90:10:20],[-90:10:20],[lons(1:10:end)],[lons(1:10:end)],mtitle,1,[0 360],xlims,0,'none',1,'Miller Cylindrical');

    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/','correlations/maps/composite2_',monthnames(tozmonth,0,0),'toz','_',var,'_detrend',...
            num2str(thigh(1)),'-',num2str(thigh(2)),'_',num2str(abs(lats2(1))),'-',...
            num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,1),'_',ENSOext,'_',ClLevel];

    export_fig(filename,'-png');
    %export_fig(filename,'-pdf');
end

%% plotting individuals
x = lons;
x(size(x)/2+1:end) = x(size(x)/2+1:end) - 360;
x = [x(size(x)/2+1:end);x(1:size(x)/2)];


plotind = 1;

    if plotind
        for i = 1:size(r.ind.highcl)

        subplotmaps(r.ind.highcl(i,:,:),lons,lats,{'div','RdBu'},1,p.ind.highcl(i,:,:),16,{['No.',num2str(i)]},'Longitude','latitude','Correlation','on',...
            [-.75 .75],22,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],'',1,[0 360],xlims,0,'none',1,'Miller Cylindrical');

        %rtoplot = circshift((abs(squeeze(r.ind.highcl(i,:,:)))'),size((abs(squeeze(r.ind.highcl(i,:,:)))'),2)/2,2);
        %m_contour(x,lats,rtoplot,[.6 .6],'color',[77,175,74]/255,'Linewidth',3,'LineStyle','-');

        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/','correlations/maps/Ind_','No.',sprintf('%02d',i),monthnames(tozmonth,0,0),'toz','_',var,'_detrend',...
            num2str(thigh(1)),'-',num2str(thigh(2)),'_',num2str(abs(lats2(1))),'-',...
            num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,1),'_',ENSOext,'_',ClLevel];
        
        export_fig(filename,'-png');
        end
    end
end
