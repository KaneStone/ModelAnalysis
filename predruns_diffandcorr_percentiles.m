function [dataVarMonthAve] = predruns_diffandcorr_percentiles(dataMonthArrange,dataMonthArrangeMean,toz_dataMonthArrange,varmonth,tozmonth,var,...
    lats,lons,thigh,lats2,removeENSO,Eachyear,ClLevel)

% This function plots all individual and ensemble mean correlations

fields = fieldnames(dataMonthArrange);

%[~,ENSAVE,~] = Read_in_netcdf('/Volumes/MyBook/work/data/predruns/TS/lowcl/ensave/TS_b.e11.BWTREFC2.f19_g16.ccmi34.LowCl.ENSAVE.cam.h0.nc');

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
             
    save('/Volumes/MyBook/work/data/predruns/output/NINO34/highcl_1995_2014','NINO34all','NINOmonth');
end

%Calculate correlations between ENSO and TCO


if ~exist(['/Volumes/MyBook/work/data/predruns/output/data/TS_ninoremoved_',num2str(thigh(1)),'-',num2str(thigh(2)),num2str(abs(lats2(1))),'-',...
    num2str(abs(lats2(2))),'.mat'],'file')

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

    save(['/Volumes/MyBook/work/data/predruns/output/data/TS_ninoremoved_',num2str(thigh(1)),'-',num2str(thigh(2)),num2str(abs(lats2(1))),'-',...
        num2str(abs(lats2(2)))],'dataVarMonthAve','dataMonthArrangeMean');
else
    load(['/Volumes/MyBook/work/data/predruns/output/data/TS_ninoremoved_',num2str(thigh(1)),'-',num2str(thigh(2)),num2str(abs(lats2(1))),'-',...
        num2str(abs(lats2(2)))]);
end



%% extracting upper and lower percentiles
clearvars varlower varupper

%for i = 1:size(Eachyear.lowindex,2)
    for j = 1:size(Eachyear.lowindex,3)
%         for k = 1:length(varmonth)
%             if varmonth(k) < 3
%                 varlower(:,k,j,:,:) = dataMonthArrange.highcl(Eachyear.lowindex(:,tozmonth,j),varmonth(k),j+1,:,:);
%                 varupper(:,k,j,:,:) = dataMonthArrange.highcl(Eachyear.highindex(:,tozmonth,j),varmonth(k),j+1,:,:);
%             else
%                 varlower(:,k,j,:,:) = dataMonthArrange.highcl(Eachyear.lowindex(:,tozmonth,j),varmonth(k),j,:,:);
%                 varupper(:,k,j,:,:) = dataMonthArrange.highcl(Eachyear.highindex(:,tozmonth,j),varmonth(k),j,:,:);
%             end
                 varlower(:,j,:,:) = dataVarMonthAve.highcl(Eachyear.lowindex(:,tozmonth,j),j,:,:);
                 varupper(:,j,:,:) = dataVarMonthAve.highcl(Eachyear.highindex(:,tozmonth,j),j,:,:);
                 tozcomp(:,j) = [squeeze(toz_dataMonthArrange.highcl(Eachyear.lowindex(:,tozmonth,j),tozmonth,j));...
                     squeeze(toz_dataMonthArrange.highcl(Eachyear.highindex(:,tozmonth,j),tozmonth,j))];                  
        %end
    end
%end
varcomp = permute(cat(1,varlower,varupper),[3,4,1,2]);
varcomp = varcomp(:,:,:);
tozcomp = tozcomp(:);

save(['/Volumes/MyBook/work/data/predruns/output/data/varcomp1995-2024_North'],'varcomp','tozcomp')

clearvars testlower testupper 
%%
clearvars difference p
testlower = permute(varlower,[3,4,2,1]);
testlower = nanmean(testlower,4);
testlower2 = nanmean(testlower,3);

testupper = permute(varupper,[3,4,2,1]);
testupper = nanmean(testupper,4);
testupper2 = nanmean(testupper,3);
difference(1,:,:) = testupper2 - testlower2;
difference = permute(difference,[1,3,2]);
for i = 1:size(testlower,1)
    for j = 1:size(testlower,2)
        p(1,i,j) = ttest2(testlower(i,j,:),testupper(i,j,:));
    end
end
p = permute(p,[1,3,2]);
p (p == 1) = -1;
p (p == 0) = 1;
p (p == -1) = 0;
%titles = {[num2str(thigh(1)),'-',num2str(thigh(2)),' (high chlorine)']};
    
if lats2 < 0
    clims = [-1 1];    
    xlim = [-90 0];
else
     clims = [-6 6];    
     xlim = [0 90];
end

mtitle = {[monthnames(tozmonth,0,0),' TCO - ',monthnames(varmonth,1,1),' ',ClLevel,' surface temperature differences']};

%mtitle = ['Esemble mean correlations of ',num2str(abs(lats2(1))),'-',num2str(abs(lats2(2))),hemext,' toz and ',var];

subplotmaps(difference,lons,lats,{'div','RdBu'},1,p,16,mtitle,'Longitude','latitude','Temperature difference (K)','on',...
    clims,22,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],{''},1,[0 360],xlim,0,'none',1,'Miller Cylindrical');
    
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/','correlations/maps/difference_',monthnames(tozmonth,0,0),'toz','_',var,'_detrend',...
    num2str(thigh(1)),'-',num2str(thigh(2)),'_',num2str(abs(lats2(1))),'-',...
    num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,1),'_',ClLevel];

export_fig(filename,'-png');
        
%% plot correlations

%% plotting

if removeENSO
    %ENSOext = 'rmENSO';
    ENSOext = 'rmENSO_lagged';
else
    ENSOext = [];
end

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

%cbrew = cbrewer('div','RdBu',16);    

for i = 1:size(varcomp,1)
    for j = 1:size(varcomp,2)
        [r(1,j,i),p(1,j,i)] = corr(squeeze(varcomp(i,j,:)),tozcomp);
    end
end

p (p <= .05) = 0;
p (p > .05) = 1;

titles = {[num2str(thigh(1)),'-',num2str(thigh(2)),' (high chlorine)']};
    
mtitle = '';
%mtitle = ['Esemble mean correlations of ',num2str(abs(lats2(1))),'-',num2str(abs(lats2(2))),hemext,' toz and ',var];

subplotmaps(r(1,:,:),lons,lats,{'div','RdBu'},1,p,16,titles,'Longitude','latitude','Correlation','on',...
    clims,22,[-90:10:90],[-90:10:90],[lons(1:10:end)],[lons(1:10:end)],mtitle,1,[0 360],xlims,0,'none',1,'Miller Cylindrical');

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/','correlations/maps/compositeperc_',monthnames(tozmonth,0,0),'toz','_',var,'_detrend',...
        num2str(thigh(1)),'-',num2str(thigh(2)),'_',num2str(abs(lats2(1))),'-',...
        num2str(abs(lats2(2))),'S_Tperiod-',monthnames(varmonth,1,1),'_',ENSOext,'highclonly'];

export_fig(filename,'-png');

end