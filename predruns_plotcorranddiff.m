function [] = predruns_plotcorranddiff(correlations,differences,inputs,lons,lats)

%% plot individual members as 9 or 10 member subplot
if inputs.removeENSO
    ensoext = 'removeENSO';
else
    ensoext = 'noremoveENSO';
end

if inputs.lats(1) < 0
        ylim = [-90 0];
    else
        ylim = [0 90];
end

if inputs.plotind
        
    
    mtitle = {''};
    
    for i = 1:size(correlations.individual,1)
        titles{i} = ['Correlations',' - No. ',num2str(i)]; 
        titlessave{i} = ['Correlations','_No.',num2str(i)]; 
        titles2{i} = ['Percentile correlations',' - No. ',num2str(i)]; 
        titlessave2{i} = ['PercentileCorrelations','_No.',num2str(i)]; 
    end
        
    %% plotting correlations individually
    clims = [-1 1];
    
    for i = 1:size(correlations.individual)
        if ~exist(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/Ind_',titlessave{i},'_',monthnames(inputs.varmonthtomean,1,1),'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'.png'])
            subplotmaps(correlations.individual(i,:,:),lons,lats,{'div','RdBu'},1,[],16,titles(i),'Longitude','Latitude','Correlation','on',...
                clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
            filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/Ind_',titlessave{i},'_',monthnames(inputs.varmonthtomean,1,1),'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
            export_fig(filename,'-png');
        end
        if ~exist(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/Ind_pct_',titlessave2{i},'_',monthnames(inputs.varmonthtomean,1,1),'.png'])
            subplotmaps(correlations.individual_pct(i,:,:),lons,lats,{'div','RdBu'},1,[],16,titles2(i),'Longitude','Latitude','Correlation','on',...
                clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
            filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/Ind_pct_',titlessave2{i},'_',monthnames(inputs.varmonthtomean,1,1),'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
            export_fig(filename,'-png');
            
        end
                
    end
        
    
    %% plotting differences individually
    
    for i = 1:size(correlations.individual,1)
        titles{i} = ['Differences',' - No. ',num2str(i)]; 
        titlessave{i} = ['Differences','_No.',num2str(i)]; 
    end
    
    clims = [-5 5];
    
    for i = 1:size(correlations.individual)
        if ~exist(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/differences/maps/Ind_',titlessave{i},'_',monthnames(inputs.varmonthtomean,1,1),'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'.png'])
            subplotmaps(differences.individual(i,:,:),lons,lats,{'div','RdBu'},1,[],16,titles(i),'Longitude','Latitude','Temperature difference','on',...
                clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
            filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/differences/maps/Ind_',titlessave{i},'_',monthnames(inputs.varmonthtomean,1,1),'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
            export_fig(filename,'-png');
        end
        
        
    end    
   
    
end

%% plot individual months correlations
if inputs.plotindmonths
    mtitle = {''};    
    
        
    %% plotting correlations individually
    clims = [-1 1];
    for j = 1:size(correlations.indmonths.individual,2)
        for k = 1:size(correlations.indmonths.individual,1)
            titlesindmonths{j,k} = ['Correlations',' - No. ',num2str(k),' ',monthnames(inputs.varmonth(j),0,0)]; 
            titlesindmonthssave{j,k} = ['Correlations','_No.',num2str(k),'_',monthnames(inputs.varmonth(j),0,0)]; 
        end
        for i = 1:size(correlations.indmonths.individual,1)
            toplot = squeeze(correlations.indmonths.individual(i,:,:,:));
            toplot = toplot(j,:,:);
            if ~exist(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/Ind_months','_',titlesindmonthssave{j,i},'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'.png'],'file')            
                subplotmaps(toplot,lons,lats,{'div','RdBu'},1,[],16,titlesindmonths(j,i),'Longitude','Latitude','Correlation','on',...
                    clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
                filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/Ind_months','_',titlesindmonthssave{j,i},'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
                export_fig(filename,'-png');
            end            
        end
    end       
        
    %%
    clims = [-5 5];
    for j = 1:size(differences.indmonths.individual,4)
        for k = 1:size(differences.indmonths.individual,1)
            titlesindmonths{j,k} = ['Differences',' - No. ',num2str(k),' ',monthnames(inputs.varmonth(j),0,0)]; 
            titlesindmonthssave{j,k} = ['Differences','_No.',num2str(k),'_',monthnames(inputs.varmonth(j),0,0)]; 
        end
        for i = 1:size(differences.indmonths.individual,1)
            toplot = differences.indmonths.individual(i,:,:,j);
            %toplot = toplot(j,:,:);
            if ~exist(['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/differences/maps/Ind_months','_',titlesindmonthssave{j,i},'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2))),'.png'],'file')            
                subplotmaps(toplot,lons,lats,{'div','RdBu'},1,[],16,titlesindmonths(j,i),'Longitude','Latitude','Correlation','on',...
                    clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');                
                filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/differences/maps/Ind_months','_',titlesindmonthssave{j,i},'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
                export_fig(filename,'-png');
            end
        end
    end       
          
    
end

%% plot ensemble mean correlations and (minus) std

if inputs.plotcorr
    
    %read in zonal wind
    files = dir(['/Volumes/ExternalOne/work/data/predruns/U/highCl/500hPa/','*.mat']);
    for i = 1:length(files)
        Uwind(i) = load(['/Volumes/ExternalOne/work/data/predruns/U/highCl/500hPa/',files(i).name]);
        for j = 1:length(inputs.varmonthtomean)
            Uwindmonthave(:,:,i,j) = nanmean(squeeze(Uwind(i).standardout(:,:,:,inputs.varmonthtomean(j):12:end)),3);                
        end
    end
    
    Uwindstd = squeeze(nanstd(Uwindmonthave(:,:,:),0,3))./max(abs(Uwindmonthave(:,:,:)),[],3);
    Uwindmean = squeeze(nanmean(Uwindmonthave(:,:,:),3));
    
    % constructing sig
    p = correlations.sig.composite_all_pct;
    %p = correlations.sig.corrmean;
    p (p > .05) = 1;
    p (p <= .05) = 0;
    p (p == 1) = -1;
    
    p2 = correlations.sig.individualcorr_pct;
    p2 (p2 > .05) = 1;
    p2 (p2 <= .05) = 0;
    p2 (p2 == 1) = -1;
    
    mtitle = {''};    
    toplot = cat(1,correlations.corrmean,correlations.indstd./abs(correlations.corrmean),...
        abs(correlations.corrmean) - correlations.indstd);%abs(correlations.composite_all_pct) - correlations.indstd
    toplot2 = cat(1,correlations.corrmean_pct,correlations.indstd_pct./abs(correlations.corrmean_pct),...
        abs(correlations.corrmean_pct) - correlations.indstd_pct);%abs(correlations.composite_all_pct) - correlations.indstd
    toplotp = cat(1,p,ones(size(p)),ones(size(p)));
    toplotp2 = cat(1,p2,ones(size(p2)),ones(size(p2)));
    
    clims = [-.75 .75];
    % all data
    titles = {'Ensemble mean correlation'}; 
    subplotmaps(toplot(1,:,:),lons,lats,{'div','RdBu'},1,toplotp(1,:,:),16,titles,'Longitude','Latitude','Correlation','on',...
        clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/EnsMean_correlations2_',...
        monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',ensoext,'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
    export_fig(filename,'-png');        
    %%
    clearvars ax2 ax
    titles = {'Standard deviation/Ensemble mean'}; 
    [fig,sh] = subplotmaps(toplot(2,:,:),lons,lats,{'seq','Blues'},0,toplotp(2,:,:),16,titles,'Longitude','Latitude','','on',...
        [0 5],11,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
%     ax = gca;
%     hold on
%     lontoplot = lons-180;
%     pos = get(sh,'position');
%     ax2 = axes;            
%     m_contour(lontoplot,lats,squeeze(Uwindstd)','LineWidth',2);
%     linkaxes([ax,ax2]);
%     ax2.Visible = 'off';
%     ax2.XTick = [];
%     ax2.YTick = [];
%     set(ax2,'position',pos);
%     set(ax2,'color','none');
%     m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[ylim(1) ylim(2)]);
%     %m_coast('color','k','LineWidth',1);
%     m_grid('ytick',0:15:90,'xtick',-180:60:180,'XaxisLocation','bottom','fontsize',20,'color','none','backcolor','none');
%     cnew = cbrewer('seq','Reds',9);
%     colormap(ax2,cnew)
%     ch = colorbar;
%     set(ch,'orientation','horizontal','position',[.105 .3 .725 .02],'fontsize',22);
           
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/EnsMean_sigtonoise_stddivmean2_withuwindat500hPa_std',...
        monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
    %export_fig(filename,'-png');
    print(filename,'-depsc');
    %%
    titles = {'Ensemble mean - standard deviation'}; 
    subplotmaps(toplot(3,:,:),lons,lats,{'div','RdBu'},1,toplotp(2,:,:),16,titles,'Longitude','Latitude','','on',...
        [-.3 .3],12,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/EnsMean_sigtonoise_meanminstd2_',...
        monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),ensoext,'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
    export_fig(filename,'-png');
   
    % Percentiles
    titles = {'Ensemble mean upper and lower percentile correlation'}; 
    subplotmaps(toplot2(1,:,:),lons,lats,{'div','RdBu'},1,toplotp2,16,titles,'Longitude','Latitude','Correlation','on',...
        clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/EnsMean_correlations2_pct_',...
        monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),ensoext,'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
    export_fig(filename,'-png');        
    
    titles = {'percentile correlation - Standard deviation/Ensemble mean'}; 
    subplotmaps(toplot2(2,:,:),lons,lats,{'seq','Blues'},0,toplotp2(2,:,:),16,titles,'Longitude','Latitude','','on',...
        [0 5],11,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/EnsMean_sigtonoise_stddivmean2_pct',...
        monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),ensoext,'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
    export_fig(filename,'-png');
    
    titles = {'percentile correlation - Ensemble mean - standard deviation'}; 
    subplotmaps(toplot2(3,:,:),lons,lats,{'div','RdBu'},1,toplotp2(2,:,:),16,titles,'Longitude','Latitude','','on',...
        [-.3 .3],12,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/correlations/maps/EnsMean_sigtonoise_meanminstd2_pct',...
        monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),ensoext,'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
    export_fig(filename,'-png');
    
end

%% plot ensemble mean differences and (minus) std
if inputs.plotdiff
    % constructing sig
    p = differences.ttest.composite;
    p (p == 1) = -2;
    p (p == 0) = -1;
    p (p == -2) = 0;
    
    mtitle = {''};
    
    toplot = cat(1,differences.composite,differences.indstd./abs(differences.composite));%abs(differences.composite) - differences.indstd
    toplotp = cat(1,p,ones(size(p)));
    
    clims = [-5 5];
    
    titles = {'Composite','Composite/standard deviation'}; 
    subplotmaps(toplot,lons,lats,{'div','RdBu'},1,toplotp,16,titles,'Longitude','Latitude','Temperature difference','on',...
            clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
        
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/differences/maps/Composite_sigtonoise_meanminstd_',...
        monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),ensoext,'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
        export_fig(filename,'-png');
        
    p = differences.ttest.eachyear;
    p (p == 1) = -2;
    p (p == 0) = -1;
    p (p == -2) = 0;
    
    mtitle = {''};
    
    toplot = cat(1,differences.eachyear,differences.indstd./abs(differences.eachyear));%abs(differences.eachyear) - differences.indstd
    toplotp = cat(1,p,ones(size(p)));
    
    clims = [-5 5];
    
    titles = {'Ensemble-Each year','Ensemble-Each year/standard deviation'}; 
    subplotmaps(toplot,lons,lats,{'div','RdBu'},1,toplotp,16,titles,'Longitude','Latitude','Temperature difference','on',...
            clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/differences/maps/EnsMean_Eachyear_sigtonoise_meanminstd_',...
        monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),ensoext,'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
        export_fig(filename,'-png');
                
    p = differences.ttest.ensmean;
    p (p == 1) = -2;
    p (p == 0) = -1;
    p (p == -2) = 0;
    
    mtitle = {''};
    
    toplot = cat(1,differences.ensmean,differences.indstd./abs(differences.ensmean));%abs(differences.eachyear) - differences.indstd
    toplotp = cat(1,p,ones(size(p)));
    
    clims = [-5 5];
    
    titles = {'Ensemble mean','Ensemble mean/standard deviation'}; 
    subplotmaps(toplot,lons,lats,{'div','RdBu'},1,toplotp,16,titles,'Longitude','Latitude','Temperature difference','on',...
            clims,22,[lons(1:24:end)]-180,[lons(1:24:end)]-180,[ylim(1):15:ylim(2)],[ylim(1):15:ylim(2)],mtitle,1,[0 360],ylim,0,'none',1,'Miller Cylindrical');
    
    filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/differences/maps/EnsMean_sigtonoise_meanminstd_',...
        monthnames(inputs.varmonthtomean,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),ensoext,'_',...
                num2str(abs(inputs.lats(1))),'-',num2str(abs(inputs.lats(2)))];
        export_fig(filename,'-png');    
        
end

end
