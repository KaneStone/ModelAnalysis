function [GSSmonthtoplot_areamean,GSSmonthtoplot_areamean_pct,GSSmonthtoplot_areamean_pctsame,...
    forerrorextract,forerrorextractpct,lattoplot,lontoplot,lattoplot_pct,lontoplot_pct] = ...
    predruns_calculateHSStoplot(same,difftoplot,corrtoplot,latitude,longitude,differences,...
    modelprediction_ind_months,actualdata_ind_months,medmonths,modelprediction_ind_pct_months2,...
    actualdata_ind_pct_months2,medmonths_pct2,GSS,areas_lons,areas_lats,r_months,allsig,pctsig)

for mon = 1:5


if same
    sameext = 'same';
else
    sameext = '';
end
if corrtoplot
    corrtoplotext = 'corr';
elseif difftoplot
    corrtoplotext = 'diff';
else
    corrtoplotext = 'GSS';
end

r_monthsmean = squeeze(nanmean(r_months,1));
for i = 1:size(areas_lons,1)
    lats = latitude > areas_lats(i,1) & latitude < areas_lats(i,2);
    lons = longitude > areas_lons(i,1) & longitude < areas_lons(i,2);
    latextract = latitude(lats);
    lonextract = longitude(lons);
    [latmesh,lonmesh] = meshgrid(latextract,lonextract);
    
    latmesh = latmesh';
    lonmesh = lonmesh';
    for k = 1:size(differences,1)           
        
        if corrtoplot
            if same
                mult = squeeze(nanmean(r_monthsmean(:,mon,lats,lons),1));        
                mult_pct = squeeze(nanmean(r_monthsmean(:,mon,lats,lons),1));        
            else
                mult = squeeze(nanmean(r_months(k,:,mon,lats,lons),2));        
                mult_pct = squeeze(nanmean(r_months(k,:,mon,lats,lons),2));        
            end
        elseif difftoplot            
            diff = permute(differences,[1,4,3,2]);           
            if same                
                mult = squeeze(nanmean(diff(:,mon,lats,lons),1));        
                mult_pct = squeeze(nanmean(diff(:,mon,lats,lons),1));        
            else
                mult = squeeze(diff(k,mon,lats,lons));        
                mult_pct = squeeze(diff(k,mon,lats,lons));        
            end
        else
            if ~same
                mult = squeeze(GSS.monthsall.ind(k,1,mon,lats,lons));
                mult_pct = squeeze(GSS.monthspct.ind(k,1,mon,lats,lons));
            else
                mult = squeeze(nanmean(GSS.monthsall.ind(:,1,mon,lats,lons),1));
                mult_pct = squeeze(nanmean(GSS.monthspct.ind(:,1,mon,lats,lons),1));
            end
        end        
        [maxval(i,k),maxind] = max(abs(mult(:)));
        [maxval_pct(i,k),maxind_pct] = max(abs(mult_pct(:)));        
        lattoplot(mon,k,i) = latmesh(maxind);
        lontoplot(mon,k,i) = lonmesh(maxind);
        
        lattoplot_pct(mon,k,i) = latmesh(maxind_pct);
        lontoplot_pct(mon,k,i) = lonmesh(maxind_pct);
        
        if k == size(differences,1) && same
            forerror = squeeze(allsig.percentilesm.ninetyfive(:,lats,lons));
            forerrorextract(mon,:,i) = forerror(:,maxind);
            
            forerrorpct = squeeze(pctsig.percentilesm.ninetyfive(:,lats,lons));
            forerrorextractpct(mon,:,i) = forerrorpct(:,maxind);  
            
            modelmonthsforerror_extract = [];
            datamonthsforerror_extract = [];
            medmonthsforerror_extract = [];
            modelmonthsforerror_pct_extract = [];
            datamonthsforerror_pct_extract = [];
            medmonthsforerror_pct_extract = [];
        elseif ~same
            
            modelmonthsforerror = squeeze(modelprediction_ind_months(k,:,:,lats,lons));
            modelmonthsforerror_extract(i,:,k,:) = modelmonthsforerror(:,:,maxind)';
            datamonthsforerror = squeeze(actualdata_ind_months(k,:,:,lats,lons));
            datamonthsforerror_extract(i,:,k,:) = datamonthsforerror(:,:,maxind)';
            medmonthsforerror = squeeze(medmonths(k,:,:,lats,lons));
            medmonthsforerror_extract(i,:,k,:) = medmonthsforerror(:,:,maxind)';
            
            modelmonthsforerror_pct = squeeze(modelprediction_ind_pct_months2(k,:,:,lats,lons));
            modelmonthsforerror_pct_extract(i,:,k,:) = modelmonthsforerror_pct(:,:,maxind)';
            datamonthsforerror_pct = squeeze(actualdata_ind_pct_months2(k,:,:,lats,lons));
            datamonthsforerror_pct_extract(i,:,k,:) = datamonthsforerror_pct(:,:,maxind)';
            medmonthsforerror_pct = squeeze(medmonths_pct2(k,:,:,lats,lons));
            medmonthsforerror_pct_extract(i,:,k,:) = medmonthsforerror_pct(:,:,maxind)';
            
        end
        
        %[maxrow(k,i),maxcol(k,i)] = find(mult == max(mult(:)));
        GSSmonthtoplot2 = squeeze(GSS.monthsall.ind(k,1,:,lats,lons));
        GSSmonthtoplot(:,k,i) = GSSmonthtoplot2(:,maxind);
        GSSmonthtoplot3 = squeeze(GSS.monthspct.ind(k,1,:,lats,lons));
        GSSmonthtoplot4(:,k,i) = GSSmonthtoplot3(:,maxind);
        if maxind <= size(mult,1)
            if maxind ~= 1
                GSSmonthtoplot_areamean(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                    GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                                
                
                GSSmonthtoplot_areamean_pctsame(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                    GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                                
                
            else
                GSSmonthtoplot_areamean(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),...
                    GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                            
                
                GSSmonthtoplot_areamean_pctsame(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),...
                    GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                            
            end
        elseif maxind == size(mult,1)+1
            GSSmonthtoplot_areamean(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
                GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        
            
            GSSmonthtoplot_areamean_pctsame(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
                GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        
            
        elseif maxind == numel(mult) - size(mult,1)
            GSSmonthtoplot_areamean(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)),...
                GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        
            
            GSSmonthtoplot_areamean_pctsame(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)),...
                GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                    
            
        elseif maxind > numel(mult) - size(mult,1)
            if maxind ~= numel(mult)
                GSSmonthtoplot_areamean(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                    GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                
                
                GSSmonthtoplot_areamean_pctsame(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                    GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                
                
            else                
                GSSmonthtoplot_areamean(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind-1),...
                    GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                
                
                GSSmonthtoplot_areamean_pctsame(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind-1),...
                    GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                
            end
        else
            GSSmonthtoplot_areamean(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
                GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        
            
            GSSmonthtoplot_areamean_pctsame(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
                GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        
            
        end
        
        %pct
        if maxind_pct <= size(mult_pct,1)
            if maxind_pct ~= 1               
                GSSmonthtoplot_areamean_pct(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                    GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
            else                
            
                GSSmonthtoplot_areamean_pct(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),...
                    GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
            end
        elseif maxind_pct == size(mult_pct,1)+1            
            
            GSSmonthtoplot_areamean_pct(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
                GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
            
        elseif maxind_pct == numel(mult_pct) - size(mult_pct,1)
            
            GSSmonthtoplot_areamean_pct(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)),...
                GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
            
        elseif maxind_pct > numel(mult_pct) - size(mult_pct,1)
            if maxind_pct ~= numel(mult_pct)                
                
                GSSmonthtoplot_areamean_pct(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                    GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
                
            else                
                
                GSSmonthtoplot_areamean_pct(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct-1),...
                    GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
            end
        else
           
            GSSmonthtoplot_areamean_pct(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
                GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);            
        end                
    end
end

% constructing bootstrap error
if ~same
    mp = modelmonthsforerror_extract(:,:,:) - medmonthsforerror_extract(:,:,:);
    ap = datamonthsforerror_extract(:,:,:) - medmonthsforerror_extract(:,:,:);
    mpp = modelmonthsforerror_pct_extract(:,:,:) - medmonthsforerror_pct_extract(:,:,:);
    app = datamonthsforerror_pct_extract(:,:,:) - medmonthsforerror_pct_extract(:,:,:);
    for j = 1:size(mp,1)       
                
        for i = 1:size(mp,2)
            
            mpe = squeeze(mp(j,i,:));
            ade = squeeze(ap(j,i,:));
            
            mpme = squeeze(mpp(j,i,:));
            adme = squeeze(app(j,i,:));
            
            
            bootstat(j,i,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap(mpe,ade),mpe);            
            forerrorextract(mon,i,j) = prctile(squeeze(bootstat(j,i,:)),95);
            
            bootstatm(j,i,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap(mpme,adme),mpme);                
            forerrorextractpct(mon,i,j) = prctile(squeeze(bootstatm(j,i,:)),95);

        end
    end
end



end
end
