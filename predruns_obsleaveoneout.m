function [GSS] = predruns_obsleaveoneout(data,waclon,waclat,inputs)

    ERAdata = ReadinERA('/Volumes/ExternalOne/work/data/ERA-Interim/TS/TS_ERA-Interim.nc');
    %data.obstemp,data.toz_zm_nodetrend
    surfacedata = double(permute(nanmean(data.obstemp(:,:,:,inputs.varmonthtomean),4),[3,2,1]));
    %surfacedata_pct = nanmean(data.obstemp(:,:,[data.obspct.upperind,data.obspct.lowerind],inputs.varmonthstomean),4);
    surfacedataind =  double(permute(data.obstemp(:,:,:,inputs.varmonth),[3,4,2,1]));    
    
    %find upper and lower percentiles
    tozautocorr = data.toz_zm_nodetrend;
    upperpctind = find(tozautocorr >= prctile(tozautocorr,80));
    lowerpctind = find(tozautocorr <= prctile(tozautocorr,20));
    
    precondtoz = data.toz_zm_nodetrend;       
    allin = 1:length(precondtoz);
    count = 1;
    for i = 1:length(allin)
        leaveout = allin;
        if i == 1
            leaveout (leaveout == count | leaveout == count+1) = []; 
        elseif i == 37
            leaveout (leaveout == count | leaveout == count-1) = []; 
        else
            leaveout (leaveout == count-1 | leaveout == count | leaveout == count+1) = [];         
        end
        %leaveout (leaveout == count) = [];     
        
        predictors_anom = [leaveout',ones(size(leaveout))'];                
        
        ozone_ind_b = predictors_anom\precondtoz(leaveout)';    
        ozone_anomaly = precondtoz(leaveout) - ozone_ind_b(1).*leaveout;                         
        
        ozone_pct(i) = prctile(ozone_anomaly,50);
        ozone_pct1(i) = prctile(ozone_anomaly,20);
        ozone_pct2(i) = prctile(ozone_anomaly,80);                        
        ozonelowerind(i,:) = find(ozone_anomaly <= ozone_pct1(i));
        ozoneupperind(i,:) = find(ozone_anomaly >= ozone_pct2(i));                      
        
        %ozone_left(i) = precondtoz(i) - ozone_ind_b(1).*i-nanmean(ozone_anomaly);
        %ozone_left(i) = precondtoz(i) - ozone_ind_b(1).*i-ozone_pct(i);
        ozone_left(i) = precondtoz(count) - ozone_ind_b(1).*i-ozone_pct(i);
        ozone_anomaly = ozone_anomaly - nanmedian(ozone_anomaly);
        % is ozone anomaly closer to the 80th percentile or the 20th
        % percentile
        oz1 = ozone_left(i) - ozone_pct1(i);
        oz2 = ozone_left(i) - ozone_pct2(i);        
        
        if oz1 > oz2
            predsign3(i) = -1;                                                
        else
            predsign3(i) = 1;                                               
        end
        
        ozone_left(i) = ozone_left(i);% - nanmean(ozone_left);        
        predsign(i) = sign(ozone_left(i));                                                        
        
        for j = 1:length(ERAdata.latitude)   
            
            if j == 15
                abc = 1;
            end
            b_anom = predictors_anom\squeeze(surfacedata(leaveout,j,:));
            temp_anomaly = squeeze(surfacedata(leaveout,j,:))-b_anom(1,:).*leaveout';% - b_anom(2,:);
            %temp_left(i,j,:) = squeeze(surfacedata(i,j,:))'-b_anom(1,:).*i;% - b_anom(2,:);            
            temp_left(i,j,:) = squeeze(surfacedata(count,j,:))'-b_anom(1,:).*i;% - b_anom(2,:);            
            temp_left(i,j,:) = squeeze(temp_left(i,j,:))' - nanmedian(temp_anomaly);
            temp_anomaly = temp_anomaly - nanmean(temp_anomaly);            
            
            temp_anomaly_upper = squeeze(nanmean(temp_anomaly(ozoneupperind(i,:),:)));
            temp_anomaly_lower = squeeze(nanmean(temp_anomaly(ozonelowerind(i,:),:)));    
            tempdiff(i,j,:) = temp_anomaly_upper - temp_anomaly_lower;
            signchange(i,j,:) = sign(temp_anomaly_upper-temp_anomaly_lower);%-nanmedian(temp_anomaly)                        
            
            datasign(i,j,:) = sign(temp_left(i,j,:));                        
            
             for l = 1:length(squeeze(datasign(i,j,:)))
                
                if signchange(i,j,l) == 1
                    pred(i,j,l) = 1;
                else
                    pred(i,j,l) = -1;
                end              
                
                if pred(i,j,l) == 1
                    predsign2(i,j,l) = predsign(i);
                    predsign4(i,j,l) = predsign3(i);
                    if predsign(i) == datasign(i,j,l)
                        correct(i,j,l) = 1;
                    else
                        correct(i,j,l) = 0;
                    end
                    
                    if predsign3(i) == datasign(i,j,l)
                        correct2(i,j,l) = 1;
                    else
                        correct2(i,j,l) = 0;
                    end
                    
                elseif pred(i,j,l) == -1
                    predsign2(i,j,l) = predsign(i).*-1;
                    predsign4(i,j,l) = predsign3(i).*-1;
                    if predsign(i) == datasign(i,j,l)
                        correct(i,j,l) = 0;
                    else
                        correct(i,j,l) = 1;
                    end
                    
                    if predsign3(i) == datasign(i,j,l)
                        correct2(i,j,l) = 0;
                    else
                        correct2(i,j,l) = 1;
                    end
                    
                end
                
             end         
            
             for m = 1:size(surfacedataind,2)
                b_months = predictors_anom\squeeze(surfacedataind(leaveout,m,j,:));
                temp_anomaly_months = squeeze(surfacedataind(leaveout,m,j,:)) - leaveout'.*b_months(1,:);% - b(2,:); %looks good
                %temp_left_months = squeeze(surfacedataind(i,m,j,:))' - i.*b_months(1,:);% - b(2,:); %looks good                        
                temp_left_months = squeeze(surfacedataind(count,m,j,:))' - i.*b_months(1,:);% - b(2,:); %looks good                        
                temp_left_months = squeeze(temp_left_months) - squeeze(nanmedian(temp_anomaly_months,1));
                temp_anomaly_months = squeeze(temp_anomaly_months)' - squeeze(nanmean(temp_anomaly_months,2))';
                temp_anomaly_upper_months = squeeze(nanmean(temp_anomaly_months(:,ozoneupperind(i,:)),2));
                temp_anomaly_lower_months = squeeze(nanmean(temp_anomaly_months(:,ozonelowerind(i,:)),2));                        
                
                signchange_months = sign(temp_anomaly_upper_months-temp_anomaly_lower_months);%-squeeze(nanmedian(temp_anomaly_months,1))';
                
                datasign_months(i,j,m,:) = sign(temp_left_months);
                
                for l = 1:length(squeeze(datasign_months))
                    if signchange_months(l) == 1 
                        pred_months(i,j,m,l) = 1;
                    elseif signchange_months(l) == -1 
                        pred_months(i,j,m,l) = -1;        
                    end

                    if pred_months(i,j,m,l) == 1
                        predsign_months2(i,j,m,l) = predsign(i).*1;
                        if predsign(i) == datasign_months(i,j,m,l)
                            correct_months(i,j,m,l) = 1;
                        else
                            correct_months(i,j,m,l) = 0;
                        end
                    elseif pred_months(i,j,m,l) == -1
                        predsign_months2(i,j,m,l) = predsign(i).*-1;
                        if predsign(i) == datasign_months(i,j,m,l)
                            correct_months(i,j,m,l) = 0;
                        else
                            correct_months(i,j,m,l) = 1;
                        end
                    end
                end  
                
             end                          
        end         
        count = count+1;
    end
    
    
    %% extract percentiles
%     correctpct = correct([data.obspct.upperind,data.obspct.lowerind],:,:);    
%     correctpct2 = correct2([data.obspct.upperind,data.obspct.lowerind],:,:);    
%     correct_months_pct = correct_months([data.obspct.upperind,data.obspct.lowerind],:,:,:);
%         
%     predsignpct = predsign2([data.obspct.upperind,data.obspct.lowerind],:,:);
%     datasignpct = datasign([data.obspct.upperind,data.obspct.lowerind],:,:);
%     
%     predsignpct_months = predsign_months2([data.obspct.upperind,data.obspct.lowerind],:,:,:);
%     datasignpct_months = datasign_months([data.obspct.upperind,data.obspct.lowerind],:,:,:);    
    
    
    correctpct = correct([upperpctind,lowerpctind],:,:);        
    correct_months_pct = correct_months([upperpctind,lowerpctind],:,:,:);
        
    predsignpct = predsign2([upperpctind,lowerpctind],:,:);
    datasignpct = datasign([upperpctind,lowerpctind],:,:);
    
    predsignpct_months = predsign_months2([upperpctind,lowerpctind],:,:,:);
    datasignpct_months = datasign_months([upperpctind,lowerpctind],:,:,:);    
    
    %% calculating error
    
    %% testing bootsrapping
%% testing bootsrapping
mp = permute(predsign2,[2,3,1]);
%mp = permute(modelprediction_ind-med,[3,4,1,2]);
mp = mp(:,:,:);
ad = permute(datasign,[2,3,1]);
%ad = permute(actualdata_ind-med,[3,4,1,2]);
ad = ad(:,:,:);

mpm = permute(predsign_months2,[3,2,4,1]);
mpm = mpm(:,:,:,:);
adm = permute(datasign_months,[3,2,4,1]);
adm = adm(:,:,:,:);


if ~exist('/Volumes/ExternalOne/work/data/predruns/output/HSS/obs_empcomp_Allperc.mat','file')
    bootstat = zeros(size(mp,1),size(mp,2),500);
    percentiles.eighty = zeros(size(mp,1),size(mp,2));
    percentiles.ninety = zeros(size(mp,1),size(mp,2));
    percentiles.ninetyfive = zeros(size(mp,1),size(mp,2));
    bootstatm = zeros(5,size(mp,1),size(mp,2),500);
    percentilesm.eighty = zeros(5,size(mp,1),size(mp,2));
    percentilesm.ninety = zeros(5,size(mp,1),size(mp,2));
    percentilesm.ninetyfive = zeros(5,size(mp,1),size(mp,2));
    for i = 1:size(mp,1)        
        tic;
        for j = 1:size(mp,2)       
            mpe = squeeze(mp(i,j,:));
            ade = squeeze(ad(i,j,:));
            bootstat(i,j,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap_compemp(mpe,ade),mpe);
            percentiles.eighty(i,j) = prctile(squeeze(bootstat(i,j,:)),80);
            percentiles.ninty(i,j) = prctile(squeeze(bootstat(i,j,:)),90);
            percentiles.ninetyfive(i,j) = prctile(squeeze(bootstat(i,j,:)),95);
            for k = 1:5
                mpme = squeeze(mpm(k,i,j,:));
                adme = squeeze(adm(k,i,j,:));
                bootstatm(k,i,j,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap_compemp(mpme,adme),mpme);
                percentilesm.eighty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),80);
                percentilesm.ninty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),90);
                percentilesm.ninetyfive(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),95);
            end
        end
        toc;
    end
    save('/Volumes/ExternalOne/work/data/predruns/output/HSS/obs_empcomp_Allperc.mat','bootstat','bootstatm','percentiles','percentilesm');
else
    allsig = load('/Volumes/ExternalOne/work/data/predruns/output/HSS/obs_empcomp_Allperc.mat');
end

%%


mp = permute(predsignpct,[2,3,1]);
%mp = permute(modelprediction_ind-med,[3,4,1,2]);
mp = mp(:,:,:);
ad = permute(datasignpct,[2,3,1]);
%ad = permute(actualdata_ind-med,[3,4,1,2]);
ad = ad(:,:,:);

mpm = permute(predsignpct_months,[3,2,4,1]);
mpm = mpm(:,:,:,:);
adm = permute(datasignpct_months,[3,2,4,1]);
adm = adm(:,:,:,:);


if ~exist('/Volumes/ExternalOne/work/data/predruns/output/HSS/obs_empcomp_Pctperc.mat','file')
    bootstat = zeros(size(mp,1),size(mp,2),500);
    percentiles.eighty = zeros(size(mp,1),size(mp,2));
    percentiles.ninety = zeros(size(mp,1),size(mp,2));
    percentiles.ninetyfive = zeros(size(mp,1),size(mp,2));
    bootstatm = zeros(5,size(mp,1),size(mp,2),500);
    percentilesm.eighty = zeros(5,size(mp,1),size(mp,2));
    percentilesm.ninety = zeros(5,size(mp,1),size(mp,2));
    percentilesm.ninetyfive = zeros(5,size(mp,1),size(mp,2));
    for i = 1:size(mp,1)        
        tic;
        for j = 1:size(mp,2)       
            mpe = squeeze(mp(i,j,:));
            ade = squeeze(ad(i,j,:));
            bootstat(i,j,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap_compemp(mpe,ade),mpe);
            percentiles.eighty(i,j) = prctile(squeeze(bootstat(i,j,:)),80);
            percentiles.ninty(i,j) = prctile(squeeze(bootstat(i,j,:)),90);
            percentiles.ninetyfive(i,j) = prctile(squeeze(bootstat(i,j,:)),95);
            for k = 1:5
                mpme = squeeze(mpm(k,i,j,:));
                adme = squeeze(adm(k,i,j,:));
                bootstatm(k,i,j,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap_compemp(mpme,adme),mpme);
                percentilesm.eighty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),80);
                percentilesm.ninty(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),90);
                percentilesm.ninetyfive(k,i,j) = prctile(squeeze(bootstatm(k,i,j,:)),95);
            end
        end
        toc;
    end
    save('/Volumes/ExternalOne/work/data/predruns/output/HSS/obs_empcomp_Pctperc.mat','bootstat','bootstatm','percentiles','percentilesm');
else
    allsigpct = load('/Volumes/ExternalOne/work/data/predruns/output/HSS/obs_empcomp_Pctperc.mat');
end
    
    
    %% calculate skill scores
    
    HA_ind = sum(predsign2 >= 0 & datasign >= 0,1);
    HB_ind = sum(predsign2 > datasign,1);
    HC_ind = sum(predsign2 < datasign,1);
    HD_ind = sum(predsign2 <= 0 & datasign <= 0,1);

    HSS.ind = 2.*(HA_ind.*HD_ind - HB_ind.*HC_ind)./(((HA_ind+HC_ind).*(HC_ind+HD_ind))+((HA_ind+HB_ind).*(HB_ind+HD_ind))).*100; %2(ad-bc)/(a+c)(c+d)+(a+b)(b+d)        
    
    GSS.all.ind = (sum(correct,1) - size(correct,1)./2)./...
        (size(correct,1)-size(correct,1)./2)*100;  

    GSS.pct.ind = (sum(correctpct,1) - size(correctpct,1)./2)./...
        (size(correctpct,1)-size(correctpct,1)./2)*100;      

%     GSS.all.ind = (sum(correct2,1) - size(correct2,1)./2)./...
%         (size(correct2,1)-size(correct2,1)./2)*100;  
% 
%     GSS.pct.ind = (sum(correctpct2,1) - size(correctpct2,1)./2)./...
%         (size(correctpct2,1)-size(correctpct2,1)./2)*100;      
    
    
    

    GSS.monthsall.ind = (sum(correct_months,1) - size(correct_months,1)./2)./...
        (size(correct_months,1)-size(correct_months,1)./2)*100;  
    GSS.monthspct.ind = (sum(correct_months_pct,1) - size(correct_months_pct,1)./2)./...
        (size(correct_months_pct,1)-size(correct_months_pct,1)./2)*100;      
        
    
%     GSS.monthsall.ind = permute(GSS.monthsall.ind,[1,3,2,4]);
%     GSS.monthspct.ind = permute(GSS.monthspct.ind,[1,3,2,4]);
    
p = zeros(size(GSS.all.ind));
p (GSS.all.ind < reshape(allsig.percentiles.ninetyfive,[1,size(allsig.percentiles.ninetyfive)])) = -1;

p2 = zeros(size(GSS.all.ind));
p2 (GSS.pct.ind < reshape(allsigpct.percentiles.ninetyfive,[1,size(allsigpct.percentiles.ninetyfive)])) = -1;

%% plot
titles = {'Observed HSSs','Observed HSSs during ozone extremes','Observed percentile general skill score (leave one out)'};
%toplotp = permute(cat(1,p,p2),[1,3,2]);
longitudes = double(ERAdata.longitude);
latitudes = double(ERAdata.latitude);

% subplotmaps(permute(cat(1,GSS.all.ind,GSS.pct.ind),[1,3,2]),longitudes,latitudes,{'seq','YlGnBu'},0,[],16,titles,'Longitude','Latitude','HSS','on',...
%         [0,80],11,[longitudes(1:24:end)]-180,[longitudes(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');         
toplotp = permute(cat(1,p,p2),[1,3,2]);    
[fig,fh] = subplotmaps(permute(cat(1,GSS.all.ind,GSS.pct.ind),[1,3,2]),longitudes,latitudes,{'seq','YlGnBu'},0,toplotp,22,titles,'Longitude','Latitude','HSS','on',...
    [0,80],21,[longitudes(1:40:end)]-180,[longitudes(1:40:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');         

% subplotmaps(permute(cat(1,GSS.all.ind,GSS.pct.ind),[1,3,2]),longitudes,latitudes,{'div','RdBu'},1,[],16,titles,'Longitude','Latitude','HSS','on',...
%     [-80,80],22,[longitudes(1:24:end)]-180,[longitudes(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');         
 
% subplotmaps(permute(cat(1,GSS.all2.mean,GSS.pct2.mean),[1,3,2]),longitude,latitude,{'seq','YlGnBu'},0,[],16,titles,'Longitude','Latitude','','on',...
%         [0,50],11,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');             
    
child = get(fig,'children');
axes1 = child(4);
axes2 = child(2);

axes(axes1);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','a','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',24,... % 
    'EdgeColor','none','fontweight','bold');    

axes(axes2);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','b','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',24,... % 
    'EdgeColor','none','fontweight','bold');    


set(gcf,'Renderer','Painters');
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/Obs_GSSothercolor_nodetrend_upto80_',monthnames(inputs.varmonth,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'.eps'];
print(filename,'-depsc');                

abc = 1
    
    %% 
% %Find areas of large difference and large GSS
% % difference above 3 degrees and GSS above 4
% % areas_lons = [290,320;30,120;240,290;30,120];%lons (Greenland,Russia,America,Asia)
% % areas_lats = [60,80;50,75;30,60;15,45];%lats
% lstyle = {':','-','-.','--','-'};        
% fsize = 16;
% xticklab = {'March','April','May','June','July'};
% 
% same = 0;
% corrtoplot = 0;
% difftoplot = 1;
% for mon = 1:5
% 
% if same
%     sameext = 'same';
% else
%     sameext = '';
% end
% if corrtoplot
%     corrtoplotext = 'corr';
% elseif difftoplot
%     corrtoplotext = 'diff';
% else
%     corrtoplotext = 'GSS';
% end
% areas_lons = [80,180;30,80;240,290;60,120];%lons (Greenland,East Russia, West Russia,America,Asia)
% areas_lats = [60,80;60,80;35,60;25,45];%lats
% 
% cbrewqual = cbrewer('qual','Set1',10);
% cbrewqual2 = cbrewqual([3,4,1,10,8],:);
% 
% %r_monthsmean = squeeze(nanmean(r_months,1));
% for i = 1:size(areas_lons,1)
%     lats = ERAdata.latitude > areas_lats(i,1) & ERAdata.latitude < areas_lats(i,2);
%     lons = ERAdata.longitude > areas_lons(i,1) & ERAdata.longitude < areas_lons(i,2);
%     latextract = ERAdata.latitude(lats);
%     lonextract = ERAdata.longitude(lons);
%     [latmesh,lonmesh] = meshgrid(latextract,lonextract);
%     
%     latmesh = latmesh';
%     lonmesh = lonmesh';
%     for k = 1:size(GSS.all.ind,1)    
%         
% %         mult = squeeze(GSS.monthsall.ind(k,1,mon,lats,lons));
% %         mult_pct = squeeze(GSS.monthspct.ind(k,1,mon,lats,lons));
%         
% %         mult = squeeze(GSS.all.mean(1,lats,lons));        
% %         mult_pct = squeeze(GSS.pct.mean(1,lats,lons));        
%         
%         if corrtoplot
%             if same
%                 mult = squeeze(nanmean(r_monthsmean(:,mon,lats,lons),1));        
%                 mult_pct = squeeze(nanmean(r_monthsmean(:,mon,lats,lons),1));        
%             else
%                 mult = squeeze(nanmean(r_months(:,mon,lats,lons),1));        
%                 mult_pct = squeeze(nanmean(r_months(:,mon,lats,lons),1));                                        
%                 
%             end
%         elseif difftoplot            
%             diff = permute(data.differences_ind(mon,:,:),[1,3,2]);
%             
%             mult = squeeze(diff(1,lats,lons));        
%             mult_pct = squeeze(diff(1,lats,lons));                    
%         else
%             if ~same
%                 mult = squeeze(GSS.monthsall.ind(1,mon,lats,lons));
%                 mult_pct = squeeze(GSS.monthspct.ind(1,mon,lats,lons));
%             else
%                 mult = squeeze(nanmean(GSS.monthsall.ind(1,mon,lats,lons),1));
%                 mult_pct = squeeze(nanmean(GSS.monthspct.ind(1,mon,lats,lons),1));
% %                 mult = squeeze(GSS.all.mean(1,lats,lons));        
% %                 mult_pct = squeeze(GSS.pct.mean(1,lats,lons));        
%             end
%         end
%         %mult = squeeze(nanmean(r_monthsmean(:,mon,lats,lons),1));        
%         %mult = squeeze(nanmean(r(k,:,lats,lons),2));        
%         %mult = squeeze(GSS.cheatpct.ind(k,1,lats,lons));
%         [maxval(i,k),maxind] = max(abs(mult(:)));
%         [maxval_pct(i,k),maxind_pct] = max(abs(mult_pct(:)));
%         %[~,maxind2] = max(mult2(:));
%         lattoplot(k,i) = latmesh(maxind);
%         lontoplot(k,i) = lonmesh(maxind);
%         
%         lattoplot_pct(k,i) = latmesh(maxind_pct);
%         lontoplot_pct(k,i) = lonmesh(maxind_pct);
%         
%         %[maxrow(k,i),maxcol(k,i)] = find(mult == max(mult(:)));
%         GSSmonthtoplot2 = squeeze(GSS.monthsall.ind(1,:,lats,lons));
%         
%         GSSmonthtoplot3 = squeeze(GSS.monthspct.ind(1,:,lats,lons));
%         
%         if maxind <= size(mult,1)
%             if maxind ~= 1
%                 GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
%                     GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                                
%                 
%                 GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
%                     GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                                
%                 
%             else
%                 GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),...
%                     GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                            
%                 
%                 GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),...
%                     GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                            
%             end
%         elseif maxind == size(mult,1)+1
%             GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
%                 GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
%                 GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        
%             
%             GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
%                 GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
%                 GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        
%             
%         elseif maxind == numel(mult) - size(mult,1)
%             GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
%                 GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)),...
%                 GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        
%             
%             GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
%                 GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)),...
%                 GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                    
%             
%         elseif maxind > numel(mult) - size(mult,1)
%             if maxind ~= numel(mult)
%                 GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
%                     GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                
%                 
%                 GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
%                     GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                
%                 
%             else                
%                 GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind-1),...
%                     GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                
%                 
%                 GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind-1),...
%                     GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                
%             end
%         else
%             GSSmonthtoplot_areamean(:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
%                 GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
%                 GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        
%             
%             GSSmonthtoplot_areamean_pctsame(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
%                 GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
%                 GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        
%             
%         end
%         
%         %pct
%         if maxind_pct <= size(mult_pct,1)
%             if maxind_pct ~= 1               
%                 GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
%                     GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
%             else                
%             
%                 GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),...
%                     GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
%             end
%         elseif maxind_pct == size(mult_pct,1)+1            
%             
%             GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
%                 GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
%                 GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
%             
%         elseif maxind_pct == numel(mult_pct) - size(mult_pct,1)
%             
%             GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
%                 GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)),...
%                 GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
%             
%         elseif maxind_pct > numel(mult_pct) - size(mult_pct,1)
%             if maxind_pct ~= numel(mult_pct)                
%                 
%                 GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
%                     GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
%                 
%             else                
%                 
%                 GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct-1),...
%                     GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
%             end
%         else
%            
%             GSSmonthtoplot_areamean_pct(:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
%                 GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
%                 GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
%             
%         end
%         
%         
%     end
% end
% 
% 
% areas_lons2 = areas_lons;
% areas_lats2 = areas_lats;
% 
% fig = figure;
% set(fig,'position',[100 100 1000 1200],'color','white','Visible','on');
% 
% subplot(3,2,1)
% %toplot = squeeze(nanmean(GSSmonthtoplot,2));
% toplot = squeeze(nanmean(GSSmonthtoplot_areamean,2));
% for i = 1:4
%     plot(toplot(:,i),'o','linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
%     hold on
% end
% set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
% xlabel('Month','fontsize',fsize+2);
% xlim([.5,5.5]);
% ylim([-40 80]);
% ylabel('HSS','fontsize',fsize+2);
% title('Ensemble mean','fontsize',fsize+4)
% %March TCO - surface temperature ensemble mean HSS
% lh = legend('Eastern Russia','Western Russia','Northern America','Asia');
% set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
% 
% subplot(3,2,2)
% if same
%     toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pctsame,2));
% else
%     toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pct,2));
% end
% for i = 1:4    
%     plot(toplot(:,i),'o','linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
%     hold on
% end
% set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
% xlabel('Month','fontsize',fsize+2);
% %ylabel('HSS','fontsize',fsize+2);
% title('Ensemble mean 20th percentile','fontsize',fsize+4)
% xlim([.5,5.5]);
% ylim([-40 80]);
% % plotting latitudes
% for i = 1:numel(lontoplot)
%     if lontoplot(i) > 180
%         lontoplot(i) = lontoplot(i) - 360;
%     end
% end
% 
% for i = 1:numel(lontoplot_pct)
%     if lontoplot_pct(i) > 180
%         lontoplot_pct(i) = lontoplot_pct(i) - 360;
%     end
% end
% 
% for i = 1:numel(areas_lons2)
%     if areas_lons2(i) > 180
%         areas_lons2(i) = areas_lons2(i) - 360;
%     end
% end
% 
% subplot(3,1,2);
% cbrewqual = cbrewer('qual','Set1',10);
% cbrewqual2 = cbrewqual([3,4,1,10,8],:);
% cbrewqual3 = cbrewqual2./2;
% % fig = figure;
% % set(fig,'position',[100 100 1000 400],'color','white');
% m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
% m_coast('color','k','LineWidth',1);
% m_grid('ytick',[0:15:90],'xtick',[-180:60:360],'XaxisLocation','bottom','fontsize',fsize);            
% ylab = ylabel('Latitude','fontsize',fsize+2);
% ylabpos = get(ylab,'position');
% set(ylab,'position',[ylabpos(1)-.2,ylabpos(2),ylabpos(3)]);
% xlabel('Longitude','fontsize',fsize+2);
% if corrtoplot
%     title(['Location of maximum ensemble average March TCO and ',monthnames(mon+2,0,0), ' TS correlation'],'fontsize',fsize+2);
% elseif difftoplot
%     title(['Location of maximum ensemble member March TCO and ',monthnames(mon+2,0,0), ' TS differences'],'fontsize',fsize+2);
% else
%     title(['Location of maximum ensemble member March TCO and ',monthnames(mon+2,0,0), ' HSS'],'fontsize',fsize+2);
% end
% 
% box on
% hold on
% for i = 1:4
%     %m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
%     m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
%     m_plot([areas_lons2(i,1),areas_lons2(i,2),areas_lons2(i,2),areas_lons2(i,1),areas_lons2(i,1)],[areas_lats(i,1),areas_lats(i,1),areas_lats(i,2),areas_lats(i,2),areas_lats(i,1)],...
%         'LineStyle','-','LineWidth',4,'color',cbrewqual2(i,:),'LineStyle',lstyle{i})
% end
% axes.SortMethod='ChildOrder';
% 
% % pct
% if ~same && ~corrtoplot && ~difftoplot
%     subplot(3,1,3);
%     m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
%     m_coast('color','k','LineWidth',1);
%     m_grid('ytick',[0:15:90],'xtick',[-180:60:360],'XaxisLocation','bottom','fontsize',fsize);            
%     ylab = ylabel('Latitude','fontsize',fsize+2);
%     ylabpos = get(ylab,'position');
%     set(ylab,'position',[ylabpos(1)-.2,ylabpos(2),ylabpos(3)]);
%     xlabel('Longitude','fontsize',fsize+2);
%     if corrtoplot
%         title(['Location of maximum ensemble member March TCO and ',monthnames(mon+2,0,0), ' TS correlation'],'fontsize',fsize+2);
%     elseif difftoplot
%         title(['Location of maximum ensemble member March TCO and ',monthnames(mon+2,0,0), ' TS differences'],'fontsize',fsize+2);
%     else
%         title(['Location of maximum ensemble member March TCO and ',monthnames(mon+2,0,0), ' 20th percentile HSS'],'fontsize',fsize+2);
%     end
%     box on
%     hold on
%     for i = 1:5
%         %m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
%         m_plot(lontoplot_pct(:,i)',lattoplot_pct(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
%         m_plot([areas_lons2(i,1),areas_lons2(i,2),areas_lons2(i,2),areas_lons2(i,1),areas_lons2(i,1)],[areas_lats(i,1),areas_lats(i,1),areas_lats(i,2),areas_lats(i,2),areas_lats(i,1)],...
%             'LineStyle','-','LineWidth',4,'color',cbrewqual2(i,:),'LineStyle',lstyle{i})
%     end
%     axes.SortMethod='ChildOrder';
% end
% 
% annotation('textbox',[.01 .98 1 0],'String','March TCO - surface temperature ensemble mean HSS','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
%         'EdgeColor','none','fontweight','bold');    
% 
% filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/Observations_Ind_',corrtoplotext,'_','mean_Locations_withpct_test','_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),monthnames(mon+2,0,0),sameext];
% %print(filename,'-depsc');
%end

% diffforplot = differences;
% diffforplot = circshift(diffforplot,[0,144/2,0,0]);

diffforplot = data.differences_ind;
diffforplot = circshift(diffforplot,[0,240/2,0]);
diffforplot (diffforplot > 5) = 5;
diffforplot (diffforplot < -5) = -5;




cbrewdiff = cbrewer('div','RdBu',21);
cbrewdiff = flipud(cbrewdiff);
cbrewdiff = [cbrewdiff(1:10,:);cbrewdiff(12:end,:)];

cbrewqual = cbrewer('qual','Set1',10);
cbrewqual2 = cbrewqual([4,1,10,8],:);
cbrewqual3 = cbrewqual2./2;

xticklab = {'March','April','May','June','July'};

lstyle = {':','-','-.','--','-'};    
fsize = 16;
Mks = {'d','s','o','v'};



%%
same = 1;
% corrtoplot = 0;
difftoplot = 1;
corrtoplot = 0;
 areas_lons = [90,180;30,90;240,290;60,120];%lons (Greenland,East Russia, West Russia,America,Asia)
areas_lats = [60,80;60,80;35,60;25,45];%lats

for j = 1:1
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
    
        for i = 1:size(areas_lons,1)
            lats = ERAdata.latitude > areas_lats(i,1) & ERAdata.latitude < areas_lats(i,2);
            lons = ERAdata.longitude > areas_lons(i,1) & ERAdata.longitude < areas_lons(i,2);
            latextract = ERAdata.latitude(lats);
            lonextract = ERAdata.longitude(lons);
            [latmesh,lonmesh] = meshgrid(latextract,lonextract);

            latmesh = latmesh';
            lonmesh = lonmesh';
            for k = 1:1           

                if corrtoplot               
                elseif difftoplot            
                    diff = permute(data.differences_ind(:,:,:),[1,3,2]);     
                    if same                
                        mult = squeeze(nanmean(diff(mon,lats,lons),1));        
                        mult_pct = squeeze(nanmean(diff(mon,lats,lons),1));        
                    else
                        mult = squeeze(diff(mon,lats,lons));        
                        mult_pct = squeeze(diff(mon,lats,lons));        
                    end
                else                
                end        
                [maxval(i,k),maxind] = max(abs(mult(:)));
                [maxval_pct(i,k),maxind_pct] = max(abs(mult_pct(:)));        
                lattoplot(j).g(mon,k,i) = latmesh(maxind);
                lontoplot(j).g(mon,k,i) = lonmesh(maxind);

                lattoplot_pct(j).g(mon,k,i) = latmesh(maxind_pct);
                lontoplot_pct(j).g(mon,k,i) = lonmesh(maxind_pct);

                if same
                    forerror = squeeze(allsig.percentilesm.ninetyfive(:,lats,lons));
                    forerrorextract(j).g(mon,:,i) = forerror(:,maxind);

                    forerrorpct = squeeze(allsigpct.percentilesm.ninetyfive(:,lats,lons));
                    forerrorextractpct(j).g(mon,:,i) = forerrorpct(:,maxind);  

                    modelmonthsforerror_extract = [];
                    datamonthsforerror_extract = [];
                    %medmonthsforerror_extract = [];
                    modelmonthsforerror_pct_extract = [];
                    datamonthsforerror_pct_extract = [];
                    %medmonthsforerror_pct_extract = [];
                elseif ~same

                    %modelmonthsforerror = squeeze(modelprediction_ind_months(k,:,:,lats,lons));
                    modelmonthsforerror = permute(squeeze(predsign_months2(k,:,lats,:,lons)),[1,3,2,4]);
                    modelmonthsforerror_extract(i,:,k,:) = modelmonthsforerror(:,:,maxind)';
                    %datamonthsforerror = squeeze(actualdata_ind_months(k,:,:,lats,lons));
                    datamonthsforerror = permute(squeeze(datasign_months(k,:,lats,:,lons)),[1,3,2,4]);
                    datamonthsforerror_extract(i,:,k,:) = datamonthsforerror(:,:,maxind)';
                    %medmonthsforerror = squeeze(medmonths(k,:,:,lats,lons));
                    %medmonthsforerror_extract(i,:,k,:) = medmonthsforerror(:,:,maxind)';

                    modelmonthsforerror_pct = permute(squeeze(predsignpct_months(k,:,lats,:,lons)),[1,3,2,4]);
                    modelmonthsforerror_pct_extract(i,:,k,:) = modelmonthsforerror_pct(:,:,maxind)';
                    datamonthsforerror_pct = permute(squeeze(datasignpct_months(k,:,lats,:,lons)),[1,3,2,4]);
                    datamonthsforerror_pct_extract(i,:,k,:) = datamonthsforerror_pct(:,:,maxind)';
                    %medmonthsforerror_pct = squeeze(medmonths_pct2(k,:,:,lats,lons));
                    %medmonthsforerror_pct_extract(i,:,k,:) = medmonthsforerror_pct(:,:,maxind)';

                end

                %[maxrow(k,i),maxcol(k,i)] = find(mult == max(mult(:)));
                GSSmonthtoplot2 = permute(squeeze(GSS.monthsall.ind(1,lats,:,lons)),[2,1,3]);
                GSSmonthtoplot(:,k,i) = GSSmonthtoplot2(:,maxind);
                GSSmonthtoplot3 = permute(squeeze(GSS.monthspct.ind(1,lats,:,lons)),[2,1,3]);
                GSSmonthtoplot4(:,k,i) = GSSmonthtoplot3(:,maxind);
                if maxind <= size(mult,1)
                    if maxind ~= 1
                        GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                            GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                                

                        GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                            GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                                

                    else
                        GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),...
                            GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1)],2);                            

                        GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),...
                            GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1)],2);                            
                    end
                elseif maxind == size(mult,1)+1
                    GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                        GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
                        GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        

                    GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                        GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
                        GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        

                elseif maxind == numel(mult) - size(mult,1)
                    GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                        GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)),...
                        GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        

                    GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                        GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)),...
                        GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                    

                elseif maxind > numel(mult) - size(mult,1)
                    if maxind ~= numel(mult)
                        GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                            GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                

                        GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                            GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                

                    else                
                        GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind-1),...
                            GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                                

                        GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind-1),...
                            GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                                
                    end
                else
                    GSSmonthtoplot_areamean(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot2(:,maxind),GSSmonthtoplot2(:,maxind+1),GSSmonthtoplot2(:,maxind-1),...
                        GSSmonthtoplot2(:,maxind+size(mult,1)-1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),GSSmonthtoplot2(:,maxind+size(mult,1)+1),...
                        GSSmonthtoplot2(:,maxind-size(mult,1)-1),GSSmonthtoplot2(:,maxind-size(mult,1)+1),GSSmonthtoplot2(:,maxind-size(mult,1))],2);                        

                    GSSmonthtoplot_areamean_pctsame(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind),GSSmonthtoplot3(:,maxind+1),GSSmonthtoplot3(:,maxind-1),...
                        GSSmonthtoplot3(:,maxind+size(mult,1)-1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),GSSmonthtoplot3(:,maxind+size(mult,1)+1),...
                        GSSmonthtoplot3(:,maxind-size(mult,1)-1),GSSmonthtoplot3(:,maxind-size(mult,1)+1),GSSmonthtoplot3(:,maxind-size(mult,1))],2);                        

                end

                %pct
                if maxind_pct <= size(mult_pct,1)
                    if maxind_pct ~= 1               
                        GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                            GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
                    else                

                        GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),...
                            GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1)],2);
                    end
                elseif maxind_pct == size(mult_pct,1)+1            

                    GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                        GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
                        GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);

                elseif maxind_pct == numel(mult_pct) - size(mult_pct,1)

                    GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                        GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)),...
                        GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);

                elseif maxind_pct > numel(mult_pct) - size(mult_pct,1)
                    if maxind_pct ~= numel(mult_pct)                

                        GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                            GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);

                    else                

                        GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct-1),...
                            GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);
                    end
                else

                    GSSmonthtoplot_areamean_pct(j).g(mon,:,k,i) = nanmean([GSSmonthtoplot3(:,maxind_pct),GSSmonthtoplot3(:,maxind_pct+1),GSSmonthtoplot3(:,maxind_pct-1),...
                        GSSmonthtoplot3(:,maxind_pct+size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct+size(mult,1)+1),...
                        GSSmonthtoplot3(:,maxind_pct-size(mult,1)-1),GSSmonthtoplot3(:,maxind_pct-size(mult,1)+1),GSSmonthtoplot3(:,maxind_pct-size(mult,1))],2);            
                end                
            end
        end
    
    

        %% constructing bootstrap error
%         if ~same
%             mp = modelmonthsforerror_extract(:,:,:);% - medmonthsforerror_extract(:,:,:);
%             ap = datamonthsforerror_extract(:,:,:); %- medmonthsforerror_extract(:,:,:);
%             mpp = modelmonthsforerror_pct_extract(:,:,:);% - medmonthsforerror_pct_extract(:,:,:);
%             app = datamonthsforerror_pct_extract(:,:,:);% - medmonthsforerror_pct_extract(:,:,:);
%             for j = 1:size(mp,1)       
% 
%                 for i = 1:size(mp,2)
% 
%                     mpe = squeeze(mp(j,i,:));
%                     ade = squeeze(ap(j,i,:));
% 
%                     mpme = squeeze(mpp(j,i,:));
%                     adme = squeeze(app(j,i,:));
% 
% 
%                     bootstatspec(j,i,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap_compemp(mpe,ade),mpe);            
%                     forerrorextract(2).g(mon,i,j) = prctile(squeeze(bootstatspec(j,i,:)),95);
% 
%                     bootstatmspec(j,i,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap_compemp(mpme,adme),mpme);                
%                     forerrorextractpct(2).g(mon,i,j) = prctile(squeeze(bootstatmspec(j,i,:)),95);
% 
%                 end
%             end
%         end
% 
%     end
%     same = 0;
end

%%

if corrtoplot
    corrtoplotext = 'corr';
elseif difftoplot
    corrtoplotext = 'diff';
else
    corrtoplotext = 'GSS';
end

if same
    sameext = 'same';
else
    sameext = '';
end

titles = {'Observations','Observations, ozone extremes','Ensemble member','Ensemble member 20^t^h percentiles'}; 
labels = {'a','b'};
for mon = 1:5
    count = 1;
    fig = figure;
    set(fig,'position',[100 1 1000 984],'color','white','Visible','on');
    for j = 1:1
        areas_lons2 = areas_lons;
        areas_lats2 = areas_lats;
               
        sp(count) = subplot(3,2,count);
        sppos = get(sp(count),'position');
        if count == 3
            set(sp(count),'position',[sppos(1)-.025,sppos(2)+.020,sppos(3:4)]);            
        else
            set(sp(count),'position',[sppos(1)-.025,sppos(2:4)]);
        end
        sppos = get(sp(count),'position');
        
        %toplot = squeeze(nanmean(GSSmonthtoplot,2));
        toplot = squeeze(nanmean(GSSmonthtoplot_areamean(j).g(mon,:,:,:),3));
        totimes = [-3,-1,1,3];
        for i = 1:size(areas_lats,1)   

            eh1(i) = errorbar([1:5]+(totimes(i)./20),toplot(:,i),forerrorextract(j).g(mon,:,i),...
                'linewidth',2,'LineStyle','none','color',cbrewqual2(i,:));
            hold on
            phl(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'MarkerSize',8,'MarkerEdgeColor',cbrewqual3(i,:),'MarkerFaceColor',cbrewqual2(i,:),'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));

            phl2(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},...
                'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'MarkerSize',10);
        end
        annotation('textbox',[sppos(1),sppos(2),sppos(3:4)],'String',labels{1},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+2,... % 
            'EdgeColor','none','fontweight','bold');    
        uistack(eh1,'bottom')
        uistack(phl2,'top')
        plot([0,6],[0,0],'linewidth',2,'LineStyle','--','color','k');
        set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
        
        xlabel('Month','fontsize',fsize+2);
        
        xlim([.5,5.5]);
        ylim([-80 130]);
        set(gca,'ytick',-80:20:130,'yticklabel',-80:20:130);
        ylabel('HSS','fontsize',fsize+2);
        title(titles{count},'fontsize',fsize+4)
        %March TCO - surface temperature ensemble mean HSS
        if count == 1
            lh = legend(phl,'Eastern Russia','Western Russia','Northern America','Southern Asia');                    
            set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
        end
        box on
        count = count+1;
        
        sp(count) = subplot(3,2,count);
        sppos = get(sp(count),'position');
        if count == 4            
            set(sp(count),'position',[sppos(1)-.075,sppos(2)+.020,sppos(3:4)]);
        else
            set(sp(count),'position',[sppos(1)-.075,sppos(2:4)]);
        end
        sppos = get(sp(count),'position');
        if same
            toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pctsame(j).g(mon,:,:,:),3));
        else
            toplot = squeeze(nanmean(GSSmonthtoplot_areamean_pct(j).g(mon,:,:,:),3));
        end
        for i = 1:size(areas_lats,1)   
            ehp1(i) = errorbar([1:5]+(totimes(i)./20),toplot(:,i),squeeze(forerrorextractpct(j).g(mon,:,i)),...
                'linewidth',2,'LineStyle','none','color',cbrewqual2(i,:));
            hold on
            php(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));    
            php2(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},...
                'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'MarkerSize',10);
        end
        uistack(ehp1,'bottom')
        uistack(php2,'top')
        plot([0,6],[0,0],'linewidth',2,'LineStyle','--','color','k');
        set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
        annotation('textbox',[sppos(1),sppos(2),sppos(3:4)],'String',labels{2},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+2,... % 
            'EdgeColor','none','fontweight','bold');    
        xlabel('Month','fontsize',fsize+2);
        
        %ylabel('HSS','fontsize',fsize+2);
        title(titles{count},'fontsize',fsize+4)
        xlim([.5,5.5]);
        ylim([-80 130]);
        set(gca,'ytick',-80:20:130,'yticklabel',-80:20:130);
        box on
        % plotting latitudes
        for i = 1:numel(lontoplot(j).g)
            if lontoplot(j).g(i) > 180
                lontoplot(j).g(i) = lontoplot(j).g(i) - 360;
            end
        end

        for i = 1:numel(lontoplot_pct(j).g)
            if lontoplot_pct(j).g(i) > 180
                lontoplot_pct(j).g(i) = lontoplot_pct(j).g(i) - 360;
            end
        end

        for i = 1:numel(areas_lons2)
            if areas_lons2(i) > 180
                areas_lons2(i) = areas_lons2(i) - 360;
            end
        end
        count = count+1;
    end

sp2 = subplot(3,1,3);

% fig = figure;
% set(fig,'position',[100 100 1000 400],'color','white');
m_proj('Equidistant Cylindrical','lon',[-180 180],'lat',[0 90]);
[~,h] = m_contourf(ERAdata.longitude-180,ERAdata.latitude,squeeze(nanmean(diffforplot(mon,:,:),1))',-5:.5:5,'LineStyle','none');
hold on
m_coast('color','k','LineWidth',1);
m_grid('ytick',[0:15:90],'xtick',[-180:60:360],'XaxisLocation','bottom','fontsize',fsize);  
caxis([-5 5]);
ch = colorbar;
set(ch,'YTick',[-5:1:5],'fontsize',fsize)%,
set(get(ch,'ylabel'),'string','Temperature','fontsize',fsize+2)
colormap(cbrewdiff);

ylab = ylabel('Latitude','fontsize',fsize+2);
ylabpos = get(ylab,'position');
set(ylab,'position',[ylabpos(1)-.25,ylabpos(2),ylabpos(3)]);
xlabel('Longitude','fontsize',fsize+2);
set(sp2,'position',[.14 .35 .65 .32]);
cbarrow
if corrtoplot
    title(['Location of maximum ensemble March TCO and ',monthnames(mon+2,0,0), ' TS differences'],'fontsize',fsize+2);
elseif difftoplot
    title(['Locations of ensemble March TCO and ',monthnames(mon+2,0,0), ' surface temperature differences'],'fontsize',fsize+2);
    
else
    title(['Location of maximum ensemble member March TCO and ',monthnames(mon+2,0,0), ' HSS'],'fontsize',fsize+2);
end
%a = 1
box on
hold on
for i = 1:size(areas_lats,1)
    %m_plot(lontoplot(:,i)',lattoplot(:,i)','LineStyle','none','Marker','o','MarkerSize',15,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',2)
    m_plot([areas_lons2(i,1),areas_lons2(i,2),areas_lons2(i,2),areas_lons2(i,1),areas_lons2(i,1)],[areas_lats(i,1),areas_lats(i,1),areas_lats(i,2),areas_lats(i,2),areas_lats(i,1)],...
        'LineStyle','-','LineWidth',4,'color',cbrewqual2(i,:),'LineStyle',lstyle{i})
    %mp(i) = m_plot(lontoplot(2).g(mon,:,i)',lattoplot(2).g(mon,:,i)','LineStyle','none','Marker',Mks{i},'MarkerSize',12,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',1);
    mp2(i) = m_plot(lontoplot(1).g(mon,:,i)',lattoplot(1).g(mon,:,i)','LineStyle','none','Marker',Mks{i},'MarkerSize',12,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor','k','LineWidth',2);
end
%axes.SortMethod='ChildOrder';

annotation('textbox',[.01 .995 .925 0],'String','March TCO - surface temperature HSSs','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');    

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/Obs_CompEmp_Ind_',corrtoplotext,'_','mean_Locations_withpct_test','_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),monthnames(mon+2,0,0),'both'];
print(filename,'-depsc');
end


end
