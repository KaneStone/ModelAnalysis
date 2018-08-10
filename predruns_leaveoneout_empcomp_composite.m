function [] = predruns_leaveoneout_empcomp_composite(surfacedata,surfacedataind,tozdata,inputs,latitude,longitude,pct,differences,diffcomp,alldata)

% over 1995-2014. Leave one year out of the regression and see if I can
% predict the remaining years anomaly

% Do this for the ensemble mean and produce a map metric of predictability
% for the ensemble average.

% Using the high signal to noise plots, extract regions near these places
% for each 

% keep leaving years out

% find average temperature anomaly of upper 50th percentile of ozone
% find average temperature anomaly of lower 50th percentile of ozone
% two cases that can be used as a binary prediction (shift in pdf)

%% rearrange data into composite

%% rearrange data into composite

surfacedata_composite = permute(surfacedata,[3,4,2,1]);
surfacedata_composite = surfacedata_composite(:,:,:); 

for i = 1:9
    precondtoz(:,i) = squeeze(tozdata(i,inputs.tozmonth,:));%-nanmean(squeeze(tozdata(i,inputs.tozmonth,:)));    
end

tozdata_composite = precondtoz(:);

surfacedata = double(surfacedata);
surfacedataind = double(surfacedataind);

surfacedatacomposite = permute(surfacedata,[3,4,2,1]);
surfacedatacomposite = surfacedatacomposite(:,:,:);
surfacedatacomposite = permute(surfacedatacomposite,[3,1,2]);
surfacedatacompositeind = permute(surfacedataind,[3,4,5,2,1]);
surfacedatacompositeind = surfacedatacompositeind(:,:,:,:);
surfacedatacompositeind = permute(surfacedatacompositeind,[4,1,2,3]);

%% each individual separately 


%%
% pred = zeros(9,28,96,144);
% pred_months = zeros(9,28,96,5,144);
% pred2 = zeros(9,28,96,144);
% correct = zeros(9,28,96,144);
% correct2 = zeros(9,28,96,144);
% correct_months = zeros(9,28,96,5,144);


pred = zeros(270,96,144);
pred_months = zeros(270,96,5,144);
pred2 = zeros(270,96,144);
correct = zeros(270,96,144);
correct2 = zeros(270,96,144);
correct_months = zeros(270,96,5,144);


    tic;
    allin = 1:270;
    
    %tozautocorr = precondtoz(2:end-1,k);
%     tozautocorr = precondtoz(:,k);
upperpctind = find(tozdata_composite >= prctile(tozdata_composite,80));
lowerpctind = find(tozdata_composite <= prctile(tozdata_composite,20));
    count = 1;
    for i = 1:length(allin)
        leaveout = allin;        
        %leaveout (leaveout == i) = [];                         
        if i == 1
            leaveout (leaveout == count | leaveout == count+1) = []; 
        elseif i == 270
            leaveout (leaveout == count | leaveout == count-1) = []; 
        else
            leaveout (leaveout == count-1 | leaveout == count | leaveout == count+1) = [];         
        end
        
        predictors_anom = [leaveout',ones(size(leaveout))'];        
        
%         ozone_ind_b = predictors_anom\precondtoz(leaveout,k);        
%         ozone_anomaly = precondtoz(leaveout,k) - ozone_ind_b(1).*leaveout';% - ozone_ind_b(2); 
        
        %ozone_anomaly = repmat(ozone_anomaly,[1,length(longitude)]);                            
        %ozone_left(k,i) = precondtoz(count,1) - ozone_ind_b(1).*count;% - ozone_ind_b(2);    
        
        % calculating 50th percentile 
        ozone_pct(i) = prctile(tozdata_composite(leaveout),50);
        ozone_pct1(i) = prctile(tozdata_composite(leaveout),20);
        ozone_pct2(i) = prctile(tozdata_composite(leaveout),80);
        
%         ozone_pct(k,i) = prctile(ozone_anomaly,50);
%         ozone_pct1(k,i) = prctile(ozone_anomaly,20);
%         ozone_pct2(k,i) = prctile(ozone_anomaly,80);
        
        %ozone_pct(k,i) = prctile(ozone_anomaly,50);
        
        ozonelowerind(i).m = find(tozdata_composite(leaveout) <= ozone_pct1(i));
        ozoneupperind(i).m = find(tozdata_composite(leaveout) >= ozone_pct2(i));                      
        
%         ozonelowerind(i).m = find(ozone_anomaly <= ozone_pct(i));
%         ozoneupperind(i).m = find(ozone_anomaly >= ozone_pct(i));                      
        
        ozone_left(i) = tozdata_composite(count) - ozone_pct(i);
        %ozone_left(i) = ozone_left(i) - ozone_pct(i);
        %ozone_left(i) = precondtoz(i,1) - ozone_ind_b(1).*i - ozone_ind_b(2) - ozone_pct(i);    
        predsign(i) = sign(ozone_left(i));        
%         ozone_left = ozone_left - nanmedian(ozone_anomaly);
%         ozone_anomaly = ozone_anomaly - nanmedian(ozone_anomaly);
        
        %ozone_anomalysign = sign(ozone_anomaly);
                
        
        % read in ENSO 
%         ENSOyearind = predruns_removeENSOforpred(alldata(:,:,leaveout,:,:),latitude,longitude,0);
%         ENSOall = predruns_removeENSOforpred(alldata,latitude,longitude,0);
%         laglength = 3;        
        
        for j = 1:length(latitude)                           
            % empirical model 
%             for k = 1:9
%                 leaveoutind = [count-1,count,count+1]
%                 if leaveoutind <=
%                 
%             end
            b = predictors_anom\squeeze(surfacedatacomposite(leaveout,j,:));
            temp_anomaly = squeeze(surfacedatacomposite(leaveout,j,:)) - leaveout'.*b(1,:);% - b(2,:); %looks good
            temp_left = squeeze(surfacedatacomposite(count,j,:))' - count.*b(1,:);% - b(2,:); %looks good            
            
            temp_left = squeeze(temp_left) - squeeze(nanmedian(temp_anomaly));
            temp_anomaly = squeeze(temp_anomaly) - squeeze(nanmedian(temp_anomaly));

            temp_anomaly_upper = squeeze(nanmean(temp_anomaly(ozoneupperind(i).m,:)));
            temp_anomaly_lower = squeeze(nanmean(temp_anomaly(ozonelowerind(i).m,:)));            
                        
            if j == 84 && i > 1 && i < 30
                templow(i,:) = temp_anomaly(ozoneupperind(i).m,29);
                temphigh(i,:) = temp_anomaly(ozonelowerind(i).m,29);
            end
            
            signchange = sign(temp_anomaly_upper-temp_anomaly_lower);
            
            datasign(i,j,:) = sign(temp_left);
                  
            if j == 88 && i == 30
                abc = 1;
            end
            
            for l = 1:length(squeeze(datasign(i,j,:)))
%                 if signlower(i,j,l) < 0 
%                     pred(i,j,l) = 1;
%                 elseif signlower(i,j,l) > 0 
%                     pred(i,j,l) = -1;        
%                 end
                
                if signchange(l) == 1
                    pred(i,j,l) = 1;
                else
                    pred(i,j,l) = -1;
                end
                
%                 if pred(i,j,l) == 1
%                     if predsign(i) == datasign(i,j,l)
%                         correct(i,j,l) = 1;
%                     else
%                         correct(i,j,l) = 0;
%                     end
%                 elseif pred(i,j,l) == -1
%                     if predsign(i) == datasign(i,j,l)
%                         correct(i,j,l) = 0;
%                     else
%                         correct(i,j,l) = 1;
%                     end
%                 end
                
                if pred(i,j,l) == 1
                    presign2(i,j,l) = predsign(i);
                    if predsign(i) == datasign(i,j,l)
                        correct(i,j,l) = 1;
                    else
                        correct(i,j,l) = 0;
                    end
                elseif pred(i,j,l) == -1
                    presign2(i,j,l) = predsign(i).*-1;
                    if predsign(i) == datasign(i,j,l)
                        correct(i,j,l) = 0;
                    else
                        correct(i,j,l) = 1;
                    end
                end
                
            end  
            
            for m = 1:size(surfacedataind,3)
                b_months = predictors_anom\squeeze(surfacedatacompositeind(leaveout,m,j,:));
                temp_anomaly_months = squeeze(surfacedatacompositeind(leaveout,m,j,:)) - leaveout'.*b_months(1,:);% - b(2,:); %looks good
                temp_left_months = squeeze(surfacedatacompositeind(count,m,j,:))' - count.*b_months(1,:);% - b(2,:); %looks good                        
                temp_left_months = squeeze(temp_left_months) - squeeze(nanmedian(temp_anomaly_months,1));
                temp_anomaly_months = squeeze(temp_anomaly_months)' - squeeze(nanmedian(temp_anomaly_months,1))';
                temp_anomaly_upper_months = squeeze(nanmean(temp_anomaly_months(:,ozoneupperind(i).m),2));
                temp_anomaly_lower_months = squeeze(nanmean(temp_anomaly_months(:,ozonelowerind(i).m),2));                        
                
                signchange_months = sign(temp_anomaly_upper_months-temp_anomaly_lower_months);
                
                datasign_months(i,j,m,:) = sign(temp_left_months);
                
                for l = 1:length(squeeze(datasign_months(i,j,m,:)))
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
    toc;

%% emp

% extracting pct
% for i = 1:9
%     correctpct(i,:,:,:) = correct(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:);   
%     correct_months_pct(i,:,:,:,:) = correct_months(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:,:);
%         
%     predsignpct(i,:,:,:) = presign2(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:);
%     datasignpct(i,:,:,:) = datasign(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:);
%     predsignpct_months(i,:,:,:,:) = predsign_months2(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:,:);
%     datasignpct_months(i,:,:,:,:) = datasign_months(i,[pct.highCl.ind.highind(i,:),pct.highCl.ind.lowind(i,:)],:,:,:);

    correctpct = correct([upperpctind,lowerpctind],:,:);   
    correct_months_pct = correct_months([upperpctind,lowerpctind],:,:,:);
        
    predsignpct = presign2([upperpctind,lowerpctind],:,:);
    datasignpct = datasign([upperpctind,lowerpctind],:,:);
    predsignpct_months = predsign_months2([upperpctind,lowerpctind],:,:,:);
    datasignpct_months = datasign_months([upperpctind,lowerpctind],:,:,:);

% end

%%

GSS.all.ind = (sum(correct) - size(correct,1)./2)./...
    (size(correct,1)-size(correct,1)./2)*100;  
GSS.pct.ind = (sum(correctpct) - size(correctpct,1)./2)./...
    (size(correctpct,1)-size(correctpct,1)./2)*100;  
GSS.all.mean(1,:,:) = GSS.all.ind;
GSS.pct.mean(1,:,:) = GSS.pct.ind;

GSS.monthsall.ind = (sum(correct_months,1) - size(correct_months,1)./2)./...
    (size(correct_months,1)-size(correct_months,1)./2)*100;  
GSS.monthspct.ind = (sum(correct_months_pct,1) - size(correct_months_pct,1)./2)./...
    (size(correct_months_pct,1)-size(correct_months_pct,1)./2)*100;  
GSS.monthsall.mean(1,:,:,:) = GSS.monthsall.ind;
GSS.monthspct.mean(1,:,:,:) = GSS.monthspct.ind;
  
%% testing bootsrapping
mp = permute(presign2,[3,4,1,2]);
%mp = permute(modelprediction_ind-med,[3,4,1,2]);
mp = mp(:,:,:);
ad = permute(datasign,[3,4,1,2]);
%ad = permute(actualdata_ind-med,[3,4,1,2]);
ad = ad(:,:,:);

mpm = permute(predsign_months2,[4,3,5,1,2]);
mpm = mpm(:,:,:,:);
adm = permute(datasign_months,[4,3,5,1,2]);
adm = adm(:,:,:,:);


if ~exist('/Volumes/MyBook/work/data/predruns/output/HSS/empcomp_Allperc.mat','file')
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
    save('/Volumes/MyBook/work/data/predruns/output/HSS/empcomp_Allperc.mat','bootstat','bootstatm','percentiles','percentilesm');
else
    allsig = load('/Volumes/MyBook/work/data/predruns/output/HSS/empcomp_Allperc.mat');
end

p = zeros(size(GSS.all.mean));
p (GSS.all.mean < reshape(allsig.percentiles.ninetyfive,[1,size(allsig.percentiles.ninetyfive)])) = -1;

%%

mp = permute(predsignpct,[3,4,1,2]);
%mp = permute(modelprediction_ind-med,[3,4,1,2]);
mp = mp(:,:,:);
ad = permute(datasignpct,[3,4,1,2]);
%ad = permute(actualdata_ind-med,[3,4,1,2]);
ad = ad(:,:,:);

mpm = permute(predsignpct_months,[4,3,5,1,2]);
mpm = mpm(:,:,:,:);
adm = permute(datasignpct_months,[4,3,5,1,2]);
adm = adm(:,:,:,:);


if ~exist('/Volumes/MyBook/work/data/predruns/output/HSS/empcomp_Pctperc.mat','file')
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
    save('/Volumes/MyBook/work/data/predruns/output/HSS/empcomp_Pctperc.mat','bootstat','bootstatm','percentiles','percentilesm');
else
    allsigpct = load('/Volumes/MyBook/work/data/predruns/output/HSS/empcomp_Pctperc.mat');
end

p2 = zeros(size(GSS.all.mean));
p2 (GSS.pct.mean < reshape(allsigpct.percentiles.ninetyfive,[1,size(allsigpct.percentiles.ninetyfive)])) = -1;


%% plot     
titles = {'Ensemble composite HSSs (leave three out)','Ensemble composite HSSs during ozone extremes (leave three out)','Observed percentile general skill score (leave one out)'};
toplotp = permute(cat(1,p,p2),[1,3,2]);
[fig,fh] = subplotmaps(permute(cat(1,GSS.all.mean,GSS.pct.mean),[1,3,2]),longitude,latitude,{'seq','YlGnBu'},0,toplotp,16,titles,'Longitude','Latitude','HSS','on',...
        [0,50],11,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');         
 
% subplotmaps(permute(cat(1,GSS.all2.mean,GSS.pct2.mean),[1,3,2]),longitude,latitude,{'seq','YlGnBu'},0,[],16,titles,'Longitude','Latitude','','on',...
%         [0,50],11,[longitude(1:24:end)]-180,[longitude(1:24:end)]-180,[0:15:90],[0:15:90],{''},1,[0 360],[0 90],0,'none',1,'Miller Cylindrical');             
    
child = get(fig,'children');
axes1 = child(4);
axes2 = child(2);

axes(axes1);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','c','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',18,... % 
    'EdgeColor','none','fontweight','bold');    

axes(axes2);
sppos = get(gca,'position');
annotation('textbox',[sppos(1),sppos(2)+.01,sppos(3:4)],'String','d','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',18,... % 
    'EdgeColor','none','fontweight','bold');    


set(gcf,'Renderer','Painters');
filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/predictions/maps/CompEmp_Ensmean_GSSothercolor_nodetrend_upto80_',monthnames(inputs.varmonth,1,1),'_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),'.eps'];
print(filename,'-depsc');            
    
%% plot skill lines    
    %% plot skill lines    

%Find areas of large difference and large GSS
% difference above 3 degrees and GSS above 4
% areas_lons = [290,320;30,120;240,290;30,120];%lons (Greenland,Russia,America,Asia)
% areas_lats = [60,80;50,75;30,60;15,45];%lats

diffforplot = differences;
diffforplot = circshift(diffforplot,[0,144/2,0,0]);

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
%  areas_lons = [80,180;30,80;240,290;60,120];%lons (Greenland,East Russia, West Russia,America,Asia)
% areas_lats = [60,80;60,80;35,60;25,45];%lats

areas_lons = [85,180;30,85;240,290;30,120]; %lons (East Russia, West Russia,America,Asia)
areas_lats = [55,80;55,80;40,65;20,45]; %lats



for j = 1:2
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
            lats = latitude > areas_lats(i,1) & latitude < areas_lats(i,2);
            lons = longitude > areas_lons(i,1) & longitude < areas_lons(i,2);
            latextract = latitude(lats);
            lonextract = longitude(lons);
            [latmesh,lonmesh] = meshgrid(latextract,lonextract);

            latmesh = latmesh';
            lonmesh = lonmesh';
            for k = 1:size(diffcomp,1)           

                if corrtoplot               
                elseif difftoplot            
                    diff = permute(differences,[1,4,3,2]);           
                    diff2 = permute(diffcomp,[1,3,2]);
                    if same                
                        %mult = squeeze(nanmean(diff(:,mon,lats,lons),1));        
                        mult = squeeze(nanmean(diff2(:,lats,lons),1));        
                        mult_pct = squeeze(nanmean(diff(:,mon,lats,lons),1));        
                    else
                        mult = squeeze(diff(k,mon,lats,lons));        
                        mult_pct = squeeze(diff(k,mon,lats,lons));        
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
                    modelmonthsforerror = permute(squeeze(predsign_months2(:,lats,:,lons)),[1,3,2,4]);
                    modelmonthsforerror_extract(i,:,k,:) = modelmonthsforerror(:,:,maxind)';
                    %datamonthsforerror = squeeze(actualdata_ind_months(k,:,:,lats,lons));
                    datamonthsforerror = permute(squeeze(datasign_months(:,lats,:,lons)),[1,3,2,4]);
                    datamonthsforerror_extract(i,:,k,:) = datamonthsforerror(:,:,maxind)';
                    %medmonthsforerror = squeeze(medmonths(k,:,:,lats,lons));
                    %medmonthsforerror_extract(i,:,k,:) = medmonthsforerror(:,:,maxind)';

                    modelmonthsforerror_pct = permute(squeeze(predsignpct_months(:,lats,:,lons)),[1,3,2,4]);
                    modelmonthsforerror_pct_extract(i,:,k,:) = modelmonthsforerror_pct(:,:,maxind)';
                    datamonthsforerror_pct = permute(squeeze(datasignpct_months(:,lats,:,lons)),[1,3,2,4]);
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
        if ~same
            mp = modelmonthsforerror_extract(:,:,:);% - medmonthsforerror_extract(:,:,:);
            ap = datamonthsforerror_extract(:,:,:); %- medmonthsforerror_extract(:,:,:);
            mpp = modelmonthsforerror_pct_extract(:,:,:);% - medmonthsforerror_pct_extract(:,:,:);
            app = datamonthsforerror_pct_extract(:,:,:);% - medmonthsforerror_pct_extract(:,:,:);
            for x = 1:size(mp,1)       

                for i = 1:size(mp,2)

                    mpe = squeeze(mp(x,i,:));
                    ade = squeeze(ap(x,i,:));

                    mpme = squeeze(mpp(x,i,:));
                    adme = squeeze(app(x,i,:));


                    bootstatspec(x,i,:) = bootstrp(500, @(mpe) predruns_calcHeidke_forbootstrap_compemp(mpe,ade),mpe);            
                    forerrorextract(2).g(mon,i,x) = prctile(squeeze(bootstatspec(x,i,:)),95);

                    bootstatmspec(x,i,:) = bootstrp(500, @(mpme) predruns_calcHeidke_forbootstrap_compemp(mpme,adme),mpme);                
                    forerrorextractpct(2).g(mon,i,x) = prctile(squeeze(bootstatmspec(x,i,:)),95);

                end
            end
        end

    end
    same = 0;
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

titles = {'Ensemble composite','Ensemble composite, ozone extremes','Ensemble members','Ensemble members, ozone extremes'}; 
labels = {'c','d','a','b'};
for mon = 1:5
    count = 1;
    fig = figure;
    set(fig,'position',[100 1 1000 984],'color','white','Visible','on');
    for j = 1:2
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
            phltemp(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'MarkerSize',8,'MarkerEdgeColor',cbrewqual3(i,:),'MarkerFaceColor',cbrewqual2(i,:),'linewidth',2,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));
            phl(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},'MarkerSize',8,'MarkerEdgeColor',cbrewqual3(i,:),'MarkerFaceColor',cbrewqual2(i,:),'linewidth',3,'LineStyle',lstyle{i},'color',cbrewqual2(i,:));

            phl2(i) = plot([1:5]+(totimes(i)./20),toplot(:,i),Mks{i},...
                'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'MarkerSize',10);
        end
        uistack(eh1,'bottom')
        uistack(phl2,'top')
        plot([0,6],[0,0],'linewidth',2,'LineStyle','--','color','k');
        set(gca,'fontsize',fsize,'xtick',1:5,'xticklabel',xticklab)
        if count == 3
            xlabel('Month','fontsize',fsize+2);
        end
        xlim([.5,5.5]);
        ylim([-30 80]);
        set(gca,'ytick',-20:10:100,'yticklabel',-20:10:100);
        ylabel('HSS','fontsize',fsize+2);
        title(titles{count},'fontsize',fsize+4)
        %March TCO - surface temperature ensemble mean HSS
        if count == 1
            lh = legend(phltemp,'eastern Russia','western Russia','North America','southern Asia');                    
            set(lh,'box','off','fontsize',fsize-2,'location','NorthEast')
        end
        annotation('textbox',sppos,'String',labels{count},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+2,... % 
            'EdgeColor','none','fontweight','bold');    
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
        if count == 4
            xlabel('Month','fontsize',fsize+2);
        end
         annotation('textbox',sppos,'String',labels{count},'FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+2,... % 
            'EdgeColor','none','fontweight','bold');    
        %ylabel('HSS','fontsize',fsize+2);
        title(titles{count},'fontsize',fsize+4)
        xlim([.5,5.5]);
        ylim([-30 80]);
        set(gca,'ytick',-20:10:100,'yticklabel',-20:10:100);
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
[~,h] = m_contourf(longitude-180,latitude,squeeze(nanmean(diffforplot(:,:,:,mon),1))',-5:.5:5,'LineStyle','none');
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
set(sp2,'position',[.14 .075 .65 .32]);
sppos = get(sp2,'position');
annotation('textbox',[sppos(1),sppos(2)-.04,sppos(3:4)],'String','c','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','left','fontsize',fsize+2,... % 
            'EdgeColor','none','fontweight','bold');    
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
    mp(i) = m_plot(lontoplot(2).g(mon,:,i)',lattoplot(2).g(mon,:,i)','LineStyle','none','Marker',Mks{i},'MarkerSize',12,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor',cbrewqual3(i,:),'LineWidth',1);
    mp2(i) = m_plot(lontoplot(1).g(mon,:,i)',lattoplot(1).g(mon,:,i)','LineStyle','none','Marker','p','MarkerSize',20,'MarkerFaceColor',cbrewqual2(i,:),'MarkerEdgeColor','k','LineWidth',2);
end
%axes.SortMethod='ChildOrder';

annotation('textbox',[.01 .995 .925 0],'String','March TCO - surface temperature ensemble composite HSSs','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+6,... % 
        'EdgeColor','none','fontweight','bold');    

filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MyPapers/Mine/OzonePred/Draft/Figures/CompEmp_Ind_',corrtoplotext,'_','mean_Locations_withpct_test','_',num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2)),monthnames(mon+2,0,0),'both'];
print(filename,'-depsc');
end
end