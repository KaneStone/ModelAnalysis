function [dataVarMonthAve,dataVarMonth,NINO34all,differences,correlations] = predruns_diffandcorr_percentiles(...
    dataMonthArrange,toz_dataMonthArrange,lats,lons,Eachyear,pct,ClLevel,inputs)

% This function calculates the ensemble mean and individual correlation and percentile
% differeces

[dataVarMonthAve,dataVarMonth,NINO34all] = predruns_removeENSO(dataMonthArrange,lats,lons,inputs,ClLevel);


%% extracting upper and lower percentiles from the lowest and highest two simulations from each year
clearvars varlower varupper

for j = 1:size(Eachyear.(ClLevel).lowindex,3)
     varlower(:,j,:,:) = dataVarMonthAve(Eachyear.(ClLevel).lowindex(:,inputs.tozmonth,j),j,:,:);
     varupper(:,j,:,:) = dataVarMonthAve(Eachyear.(ClLevel).highindex(:,inputs.tozmonth,j),j,:,:);
     
     varlower_indmonths(:,j,:,:,:) = dataVarMonth(Eachyear.(ClLevel).lowindex(:,inputs.tozmonth,j),j,:,:,:);
     varupper_indmonths(:,j,:,:,:) = dataVarMonth(Eachyear.(ClLevel).highindex(:,inputs.tozmonth,j),j,:,:,:);
     
     toztemp(:,:,j) = [squeeze(toz_dataMonthArrange(Eachyear.(ClLevel).highindex(:,inputs.tozmonth,j),inputs.tozmonth,j)),...
         squeeze(toz_dataMonthArrange(Eachyear.(ClLevel).lowindex(:,inputs.tozmonth,j),inputs.tozmonth,j))];                          
end

lower_ra = permute(varlower,[3,4,2,1]);
lower_ra = lower_ra(:,:,:);
upper_ra = permute(varupper,[3,4,2,1]);
upper_ra = upper_ra(:,:,:);

lower_ra_indmonths = permute(varlower_indmonths,[4,5,3,2,1]);
lower_ra_indmonths = lower_ra_indmonths(:,:,:,:);
upper_ra_indmonths = permute(varupper_indmonths,[4,5,3,2,1]);
upper_ra_indmonths = upper_ra_indmonths(:,:,:,:);

varcomp = cat(3,upper_ra,lower_ra);
varcomp_indmonths = cat(4,upper_ra_indmonths,lower_ra_indmonths);
tozcomp = cat(3,squeeze(toztemp(:,1,:)),squeeze(toztemp(:,2,:)));
tozcomp = tozcomp(:);

%% extracting upper and lower percentiles from the 20th percentile of the composite
dataVarMonthAve_temp = permute(dataVarMonthAve,[4,3,2,1]);
dataVarMonthAve_temp = dataVarMonthAve_temp(:,:,:,:);

dataVarMonth_temp = permute(dataVarMonth,[5,4,3,2,1]);
%dataVarMonth_temp = dataVarMonthAve_temp(:,:,:,:);

var_composite_lower = dataVarMonthAve_temp(:,:,pct.(ClLevel).ens.lowind.a);
var_composite_upper = dataVarMonthAve_temp(:,:,pct.(ClLevel).ens.highind.a);

var_composite_ind_lower = dataVarMonth_temp(:,:,:,pct.(ClLevel).ens.lowind.a);
var_composite_ind_upper = dataVarMonth_temp(:,:,:,pct.(ClLevel).ens.highind.a);

toz_dataMonthArrange_ra = squeeze(toz_dataMonthArrange(:,inputs.tozmonth,:))';
toz_dataMonthArrange_ra = toz_dataMonthArrange_ra(:); 
toz_composite_lower = toz_dataMonthArrange_ra(pct.(ClLevel).ens.lowind.a);
toz_composite_upper = toz_dataMonthArrange_ra(pct.(ClLevel).ens.highind.a);
toz_composite = [toz_composite_upper;toz_composite_lower];
%% extracing upper and lower percentiles for each simulation separately

for i = 1:size(dataVarMonthAve,1)
    var_individual_lower(i,:,:,:) = dataVarMonthAve(i,pct.(ClLevel).ind.lowind(i,:),:,:);
    var_individual_upper(i,:,:,:) = dataVarMonthAve(i,pct.(ClLevel).ind.highind(i,:),:,:);    
    toz_individual_lower(:,i) = toz_dataMonthArrange(i,inputs.tozmonth,pct.(ClLevel).ind.lowind(i,:));
    toz_individual_upper(:,i) = toz_dataMonthArrange(i,inputs.tozmonth,pct.(ClLevel).ind.highind(i,:));
    var_individual_lower_indmonths(i,:,:,:,:) = dataVarMonth(i,pct.(ClLevel).ind.lowind(i,:),:,:,:);
    var_individual_upper_indmonths(i,:,:,:,:) = dataVarMonth(i,pct.(ClLevel).ind.highind(i,:),:,:,:);    
end

%% Calculating differences

differences.eachyear(1,:,:) = [nanmean(upper_ra,3) - nanmean(lower_ra,3)]';
differences.indmonths.eachyear(:,:,:) = permute(nanmean(upper_ra_indmonths,4) - nanmean(lower_ra_indmonths,4),[3,2,1]);
differences.composite(1,:,:) = nanmean(var_composite_upper,3) - nanmean(var_composite_lower,3);
differences.indcomposite(1,:,:,:) = nanmean(var_composite_ind_upper,4) - nanmean(var_composite_ind_lower,4);
differences.individual = permute(squeeze(nanmean(var_individual_upper,2)) - ...
    squeeze(nanmean(var_individual_lower,2)),[1,3,2]);
differences.ensmean = nanmean(permute(squeeze(nanmean(var_individual_upper,2)) - ...
    squeeze(nanmean(var_individual_lower,2)),[1,3,2]),1);
differences.indmonths.individual = permute(squeeze(nanmean(var_individual_upper_indmonths,2)) - ...
    squeeze(nanmean(var_individual_lower_indmonths,2)),[1,4,3,2]);


forsig_ind_upper = permute(var_individual_upper,[3,4,1,2]);
forsig_ind_lower = permute(var_individual_lower,[3,4,1,2]);
forsig_ind_upper = forsig_ind_upper(:,:,:); 
forsig_ind_lower = forsig_ind_lower(:,:,:);

for i = 1:size(upper_ra,1)
    for j = 1:size(upper_ra,2)
        differences.ttest.eachyear(1,j,i) = ttest2(squeeze(upper_ra(i,j,:)),...
            squeeze(lower_ra(i,j,:)));
        differences.ttest.composite(1,j,i) = ttest2(squeeze(var_composite_upper(j,i,:)),...
            squeeze(var_composite_lower(j,i,:)));        
        differences.ttest.ensmean(1,j,i) = ttest2(squeeze(forsig_ind_upper(i,j,:)),...
            squeeze(forsig_ind_lower(i,j,:)));        
    end
end

%% calculating correlations

dataVarMonthAve_ra = permute(dataVarMonthAve,[3,4,2,1]);
dataVarMonthAve_ra = dataVarMonthAve_ra(:,:,:);

%construct correlation percentile
lowerind_comp = permute(var_individual_lower,[4,3,2,1]);
lowerind_comp = lowerind_comp(:,:,:);
upperind_comp = permute(var_individual_upper,[4,3,2,1]);
upperind_comp = upperind_comp(:,:,:);
comp_pct = cat(3,upperind_comp,lowerind_comp);

toz_comp_pct = [toz_individual_upper(:);toz_individual_lower(:)];
for i = 1:size(varcomp,1)
    i
    for j = 1:size(varcomp,2)
        %composite percentile
        [correlations.composite_pct(1,j,i),correlations.sig.composite_pct(1,j,i)] = ...
            corr([squeeze(var_composite_upper(j,i,:));squeeze(var_composite_lower(j,i,:))],toz_composite);
        
        %individual percentile mean
        [correlations.individualcorr_pct(1,j,i),correlations.sig.individualcorr_pct(1,j,i)] = ...
            corr(squeeze(comp_pct(j,i,:)),toz_comp_pct);
        
        %mean for upper and lower 2 runs from each
        %year
        [correlations.eachyear_pct(1,j,i),correlations.sig.eachyear_pct(1,j,i)] = ...
            corr(squeeze(varcomp(i,j,:)),tozcomp);
        
        %all data in correlation     
        [correlations.composite_all_pct(1,j,i),correlations.sig.composite_all_pct(1,j,i)] = ...
            corr(squeeze(dataVarMonthAve_ra(i,j,:)),toz_dataMonthArrange_ra);
        
%         if length(inputs.varmonth) > 1
%             [correlations.indmonths.composite_pct(1,j,i),correlations.indmonths.sig.composite_pct(1,j,i)] = ...
%                 corr([squeeze(var_composite_upper(j,i,:));squeeze(var_composite_lower(j,i,:))],toz_composite);
%             [correlations.indmonths.eachyear_pct(1,j,i),correlations.indmonths.sig.eachyear_pct(1,j,i)] = ...
%                 corr(squeeze(varcomp(i,j,:)),tozcomp);
%             [correlations.indmonths.composite_all_pct(1,j,i),correlations.indmonths.sig.composite_all_pct(1,j,i)] = ...
%                 corr(squeeze(dataVarMonthAve_ra(i,j,:)),toz_dataMonthArrange_ra);
%         end
        
        
        for k = 1:size(dataVarMonthAve,1)
            %individual correlations
            [correlations.individual(k,j,i),correlations.sig.individual(k,j,i)] = ...
                corr(squeeze(dataVarMonthAve(k,:,i,j))',squeeze(toz_dataMonthArrange(k,inputs.tozmonth,:)));
            %individual correlations (percentiles)
            [correlations.individual_pct(k,j,i),correlations.sig.individual_pct(k,j,i)] = ...
                corr(squeeze(dataVarMonthAve(k,[pct.(ClLevel).ind.highind(k,:),pct.(ClLevel).ind.lowind(k,:)],i,j))',...
                squeeze(toz_dataMonthArrange(k,inputs.tozmonth,[pct.(ClLevel).ind.highind(k,:),pct.(ClLevel).ind.lowind(k,:)])));
            for l = 1:size(dataVarMonth,3)
                %individual correlations and months                
                [correlations.indmonths.individual(k,l,j,i),correlations.indmonths.sig.individual(k,l,j,i)] = ...
                    corr(squeeze(dataVarMonth(k,:,l,i,j))',squeeze(toz_dataMonthArrange(k,inputs.tozmonth,:)));
                %individual correlations and months (percentiles)
                [correlations.indmonths.individual_pct(k,l,j,i),correlations.indmonths.sig.individual_pct(k,l,j,i)] = ...
                    corr(squeeze(dataVarMonth(k,[pct.(ClLevel).ind.highind(k,:),pct.(ClLevel).ind.lowind(k,:)],l,i,j))',...
                    squeeze(toz_dataMonthArrange(k,inputs.tozmonth,[pct.(ClLevel).ind.highind(k,:),pct.(ClLevel).ind.lowind(k,:)])));
            end
        end
    end
end

correlations.corrmean = nanmean(correlations.individual,1);
correlations.corrmean_pct = nanmean(correlations.individual_pct,1);
correlations.indmonths.corrmean = nanmean(correlations.indmonths.individual,1);
correlations.indmonths.corrmean_pct = nanmean(correlations.indmonths.individual_pct,1);
correlations.sig.corrmean = nanmean(correlations.sig.individual,1);
correlations.sig.corrmean_pct = nanmean(correlations.sig.individual_pct,1);
correlations.indmonths.sig.corrmean = nanmean(correlations.indmonths.sig.individual,1);
correlations.indmonths.sig.corrmean_pct = nanmean(correlations.indmonths.sig.individual_pct,1);

%% calculate individual standard deviations of correlations and percentile differences

correlations.indstd = std(correlations.individual,0,1);
correlations.indstd_pct = std(correlations.individual_pct,0,1);
differences.indstd = std(differences.individual,0,1);

correlations.indmonths.indstd = std(correlations.indmonths.individual,0,1);
correlations.indmonths.indstd_pct = std(correlations.indmonths.individual_pct,0,1);
differences.indmonths.indstd = std(differences.indmonths.individual,0,1);

end
