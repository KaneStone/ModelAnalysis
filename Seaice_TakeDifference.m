function [output] = Seaice_TakeDifference(inputs,mons,surfacedata,nofiles,pct,observations,SSWSPVcomposite,sig)

% This extracts the low and high values for the input months depending on
% whether we are using wind or ozone extremes
% - coincidence is coincident ozone and wind extremes

if inputs.nocoincidence || inputs.coincidence
    load(['/Volumes/ExternalOne/work/data/predruns/output/TCOwindextremes/coincidence_wind','40','and','0','.mat'])
    ozoneind.upper.noSSW(1) = [];
    ozoneind.lower.noSPV(1) = [];
    ozoneind.upper.SSW(1) = [];
    ozoneind.lower.SPV(1) = [];
end

if inputs.removeENSO
    data = surfacedata.highCl.dataMonthArrange_re;
else
    data = surfacedata.highCl.dataMonthArrange;
end
%mons = mons - (inputs.varmonth(1) -1);
for l = 1:length(mons)
    
    for i = 1:nofiles


        if inputs.compareSSW
            extremeslow(i).w(:,:,:) = permute(squeeze(nanmean(...
                data(i,mons(l,:),...
                SSWSPVcomposite.model.SPV(i+1,:),:,:),2)),[2,3,1]); %+1 to account for only 9 toz files
            extremeshigh(i).w(:,:,:) = permute(squeeze(nanmean(...
                data(i,mons(l,:),SSWSPVcomposite.model.SSW(i+1,:),:,:),2)),[2,3,1]);        
        end

        if ~inputs.coincidence && ~inputs.nocoincidence
            % ozone only
            extremeslow(:,:,:,l,i) = permute(squeeze(nanmean(...
                data(i,mons(l,:),pct.highCl.ind.lowind(i,:),:,:),2)),[2,3,1]);
            extremeshigh(:,:,:,l,i) = permute(squeeze(nanmean(...
                data(i,mons(l,:),pct.highCl.ind.highind(i,:),:,:),2)),[2,3,1]);        
        elseif inputs.nocoincidence
            if length(ozoneind.lower.noSPV(i).a) == 1
                extremeslow(i).a(:,:,:,l) = squeeze(nanmean(...
                    data(i,mons(l,:),ozoneind.lower.noSPV(i).a,:,:),2));    
            elseif length(ozoneind.lower.noSPV(i).a) > 1
                extremeslow(i).a(:,:,:,l) = permute(squeeze(nanmean(...
                    data(i,mons(l,:),ozoneind.lower.noSPV(i).a,:,:),2)),[2,3,1]);   
            end
            if ~isempty(ozoneind.upper.noSSW(i).a)
                if length(ozoneind.upper.noSSW(i).a) == 1
                    extremeshigh(i).a(:,:,:,l) = squeeze(nanmean(...
                        data(i,mons(l,:),ozoneind.upper.noSSW(i).a,:,:),2));         
                elseif length(ozoneind.upper.noSSW(i).a) > 1
                    extremeshigh(i).a(:,:,:,l) = permute(squeeze(nanmean(...
                        data(i,mons(l,:),ozoneind.upper.noSSW(i).a,:,:),2)),[2,3,1]);   
                end                                            
            else
                extremeshigh(i).a = [];
            end
        elseif inputs.coincidence
            if length(ozoneind.lower.SPV(i).a) == 1
                extremeslow(i).a(:,:,:,l) = squeeze(nanmean(...
                    data(i,mons(l,:),ozoneind.lower.SPV(i).a,:,:),2));    
            elseif length(ozoneind.lower.SPV(i).a) > 1
                extremeslow(i).a(:,:,:,l) = permute(squeeze(nanmean(...
                    data(i,mons(l,:),ozoneind.lower.SPV(i).a,:,:),2)),[2,3,1]);   
            end
            if ~isempty(ozoneind.upper.SSW(i).a)
                if length(ozoneind.upper.SSW(i).a) == 1
                    extremeshigh(i).a(:,:,:,l) = squeeze(nanmean(...
                        data(i,mons(l,:),ozoneind.upper.SSW(i).a,:,:),2));         
                elseif length(ozoneind.upper.SSW(i).a) > 1
                    extremeshigh(i).a(:,:,:,l) = permute(squeeze(nanmean(...
                        data(i,mons(l,:),ozoneind.upper.SSW(i).a,:,:),2)),[2,3,1]);   
                end                                            
            else
                extremeshigh(i).a = [];
            end
        end

    end
%     %%
%     if inputs.coincidence
%         output.extremeslow = permute(cat(1,extremeslowtemp(:).a),[2,3,1]); 
%         output.extremeshigh = permute(cat(1,extremeshightemp(:).a),[2,3,1]);
%     end

end

% there could be an issue with detrending, so I might need to implement
% this here.
%lowmean (lowmean < 0) = 0;
%highmean (highmean < 0) = 0;

% take difference as high - low
if inputs.nocoincidence || inputs.coincidence
    
%     output.extremesLow.ind = squeeze(nanmean(extremeslow,3)); % lat lon months members
%     output.extremesHigh.ind = squeeze(nanmean(extremeshigh,3));
    output.extremesLow.ens = squeeze(nanmean(cat(3,extremeslow(:).a),3));
    output.extremesHigh.ens = squeeze(nanmean(cat(3,extremeshigh(:).a),3));
%     output.difference.ind = output.extremesHigh.ind - output.extremesLow.ind;
    output.difference.ens = output.extremesHigh.ens - output.extremesLow.ens;
else
    output.extremesLow.ind = squeeze(nanmean(extremeslow,3)); % lat lon months members
    output.extremesHigh.ind = squeeze(nanmean(extremeshigh,3));
    output.extremesLow.ens = nanmean(output.extremesLow.ind,4);
    output.extremesHigh.ens = nanmean(output.extremesHigh.ind,4);
    output.difference.ind = output.extremesHigh.ind - output.extremesLow.ind;
    output.difference.ens = output.extremesHigh.ens - output.extremesLow.ens;
    %output.difference.ens = nanmean(output.extremesHigh.ind - output.extremesLow.ind,4);
end
%output.difference.ens (isnan(output.difference.ens)) = 0;
%output.difference.ind (isnan(output.difference.ind)) = 0;

% calculate significance
if inputs.nocoincidence || inputs.coincidence
    extremeslowcomp = cat(3,extremeslow(:).a);
    extremeshighcomp = cat(3,extremeshigh(:).a);
    output.difference.ttest = squeeze(ttest2(permute(extremeslowcomp,[3,1,2,4]),...
    permute(extremeshighcomp,[3,1,2,4]),'Alpha',sig));
else
    for i = 1:length(mons)                
        extremeslowcomptemp = squeeze(extremeslow(:,:,:,i,:));
        extremeslowcomp(:,:,:,i) = extremeslowcomptemp(:,:,:);
        extremeshighcomptemp = squeeze(extremeshigh(:,:,:,i,:));
        extremeshighcomp(:,:,:,i) = extremeshighcomptemp(:,:,:);
    end    
    
    % interpolate significance
%     [X,Y] = meshgrid(surfacedata.highCl.data(1).lon,surfacedata.highCl.data(1).lat);
%     [Xq,Yq] = meshgrid(surfacedata.highCl.data(1).lon(1):.25:surfacedata.highCl.data(1).lon(end),...
%         surfacedata.highCl.data(1).lat(1):.25:surfacedata.highCl.data(1).lat(end));
%     output.interplon = Xq(1,:)';
%     output.interplat = Yq(:,1);
%     for i = 1:size(extremeslowcomp,3)
%         for j = 1:size(extremeslowcomp,4)
%             extremeslowcomp_int(:,:,i,j) = interp2(X,Y,extremeslowcomp(:,:,i,j),Xq,Yq,'nearest');
%             extremeshighcomp_int(:,:,i,j) = interp2(X,Y,extremeshighcomp(:,:,i,j),Xq,Yq,'nearest');
%         end
%     end
    output.difference.ttest = squeeze(ttest2(permute(extremeslowcomp,[3,1,2,4]),...
    permute(extremeshighcomp,[3,1,2,4]),'Alpha',sig));



end



%% observations

%% extract same time period as toz or wind
%observationsextract = observations.montharrange(:,observations.timeperiod >= observations.toz.timeperiod(1) & observations.timeperiod <= observations.toz.timeperiod(2),:,:);
if inputs.compareERA
    if inputs.removeENSO
        observationsextract = observations.montharrange_re;
    else
        observationsextract = observations.montharrange;
    end
    %observationsextract_modeldim = observations.montharrange_modeldim(:,observations.timeperiod >= observations.toz.timeperiod(1) & observations.timeperiod <= observations.toz.timeperiod(2),:,:);
    if inputs.seasons
        for l = 1:length(mons)
            output.observations.lowdata(l,:,:,:) = nanmean(observationsextract(mons(l,:),observations.toz.lowind,:,:));
            output.observations.highdata(l,:,:,:) = nanmean(observationsextract(mons(l,:),observations.toz.highind,:,:));    
        end    
        output.observations.lowdatamean  = squeeze(nanmean(output.observations.lowdata,2));
        output.observations.highdatamean  = squeeze(nanmean(output.observations.highdata,2));
        output.observations.difference = output.observations.highdatamean - output.observations.lowdatamean;
    else
        for l = 1:length(mons)
            output.observations.lowdata(l,:,:,:) = nanmean(observationsextract(mons(l,:),observations.toz.lowind,:,:));
            output.observations.highdata(l,:,:,:) = nanmean(observationsextract(mons(l,:),observations.toz.highind,:,:));    
            output.observations.difference = output.observations.highdata - output.observations.lowdata;
        end    
    end

%     [X,Y] = meshgrid(observations.latitude,observations.longitude);
%     [Xq,Yq] = meshgrid(observations.latitude(1):-.25:observations.latitude(end),...
%         observations.longitude(1):.25:observations.longitude(end));
%     output.observations.interplat = Xq(1,:)';
%     output.observations.interplon = Yq(:,1);
%     for i = 1:size(output.observations.lowdata,1)
%         for j = 1:size(output.observations.lowdata,2)
%             obslowcomp_int(i,j,:,:) = interp2(X,Y,squeeze(output.observations.lowdata(i,j,:,:)),Xq,Yq,'nearest');
%             obshighcomp_int(i,j,:,:) = interp2(X,Y,squeeze(output.observations.highdata(i,j,:,:)),Xq,Yq,'nearest');
%         end
%     end

    output.observations.ttest = squeeze(ttest2(permute(output.observations.lowdata,[2,1,3,4]),...
     permute(output.observations.highdata,[2,1,3,4]),'Alpha',sig));
else
    output.observations.ttest = zeros(3,240,121);
    output.observations.difference = zeros(3,240,121);
    output.observations.latitude = 90:-1.5:-90;
    output.observations.longitude = 0:1.5:360-1.5;
end

% output.observations.ttest = squeeze(ttest2(permute(obslowcomp_int,[2,1,3,4]),...
%  permute(obshighcomp_int,[2,1,3,4]),'Alpha',sig));

% save observations and model differences

if inputs.seasons && strcmp(inputs.var,'SNOWHICE')
    observeddifferences = output.observations.difference;
    modelensdifferences = output.difference.ens;
    if ~exist('/Volumes/ExternalOne/work/data/predruns/output/data/detrendedENSOremoved_differences.mat')
        save('/Volumes/ExternalOne/work/data/predruns/output/data/detrendedENSOremoved_differences.mat','observeddifferences','modelensdifferences');
    end
end

end