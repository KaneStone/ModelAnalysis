function [deviation_final, deviation_final_percent, deviationnorm,stdall] = constructdeviations(data,allyears,...
    yearstoexclude,deviationyear,vert,sa)

deviation_final = [];
deviation_final_percent = [];
% data must be in format of [month,year,latitude,pressure];

deviationindex = find(allyears == deviationyear);
otheryearsindex = 1:length(allyears);
for i = 1:length(yearstoexclude)    
    otheryearsindex(allyears == yearstoexclude(i)) = NaN;    
end
otheryearsindex(isnan(otheryearsindex)) = [];
if vert
    Meanall = nanmean(data,2);
    %Meanall = nanmean(data(:,1:end-2,:,:),2);
    Meanall = repmat(Meanall,[1, size(data,2),1,1]);
%     Mean1 = squeeze(nanmean(data(:,otheryearsindex(1:end-1),:,:),2));
%     Mean2 = squeeze(nanmean(data(:,otheryearsindex(1:end-1)+1,:,:),2));
    
    if sa
        stdall = nanstd(data,0,2);
        %stdall = nanstd(data(:,1:end-2,:,:),0,2);
        stdall = repmat(stdall,[1,size(data,2),1,1]);
%         std1 = squeeze(nanstd(data(:,otheryearsindex(1:end-1),:,:),0,2));
%         std2 = squeeze(std(data(:,otheryearsindex(1:end-1)+1,:,:),0,2));
    else
        stdall = 1;
    end
    
    deviationnorm = squeeze(data-Meanall)./stdall;
    %deviationnorm = squeeze(data./Meanall)-1;
    %deviationnorm2 = (squeeze(data(:,deviationindex+1,:,:))-Mean2)./std2;
    
%     deviation = (squeeze(data(:,deviationindex,:,:))-Mean1)./Mean1;
%     deviation2 = (squeeze(data(:,deviationindex+1,:,:))-Mean2)./Mean2;
%     deviation_percent = (squeeze(data(:,deviationindex,:,:))-Mean1)./Mean1*100;
%     deviation2_percent = (squeeze(data(:,deviationindex+1,:,:))-Mean2)./Mean2*100;
%     deviation_final = [deviation(1:12,:,:); deviation2(1:12,:,:)];
%     deviation_final_percent = [deviation_percent(1:12,:,:); deviation2_percent(1:12,:,:)];
    
else
%     Mean1 = squeeze(nanmean(data(:,otheryearsindex(1:end-1),:),2));
%     Mean2 = squeeze(nanmean(data(:,otheryearsindex(1:end-1)+1,:),2));
%     deviation = squeeze(data(:,deviationindex,:))-Mean1;
%     deviation2 = squeeze(data(:,deviationindex+1,:))-Mean2;
%     deviation_percent = (squeeze(data(:,deviationindex,:))-Mean1)./Mean1*100;
%     deviation2_percent = (squeeze(data(:,deviationindex+1,:))-Mean2)./Mean2*100;
%     deviation_final = [deviation(1:12,:); deviation2(1:12,:)];
%     deviation_final_percent = [deviation_percent(1:12,:); deviation2_percent(1:12,:)];

    Meanall = nanmean(data,2);
    Meanall = repmat(Meanall,[1, size(data,2),1]);
    stdall = nanstd(data,0,2);    
    stdall = repmat(stdall,[1,size(data,2),1]);    
    deviationnorm = squeeze(data-Meanall)./stdall;

end