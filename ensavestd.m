function [] = ensavestd(data)
clearvars anom
for i = 1:length(data)-2
    for j = 1:12
        anom(:,:,:,i) = data(10).percent - data(i).percent;
        
    end
end

anomstd_ensave = std(data(11).percent_residuals_months,0,1);

anom = permute(anom,[2,3,1,4]);
anom = permute(anom(:,:,:),[3,1,2]);
anomstd = std(anom,0,1);
plot(anom(:,27,85));
%%
% for i = 1:10
%     plot(data(i).percent(12:12:end,27,85),'LineWidth',2);
%     hold on
% end

end
