function [years] = CCMI_years(date,ifwac)

% function to extract only year vector from CESM CCMI date.

temp = num2str(date);
for j = 1:length(temp)
    years(j) = str2double(temp(j,1:4));
end
clearvars temp
if ifwac
    years = circshift(years,[0,1]);
    years(1) = years(2);
end

end
