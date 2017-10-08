function [monthout] = monthnames(monthinput,stringtogether,shortnames)

if shortnames
    monthtitles = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',...
        'Sep','Oct','Nov','Dec'};
else
    monthtitles = {'January','February','March','April','May','June','July','August',...
        'September','October','November','December'};
end

for i = 1:length(monthinput)
    monthout{i} = monthtitles{monthinput(i)};
end
if length(monthinput) == 1
    monthout = [monthout{:}];
else
    if stringtogether
        monthout = [monthout{:}];
    end
end
end