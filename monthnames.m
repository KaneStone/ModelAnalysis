function [monthout] = monthnames(monthinput,stringtogether,type)

ext = [];
for i = 1:length(monthinput)
    if monthinput(i) > 12
        monthinput(i) = monthinput(i)-12; 
        ext = '+';    
    end    
end

if strcmp(type,'short')
    monthtitles = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',...
        'Sep','Oct','Nov','Dec'};
elseif strcmp(type,'long')
    monthtitles = {'January','February','March','April','May','June','July','August',...
        'September','October','November','December'};
elseif strcmp(type,'single')
    monthtitles = {'J','F','M','A','M','J','J','A',...
        'S','O','N','D'};
end

for i = 1:length(monthinput)
    monthout{i} = monthtitles{monthinput(i)};
end
if length(monthinput) == 1
    monthout = [monthout{:},ext];
else
    if stringtogether
        monthout = [monthout{:},ext];
    end
end
end
