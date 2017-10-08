function [data] = ReadinERA(directory)
% directory including file
info = ncinfo(directory);
for i = 1:length(info.Variables)
    data.(info.Variables(i).Name) = ncread(directory,info.Variables(i).Name);
end

end