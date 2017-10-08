function [data, monthyear, SBUVTCO, yearmean, yearmeananomaly] = ReadSBUV(filename,lats,anomyear,months)

fid = fopen(filename);
keepgoing  = 1;
while keepgoing
    line = fgetl(fid); 
    if strcmp(line,'Jan 1970');
        monthyear{1} = line;
        keepgoing = 0;
        headers = fgetl(fid);
        formatspec = '%f %f %f %f %f %f %f %f %f %f %f %f %f';
        for i = 1:(2015-1970)*12+12
            data(:,:,i) = fscanf(fid,formatspec,[13,36])';
            temp = fgetl(fid);
            monthyear{i+1} = fgetl(fid);
            temp = fgetl(fid);            
        end            
    else keepgoing = 1;
    end
end

SBUVrefstart = char(monthyear{1, 1});
SBUVrefend = char(monthyear{1, end-1});
SBUVref(1) = str2double(SBUVrefstart(5:8));
SBUVref(2) = str2double(SBUVrefend(5:8));

SBUVTCO = squeeze(data(:,13,:));
Latitudes = data(:,1:2,1);
SBUVTCO (SBUVTCO == 0) = NaN;
Latitudes_middle = (Latitudes(:,1) + Latitudes(:,2))/2;

[yearmeananomaly,yearmean] = TCOanomaly(SBUVTCO,lats,Latitudes_middle,anomyear,SBUVref,months);

% figure;
% plot(yearmean)
% hold on
% plot(1:50,zeros(1,50));
% ylim([-5 11]);
% set(gca,'ytick',-5:1:20,'yticklabel',-5:1:20,'xtick',1:1:50,'xticklabel',1970:1:2030);

end


