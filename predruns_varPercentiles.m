function [pct,dataout,Eachyear] = predruns_varPercentiles(data,datamontharrange,monthin,percentile,nofiles)

% Using the percentile function - calulates the percentiles and associated
% index locations of the upper and lower user input: percentile

%%

yearlength = size(data,2)/nofiles;

% take for ensmeble composite
count = 1;
for i = monthin
    pct.ens.lowpercentile(count) = prctile(data(monthin,:),percentile);
    pct.ens.highpercentile(count) = prctile(data(monthin,:),100-percentile);
    pct.ens.lowind(count).a = find(data(monthin,:) <= pct.ens.lowpercentile(count));
    pct.ens.highind(count).a = find(data(monthin,:) >= pct.ens.highpercentile(count));  
    count = count+1;
end

% take for each simulation seperately
for i = 1:nofiles
    pct.ind.lowpercentile(i) = prctile(datamontharrange(i,monthin,:),percentile);
    pct.ind.highpercentile(i) = prctile(datamontharrange(i,monthin,:),100-percentile);
    pct.ind.lowind(i,:) = find(datamontharrange(i,monthin,:) <= pct.ind.lowpercentile(i));
    pct.ind.highind(i,:) = find(datamontharrange(i,monthin,:) >= pct.ind.highpercentile(i));  
end
    

for i = 1:nofiles
    if i == 1        
        pct.ens.lowindrestruct(i).a = pct.ens.lowind.a(pct.ens.lowind.a <= yearlength*i & pct.ens.lowind.a >= 0)-(i-1)*yearlength;
        pct.ens.highindrestruct(i).a = pct.ens.highind.a(pct.ens.highind.a <= yearlength*i & pct.ens.highind.a >= 0)-(i-1)*yearlength;
    else
        pct.ens.lowindrestruct(i).a = pct.ens.lowind.a(pct.ens.lowind.a <= yearlength*i & pct.ens.lowind.a > yearlength*(i-1))-(i-1)*yearlength;
        pct.ens.highindrestruct(i).a = pct.ens.highind.a(pct.ens.highind.a <= yearlength*i & pct.ens.highind.a > yearlength*(i-1))-(i-1)*yearlength;
    end
    temp.lowextract(i).a = squeeze(datamontharrange(i,monthin,pct.ens.lowindrestruct(i).a));
    temp.highextract(i).a = squeeze(datamontharrange(i,monthin,pct.ens.highindrestruct(i).a));
        
end
dataout.lowextract = vertcat(temp.lowextract(:).a);
dataout.highextract = vertcat(temp.highextract(:).a);

%% taking upper and lower for each year
[out,index] = sort(datamontharrange,1);
Eachyear.lowextract = out(1:2,:,:);
Eachyear.lowindex = index(1:2,:,:);
Eachyear.highextract = out(end-1:end,:,:);
Eachyear.highindex = index(end-1:end,:,:);
end