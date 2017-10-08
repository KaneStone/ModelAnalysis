function [pct] = predruns_varPercentiles(data,monthin,percentile,nofiles)

yearlength = size(data,2)/nofiles;

count = 1;
for i = monthin
    pct.lowpercentile(count) = prctile(data(monthin,:),percentile);
    pct.highpercentile(count) = prctile(data(monthin,:),100-percentile);
    pct.lowind(count).a = find(data(monthin,:) <= pct.lowpercentile(count));
    pct.highind(count).a = find(data(monthin,:) >= pct.highpercentile(count));  
    count = count+1;
end

for i = 1:nofiles
    if i == 1        
        pct.lowindrestruct(i).a = pct.lowind.a(pct.lowind.a <= yearlength*i & pct.lowind.a >= 0)-(i-1)*yearlength;
        pct.highindrestruct(i).a = pct.highind.a(pct.highind.a <= yearlength*i & pct.highind.a >= 0)-(i-1)*yearlength;
    else
        pct.lowindrestruct(i).a = pct.lowind.a(pct.lowind.a <= yearlength*i & pct.lowind.a > yearlength*(i-1))-(i-1)*yearlength;
        pct.highindrestruct(i).a = pct.highind.a(pct.highind.a <= yearlength*i & pct.highind.a > yearlength*(i-1))-(i-1)*yearlength;
    end
        
end

end