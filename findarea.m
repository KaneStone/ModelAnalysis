function [difference] = findarea(cindhigh,cindlow,cindhighnanind,cindlownanind,lonbounds,latbounds,count,toplot,tactogetherind,tactogetherindlow)
% This function finds the sea ice area within a lat-lon boundary.

% There are different cases

% Case 1. Area is terminated along south boundary (tack together segments)
% Case 2. Area is enclosed
% Case 3. Area is bounded by west boundary only
% Case 4. Area is bounded by east boundary only
% Case 5. Area is bounded by north boundary only
% Case 6. Area is bounded by south boundary only
% Case 7. Area is bounded by west and north boundaries
% Case 8. Area is bounded by east and north boundaries
% Case 9. Area is bounded by west and south boundaries
% Case 10. Area is bounded by west and east boundaries


% first thing is sort from west to east.

if ~isempty(cindhigh)
    for l = 1:length(cindhighnanind)-1
        sortindexcind(l).s = cindhigh(:,cindhighnanind(l)+1:cindhighnanind(l+1)-1);
        minhighcl(l) = min(sortindexcind(l).s(1,:));
    end

    [~,ind] = sort(minhighcl);

    for l = 1:length(ind)
        cindhigh2(l).s = sortindexcind(ind(l)).s;
    end
end

if ~isempty(cindlow)
    for l = 1:length(cindlownanind)-1
        sortindexcind(l).s = cindlow(:,cindlownanind(l)+1:cindlownanind(l+1)-1);
        minlowcl(l) = min(sortindexcind(l).s(1,:));
    end

    [~,ind] = sort(minlowcl);

    for l = 1:length(ind)
        cindlow2(l).s = sortindexcind(ind(l)).s;
    end
end

if toplot
    figure;
end

%%
cindhighnanind = [find(isnan(cindhigh(1,:))),length(cindhigh(1,:))+1];
cindlownanind = [find(isnan(cindlow(1,:))),length(cindlow(1,:))+1];
if ~isempty(cindhigh)
    for l = 1:length(cindhighnanind)-1

        cindhightemplat(l).c = cindhigh(2,cindhighnanind(l)+1:cindhighnanind(l+1)-1);
        cindhightemplon(l).c = cindhigh(1,cindhighnanind(l)+1:cindhighnanind(l+1)-1);        

        % case 1.
        if ismember(count,tactogetherind)
            %length(cindhigh2) == 2 && cindhigh2(l).s(2,end) == latbounds(1) && cindhigh2(l+1).s(2,1) == latbounds(1)
            segnum = length(cindhigh2);
            cindhightack = [cindhigh2(:).s];
            cindhightacklon = cindhightack(1,:);
            cindhightacklat = cindhightack(2,:);
            cindhightacklon(1,end+1:end+4) = [lonbounds(2),lonbounds(2),lonbounds(1),lonbounds(1)];
            cindhightacklat(1,end+1:end+4) = [cindhightacklat(end),latbounds(2),latbounds(2),cindhightacklat(1)];        
            areahigh(l) = areaint(cindhightacklat,cindhightacklon,almanac('earth','geoid','kilometers'));        
            cindhightemplon(l).c = cindhightacklon;
            cindhightemplat(l).c = cindhightacklat;
            break
        elseif cindhigh(1,cindhighnanind(l)+1) == cindhigh(1,cindhighnanind(l+1)-1) && cindhigh(2,cindhighnanind(l)+1) == cindhigh(2,cindhighnanind(l+1)-1)
            areahigh(l) = areaint(cindhightemplat(l).c,cindhightemplon(l).c,almanac('earth','geoid','kilometers'));        
        % case 2. and case 3.
        elseif cindhigh(1,cindhighnanind(l)+1) == cindhigh(1,cindhighnanind(l+1)-1) && cindhigh(2,cindhighnanind(l)+1) ~= cindhigh(2,cindhighnanind(l+1)-1) 
            cindhightemplat(l).c(end+1) = cindhightemplat(l).c(1);
            cindhightemplon(l).c(end+1) = cindhightemplon(l).c(1);        
            areahigh(l) = areaint(cindhightemplat(l).c,cindhightemplon(l).c,almanac('earth','geoid','kilometers'));                
        % case 4. and case 5.
        elseif cindhigh(1,cindhighnanind(l)+1) ~= cindhigh(1,cindhighnanind(l+1)-1) && cindhigh(2,cindhighnanind(l)+1) == cindhigh(2,cindhighnanind(l+1)-1)        
            cindhightemplat(l).c(end+1) = cindhightemplat(l).c(1);
            cindhightemplon(l).c(end+1) = cindhightemplon(l).c(1);        
            areahigh(l) = areaint(cindhightemplat(l).c,cindhightemplon(l).c,almanac('earth','geoid','kilometers'));        
        % case 6.    
        elseif cindhigh(1,cindhighnanind(l)+1) > cindhigh(1,cindhighnanind(l+1)-1) && cindhigh(2,cindhighnanind(l)+1) < cindhigh(2,cindhighnanind(l+1)-1)        
            cindhightemplat(l).c(end+1:end+2) = [cindhightemplat(l).c(end),cindhightemplat(l).c(1)];
            cindhightemplon(l).c(end+1:end+2) = [cindhightemplon(l).c(1),cindhightemplon(l).c(1)];        
            areahigh(l) = areaint(cindhightemplat(l).c,cindhightemplon(l).c,almanac('earth','geoid','kilometers'));        
        % case 7.
        elseif cindhigh(1,cindhighnanind(l)+1) ~= lonbounds(1) && cindhigh(1,cindhighnanind(l)+1) < cindhigh(1,cindhighnanind(l+1)-1) && cindhigh(2,cindhighnanind(l)+1) > cindhigh(2,cindhighnanind(l+1)-1)        
            cindhightemplat(l).c(end+1:end+2) = [cindhightemplat(l).c(1),cindhightemplat(l).c(1)];
            cindhightemplon(l).c(end+1:end+2) = [cindhightemplon(l).c(end),cindhightemplon(l).c(1)];        
            areahigh(l) = areaint(cindhightemplat(l).c,cindhightemplon(l).c,almanac('earth','geoid','kilometers'));        
            figure;
            %plot(cindhightemplon(l).c,cindhightemplat(l).c);
        % case 8.
        elseif cindhigh(1,cindhighnanind(l)+1) == lonbounds(2) && cindhigh(2,cindhighnanind(l+1)-1) == latbounds(2)
            %cindhigh(1,cindhighnanind(l)+1) < cindhigh(1,cindhighnanind(l+1)-1) && cindhigh(2,cindhighnanind(l)+1) > cindhigh(2,cindhighnanind(l+1)-1) && cindhigh(1,2) ~= lonbounds(1) && cindhigh(1,end) ~= lonbounds(2)        
            cindhightemplat(l).c(end+1:end+2) = [cindhightemplat(l).c(end),cindhightemplat(l).c(1)];
            cindhightemplon(l).c(end+1:end+2) = [cindhightemplon(l).c(1),cindhightemplon(l).c(1)];        
            areahigh(l) = areaint(cindhightemplat(l).c,cindhightemplon(l).c,almanac('earth','geoid','kilometers'));        
        % case 9   
        elseif cindhigh(1,cindhighnanind(l)+1) == lonbounds(1) && cindhigh(1,cindhighnanind(l+1)-1) == lonbounds(2)        
            cindhightemplat(l).c(end+1:end+3) = [latbounds(2),latbounds(2),cindhightemplat(l).c(1)];
            cindhightemplon(l).c(end+1:end+3) = [lonbounds(2),lonbounds(1),lonbounds(1)];                        
            areahigh(l) = areaint(cindhightemplat(l).c,cindhightemplon(l).c,almanac('earth','geoid','kilometers'));    
        end
    end    
    if toplot
        for i = 1:length(cindhightemplon)
            hold on
            plot(cindhightemplon(i).c,cindhightemplat(i).c,'b');        
        end
    end
else
    areahigh = 0;
end

if ~isempty(cindlow)
    for l = 1:length(cindlownanind)-1        

        cindlowtemplat(l).c = cindlow(2,cindlownanind(l)+1:cindlownanind(l+1)-1);
        cindlowtemplon(l).c = cindlow(1,cindlownanind(l)+1:cindlownanind(l+1)-1);

        % case 1.
        if ismember(count,tactogetherindlow)
            cindlowtack = [cindlow2(:).s];
            cindlowtacklon = cindlowtack(1,:);
            cindlowtacklat = cindlowtack(2,:);
            cindlowtacklon(1,end+1:end+4) = [lonbounds(2),lonbounds(2),lonbounds(1),lonbounds(1)];
            cindlowtacklat(1,end+1:end+4) = [cindlowtacklat(end),latbounds(2),latbounds(2),cindlowtacklat(1)];

            arealow(l) = areaint(cindlowtacklat,cindlowtacklon,almanac('earth','geoid','kilometers'));
            %arealow(l) = areaint(cindlowtemplat(l).c,cindlowtemplon(l).c,almanac('earth','geoid','kilometers'));

            cindlowtemplon(l).c = cindlowtacklon;
            cindlowtemplat(l).c = cindlowtacklat;

        break
        elseif cindlow(1,cindlownanind(l)+1) == cindlow(1,cindlownanind(l+1)-1) && cindlow(2,cindlownanind(l)+1) == cindlow(2,cindlownanind(l+1)-1)
            arealow(l) = areaint(cindlowtemplat(l).c,cindlowtemplon(l).c,almanac('earth','geoid','kilometers'));
        % case 2. and case 3.
        elseif cindlow(1,cindlownanind(l)+1) == cindlow(1,cindlownanind(l+1)-1) && cindlow(2,cindlownanind(l)+1) ~= cindlow(2,cindlownanind(l+1)-1)        
            cindlowtemplat(l).c(end+1) = cindlowtemplat(l).c(1);
            cindlowtemplon(l).c(end+1) = cindlowtemplon(l).c(1);                
            arealow(l) = areaint(cindlowtemplat(l).c,cindlowtemplon(l).c,almanac('earth','geoid','kilometers'));

        % case 4. and case 5.
        elseif cindlow(1,cindlownanind(l)+1) ~= cindlow(1,cindlownanind(l+1)-1) && cindlow(2,cindlownanind(l)+1) == cindlow(2,cindlownanind(l+1)-1)        
            cindlowtemplat(l).c(end+1) = cindlowtemplat(l).c(1);
            cindlowtemplon(l).c(end+1) = cindlowtemplon(l).c(1);                
            arealow(l) = areaint(cindlowtemplat(l).c,cindlowtemplon(l).c,almanac('earth','geoid','kilometers'));
        % case 6.    
        elseif cindlow(1,cindlownanind(l)+1) > cindlow(1,cindlownanind(l+1)-1) && cindlow(2,cindlownanind(l)+1) < cindlow(2,cindlownanind(l+1)-1)        
            cindlowtemplat(l).c(end+1:end+2) = [cindlowtemplat(l).c(end),cindlowtemplat(l).c(1)];
            cindlowtemplon(l).c(end+1:end+2) = [cindlowtemplon(l).c(1),cindlowtemplon(l).c(1)];                
            arealow(l) = areaint(cindlowtemplat(l).c,cindlowtemplon(l).c,almanac('earth','geoid','kilometers'));
        % case 7.
        elseif cindlow(1,cindlownanind(l)+1) ~= lonbounds(1) && cindlow(1,cindlownanind(l)+1) < cindlow(1,cindlownanind(l+1)-1) && cindlow(2,cindlownanind(l)+1) > cindlow(2,cindlownanind(l+1)-1)        
            cindlowtemplat(l).c(end+1:end+2) = [cindlowtemplat(l).c(1),cindlowtemplat(l).c(1)];
            cindlowtemplon(l).c(end+1:end+2) = [cindlowtemplon(l).c(end),cindlowtemplon(l).c(1)];                
            arealow(l) = areaint(cindlowtemplat(l).c,cindlowtemplon(l).c,almanac('earth','geoid','kilometers'));
        % case 8.
        elseif cindlow(1,cindlownanind(l+1)-1) == lonbounds(2) && cindlow(2,cindlownanind(l)+1) == latbounds(2)
            cindlowtemplat(l).c(end+1:end+2) = [cindlowtemplat(l).c(end),cindlowtemplat(l).c(1)];
            cindlowtemplon(l).c(end+1:end+2) = [cindlowtemplon(l).c(1),cindlowtemplon(l).c(1)];                
            arealow(l) = areaint(cindlowtemplat(l).c,cindlowtemplon(l).c,almanac('earth','geoid','kilometers'));
        % case 9   
        elseif cindlow(1,cindlownanind(l)+1) == lonbounds(1) && cindlow(1,cindlownanind(l+1)-1) == lonbounds(2)               
            cindlowtemplat(l).c(end+1:end+3) = [latbounds(2),latbounds(2),cindlowtemplat(l).c(1)];
            cindlowtemplon(l).c(end+1:end+3) = [lonbounds(2),lonbounds(1),lonbounds(1)];                
            arealow(l) = areaint(cindlowtemplat(l).c,cindlowtemplon(l).c,almanac('earth','geoid','kilometers'));        
        elseif cindlow(1,cindlownanind(l)+1) == lonbounds(1) && cindlow(2,cindlownanind(l+1)-1) == latbounds(2)               
            cindlowtemplat(l).c(end+1:end+2) = [latbounds(2),cindlowtemplat(l).c(1)];
            cindlowtemplon(l).c(end+1:end+2) = [lonbounds(1),lonbounds(1)];                
            arealow(l) = areaint(cindlowtemplat(l).c,cindlowtemplon(l).c,almanac('earth','geoid','kilometers'));        

        end
    end
    if toplot
        for i = 1:length(cindlowtemplon)
            hold on
            plot(cindlowtemplon(i).c,cindlowtemplat(i).c,'r');
        end
    end
else
    arealow = 0;
end

indtosubtracthigh = [3,4,7,13];
indtosubtractlow = [3,4,7,13];

if ismember(count,indtosubtracthigh)
    if count == 4
        difference = (areahigh(2) - sum(areahigh([1,3])))-(arealow(2) - sum(arealow([1,3])));
    else
        difference = diff(areahigh)-diff(arealow);
    end
else
    difference = sum(areahigh)-sum(arealow);
end

end