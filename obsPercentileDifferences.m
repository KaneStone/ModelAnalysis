function [difference,p,pct,TSmonths] = obsPercentileDifferences(surfacetemperature,ozone,percentile,mons,tozmon,longitude,latitude,hem,plotting)

pct.lowerpercentile = prctile(ozone,percentile);
pct.upperpercentile = prctile(ozone,100-percentile);
pct.lowerind = find(ozone <= pct.lowerpercentile);
pct.upperind = find(ozone >= pct.upperpercentile);

if ~hem
    %surfacetemperature = cat(3,surfacetemperature,zeros(size(surfacetemperature,1),size(surfacetemperature,2),2));
    surfacetemperature = cat(3,surfacetemperature,surfacetemperature(:,:,end-11:end-10));
    surfacetemperature (surfacetemperature == 0) = NaN;
end

for i = 1:12    
    for j = 1:size(surfacetemperature,1)
        for k = 1:size(surfacetemperature,2)
            if hem
                TSmonths(j,k,:,i) = detrend(squeeze(surfacetemperature(j,k,i:12:end)))+...
                    nanmean(squeeze(surfacetemperature(j,k,12+i:12:end-2)));
            else
                if i < 5
                    TSmonths(j,k,:,i) = detrend(squeeze(surfacetemperature(j,k,12+i:12:end)));                                
                else
                    TSmonths(j,k,:,i) = detrend(squeeze(surfacetemperature(j,k,i:12:end)));
                end        
            end
        end
    end
end

upper = TSmonths(:,:,pct.upperind,mons);
lower = TSmonths(:,:,pct.lowerind,mons);
difference = nanmean(upper(:,:,:),3) - nanmean(lower(:,:,:),3);
difference = permute(reshape(difference,[size(difference),1]),[3,1,2]);
for i = 1:size(upper,1)
    for j = 1:size(upper,2)
       p(1,i,j) = ttest2(squeeze(upper(i,j,:)),squeeze(lower(i,j,:)));
    end
end


p (p == 0) = -1;
p (p == 1) = 0;

% Get areas of south east Australia

%%
% Ausstart

%% plotting

% sig = .05;
% polar_p_toplot = reshape(polar_p,[1,size(polar_p)]);
% polar_p_toplot (polar_p_toplot <= sig) = 0;
% polar_p_toplot (polar_p_toplot > sig) = 1;
%

if plotting

    areaslon = [70,70,90,90,70;...
        [300,300,320,320,300]-360;...
        [245,245,270,270,245]-360;...
        70,70,90,90,70];
    areaslat = [65,75,75,65,65;...
        65,78,78,65,65;...
        43,65,65,43,43;...
        25,35,35,25,25];

    if ~hem
        ylims = [-90 0];
        clims = [-3 3];
        contourtitle = {[monthnames(tozmon,0,0),' TCO - ',monthnames(mons,1,1),' ERA-interim surface air temperature differences']};
        hemext = 'south';
    else
        ylims = [0 90];
        clims = [-6 6];
        contourtitle = {[monthnames(tozmon,0,0),' TCO - ',monthnames(mons,1,1),' ERA-interim surface air temperature differences']};
        hemext = 'north';

    end

    subplotmaps(difference,longitude,latitude,{'div','RdBu'},1,p,16,contourtitle,'Longitude','','Temperature difference (K)','on',...
        clims,22,[-90:10:20],[-90:10:20],[longitude(2:10:end)],[longitude(2:10:end)],'',1,[0 360],ylims,0,'none',1,'Miller Cylindrical');

    for i = 1:size(areaslon,1)
        m_line(areaslon(i,:),areaslat(i,:),'linewi',4,'color','k'); 
    end

    if ~hem
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/compositedifferences/observations/South_',num2str(percentile),'th_percentile_',monthnames(mons,1,1)];
    else
        filename = ['/Users/kanestone/Dropbox (MIT)/Work_Share/MITwork/predruns/compositedifferences/observations/North_',num2str(percentile),'th_percentile',monthnames(mons,1,1)];
    end

    export_fig(filename,'-png');
end

end
