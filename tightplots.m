function [axpos_reorder] = tightplots(gcf,percentwidth,percentheight,cbaradjust,adjusttitles)

    % get subplot axes and position
    ax = get(gcf,'children');
    count = 1;
    colorbarexists = 0;
    for i = 1:length(ax)
        if strcmp(get(ax(i),'type'),'axes') 
            axpos(count,:) = get(ax(i),'position');
            if axpos(count,:) == [0 0 1 1]
                axpos(count,:) = [];
                continue;
            end
            axkeep(count) = ax(i);
            if i == length(ax)
                axwidth = axpos(count,3);
                axheight = axpos(count,4);
            end   
            count = count+1;
        elseif strcmp(get(ax(i),'type'),'colorbar') 
            colorbarexists = 1;
            chpos = get(ax(i),'position');   
            cind = i;
        end
    end         
          
    % order as bottom row upwards
    [~,axorder] = sort(sum(axpos(:,1:2),2));
    axpos_reorder = axpos(axorder,:);    
    axkeep = axkeep(axorder);        
    
    wantedwidthadjustment = (percentwidth./100).*axpos_reorder(1,3);        
    wantedheightadjustment = (percentheight./100).*axpos_reorder(1,4);        
    
    for i = 1:size(axpos,1)
        % find percentage of sub plot width
        rownumber = round((axpos_reorder(i,1)+axpos_reorder(i,3))./(axpos_reorder(1,1)+axpos_reorder(1,3)));
        colnumber = round((axpos_reorder(i,2)+axpos_reorder(i,4))./(axpos_reorder(1,2)+axpos_reorder(1,4)));
        if axpos_reorder(i,1) ~= min(axpos_reorder(:,1))        
            set(axkeep(i),'position',[axpos_reorder(i,1)-wantedwidthadjustment*(rownumber-1),...
                axpos_reorder(i,2),axpos_reorder(i,3),axpos_reorder(i,4)]);            
            axpos_reorder(i,:) = get(axkeep(i),'position');
        end
        if axpos_reorder(i,2) ~= min(axpos_reorder(:,2))        
            set(axkeep(i),'position',[axpos_reorder(i,1),axpos_reorder(i,2)-...
                wantedheightadjustment*(colnumber-1),axpos_reorder(i,3),axpos_reorder(i,4)]);
            axpos_reorder(i,:) = get(axkeep(i),'position');
        end
        
    end
    
    % adjust colorbar
    if colorbarexists        
        set(ax(cind),'position',[max(axpos_reorder(:,1))+axpos_reorder(1,3)+cbaradjust,...
            min(axpos_reorder(:,2)),chpos(3),max(axpos_reorder(:,2))+axpos_reorder(1,4)-axpos_reorder(1,2)]);        
        cbarrow
    end
    
    
    %adjust titles
    if adjusttitles
        for i = 1:size(axpos,1)
            titpos = get(get(axkeep(i),'title'),'position');
            set(get(axkeep(i),'title'),'position',[titpos(1),titpos(2)+adjusttitles,titpos(3)]);
        end
    end
end