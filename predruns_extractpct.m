function [varextract,varextractmean,vardifference] = predruns_extractpct(dataMonthArrange,pct_highcl,...
    pct_lowcl,pct_lowGHG,varmonth,highclnofiles,lowclnofiles,lowGHGnofiles)

% highcl
count = 1;
count2 = 1;
for i = 1:highclnofiles    
    if length(varmonth) > 1                
        for j = 1:length(varmonth)
            if varmonth(j) <= 4 
                yearplus = 1;
            else yearplus = 0;
            end
            temp.highcl.lowind(i).a(j,:,:,:) = squeeze(dataMonthArrange.highcl(i,varmonth(j),...
                pct_highcl.lowindrestruct(i).a+yearplus,:,:));
            temp.highcl.highind(i).a(j,:,:,:) = squeeze(dataMonthArrange.highcl(i,varmonth(j),...
                pct_highcl.highindrestruct(i).a+yearplus,:,:));
             %[member,month,ind,lat,lon]
        end        
        varextract.highcl.lowind(count:count-1+length(pct_highcl.lowindrestruct(i).a),:,:) = ...
            squeeze(nanmean(temp.highcl.lowind(i).a,1));
        varextract.highcl.highind(count2:count2-1+length(pct_highcl.highindrestruct(i).a),:,:) = ...
            squeeze(nanmean(temp.highcl.highind(i).a,1));                                     
    else
        if varmonth <= 2 
            yearplus = 1;
        else yearplus = 0;
        end
        varextract.highcl.lowind(count:count-1+length(pct_highcl.lowindrestruct(i).a),:,:) = ...
            squeeze(dataMonthArrange.highcl(i,varmonth,pct_highcl.lowindrestruct(i).a + yearplus,:,:));
        varextract.highcl.highind(count2:count2-1+length(pct_highcl.highindrestruct(i).a),:,:) = ...
            squeeze(dataMonthArrange.highcl(i,varmonth,pct_highcl.highindrestruct(i).a + yearplus,:,:));                                
    end
    count = count + length(pct_highcl.lowindrestruct(i).a);
    count2 = count2 + length(pct_highcl.highindrestruct(i).a);
end

% lowcl
count = 1;
count2 = 1;
 for i = 1:lowclnofiles
     if length(varmonth) > 1                
        for j = 1:length(varmonth)
            if varmonth(j) <= 4 
                yearplus = 1;
            else yearplus = 0;
            end
            temp.lowcl.lowind(i).a(j,:,:,:) = squeeze(dataMonthArrange.lowcl(i,varmonth(j),...
                pct_lowcl.lowindrestruct(i).a+yearplus,:,:));
            temp.lowcl.highind(i).a(j,:,:,:) = squeeze(dataMonthArrange.lowcl(i,varmonth(j),...
                pct_lowcl.highindrestruct(i).a+yearplus,:,:));
            %[member,month,ind,lat,lon]
        end
        varextract.lowcl.lowind(count:count-1+length(pct_lowcl.lowindrestruct(i).a),:,:) = ...
            squeeze(nanmean(temp.lowcl.lowind(i).a,1));
        varextract.lowcl.highind(count2:count2-1+length(pct_lowcl.highindrestruct(i).a),:,:) = ...
            squeeze(nanmean(temp.lowcl.highind(i).a,1));                                               
     else
        if varmonth <= 4 
            yearplus = 1;
        else yearplus = 0;
        end
        varextract.lowcl.lowind(count:count-1+length(pct_lowcl.lowindrestruct(i).a),:,:) = ...
            squeeze(dataMonthArrange.lowcl(i,varmonth,pct_lowcl.lowindrestruct(i).a + yearplus,:,:));
        varextract.lowcl.highind(count2:count2-1+length(pct_lowcl.highindrestruct(i).a),:,:) = ...
            squeeze(dataMonthArrange.lowcl(i,varmonth,pct_lowcl.highindrestruct(i).a + yearplus,:,:));                                
    end
    count = count + length(pct_lowcl.lowindrestruct(i).a);
    count2 = count2 + length(pct_lowcl.highindrestruct(i).a);
 end
   
% lowGHG
count = 1;
count2 = 1;
for i = 1:lowGHGnofiles    
    if length(varmonth) > 1                
        for j = 1:length(varmonth)
            if varmonth(j) <= 4 
                yearplus = 1;
            else yearplus = 0;
            end
            temp.lowGHG.lowind(i).a(j,:,:,:) = squeeze(dataMonthArrange.lowGHG(i,varmonth(j),...
                pct_lowGHG.lowindrestruct(i).a+yearplus,:,:));
            temp.lowGHG.highind(i).a(j,:,:,:) = squeeze(dataMonthArrange.lowGHG(i,varmonth(j),...
                pct_lowGHG.highindrestruct(i).a+yearplus,:,:));
             %[member,month,ind,lat,lon]
        end        
        varextract.lowGHG.lowind(count:count-1+length(pct_lowGHG.lowindrestruct(i).a),:,:) = ...
            squeeze(nanmean(temp.lowGHG.lowind(i).a,1));
        varextract.lowGHG.highind(count2:count2-1+length(pct_lowGHG.highindrestruct(i).a),:,:) = ...
            squeeze(nanmean(temp.lowGHG.highind(i).a,1));                                     
    else
        if varmonth <= 2 
            yearplus = 1;
        else yearplus = 0;
        end
        varextract.lowGHG.lowind(count:count-1+length(pct_lowGHG.lowindrestruct(i).a),:,:) = ...
            squeeze(dataMonthArrange.lowGHG(i,varmonth,pct_lowGHG.lowindrestruct(i).a + yearplus,:,:));
        varextract.lowGHG.highind(count2:count2-1+length(pct_lowGHG.highindrestruct(i).a),:,:) = ...
            squeeze(dataMonthArrange.lowGHG(i,varmonth,pct_lowGHG.highindrestruct(i).a + yearplus,:,:));                                
    end
    count = count + length(pct_lowGHG.lowindrestruct(i).a);
    count2 = count2 + length(pct_lowGHG.highindrestruct(i).a);
end
 
 
varextractmean.lowcl.lowind = squeeze(nanmean(varextract.lowcl.lowind));
varextractmean.lowcl.highind = squeeze(nanmean(varextract.lowcl.highind));
varextractmean.highcl.lowind = squeeze(nanmean(varextract.highcl.lowind));
varextractmean.highcl.highind = squeeze(nanmean(varextract.highcl.highind));
varextractmean.lowGHG.lowind = squeeze(nanmean(varextract.lowGHG.lowind));
varextractmean.lowGHG.highind = squeeze(nanmean(varextract.lowGHG.highind));

vardifference.lowcl = varextractmean.lowcl.highind - varextractmean.lowcl.lowind; 
vardifference.highcl = varextractmean.highcl.highind - varextractmean.highcl.lowind;
vardifference.lowGHG = varextractmean.lowGHG.highind - varextractmean.lowGHG.lowind;

end
