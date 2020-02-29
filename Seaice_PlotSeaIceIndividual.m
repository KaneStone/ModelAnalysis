function [] = Seaice_PlotSeaIceIndividual(inputs,difference,latitude,longitude,mons,season,combine)

if season
    montype = 'single';
else
    montype = 'long';
end

if strcmp(inputs.var,'ICEFRAC')
    varname = 'sea ice fraction'; 
    clims = [-.4 .4];
    ctitle = {'Sea ice fraction'};
    patchcoast = 1;
else
    varname = 'snow cover';
    clims = [-.1 .1];
    ctitle = {'Meters'};   
    patchcoast = 0;
end

for i = 1:size(difference,3)
    for j = 1:size(difference,4)
        longtouse = [longitude;longitude(1)+360];

        title = {['No. ',sprintf('%02d',j),', Change in ',...
            monthnames(mons(i),1,montype),' ',...
            varname,' due to ',monthnames(inputs.tozmonth,1,'long'), ' ozone extremes']};        
        if combine
             toplotdifferencefinal = permute(difference(:,:,:,j),[3,2,1]);
        else
            toplotdifferencefinal = permute(difference(:,:,i,j),[3,2,1]);
        end
        toplotdifferencefinal = cat(2,toplotdifferencefinal,toplotdifferencefinal(:,1,:));

        filedir = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/predruns/ozoneseaice/Individual/';
        
        ch = subplotmaps(toplotdifferencefinal,longtouse,latitude,{'div','RdBu'},1,[],16,title,[],[],ctitle,'on',...
            clims,18,[longtouse(1:24:end-1)],[longtouse(1:24:end-1)],[0:15:90],[0:15:90],{''},1,[0 360],[45 90],0,'none',1,'Stereographic',patchcoast);
        gh = get(gca,'title');
        ghpos = get(gh,'position');
        set(gh,'Position',[ghpos(1) ghpos(2)+ghpos(2)./50  ghpos(3)])     
        filename = [filedir,'Ind_','No.',sprintf('%02d',j),'_',sprintf('%02d',mons(i)),'_',inputs.var,'_',...
            monthnames(inputs.tozmonth,1,'long'),'_Arcticozoneextremes_over',...
            num2str(inputs.timeperiodvar(1)),'-',num2str(inputs.timeperiodvar(2))];
        %print(filename,'-depsc');
        export_fig(filename,'-png');
        close 1
    end
    if combine
        break;
    end
end
end