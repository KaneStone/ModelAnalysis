function [] = layer_deviation(MLSdeviation,WACCMdeviation,deviation_year,pressureplot,...
    OMIdeviationnorm,MAMTCO_deviationnorm,VCMAMTCO_deviationnorm,...
    Extinction_vertical_deviation_final,CALIPSO,MLSpressure,Latitudes,omilats,include_sulf,sa)

    
%% plotting deviations
montitle = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',...
    'Oct','Nov','Dec'};

    
fsize = 18;
montoplot = 5:12;
createfig('medium','on');

%minmax = [-67.5 12.5];
if sa
    minmax = [-2.5 2.5];
    contourinterval = .25;    
else
%    minmax = [-75 75];
%    contourinterval = 15;    
     minmax = [-1.5e12 1.5e12];
     contourinterval = .15e12;    
end


sulfurtoplot = Extinction_vertical_deviation_final;
sulcbrew = cbrewer('seq','YlOrBr',13);
sulcbrew = sulcbrew(3:12,:);  
sulf_minmax = [0 4e-6];
sulf_interval = .4e-6;
exttimes = 1e6;

cbrew = cbrewer('div','RdBu',21);
cbrew1 = flipud([cbrew(1:10,:);cbrew(12:21,:)]);
     
[~,presind] = min(abs(MLSpressure-pressureplot));

MLS_deviation_toplot = squeeze(MLSdeviation(montoplot,presind,:));
WACCM_deviation_toplot = squeeze(WACCMdeviation(:,montoplot,presind,:));
WACCM_Extinction = squeeze(Extinction_vertical_deviation_final(:,montoplot,presind,:)); 
CALIPSO_backscatter = squeeze(CALIPSO(:,presind,:));
CALIPSO_backscatter = [CALIPSO_backscatter;ones(1,34)*-9999];
CALIPSO_backscatter (CALIPSO_backscatter == -9999) = NaN;
latdevtitle = {['MLS at ',num2str(pressureplot),' hPa'],'OMI TCO',['MAM at ',num2str(pressureplot),' hPa'],...
    'MAM TCO',['VC-MAM at ',num2str(pressureplot),' hPa'],'VC-MAM TCO'};

for i = 1:6
    sp = subplot(3,2,i);
    
    if i == 1            
        [~,h(i)] = contourf(1:size(MLS_deviation_toplot,1),Latitudes(1:size(MLS_deviation_toplot,2))...
            ,MLS_deviation_toplot',minmax(1):contourinterval:minmax(2),'LineStyle','-');   
    elseif i == 2
        [~,h(i)] = contourf(1:8,omilats(1:63),squeeze(OMIdeviationnorm(end,5:12,1:63))',minmax(1):contourinterval:minmax(2),'LineStyle','-');            
    elseif i == 3
        [~,h(i)] = contourf(1:size(squeeze(WACCM_deviation_toplot(2,:,:)),1),Latitudes(1:size(MLS_deviation_toplot,2))...
            ,squeeze(WACCM_deviation_toplot(2,:,:))',minmax(1):contourinterval:minmax(2),'LineStyle','-');           
    elseif i == 4
        [~,h(i)] = contourf(1:8,Latitudes(1:34),squeeze(MAMTCO_deviationnorm(5:12,end,1:34))',minmax(1):contourinterval:minmax(2),'LineStyle','-');        
        %title('OMI','fontsize',22)
    elseif i == 5
        [~,h(i)] = contourf(1:size(squeeze(WACCM_deviation_toplot(3,:,:)),1),Latitudes(1:size(MLS_deviation_toplot,2)),...
            squeeze(WACCM_deviation_toplot(3,:,:))',minmax(1):contourinterval:minmax(2),'LineStyle','-');               
        %title('MAM','fontsize',22)
    elseif i == 6
        [~,h(i)] = contourf(1:8,Latitudes(1:34),squeeze(VCMAMTCO_deviationnorm(5:12,end,1:34))',minmax(1):contourinterval:minmax(2),'LineStyle','-');
        %title('VC-MAM','fontsize',22);
         
    %colorbar         
    end        
    %ylim([-90 -30])
    if i == 5 || i == 6
        xlabel('Month','fontsize',fsize+2);
    end
    if i == 1 || i == 3 || i == 5
        ylabel('Latitude','fontsize',fsize+2);
    end
    
    sp_pos(i,:) = get(sp,'position');
    if i == 1 || i == 3 || i == 5
        set(sp,'position',[sp_pos(i,1)-.02,sp_pos(i,2)-.01,sp_pos(i,3),sp_pos(i,4)-sp_pos(i,4)./10]);        
    else 
         set(sp,'position',[sp_pos(i,1)-.04,sp_pos(i,2)-.01,sp_pos(i,3),sp_pos(i,4)-sp_pos(i,4)./10]);        
    end
%         set(sp,'position',[sp_pos(i,1)-.02,sp_pos(i,2)-.01,sp_pos(i,3),sp_pos(i,4)-sp_pos(i,4)./10]);        
%     elseif i == 3 
%         set(sp,'position',[sp_pos(i,1)-.02,sp_pos(i,2)-.01,sp_pos(i,3),sp_pos(i,4)-sp_pos(i,4)./10]);        
%     end
    sp_pos(i,:) = get(sp,'Position');
    
    colormap(cbrew1);  
    caxis([minmax(1) minmax(2)]);
    set(gca,'color',[.8 .8 .8]);
    set(gca,'xtick',1:1:8,'xticklabel',montitle([montoplot]),'fontsize',fsize-2);
    set(gca,'ytick',-85:10:25,'yticklabel',-85:10:25,'fontsize',fsize-2); 
    if i == 6;
        ch = colorbar;
        cbaxloc = get(ch,'Position');
        suboutpos = get(sp,'OuterPosition');

        set(ch,'Position',[cbaxloc(1)+cbaxloc(1)/12, cbaxloc(2),...
            cbaxloc(3), sp_pos(1,2)+sp_pos(1,4)-sp_pos(6,2)],...
            'Box','on','YAxisLocation','right','fontsize',fsize-2); %
        
        set(get(ch,'ylabel'),'string','Normalized anomaly','fontsize',fsize)            
        
        set(ch,'Ytick',minmax(1):contourinterval*2:minmax(2))                         
    end           
    title(latdevtitle{i},'fontsize',fsize+2)        
 
if include_sulf                    
    haxes1 = gca;
    haxes1_pos = get(haxes1,'Position');
    haxes2 = axes('Position',haxes1_pos,...
          'XAxisLocation','bottom',...
          'YAxisLocation','left',...
          'Color','none');
    hold on   

    if i == 1 
        [~,h2(i)] = contour(1:size(squeeze(CALIPSO_backscatter),1)...
            ,Latitudes(1:size(CALIPSO_backscatter,2)),squeeze(CALIPSO_backscatter(:,:))',...
            sulf_minmax(1):sulf_interval:sulf_minmax(2),'LineWidth',2);   
    elseif i == 2
        [~,h(i)] = contour(1:size(squeeze(WACCM_Extinction(2,:,:)),1)...
            ,Latitudes(1:size(WACCM_Extinction,3)),squeeze(WACCM_Extinction(2,:,:))',...
            sulf_minmax(1):sulf_interval:sulf_minmax(2),'LineWidth',2);   
    elseif i == 3;
        [~,h(i)] = contour(1:size(squeeze(WACCM_Extinction(2,:,:)),1)...
            ,Latitudes(1:size(WACCM_Extinction,3)),squeeze(WACCM_Extinction(3,:,:))',...
            sulf_minmax(1):sulf_interval:sulf_minmax(2),'LineWidth',2);   
    end
         %ylim([-90 -30])
         colormap(haxes2,sulcbrew)            
         caxis([sulf_minmax(1) sulf_minmax(2)]);
         set(gca,'fontsize',fsize-2,'xtick',1:1:8,'xticklabel',montitle([montoplot]));
            set(gca,'ytick',-85:10:25,'yticklabel',-85:10:25,'fontsize',fsize-2); 
end
%          ch1 = colorbar('south');
%          set(ch1,'AxisLocation','out')
%          cbaxloc1 = get(ch1,'Position');
% %             if percent && ~ext
% %                 set(get(ch1,'ylabel'),'string','WACCM MAM chemical sulfate aerosol (percent)','fontsize',fsize)
% %             elseif percent && ext
% %                 set(get(ch1,'ylabel'),'string','WACCM absolute 2015 550 nm extinction (\times 10^{-5} m^-^1)','fontsize',fsize)                
% %             else
% %                 set(get(ch1,'ylabel'),'string','WACCM MAM chemical sulfate aerosol (cm^2/cm^3)','fontsize',fsize)                
% %             end
%             set(ch1,'Position',[sp_pos(1,1), .075,...
%             sp_pos(4,1)+sp_pos(4,3)-sp_pos(3,1), cbaxloc1(4)],...
%             'Box','on','YAxisLocation','right','fontsize',fsize-2);                                                                
%             set(ch1,'Ticks',minmax(1):contourinterval:minmax(2))             
%         end                               
%     end        
% end    
% 
% mtit([num2str(deviation_year),'-',num2str(deviation_year+1),...
%         ' deviation from the mean',' at ',sprintf('%.2f',MLSPressure(presind)),' hPa'],'xoff',-.01,'yoff',.08,'fontsize',fsize+8);    
% 
% if percent
%     filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/MLS/TimeLatitude/',...
%         'O3_Conc_',num2str(deviation_year),'_deviation_at_',sprintf('%.2f',MLSPressure(presind)),'_','percent','.png'];            
% else
%     filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/MLS/TimeLatitude/',...
%         'O3_Conc_',num2str(deviation_year),'_deviation_at_',sprintf('%.2f',MLSPressure(presind)),'_','absolute','.png'];            
%     filename2 = ['/Users/kanestone/Dropbox/Work_Share/MITwork/MLS/TimeLatitude/',...
%         'O3_Conc_',num2str(deviation_year),'_deviation_at_',sprintf('%.2f',MLSPressure(presind)),'_','absolute','.pdf'];            
% end
% 
% export_fig(filename,'-png');        
% export_fig(filename,'-pdf');        
end

annotation('textbox',[.31 .9 .1 .1],'String',...
    ['2015 normalized anomalies'],'LineStyle','none','fontsize',fsize+10);

if include_sulf
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatTime/2015anomaly/',...
        num2str(deviation_year),'_',num2str(pressureplot),'_NormAnom','_','.pdf'];
else
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatTime/2015anomaly/',...
        num2str(deviation_year),'_','NormAnom_',num2str(pressureplot),'_withaerosols','.pdf'];
end

export_fig(filename,'-pdf');

end






