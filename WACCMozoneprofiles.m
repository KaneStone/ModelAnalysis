%WACCM vertical ozone
clear all
close all
variable = 'O3'; %O3, NO2, NO, NOX, CLO

%file locations
directory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/',variable,'/'];
files = dir([directory,variable,'*']);

TAdirectory = '/Users/kanestone/work/projects/WACCM/netcdffiles/Temperature/';
TAfiles = dir([TAdirectory,'TA*']);

HCdirectory = '/Users/kanestone/work/projects/WACCM/netcdffiles/hybrid_ht_conversion/';
HCfiles = dir([HCdirectory,'hybrid*']);

datedir = '/Users/kanestone/work/projects/WACCM/netcdffiles/date/';
datefiles = dir([datedir,'date*']);

%colormaps
cbrew = cbrewer('qual','Set1',10);
cbrew([5,7],:) = [];
cbrewfade = cbrewer('qual','Pastel1',10);
cmap = flipud(cbrewer('div','RdBu',10));

names = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
namestitle = {'CCMI','MAM','VC-MAM','Chem-only','Chem-only-fixedSSTs','Chem-only-noleap'};
%names = {'CCMI','ChemonlyfixedSSTs','Chemonlynoleap'};
type = 'trend'; %trend, difference, 'raw','latitude_pressure
yearmin = 20000201;
yearmax = 20160101;
remove2002 = 'Southern';

latitudes = [[-60 60];[-60 60]];
szlat = size(latitudes);
percent = 1;
convert_to_ND = 1;
volcanic = 0;
linepercent = 0;
check2010 = 1;

zonalmean = zeros(size(latitudes,1),12,6,88,16);
zonalmean (zonalmean == 0) = NaN;

[NumberDensity, MolConc, TAdata, Pressure, Altitude, GMH, Latitudes, Longitudes, datedata]  = ReadWACCMvertical(variable,'monthly','/Users/kanestone/work/projects/WACCM/netcdffiles/',0);


%     %read in O3
%     [info.(names{i}), data.(names{i}), attributes.(names{i})] = ...
%         Read_in_netcdf([directory,files(i).name]);
%     
%     if strcmp(variable,'NONO2')
%         %data.(names{i}).NONO2 = permute(data.(names{i}).NONO2,[3,2,1]);
%     end
%     
%     %read in temperature
%     [TAinfo.(names{i}), TAdata.(names{i}), TAattributes.(names{i})] = ...
%         Read_in_netcdf([TAdirectory,TAfiles(i).name]);
%     
%     %read in hybrid_ht conversion parameters
%     [HCinfo.(names{i}), HCdata.(names{i}), HCattributes.(names{i})] = ...
%         Read_in_netcdf([HCdirectory,HCfiles(i).name]);
%     
%     %read in dates
%     [dateinfo.(names{i}), datedata.(names{i}), dateattributes.(names{i})] = ...
%         Read_in_netcdf([datedir,datefiles(i).name]);
%     
%     %finding latitudes
%     
%     
%     
%     %calculating to pressure
%     Pressure.(names{i}) = zeros(size(HCdata.(names{i}).PS,1),...
%         size(HCdata.(names{i}).hyam,1),size(HCdata.(names{i}).PS,2));    
%     hyam = permute(repmat(HCdata.(names{i}).hyam,[1,size(HCdata.(names{i}).PS)]),[2,1,3]);
%     hybm = permute(repmat(HCdata.(names{i}).hybm,[1,size(HCdata.(names{i}).PS)]),[2,1,3]);
%     P0 = permute(repmat(HCdata.(names{i}).P0,[size(HCdata.(names{i}).hyam,1),size(HCdata.(names{i}).PS)]),[2,1,3]);
%     PS = double(permute(repmat(HCdata.(names{i}).PS,[1,1,size(HCdata.(names{i}).hyam,1)]),[1,3,2]));
%     Pressure.(names{i}) = hyam .* P0 + hybm .* PS;
%     
%     %Altitude.(names{i}) = -(log(Pressure.(names{i})./PS)*1.38066e-23.*TAdata.(names{i}).T)./((28.97*1.66e-27)*9.80616);        
%     
%     Altitude.(names{i}) = -(log(Pressure.(names{i})./101300.25)*8.31.*273.15)./(9.806*.0289);        

%     %converting to number density
%     if convert_to_ND
%         convert.(names{i}) = vmr2conc(data.(names{i}).(variable),TAdata.(names{i}).T,Pressure.(names{i}),variable,'conc');
%     else 
%         convert.(names{i}) = data.(names{i}).(variable);
%     end
    
for i = 1:length(names)  
    
    for p = 1:szlat(1)
         latindex(p).p = find(Latitudes >= latitudes(p,1) & ...
             Latitudes <= latitudes(p,2));
    end
    
    %extracting dates
    NumberDensity_dates.(names{i}) = NumberDensity.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date < yearmax);
    datetemp = datedata.(names{i}).date(datedata.(names{i}).date >= yearmin & ...
        datedata.(names{i}).date < yearmax);
    
    Pressure_dates.(names{i}) = Pressure.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date < yearmax);
    
    GMH_dates.(names{i}) = GMH.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date < yearmax);
    
    switch remove2002
        case 'all_lats'
            NumberDensity_dates.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            Pressure_dates.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            GMH_dates.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        case 'Southern'        
             NumberDensity_dates.(names{i})(Latitudes < 0,...
                datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            %GMH_dates.(names{i})(Latitudes < 0,...
            %    datetemp >= 20021001 & datetemp < 20031001) = NaN; 
    end
    NumberDensity_dates.(names{i}) = cat(3,NumberDensity_dates.(names{i}),NumberDensity_dates.(names{i})(:,:,end-1:end));
    Pressure_dates.(names{i}) = cat(3,Pressure_dates.(names{i}),Pressure_dates.(names{i})(:,:,end-1:end));
    GMH_dates.(names{i}) = cat(3,GMH_dates.(names{i}),GMH_dates.(names{i})(:,:,end-1:end));
    %calculating trends b = [lats,pressure,month,run,beta]
    for k = 1:12        
        temp = size(nanmean(NumberDensity_dates.(names{i})(latindex(1).p,:,k:12:end)),3);
        for p = 1:szlat(1)
            zonalmean(p,k,i,:,1:temp) = squeeze(nanmean(NumberDensity_dates.(names{i})(latindex(p).p,:,k:12:end)));          
            zonalmeanpressure(p,k,i,:,1:temp) = squeeze(nanmean(Pressure_dates.(names{i})(latindex(p).p,:,k:12:end)));          
            zonaltimemeanpressure(p,k,i,:) = nanmean(zonalmeanpressure(p,k,i,:,:),5);
            
            zonalmeanalt(p,k,i,:,1:temp) = squeeze(nanmean(GMH_dates.(names{i})(latindex(p).p,:,k:12:end)));          
            zonaltimemeanalt(p,k,i,:) = nanmean(zonalmeanalt(p,k,i,:,:),5);
            
            DU_coeff = 1e5*1.38e-21*1e3*(273.1/10.13);
            
            for o = 1:size(zonalmean,4)
                if o ~= 1
                    zonalmeanDU(p,k,i,o,1:temp) = DU_coeff*zonalmean(p,k,i,o,1:temp).*...
                        (zonalmeanalt(p,k,i,o-1,1:temp)-zonalmeanalt(p,k,i,o,1:temp))./1000;
                else
                    zonalmeanDU(p,k,i,o,1:temp) = DU_coeff*zonalmean(p,k,i,o,1:temp).*...
                        zonalmeanalt(p,k,i,o,1:temp)./1000;
                end
            end                 
            
            for l = 1:size(NumberDensity_dates.(names{i}),2)
                
                stats(p,l,k,i) = regstats(squeeze(zonalmean(p,k,i,l,:)),...
                    1:size(zonalmean,5),'linear',{'tstat','rsquare','yhat','r'});
                [b(p,l,k,i,:),bint(p,l,k,i,:,:),~,~] = regress(squeeze(zonalmean(p,k,i,l,:)),...
                    [ones(1,size(zonalmean,5));1:size(zonalmean,5)]',.8);
                
                
                conf(p,l,k,i,:) = bint(p,l,k,i,:,2) - b(p,l,k,i,:);
            end   
        end
    end 
    
    
end

%%
if linepercent
    %Pressure_dates.ChemonlyfixedSSTs(48,:,1)/100
    cbrewqual = zeros(10,3);
    cbrewqual(3:10,:) = cbrewer('qual','Set1',8);

    %cbrewqual(6,1:2) = .6
    %sollats = [-5 5;-5 5; -60 60; -60 60];
    sollats = [-5 5; -5 5; -60 60];
    vartouse = 'ChemonlyfixedSSTs';
    vartousetitle = 'Chem-only-fixedSSTs';

    figure;
    fi = gcf;
    set(fi,'color','white','position',[100 100 1200 1200]);

    if strcmp(variable,'O3')
        yinterval = 1;
    else
        yinterval = 5;
    end

    for k = 1:4
        if k == 1
            solpres = [[.0005 .005]; [2.5 5]; [25 50]; [25 50]]*100;
        elseif k == 2
            solpres = [[.0005 .005]; [2.5 5]; [2.5 5]; [2.5 5]]*100;
        elseif k == 3
            solpres = [[.0005 .005]; [2.5 5]; [.25 .5]; [.25 .5]]*100;
        elseif k == 4
            solpres = [[.0005 .005]; [2.5 5]; [.025 .05]; [.025 .05]]*100;
        end

        clearvars solarlinenames
        for i = 1:length(sollats);
            %Find pressures    
            if i == 1
                solarlinenames{i} = [num2str(sollats(i,1)),' to ',num2str(sollats(i,2)),'N',', ',...
                    num2str(solpres(i,1)/100),' to ',num2str(solpres(i,2)/100),' hPa (scaled)'];
            elseif i == 2
                solarlinenames{i} = [num2str(sollats(i,1)),' to ',num2str(sollats(i,2)),'N',', ',...
                    num2str(solpres(i,1)/100),' to ',num2str(solpres(i,2)/100),' hPa'];
            else
                solarlinenames{i} = [num2str(sollats(i,1)),' to ',num2str(sollats(i,2)),'N'];
            end

            latindex = find(data.(vartouse).lat >= sollats(i,1) & ...
                data.ChemonlyfixedSSTs.lat <= sollats(i,2));

            %mean pressures first

            %W = cosd(latitudes);
            %W (latitudes <latboundaries(1) | latitudes > latboundaries(2)) = [];
            %yearmeananomaly(i) = nansum(temp(:).*W1(:))./nansum(W1(:));            

            solpreslatmean(i,:,:) = nanmean(Pressure_dates.(vartouse)(latindex,:,:),1);

            solpindex(i).t = find(nanmean(solpreslatmean(i,:,:),3) > solpres(i,1) & ...
                nanmean(solpreslatmean(i,:,:),3) < solpres(i,2));

            solpressure(i).t = squeeze(nanmean(NumberDensity_dates.(vartouse)(:,solpindex(i).t,:),2));

            [yearmeananomaly(i,:), yearmean(i,:)] = TCOanomaly(solpressure(i).t, sollats(i,:), data.(vartouse).lat, [2000 2014], [2000 2015], 1:12);

            sollines(i).t = squeeze(nanmean(NumberDensity_dates.(vartouse)(latindex,solpindex(i).t,:)));
            sollines1(i,:) = squeeze(nanmean(sollines(i).t,1));
            if i == 1
                sollinespercent(i,:) = (sollines1(i,:) - nanmean(sollines1(i,:)))./nanmean(sollines1(i,:))*10;    
            else
                sollinespercent(i,:) = (sollines1(i,:) - nanmean(sollines1(i,:)))./nanmean(sollines1(i,:))*100;    
            end
        end


        for j = 1:length(sollats);
            count = 1;
            for i = 1:180/12
                yearly_averages(j,i) = nanmean(sollinespercent(j,count:count+11),2);
                count = count+12;
            end
        end

        sp = subplot(2,2,k);
        sp_pos = get(sp,'position');    
        for i = 1:length(sollats)
            if i == 1;
                ph1(i) = plot(yearmeananomaly(i,:)/10,'color',cbrewqual(i,:),'LineWidth',2);
            elseif i == 2
                %ph1(i) = plot(yearly_averages(i,:),'LineStyle','--','color',cbrewqual(i,:),'LineWidth',2);
                ph1(i) = plot(yearmeananomaly(i,:),'LineStyle','--','color',cbrewqual(i,:),'LineWidth',2);
            else
                %ph1(i) = plot(yearly_averages(i,:),'color',cbrewqual(i,:),'LineWidth',2);
                ph1(i) = plot(yearmeananomaly(i,:),'color',cbrewqual(i,:),'LineWidth',2);
            end
            hold on
        end
        if k == 1
            set(sp,'position',[sp_pos(1)+.025 sp_pos(2)-.025 sp_pos(3:4)]);        
            ylabel(['% ',variable],'fontsize',16);
        elseif k == 2
            set(sp,'position',[sp_pos(1)-.025 sp_pos(2)-.025 sp_pos(3:4)]);        
        elseif k == 3
            set(sp,'position',[sp_pos(1)+.025 sp_pos(2)+.025 sp_pos(3:4)]);        
            xlabel('year','fontsize',16);
            ylabel(['% ',variable],'fontsize',16);
        elseif k == 4
            set(sp,'position',[sp_pos(1)-.025 sp_pos(2)+.025 sp_pos(3:4)]);        
            xlabel('year','fontsize',16);

        end
        set(gca,'xtick',1:2:16,'xticklabel',2000:2:2050,'fontsize',14);
        set(gca,'ytick',-100:yinterval:100)


        xlim([0 16]);
        title([num2str(solpres(i,1)/100),' to ',num2str(solpres(i,2)/100),' hPa'],'fontsize',18);
        if k == 4
            %lh = columnlegend(2,solarlinenames,'fontsize',18,'boxoff','location','user');   
        end


    end
    pmtit = mtit([vartousetitle,' - Yearly average ', variable, ' percent differences from 2000-2014 mean'],...
             'fontsize',22,'color','k',...
             'xoff',0,'yoff',.075);

    % h_text = findobj(lh,'type','text');
    % 
    % set(h_text,'FontUnits','points','FontSize',20)


    %creating axes for title
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
        'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

    lh = legend(ph1(1:2),solarlinenames{1:2});        
    set(lh,'fontsize',14,'box','off','Orientation','vertical','position',[0.45, .05, .01, .01]);



    ph2(1) = plot(1,'color',cbrewqual(3,:),'LineWidth',2);
    hold on
    ph2(2) = plot(1,'color',cbrewqual(4,:),'LineWidth',2);
    set(ha,'Visible','off');
    %set(ph2,'Visible','off');

    %lh2 = legend(ph2,solarlinenames{3:4});        
    lh2 = legend(ph2,solarlinenames{3});        
    set(lh2,'fontsize',14,'box','off','Orientation','vertical','position',[0.65, .05, .01, .01]);

    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/Solarcycleanalysis/',vartouse,'_',variable,'_MidLatitudeaDifferentAreas.pdf'];
    export_fig(fi,filename);
end

%h_text = findobj(lh2,'type','text');
 
%set(h_text,'TextColor',[.5 .5 .5],'EdgeColor',[0 0 0]);
     
%%
%interp pressure values
% pres = [1000:-50:150 ...
%     100:-5:15 ...
%     10:-.5:1.5 ...
%     1:-.05:.15 ... 
%     .1:-.005:.015 ...
%     .01:-.0005:.0015 ...
%     .001:-.00005:.00015];


pres = [1000:-50:150 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.05:.1];

preslog = log(pres);

%prestick = [1000 500 250 100 50 25 10 5 2.5 1 .5 .25 .1 .05 .025 .01 .005 .0025 ...
%    .001];
prestick = [1000 500 250 100 50 25 10 5 2.5 1 .5 .25 .1 .05 .025 .01 .005 .0025 ...
    .001 .0005 .00025 .0001];

presticklog = log(prestick);

mon = {'January','February','March','April','May','June','July','August','September',...
    'October','November','December'};

%plotting
if convert_to_ND
    xname = 'molec/cm^{3}';
else
    xname = 'mol/mol';
end
fsize = 20;
lwidth = 2;


switch type
    %%
    case 'raw'
        pressuretimeseries = 1;
        count = 1;
        for rawmonth = 10;
            rawlatitude = [75];
            rawpressures = [];
            if length(rawlatitude) == 1            
                [~, rawlatindex]  = min(abs(data.CCMI.lat - rawlatitude));
            else
                rawlatindex = find(Latitudes >= rawlatitude(1) & ...
                Latitudes <= rawlatitude(2));
            end

            for i = 1:length(names) 
                for l = 1:length(NumberDensity_dates.(names{i}))
                    for k = 1:length(rawlatitude)
                        rawozone(i,l,k,:) = interp1(log(squeeze(Pressure.(names{i})(rawlatindex(k),:,l))/100),...
                            squeeze(NumberDensity_dates.(names{i})(rawlatindex(k),:,l)),log(pres));
                    end        
                    %rawdata(i,:,:) = squeeze(NumberDensity_dates.(names{i})(rawlatindex,:,rawmonth:12:end)); 
                end            
            end
            rawozone = squeeze(rawozone(:,rawmonth:12:end,:,:));        
            rawmamvcman = (squeeze(rawozone(2,10,:))-squeeze(rawozone(3,10,:)))./squeeze(rawozone(3,10,:))*100;

            mon = {'January','February','March','April','May','June','July','August','September','October','November','December'};
            
            if ~pressuretimeseries
                if count == 1;
                    figure;
                    fig = gcf;                            
                    set(fig,'color','white','position',[100 100 700 1000]);
                end                 
                subplot(3,1,count)
                ph = plot(rawmamvcman,preslog,'LineWidth',lwidth);                    
                set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),...
                    'fontsize',fsize);
                ylabel('Pressure (hPa)','fontsize',fsize+2);
                xlabel('% ozone','fontsize',fsize+2);
                title(strcat(mon(rawmonth),' MAM - VC-MAM at','{ }',num2str(rawlatitude),'N'),'fontsize',fsize+4);
                ylim([preslog(55) preslog(1)]);
            elseif pressuretimeseries
                if count == 1;
                    figure;
                    fig = gcf;
                    set(fig,'color','white','position',[100 100 700 1000]);
                end           
                subplot(3,1,count);
                ph2forleg = plot(squeeze(rawozone(2,:,16)),'k','LineWidth',lwidth);
                hold on       
                ph3forleg = plot(squeeze(rawozone(3,:,16)),'--k','LineWidth',lwidth);
                ph2 = plot(squeeze(rawozone(2,:,16:18)),'LineWidth',lwidth);        
                ax = gca;
                ax.ColorOrderIndex = 1;        
                ph3 = plot(squeeze(rawozone(3,:,16:18)),'--','LineWidth',lwidth);
                ylabel('molecules/cm^3','fontsize',fsize+2);
                xlabel('Year','fontsize',fsize+2);
                title(strcat(mon(rawmonth),' ozone at','{ }',num2str(rawlatitude),'N'),'fontsize',fsize+4);
                set(gca,'xtick',1:2:15,'xticklabel',2000:2:2015,'fontsize',fsize);
                xlim([0 16])
            %lh = legend([ph2' ph2_forleg ph3_forleg],'250 hPa', '200 hPa', '150 hPa', 'MAM' , 'VC-MAM');            
            end
            count = count+1;
            clearvars rawozone
        end
        if pressuretimeseries
            lh = legend([ph2forleg ph3forleg ph2'],'MAM','VC-MAM','250 hPa', '200 hPa', '150 hPa');
            set(lh,'orientation','horizontal','position',[.5 .02 .01 .01],'fontsize',fsize+2,'box','off')
            filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/70North/',sprintf('%02d',rawmonth),'_',num2str(rawlatitude),'_Verticalvolcanic.pdf'];
            export_fig(fig,filename);
        else            
            export_fig(fig,['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/70North/',sprintf('%02d',rawmonth),'_',num2str(rawlatitude),'_TimeSeries.pdf']);
        end
    %%
    case 'difference'   
        lnames = {[num2str(latitudes(1,1)),' to ',num2str(latitudes(1,2)),'{\circ}N ','2006 - 2012 average'],...
            [num2str(latitudes(2,1)),' to ',num2str(latitudes(2,2)),'{\circ}N ','2006 - 2012 average'],...
            [num2str(latitudes(1,1)),' to ',num2str(latitudes(1,2)),'{\circ}N ','2000 - 2005 ; 2013-2014 average'],...
            [num2str(latitudes(2,1)),' to ',num2str(latitudes(2,2)),'{\circ}N ','2000 - 2005 ; 2013-2014 average']};
        for j = 1:12;    
            figure;
            fig = gcf;
            set(fig,'color','white','position',[100 100 840 560]);    
            for i = 1:2
                
                for k = 1:length(data.CCMI.lev)
                    [hvc(i,j,k),pvc(i,j,k)] = ttest2(flipud(squeeze(zonalmean(i,j,2,k,[1,2,3,4,6,14,15]))),...
                        flipud(squeeze(zonalmean(i,8,3,k,[1,2,3,4,6,14,15]))),'Alpha',.05);

                    [h(i,j,k),p(i,j,k)] = ttest2(flipud(squeeze(zonalmean(i,j,2,k,7:13))),...
                        flipud(squeeze(zonalmean(i,8,3,k,7:13))),'Alpha',.05);
                end

                ozonemamvc = nanmean(flipud(squeeze(zonalmean(i,j,2,:,[1,2,3,4,5,6,14,15]))),2);
                ozonevcmamvc = nanmean(flipud(squeeze(zonalmean(i,j,3,:,[1,2,3,4,5,6,14,15]))),2);
                pressuremamvc = nanmean(flipud(squeeze(zonalmeanpressure(i,j,2,:,[1,2,3,4,5,6,14,15]))),2);
                pressurevcmamvc = nanmean(flipud(squeeze(zonalmeanpressure(i,j,3,:,[1,2,3,4,5,6,14,15]))),2);

                ozonemamvcpres = interp1(log(pressuremamvc/100),ozonemamvc,log(pres));
                ozonevcmamvcpres = interp1(log(pressuremamvc/100),ozonevcmamvc,log(pres));

                ozonemam = nanmean(flipud(squeeze(zonalmean(i,j,2,:,7:13))),2);
                ozonevcmam = nanmean(flipud(squeeze(zonalmean(i,j,3,:,7:13))),2);
                pressuremam = nanmean(flipud(squeeze(zonalmeanpressure(i,j,2,:,7:13))),2);
                pressurevcmam = nanmean(flipud(squeeze(zonalmeanpressure(i,j,3,:,7:13))),2);

                ozonemampres = interp1(log(pressuremamvc/100),ozonemam,log(pres));
                ozonevcmampres = interp1(log(pressuremamvc/100),ozonevcmam,log(pres));

                if percent
                    volclean = (ozonemamvcpres - ozonevcmamvcpres)./ozonevcmamvcpres*100;    
                    volc = (ozonemampres - ozonevcmampres)./ozonevcmampres*100;    
                else
                    volclean = (ozonemamvcpres - ozonevcmamvcpres);
                    volc = (ozonemampres - ozonevcmampres);
                end

                ph(i) = plot(volc,preslog,'LineWidth',lwidth);
                hold on
                ph1(i) = plot(volclean,preslog,'LineWidth',lwidth);        
            end
            set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),...
                'fontsize',fsize);
            ylabel('Pressure (hPa)','fontsize',fsize+2);
            if percent
                xlabel('% O3','fontsize',fsize+2);
                title([mon{j},' MAM - VC-MAM ','% O3 difference'],'fontsize',fsize+4)
            else
                xlabel(xname,'fontsize',fsize+2);
                title([mon{j},' MAM - VC-MAM ','O3 difference'],'fontsize',fsize+4)
            end

            lh = legend([ph ph1],lnames,'location','NorthWest','fontsize',fsize,'box','off');
            if percent
                export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/percent/',...
                    variable,'_',sprintf('%02d',j),mon{j},'_percent.pdf']);
            else
                export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/',...
                    variable,'_',sprintf('%02d',j),mon{j},'.pdf']);
            end
        end
    case 'trend'
                                   
        %% Interp onto regular km levels
       latboundaries = [30 60];
        alts = 1000:1000:60000;
        for k = 1:length(names)
            for i = 1:size(NumberDensity_dates.(names{k}),1)
                for j = 1:size(NumberDensity_dates.(names{k}),3)                                        
                    datainterpalt(k,i,:,j) = interp1(squeeze(GMH_dates.(names{k})(i,:,j)),squeeze(NumberDensity_dates.(names{k})(i,:,j)),alts,'linear');                                        
                end
            end
            solnorm(k,:,:) = squeeze(nanmean(NumberDensity_dates.(names{k})(:,11:14,:),2));
        end
        DU_coeff = 1e5*1.38e-21*1e3*(273.1/10.13);
        solnormtemp = repmat(solnorm,1,1,length(alts),1);
        
        %solnormscale = (solnorm - min(solnorm(:)))./(max(solnorm(:)) - min(solnorm(:)));        
        
        for i = 1:length(alts)            
%             datainterpaltnorm(:,:,i,:) = squeeze(datainterpalt(:,:,i,:)) ./ ...
%                 squeeze(repmat(nanmean(datainterpalt(:,:,i,:),4),1,1,size(datainterpalt,4))) - ...
%                 ((solnorm) ./ repmat(nanmean(solnorm(:,:,:),3),1,1,size(datainterpalt,4))) .* ...
%                 squeeze(repmat(nanmean(datainterpalt(:,:,i,:),4),1,1,size(datainterpalt,4))) ...
%                 + squeeze(datainterpalt(:,:,i,:));                  
            datatemp = squeeze(datainterpalt(:,:,i,:));             
            slt(:,:,i,:) = min(datatemp(:)) + (solnorm - min(solnorm(:))).*(max(datatemp(:))-min(datatemp(:)))./ ...
                (max(solnorm(:))-min(solnorm(:)));
            datainterpaltnorm(:,:,i,:) = datainterpalt(:,:,i,:) - slt(:,:,i,:);
            
             
        end
        datainterpaltDU = datainterpalt*DU_coeff;
        
        %normalising to account for solar cycle (xs = xt-x0[xt]/[x0])
        
            
        
        
        
        %% Calculate vertical anomalies        
        mons = 1:2;        
        if length(mons) ~= 12
            fortitlemons = sprintf('%02d',mons);
        else
            fortitlemons = 'yearly average';
        end
        
        diffnames = {'Volcanoes','Dynamics','Chemistry','ssts'};
        vertdiff(1,:,:,:) = datainterpaltDU(2,:,:,:)-datainterpaltDU(3,:,:,:);
        vertdiff(2,:,:,:) = datainterpaltDU(3,:,:,:)-datainterpaltDU(4,:,:,:);
        vertdiff(3,:,:,:) = datainterpaltDU(5,:,:,:);
        vertdiff(4,:,:,:) = datainterpaltDU(4,:,:,:)-datainterpaltDU(5,:,:,:);
        
        for i = 1:length(names)            
            [yearmeananomaly(i,:,:), yearmean(i,:,:)] = Verticalanomaly(squeeze(datainterpaltDU(i,:,:,:)),...
                alts, latboundaries, Latitudes, [2000 2014], [2000 2015], mons,1);
        end        
        
        for i = 1:length(diffnames)            
            [yearmeananomalydiff(i,:,:), yearmeandiff(i,:,:)] = Verticalanomaly(squeeze(vertdiff(i,:,:,:)),...
                alts, latboundaries, Latitudes, [2000 2014], [2000 2015], mons,1);
        end        
        
        %% Calculate vertical percent trends
        for k = 1:length(names)
            for i = 1:length(alts)
                [bvert(k,i,:), bvertint(k,i,:,:),~,~] = regress(squeeze(yearmean(k,:,i))',...
                    [ones(1,size(yearmeananomaly,2)); 1:size(yearmeananomaly,2)]',.2);
                confvert(k,i,:) = bvertint(k,i,:,1) - bvert(k,i,:);
            end
        end
        
        % calculating difference percent trends
        for k = 1:length(diffnames)
            for i = 1:length(alts)
                [bvertdiff(k,i,:), bvertintdiff(k,i,:,:),~,~] = regress(squeeze(yearmeandiff(k,:,i))',...
                    [ones(1,size(yearmeandiff,2)); 1:size(yearmeandiff,2)]',.2);
                confvertdiff(k,i,:) = bvertintdiff(k,i,:,1) - bvertdiff(k,i,:);
            end
        end
        
        runnum = 2;        
        colours = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560];
        
        ph1 = createFigure('small');
        for i = 1:3
           hph(i,:) = herrorbar(squeeze(bvert(runnum,:,2))',1:60,squeeze(confvert(runnum,:,2))');
           set(hph(i,2),'LineWidth',2,'color',colours(i,:))
           set(hph(i,1),'LineWidth',1,'color',colours(i,:)) 
           hold on           
           
           if i == 1
               runnum = 3;
           elseif i == 2;
               runnum = 5;
           elseif i == 3;
               plot(zeros(1,60),1:60,'--k');
               xlabel('DU/km','fontsize',20);
               ylabel('Altitude (km)','fontsize',20);
               title(['WACCM vertical trends during ',fortitlemons],'fontsize',22);
               lh2 = legend(hph(:,2),names{[2,3,5]},'location','NorthEast');
               set(lh2,'fontsize',20,'box','off');
               
               filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/Trends/',...
                   variable,'_',num2str(latboundaries(1)),'-',num2str(latboundaries(2)),'_',fortitlemons,'_','verticaltrends.pdf'];
               export_fig(ph1,filename,'-pdf');
               
           end           
        end
        ph2 = createFigure('small');
        for i = 1:4
           hph(i,:) = herrorbar(squeeze(bvertdiff(i,:,2))',1:60,squeeze(confvertdiff(i,:,2))');
           set(hph(i,2),'LineWidth',2,'color',colours(i,:))
           set(hph(i,1),'LineWidth',1,'color',colours(i,:)) 
           hold on           
           
           if i == 1               
           elseif i == 2;               
           elseif i == 4;
               plot(zeros(1,60),1:60,'--k');
               xlabel('molecules/cm^3','fontsize',20);
               ylabel('DU/km','fontsize',20);
               title(['WACCM vertical trend contributions during ',fortitlemons],'fontsize',22);
               lh2 = legend(hph(:,2),diffnames,'location','NorthEast');
               set(lh2,'fontsize',20,'box','off');
               
               filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/Trends/',...
                   variable,'_',num2str(latboundaries(1)),'-',num2str(latboundaries(2)),'_',fortitlemons,'_','verticaltrends_diff.pdf'];
               export_fig(ph2,filename,'-pdf');
               
           end
           
        end
        
        
        %xlim([-1 1])
    case 'latitude_pressure'
        %%
        for i = 1:length(names);            
            %August
%             for l = 1:15
%                 for j = 1:12
%                     a20002002(i,l,j,:,:) = nanmean(NumberDensity_dates.(names{i})(:,:,[12+j,24+j]),3);
%                     apressure20002002(i,l,j,:,:) = nanmean(Pressure_dates.(names{i})(:,:,[12+j,24+j]),3);
% 
%                     a2008(i,l,j,:,:) = NumberDensity_dates.(names{i})(:,:,(12*9)+j);
%                     apressure2008(i,l,j,:,:) = Pressure_dates.(names{i})(:,:,(12*9)+j);
% 
%                     a2013(i,l,j,:,:) = NumberDensity_dates.(names{i})(:,:,(13*12)+j);
%                     apressure2013(i,l,j,:,:) = Pressure_dates.(names{i})(:,:,(13*12)+j);

%                   for k = 1:96
%                         a20002002int(i,l,j,k,:) = interp1(log(squeeze(apressure20002002(i,l,j,k,:))/100),squeeze(a20002002(i,l,j,k,:)),log(pres));
%                         a2008int(i,l,j,k,:) = interp1(log(squeeze(apressure2008(i,l,j,k,:))/100),squeeze(a2008(i,l,j,k,:)),log(pres));
%                         a2013int(i,l,j,k,:) = interp1(log(squeeze(apressure2013(i,l,j,k,:))/100),squeeze(a2013(i,l,j,k,:)),log(pres));
%                   end

            for l = 1:length(NumberDensity_dates.(names{i}))                  
                for k = 1:96
                    aozone(i,l,k,:) = interp1(log(squeeze(Pressure.(names{i})(k,:,l))/100),squeeze(NumberDensity_dates.(names{i})(k,:,l)),log(pres));
                end
            end
%                 end
%             end                 
        end
        
        %% Plotting
        %minmax = [-19.5 10.5];
        %minmax = [-18.5 8.5];
        minmax = [-100 10];        
        cmap = flipud(cbrewer('seq','Blues',32));
        %cmap = flipud(cbrewer('div','RdBu',37));        
        
        cmap(23:end,:) = [];
        cmap(1,:) = [.3 .3 .5];
        cmap(end,:) = [.5 .3 .3];
%         cmap(28:end,:) = [];
%         cmap(18,:) = [.9 .9 .9];
%         cmap(19,:) = [.9 .9 .9];
%         cmap(20,:) = [.9 .9 .9];
        
        cmap(20,:) = [.9 .9 .9];
        cmap(21,:) = [.9 .9 .9];
        cmap(22,:) = [.9 .9 .9];
        
        %contourstep = 1;
        contourstep = 5;
        %volcanic influence 
%         for i = 1:12
%             vola20002002(i,l,:,:) = squeeze((a20002002int(2,l,i,:,:) - a20002002int(3,l,i,:,:))./a20002002int(3,l,i,:,:)*100);
%             vola2008(i,l,:,:) = squeeze((a2008int(2,l,i,:,:) - a2008int(3,l,i,:,:))./a2008int(3,l,i,:,:)*100);
%             vola2013(i,l,:,:) = squeeze((a2013int(2,l,i,:,:) - a2013int(3,l,i,:,:))./a2013int(3,l,i,:,:)*100);             
%         end
        mon = repmat({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'},[1, 16]);
        years = 2000:2015;
        
        filename2 = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/All/',variable,'_LatPres.pdf'];

        lats = Latitudes;
         
        if volcanic
            
            figure;        
            fig = gcf;
            set(fig,'color','white','position',[100 100 840 560]); 
            count = 1;
            for i = 1:size(aozone,2)/12
                %volcanoes
                for j = 1:12
                    volc = squeeze((aozone(2,count,:,:) - aozone(3,count,:,:))./aozone(3,count,:,:)*100);
                    %volc (volc < -17.5) = -18.5;                                   
                    contourf(lats,preslog,squeeze(volc)',minmax(1):contourstep:minmax(2),'LineStyle','none');                
                    colormap(cmap);
                    ch = colorbar;
                    caxis([minmax(1) minmax(2)])
                    set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
                    set(get(ch,'ylabel'),'string','% ozone','fontsize',fsize+4)                                
                    set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),...
                            'fontsize',fsize-2);
                    ylabel('Pressure (hPa)','fontsize',fsize+2);   
                    xlabel('Latitude ({\circ}N)','fontsize',fsize+2);   
                    title(strcat(mon(count),'{ }',num2str(years(i)),' (MAM - VC-MAM)'),'fontsize',fsize+4);                
                    F(count) = getframe;                
                    %filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/',sprintf('%02d',i),'LatPres.png'];
                    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/All/',variable,'_',sprintf('%03d',count),'LatPres.png'];

                    export_fig(filename,'-png')
    %                 if count == 1;
    %                     print(filename2,'-dps');
    %                 else
    %                     print(filename2,'-dps','-append');
    %                 end
                    count = count+1;
                end
            end     
        end
        %% Checking 2010 against other years
        if check2010
            for i = 1:length(names)
                average2005(i,:,:) = squeeze(nanmean(aozone(i,10*12+1:11*12,:,:),2));
                averagerest(i,:,:) = squeeze(nanmean(cat(2,aozone(i,1:10*12,:,:),aozone(i,11*12+1:15*12,:,:)),2));
            end
        end       
        %plotting
        figure;        
        fig = gcf;
        set(fig,'color','white','position',[100 100 840 1000]);    
        
%         min2005 = floor(min(average2005(5,:,:)));
%         amx2005 = floor(max(average2005(5,:,:)));
        
        cmap2 = flipud(cbrewer('div','RdBu',11));
        cmap2(1,:) = [.3 .3 .5];
        cmap2(end,:) = [.5 .3 .3];
        cmap2(6,:) = [.9 .9 .9];
        contourstep = 5e11;
        minmax = [5e11 6e12];

        subplot(3,1,1)
        contourf(lats,preslog,squeeze(average2005(2,:,:))',minmax(1):contourstep:minmax(2),'LineStyle','none');                
        colormap(cmap2);
        ch = colorbar;                  
        caxis([minmax(1) minmax(2)])
        set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        %set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        set(get(ch,'ylabel'),'string',[variable, ' molecules/cm^3'],'fontsize',fsize)                                
        set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),...
                'fontsize',fsize-2);
        ylim([preslog(55) preslog(1)]);
        ylabel('Pressure (hPa)','fontsize',fsize+2);   
        xlabel('Latitude ({\circ}N)','fontsize',fsize+2);   
        title([variable,' 2010 MAM'],'fontsize',fsize+4);   
        
        subplot(3,1,2)
        contourf(lats,preslog,squeeze(averagerest(2,:,:))',minmax(1):contourstep:minmax(2),'LineStyle','none');                
        colormap(cmap2);
        ch = colorbar;                  
        caxis([minmax(1) minmax(2)])
        set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        %set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        set(get(ch,'ylabel'),'string',[variable, ' molecules/cm^3'],'fontsize',fsize)                                
        set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),...
                'fontsize',fsize-2);
        ylim([preslog(55) preslog(1)]);
        ylabel('Pressure (hPa)','fontsize',fsize+2);   
        xlabel('Latitude ({\circ}N)','fontsize',fsize+2);   
        title([variable,' All other years MAM'],'fontsize',fsize+4);   
        
        subplot(3,1,3)
        contourf(lats,preslog,(squeeze(average2005(2,:,:))' - squeeze(averagerest(2,:,:))')./...
            (squeeze(averagerest(2,:,:))')*100,-11:2:11,'LineStyle','none');                
        colormap(cmap2);
        ch = colorbar;                  
        caxis([-11 11])
        set(ch,'Ticks',-9:2:9)
        %set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        set(get(ch,'ylabel'),'string',['% ', variable],'fontsize',fsize)                                
        set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),...
                'fontsize',fsize-2);
        ylim([preslog(55) preslog(1)]);
        ylabel('Pressure (hPa)','fontsize',fsize+2);   
        xlabel('Latitude ({\circ}N)','fontsize',fsize+2);   
        title([variable,' 2010 percent difference from other years MAM'],'fontsize',fsize+4);   
        
        
         filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/Differences/',variable,'_','check2010','LatPres.png'];
         export_fig(fig,filename);
        
        %%        
        chem200608 = squeeze(nanmean(aozone(5,8*12+1:11*12,:,:),2));
        chem201315 = squeeze(nanmean(aozone(5,12*12+1:15*12,:,:),2));
        chem200002 = squeeze(nanmean(aozone(5,0*12+1:3*12,:,:),2));        
        
        if strcmp(variable,'O3') || strcmp(variable,'CLO')
            cmap2 = flipud(cbrewer('div','RdBu',11));
            cmap2(1,:) = [.3 .3 .5];
            cmap2(end,:) = [.5 .3 .3];
            cmap2(6,:) = [.9 .9 .9];
            contourstep = 1;
            minmax = [-5.5 5.5];
        elseif strcmp(variable,'NO2') || strcmp(variable,'NOX') || strcmp(variable,'NONO2')
            cmap2 = flipud(cbrewer('div','RdBu',11));
            cmap2(1,:) = [.3 .3 .5];
            cmap2(end,:) = [.5 .3 .3];
            cmap2(6,:) = [.9 .9 .9];
            minmax = [-88 88];
            contourstep = 16;
        end
        
        figure;        
        fig = gcf;
        set(fig,'color','white','position',[100 100 840 560]);         
        
        %Chemistry     
        
        diff13150608 = (chem201315'-chem200608')./chem200608'*100;      
        diff13150608 (diff13150608 < minmax(1)) = minmax(1);
        diff13150608 (diff13150608 > minmax(2)) = minmax(2);
        
        contourf(lats,preslog,diff13150608,minmax(1):contourstep:minmax(2),'LineStyle','none');                
        colormap(cmap2);
        ch = colorbar;        
        caxis([minmax(1) minmax(2)])
        set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        %set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        set(get(ch,'ylabel'),'string',['% ', variable],'fontsize',fsize+4)                                
        set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),...
                'fontsize',fsize-2);
        ylabel('Pressure (hPa)','fontsize',fsize+2);   
        xlabel('Latitude ({\circ}N)','fontsize',fsize+2);   
        title([variable,' 2012-2014',' - ','2008-2010',' Chem-only-fSSTs'],'fontsize',fsize+4);      
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/Differences/',variable,'_','1214_0810','LatPres.png'];
        export_fig(fig,filename);
        
        figure;        
        fig = gcf;
        set(fig,'color','white','position',[100 100 840 560]); 
        %minmax = [-5.5 5.5];                 
        
        diff00020608 = (chem200002'-chem200608')./chem200608'*100;
        diff00020608 (diff00020608 < minmax(1)) = minmax(1);
        diff00020608 (diff00020608 > minmax(2)) = minmax(2);
        
               
        contourf(lats,preslog,diff00020608,minmax(1):contourstep:minmax(2),'LineStyle','none');                
        colormap(cmap2);
        ch = colorbar;
        caxis([minmax(1) minmax(2)])
        set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        %set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        set(get(ch,'ylabel'),'string',['% ', variable],'fontsize',fsize+4)                                
        set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),...
                'fontsize',fsize-2);
        ylabel('Pressure (hPa)','fontsize',fsize+2);   
        xlabel('Latitude ({\circ}N)','fontsize',fsize+2);   
        title([variable,' 2000-2002',' - ','2008-2010',' Chem-only-fSSTs'],'fontsize',fsize+4); 
        
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/Differences/',variable,'_0002_0810','LatPres.png'];
        export_fig(fig,filename);
        
        figure;        
        fig = gcf;
        set(fig,'color','white','position',[100 100 840 560]); 
        %minmax = [-5.5 5.5];                 
                
        diff13150002 = (chem201315'-chem200002')./chem200002'*100;
        diff13150002 (diff13150002 < minmax(1)) = minmax(1);
        diff13150002 (diff13150002 > minmax(2)) = minmax(2);
        
        
        contourf(lats,preslog,diff13150002,minmax(1):contourstep:minmax(2),'LineStyle','none');                
        colormap(cmap2);
        ch = colorbar;
        caxis([minmax(1) minmax(2)])
        set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        %set(ch,'Ticks',minmax(1)+contourstep:contourstep*1:minmax(2)-contourstep)
        set(get(ch,'ylabel'),'string',['% ', variable],'fontsize',fsize+4)                                
        set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),...
                'fontsize',fsize-2);
        ylabel('Pressure (hPa)','fontsize',fsize+2);   
        xlabel('Latitude ({\circ}N)','fontsize',fsize+2);   
        title([variable,' 2012-2014',' - ','2000-2002',' Chem-only-fSSTs'],'fontsize',fsize+4);    
        
        filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/maps/LatPres/Differences/',variable,'_1214_0002','LatPres.png'];
        export_fig(fig,filename);
                                  
end


