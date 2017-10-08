% Surface temperature - O3 comparison

%% Read in netcdf

clear all
constantyears = 1;
plotting_detrend = 0;
plotting_correlations = 0;
hemisphere = 'north';
directory = ['/Volumes/MyBook/work/data/CESM-CCMI/O3/43_73hPa/',hemisphere,'/'];
files = dir([directory,'*.nc*']);
surface = 0;
[years,PSdata,yearfiles,PSfiles] = ReadinCESMDateandPS(1,1);

%% matching existing files
x = zeros(1,length(PSfiles));
for i = 1:length(PSfiles)
    for j = 1:length(files)
        if regexp(files(j).name,PSfiles(i).name)
            x(i) = 1;        
        end
    end
end
PSdata(~x) = [];

x1 = zeros(1,length(years));
for i = 1:length(years)
    for j = 1:length(files)
        if regexp(files(j).name,yearfiles(i).name)
            x1(i) = 1;        
        end
    end
end

years(~x1) = [];
%% read in hybrid_height coordinates
[~,hhc,~] = Read_in_netcdf('/Volumes/MyBook/work/data/CESM-CCMI/hybrid_coords/hhc.nc');

%% read in O3 data
for i = 1:length(files)-1;    
    [data(i).attributes, data(i).data,data(i).info] = Read_in_netcdf([directory,files(i).name]);
    data(i).data.O3 = squeeze(data(i).data.O3);  
    if i == 1
        for j = 1:length(data(i).data.lev)
            levind(j) = find(hhc.lev == data(i).data.lev(j));            
        end
        for j = 1:length(data(i).data.lat)
            latind(j) = find(PSdata(i).p.lat == data(i).data.lat(j));
        end
    end

    for j = 1:size(data(i).data.O3,2)
        O3weighted(i).wa(j,:) = weightedaverage(squeeze(data(i).data.O3(:,j,:)),data(i).data.lat);
    end
        
    ap = permute(repmat(hhc.ap(levind).*100000,1,size(data(i).data.O3,1),size(data(i).data.O3,3)),[2,1,3]);
    b = permute(repmat(hhc.b(levind),1,size(data(i).data.O3,1),size(data(i).data.O3,3)),[2,1,3]);
    PS = permute(repmat(squeeze(nanmean(PSdata(i).p.PS(:,latind,:))),1,1,length(data(i).data.lev)),[1,3,2]);        
    Pressure(i).p =  ap + b.*PS;
        
end

%% Read in surface temperature
if surface
    TSdirectory = '/Volumes/MyBook/work/data/CESM-CCMI/TS/';
    var = 'TS';
else
    TSdirectory = '/Volumes/MyBook/work/data/CESM-CCMI/T/200hPa/';
    var = 'T';
end
TSfiles = dir([TSdirectory,'*.nc*']);
if strcmp(hemisphere,'south')
    meanmonths = [12,1,2];
else
    meanmonths = [3,4];
end

for i = 1:length(TSfiles)-1;    
    [TSdata(i).attributes, TSdata(i).data,TSdata(i).info] = Read_in_netcdf([TSdirectory,TSfiles(i).name]);
    for j = 1:12
        TSdata_rearr(i).(var)(:,:,:,j) = TSdata(i).data.(var)(:,:,j:12:end);
    end
    if strcmp(hemisphere,'south')
        TSdata_MAM(i).(var) = nanmean(cat(4,TSdata_rearr(i).(var)(:,:,1:end-1,meanmonths(1)),...
            TSdata_rearr(i).(var)(:,:,2:end,meanmonths(2)),TSdata_rearr(i).(var)(:,:,2:end,meanmonths(3))),4);
    else
        TSdata_MAM(i).(var) = nanmean(cat(4,TSdata_rearr(i).(var)(:,:,:,meanmonths(1)),...
            TSdata_rearr(i).(var)(:,:,:,meanmonths(2))),4);
    end
end

%% ensemble averages
count = 1;
for i = 1:3
    O3weighted(length(files)-1+i).wa = nanmean(cat(3,O3weighted(count).wa,O3weighted(count+1).wa,O3weighted(count+2).wa),3);
    TSdata_MAM(length(files)-1+i).(var) = nanmean(cat(4,TSdata_MAM(count).(var),TSdata_MAM(count+1).(var),TSdata_MAM(count+2).(var)),4);
    years(length(files)-1+i).y = years(count).y;
    count = count+3;
end

O3weighted(length(files)-1+4).wa = nanmean(cat(3,O3weighted(count).wa,O3weighted(count+1).wa,...
    O3weighted(count+2).wa,O3weighted(count+3).wa,O3weighted(count+4).wa),3);
TSdata_MAM(length(files)-1+4).(var) = nanmean(cat(4,TSdata_MAM(count).(var),TSdata_MAM(count+1).(var),...
    TSdata_MAM(count+2).(var),TSdata_MAM(count+3).(var),TSdata_MAM(count+4).(var)),4);
years(length(files)-1+4).y = years(count).y;

runnames = {'REF-C2 no. 1','REF-C2 no. 2','REF-C2 no. 3',...
    'SEN-C2-fGHG no. 1','SEN-C2-fGHG no. 3','SEN-C2-fGHG no. 3',...
    'SEN-C2-fODS no. 1','SEN-C2-fODS no. 2','SEN-C2-fODS no. 3',...
    'REF-C1 no. 1','REF-C1 no. 2','REF-C1 no. 3','REF-C1 no. 4','REF-C1 no. 5',...
    'REF-C1SD no .1','REF-C2 ens','SEN-C2-fGHG ens','SEN-C2-fODS ens','REF-C1 ens'};
polyswitch = [1,1,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1,0,1];
%% correlations\
pastonly = 1;
futureonly = 0;
if pastonly
    TSdata_MAM(15) = [];
    O3weighted(15) = [];
    years(15) = [];
    corryears = [1960 2013];    
elseif futureonly
    TSdata_MAM([10:15,19]) = [];
    O3weighted([10:15,19]) = [];
    years([10:15,19]) = [];
    corryears = [2014 2098];
else
    corryears = [1979 2013];
end
    
if strcmp(hemisphere,'north')
    offset = 0;
else
    offset = 12;
end


if strcmp(hemisphere,'south')
    mon = 11;
else
    mon = 3;
end



count = 1
presind = 4;
ws = warning('off','all');  % Turn off warning
for i = 1:length(TSdata_MAM);
    if ~futureonly
        bp(i) = max(find(years(i).y == 1999))./12;
    else
        bp(i) = max(find(years(i).y == 2070))./12;
    end
        if constantyears
            yearstart(i) = max(find(years(i).y == corryears(1)))./12;
            yearend(i) = max(find(years(i).y == corryears(2)))./12;   
            bp(i) = bp(i) - yearstart(i);
        end
    for l = 1   
        rpoly(i).r = zeros(144,96);
        r(i).r = zeros(144,96);
        for j = 1:size(TSdata_MAM(i).(var),1)
            for k = 1:size(TSdata_MAM(i).(var),2)
                if constantyears                    
                    r(i).r(j,k) = corr(detrend(squeeze(TSdata_MAM(i).(var)(j,k,yearstart(i):yearend(i))),'linear',bp(i)),...
                        detrend(squeeze(O3weighted(i).wa(presind,yearstart(i)*12-12+mon:12:yearend(i)*12))','linear',bp(i)));
                    p = polyfit(1:numel(squeeze(TSdata_MAM(i).(var)(j,k,yearstart(i):yearend(i)))),...
                        squeeze(TSdata_MAM(i).(var)(j,k,yearstart(i):yearend(i)))', 4);
                    y = squeeze(TSdata_MAM(i).(var)(j,k,yearstart(i):yearend(i)))'...
                        - polyval(p, 1:numel(squeeze(TSdata_MAM(i).(var)(j,k,yearstart(i):yearend(i)))'));

                    if j == 1 && k == 1
                        p2 = polyfit(1:numel(squeeze(O3weighted(i).wa(presind,yearstart(i)*12-12+mon:12:yearend(i)*12))'),...
                            squeeze(O3weighted(i).wa(presind,yearstart(i)*12-12+mon:12:yearend(i)*12)),4);
                        y2 = squeeze(O3weighted(i).wa(presind,yearstart(i)*12-12+mon:12:yearend(i)*12))...
                            - polyval(p2,1:numel(squeeze(O3weighted(i).wa(presind,yearstart(i)*12-12+mon:12:yearend(i)*12))'));                            
                    end
                    rpoly(i).r(j,k) = corr(y',y2');
                    count = count + 1
%                         if i == 15 
%                             fig = createfig('small','off');
%                             plot(detrend(squeeze(TSdata_MAM(i).(var)(j,k,yearstart(i):yearend(i))),'linear',bp(i)))
%                             hold on
%                             plot(detrend(squeeze(O3weighted(i).wa(presind,yearstart(i)*12-12+mon:12:yearend(i)*12))','linear',bp(i))/1e-7)
%                             export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/detrendedlineplots/',...
%                                             var,'_',runnames{i},'_',num2str(TSdata(i).data.lon(j)),'E_',num2str(TSdata(i).data.lat(k)),'S'],'-png');  
%                             close(fig);
%                         end
                    
                    %r(i).r(j,k) = corr(detrend(squeeze(TSdata_MAM(i).(var)(j,k,yearstart(i):yearend(i)))),...
                    %    detrend(squeeze(O3weighted(i).wa(presind,yearstart(i)*12-12+mon:12:yearend(i)*12-offset))','linear',bp(i)+1));
                    
                else
                    if strcmp(hemisphere,'south')
                        %regexp(runnames{i},'REF-C1');
                        if polyswitch(i)
                            % detrending with a polynomial    
                            if l == 1
                                p1 = polyfit(years(i).y(1+offset:12:bp(i)*12),squeeze(TSdata_MAM(i).(var)(j,k,1:bp(i)-offset/12)),3);
                                x1 = years(i).y(1+offset:12:bp(i)*12);
                                y1 = polyval(p1,x1);                                                

                                p2 = polyfit(years(i).y(bp(i)*12+12:12:end),squeeze(TSdata_MAM(i).(var)(j,k,bp(i):end)),2);                        
                                x2 = years(i).y(bp(i)*12+12:12:end);
                                y2 = polyval(p2,x2);

                                if plotting_detrend
                                    fig = createfig('medium','off');
                                    plot(years(i).y(1+offset:12:end),squeeze(TSdata_MAM(i).(var)(j,k,1:end)));
                                    hold on
                                    plot(x1,y1);                        
                                    %plot(years(i).y(bp(i)*12+12:12:end-offset),squeeze(TSdata_MAM(i).(var)(j,k,bp(i)+1:end)));
                                    plot(x2,y2);

                                    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/polylineplots/',...
                                        var,'_',runnames{i},'_',num2str(TSdata(i).data.lon(j)),'E_',num2str(TSdata(i).data.lat(k)),'S'],'-png');  
                                    close(fig);
                                end

                                Temp_detrend(i).a(j,k,:) = [squeeze(TSdata_MAM(i).(var)(j,k,1:bp(i)-offset/12))-y1;squeeze(TSdata_MAM(i).(var)(j,k,bp(i):end))-y2];
                            end

                            %ozone
                            if j == 1 && k == 1
                                p1 = polyfit(years(i).y(1:12:bp(i)*12-12),squeeze(O3weighted(i).wa(l,mon:12:bp(i)*12-12))',3);
                                x1 = years(i).y(1:12:bp(i)*12-12);
                                y1 = polyval(p1,x1);                                                

                                p2 = polyfit(years(i).y(bp(i)*12:12:end-12),squeeze(O3weighted(i).wa(l,bp(i)*12-12+mon:12:end-12))',2);                        
                                x2 = years(i).y(bp(i)*12:12:end-12);
                                y2 = polyval(p2,x2);
                                %if plotting_detrend
                                    fig = createfig('medium','off');
                                    plot(years(i).y(1+offset:12:end),O3weighted(i).wa(l,mon:12:end-12)');
                                    hold on
                                    plot(x1,y1);                        
                                    %plot(years(i).y(bp(i)*12+12:12:end-offset),squeeze(TSdata_MAM(i).(var)(j,k,bp(i)+1:end)));
                                    plot(x2,y2);

                                    export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/polylineplots/',...
                                        hemisphere,'O3_',runnames{i},'_',num2str(squeeze(Pressure(1).p(1,l,1))/100),'hPa'],'-png');  
                                    close(fig);
                                %end

                                O3_detrend(l,:) = [squeeze(O3weighted(i).wa(l,mon:12:bp(i)*12-12))'-y1;squeeze(O3weighted(i).wa(l,bp(i)*12-12+mon:12:end-12))'-y2];
                            end

                            r(i).r(l,j,k) = corr(squeeze(Temp_detrend(i).a(j,k,:)),O3_detrend(l,:)');
                            count = count+1
                            if plotting_correlations
                                fig = createfig('medium','off');
                                plot(years(i).y(1:12:end-12),squeeze(Temp_detrend(i).a(j,k,:)));
                                hold on
                                plot(years(i).y(1:12:end-12),O3_detrend(l,:)./1e-7);
                                legend('T',[sprintf('%1.0f',squeeze(Pressure(1).p(1,l,1)/100)),' O3'])
                                export_fig(['/Users/kanestone/Dropbox/Work_Share/MITwork/cesm/correlations/lineplots/',...
                                    hemisphere,'_',var,'_',runnames{i},'_',num2str(TSdata(i).data.lon(j)),'E_',num2str(TSdata(i).data.lat(k)),'S_',num2str(squeeze(Pressure(1).p(1,l,1))/100),'hPa'],'-png')  
                                close(fig);
                            end

                        else
                            r(i).r(l,j,k) = corr(detrend(squeeze(TSdata_MAM(i).(var)(j,k,:)),'linear',bp(i)+1),...
                                detrend(squeeze(O3weighted(i).wa(l,mon:12:end-12))','linear',bp(i)+1));
                        end
                    else
                        r(i).r(l,j,k) = corr(detrend(squeeze(TSdata_MAM(i).(var)(j,k,:))),...
                            detrend(squeeze(O3weighted(i).wa(l,mon:12:end))','linear',bp(i)+1));
                    end
                end
            end
        end    
    end
    clearvars O3_detrend
end
    
%% saving data
latitude = TSdata(1).data.lat;
longitude = TSdata(1).data.lon;
if surface
    ext = 'TS';
else
    ext = '200hPa';
end
if constantyears
    if pastonly
        save(['/Volumes/MyBook/work/data/CESM-CCMI/O3/output/',hemisphere,'_',ext,'_O3_correlations_constantyears_past'],'latitude','longitude','r','rpoly','years','O3weighted','TSdata_MAM');
    elseif futureonly
        save(['/Volumes/MyBook/work/data/CESM-CCMI/O3/output/',hemisphere,'_',ext,'_O3_correlations_constantyears_future'],'latitude','longitude','r','rpoly','years','O3weighted','TSdata_MAM');
    else
        save(['/Volumes/MyBook/work/data/CESM-CCMI/O3/output/',hemisphere,'_',ext,'_O3_correlations_constantyears_SD'],'latitude','longitude','r','years','O3weighted','TSdata_MAM');
    end
else
    save(['/Volumes/MyBook/work/data/CESM-CCMI/O3/output/',hemisphere,'_',ext,'_O3_correlations_poly'],'latitude','longitude','r','years','O3weighted','TSdata_MAM');
end

%% test line plot
figure
for i = 1:length(TSfiles)-1;     
     plot(detrend(squeeze(TSdata_MAM(i).(var)(5,75,:)),'linear',bp(i)+1))
     hold on
end

figure
for i = 1:length(TSfiles)-1;     
     plot(squeeze(TSdata_MAM(i).(var)(5,75,:)))
     hold on
end

figure
for i = 1:length(TSfiles)-1;     
     plot(detrend(squeeze(O3weighted(i).wa(4,mon:12:end))','linear',bp(i)+1))
     hold on
end

figure
for i = 1:length(TSfiles)-1;     
     plot(squeeze(O3weighted(i).wa(4,mon:12:end))');
     hold on
end


%% test contourf plot
for i = 1:4
    figure
    contourfm(TSdata(1).data.lat,TSdata(1).data.lon,double(squeeze(r(15+i).r(4,:,:))'))
    colorbar
    load coastlines
    geoshow(coastlat, coastlon, 'Color', 'black','LineWidth',2)
end


%% test line plot
plot(O3weighted(10).wa(4,mon:12:end)-nanmean(O3weighted(10).wa(4,mon:12:end),2)+nanmean(squeeze(TSdata_MAM(10).(var)(135,75,:))));
hold on
plot(squeeze(TSdata_MAM(10).(var)(135,75,:)))
    