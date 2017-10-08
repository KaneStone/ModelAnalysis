% Plot WACCM individual vertical profiles and compare to ozonesondes
clear all
close all
variable = 'O3';
lat = -90;
comparetoozonesonde = 1;
WACCM_intercompare = 0;

month = {'9','10','11','12'}; %Needs to be four separate months
%month = {'1','2','3','4'}; %Needs to be four separate months
year = '2015';

montitles = {'January','February','March','April','May','June','July','August',...
    'September','October','November','December'};

for i = 1:length(month)
    if ~strcmp(month{i},'12')
        montoplot(i) = str2double([year,sprintf('%02d',str2double(month{i})+1),'01']);
    else
        montoplot(i) = str2double([num2str(str2double(year)+1),'01','01']);
    end

    fortitletemp = [montitles{str2double(month{i})},' ',year];    
    fortitle{i} = fortitletemp;
    forfilenametemp = montitles{str2double(month{i})};   
    forfilenametemp1(i,1:3) = forfilenametemp(1:3);        
end
forfilename = [forfilenametemp1(1,:),forfilenametemp1(end,:),year];

runs = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};

%% Reading in WACCM data
commondir = '/Users/kanestone/work/projects/WACCM/netcdffiles/';
[NumberDensity, MolConc, Temperature, Pressure, Altitude, GMH, Latitudes, Longitudes, datedata] = ...
    ReadWACCMvertical(variable,'monthly',commondir,0);

[~,latind] = min(abs(Latitudes-lat));
for j = 1:length(month)
    for i = 1:length(runs)
        temp = find(datedata.(runs{i}).date == montoplot(j)); 
        if isempty(temp)
            monind(j,i) = NaN;
        else
            monind(j,i) = temp;
        end
    end
end

%% 

if comparetoozonesonde

    station = 'South_Pole';
    addpath('/Users/kanestone/work/projects/ozonesonde/code');
    directory = ['/Users/kanestone/work/projects/ozonesonde/',station,'/ecc/'];    

    [ozonesonde_data,date, next_year, station_altitude, station_latitude, station_longitude] = ...
        read_ozonesonde_raw(station,'ECC',directory,year);        
      
    %remove data point that have same pressure
    for i = 1:length(ozonesonde_data)
        for j = 1:length(ozonesonde_data(i).profiles)-1
            if ozonesonde_data(i).profiles(j,1) == ozonesonde_data(i).profiles(j+1,1)
                ozonesonde_data(i).profiles(j+1,:) = NaN;
            end            
        end
        %ozonewacvert(:,i) = lsq_lut_piecewise(log(ozonesonde_data(i).profiles(:,1)),ozonesonde_data(i).profiles(:,2),log(squeeze(Pressure.MAM(1,:,202))./100));
        for j = length(ozonesonde_data(i).profiles):-1:1
            if isnan(ozonesonde_data(i).profiles(j,:))
                ozonesonde_data(i).profiles(j,:) = [];
            end
        end
    end
    
    % piecewise least squares interpolation onto WACCM grid
    
    
    [ozonepressureND_monthmean,ozonepressurevmr_monthmean,ozoneND_monthmean,...
        ozonevmr_monthmean,standardpressure] = ...
        ozonesonde_monthlyaverage(station,ozonesonde_data,date,station_altitude);
   for i = 1:12
       ozonewacvert(:,i) = lsq_lut_piecewise(log(standardpressure(34:170))',ozonepressureND_monthmean(34:170,i),log(squeeze(Pressure.MAM(1,:,202))./100));
   end
end

ozoneinterptest = interp1(log(standardpressure),ozonepressureND_monthmean(:,10),log(squeeze(Pressure.MAM(1,:,202))./100));

plot(ozonewacvert(:,10),log(squeeze(Pressure.MAM(1,:,202))./100))
hold on
plot(ozonepressureND_monthmean(:,10),log(standardpressure))
plot(ozoneinterptest,log(squeeze(Pressure.MAM(1,:,202))./100))

station_title = station;
station_title (station_title == '_') = ' ';

%% plotting
fsize = 18;
lwidth = 2;
fig = createfig('medium');

pres = [1000:-50:150 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.05:.1];

preslog = log(pres);
prestick = [1000 500 200 150 100 50 20 10 5 2 1 .5 .2 .1];%.05 .025 .015 .01 .005 .0025 ...
    %.001 .0005 .00025 .0001];
presticklog = log(prestick);

for j = 1:length(month)   
    for i = 1:length(runs) 
        if isnan(monind(j,i))
            ozoneprofileND(:,i) = zeros(length(pres),1);
            ozoneprofileND (ozoneprofileND == 0) = NaN;
            ozoneprofilevmr(:,i) = zeros(length(pres),1);
            ozoneprofilevmr (ozoneprofilevmr == 0) = NaN;
        else
            ozoneprofileND(:,i) = interp1(log(squeeze(Pressure.(runs{i})(latind,:,monind(j,i)))/100),...
                squeeze(NumberDensity.(runs{i})(latind,:,monind(j,i))),log(pres));                                               
            ozoneprofilevmr(:,i) = interp1(log(squeeze(Pressure.(runs{i})(latind,:,monind(j,i)))/100),...
                squeeze(MolConc.(runs{i})(latind,:,monind(j,i))),log(pres));                                               
        end
    end
    sp = subplot(2,2,j);
    hold on
    if WACCM_intercompare
        ph(1) = plot(ozoneprofileND(:,2),preslog,'LineWidth',lwidth);
        ph(2) = plot(ozoneprofileND(:,3),preslog,'LineWidth',lwidth);
        ph(3) = plot(ozoneprofileND(:,5),preslog,'LineWidth',lwidth);
    elseif comparetoozonesonde
        phoz = plot(ozonepressureND_monthmean(:,str2double(month(j))),log(standardpressure),'-','LineWidth',lwidth);
        ph = plot(ozoneprofileND(:,2),preslog,'-','LineWidth',lwidth);
        ph2 = plot(ozoneprofileND(:,3),preslog,'-','LineWidth',lwidth);        
    end

    set(gca,'YDir','reverse','ytick',fliplr(presticklog(1:1:end)),'yticklabel',fliplr(prestick(1:1:end)),...
                        'fontsize',fsize,'fontsize',fsize);
    ylim([log(.9) log(1050)]);
    if j == 1 || j == 3
        ylabel('Pressure (hPa)','fontsize',fsize+2)
    end
    if j == 3 || j == 4
        xlabel('molecules/cm^3','fontsize',fsize+2)
    end
    title([fortitle{j},' at ',sprintf('%.1f',Latitudes(latind)),' ',char(176),'S'],'fontsize',fsize+4);

    if j == length(month)
        if WACCM_intercompare 
            lh = legend('MAM','VC-MAM','Chem-only-fSSTs');
        else
            lh = legend(['Ozonesondes',' at ',station_title],'MAM','VC-MAM');
        end
        set(lh,'fontsize',fsize,'location','NorthEast','box','off');
    end
    box on
end
if WACCM_intercompare
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/Profiles/',...
        'WACCM_O3_',forfilename,'_',sprintf('%.1f',Latitudes(latind)),'.pdf'];
elseif comparetoozonesonde
    filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/Profiles/',...
        'WACCM_O3_',forfilename,'_',sprintf('%.1f',Latitudes(latind)),'_',station,'.pdf'];
end

export_fig(filename,'-pdf','-nofontswap');   
