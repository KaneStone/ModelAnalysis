% Plot WACCM time series for particular month and pressure level.

%% Define inputs

variable = 'O3'; %O3, NO2, NO, NOX, CLO

%colormaps
cbrew = cbrewer('qual','Set1',10);
cbrew([5,7],:) = [];
cbrewfade = cbrewer('qual','Pastel1',10);
cmap = flipud(cbrewer('div','RdBu',10));

names = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
namestitle = {'CCMI','MAM','VC-MAM','Chem-only','Chem-only-fixedSSTs','Chem-only-noleap'};

yearmin = 20000201;
yearmax = 20151201;
remove2002 = 'no';

%lats = -75;
lats = [-90,-85];
szlat = size(lats);
convert_to_ND = 1;

months = 10;

%preslev = [10,30,50];
preslev = [50,70,100];

%% Reading in WACCM data
[NumberDensity, MolConc, Pressure, Altitude, GMH, Latitudes, datedata] = ReadWACCMvertical(variable);

% extracting latitude ranges and dates 

%extracting dates
for i = 1:length(names)
    
    %finding latitudes
    if numel(lats) > 1
        for p = 1:szlat(1)
            latindex(p).p = find(Latitudes >= lats(p,1) & Latitudes <= lats(p,2));
        end
    else
        [~,latindex] = min(abs(Latitudes - lats));
    end
    
    %extracting dates
    NDdates.(names{i}) = NumberDensity.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date <= yearmax);
    MolConcdates.(names{i}) = MolConc.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date <= yearmax);
    datetemp = datedata.(names{i}).date(datedata.(names{i}).date >= yearmin & ...
        datedata.(names{i}).date <= yearmax);
    Pressuredates.(names{i}) = Pressure.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date <= yearmax);
    Altitudedates.(names{i}) = Altitude.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date <= yearmax);

    switch remove2002
        case 'all_lats'
            NDdates.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            MolConcdates.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            Pressuredates.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            Altitudedates.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        case 'Southern'        
            NDdates.(names{i})(Latitudes < 0,...
                datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            MolConcdates.(names{i})(Latitudes < 0,...
                datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            Altitudedates.(names{i})(Latitudes < 0,...
                datetemp >= 20021001 & datetemp < 20031001) = NaN; 
    end
    NDdates.(names{i}) = cat(3,NDdates.(names{i}),NDdates.(names{i})(:,:,end-1:end));
    MolConcdates.(names{i}) = cat(3,MolConcdates.(names{i}),MolConcdates.(names{i})(:,:,end-1:end));
    Pressuredates.(names{i}) = cat(3,Pressuredates.(names{i}),Pressuredates.(names{i})(:,:,end-1:end))./100;
    Altitudedates.(names{i}) = cat(3,Altitudedates.(names{i}),Altitudedates.(names{i})(:,:,end-1:end));
    
    [NDdatesint.(names{i}), ~] = presinterpWACCM(Pressuredates.(names{i}),NDdates.(names{i}));        
    [MolConcdatesint.(names{i}), presint] = presinterpWACCM(Pressuredates.(names{i}),MolConcdates.(names{i}));        
    
    if latindex(1).p > 1
        NDdateslatmean.(names{i}) = squeeze(nanmean(NDdates.(names{i})(:,latindex(1).p,:),2));
    end
    
end


    
%% finding presindex
for i = 1:length(preslev)
    [~,presindex(i)] = min(abs(presint - preslev(i)));
end

preslog = log(presint);

prestick = [1000 500 250 100 50 25 10 5 2.5 1 .5 .25 .1 .05 .025 .01 .005 .0025 ...
    .001 .0005 .00025 .0001];

presticklog = log(prestick);

%% plotting
monthstitle = {'January','February','March','April','May','June','July','August','September','October','November','December'};

fig = createfig('large');
lwidth = 2;
fsize = 20;
for i =1:length(presindex)
    subplot(length(presindex),1,i)
    ph = plot(squeeze(nanmean(NDdatesint.(names{2})(latindex.p,presindex(i),months:12:end),1)),'LineWidth',lwidth);
    hold on
    ph1 = plot(squeeze(nanmean(NDdatesint.(names{3})(latindex.p,presindex(i),months:12:end),1)),'LineWidth',lwidth);
    %ph = plot(squeeze(NDdatesint.(names{5})(latindex,presindex(i),:)),'LineWidth',lwidth);
    if i == length(presindex)
        xlabel('Year','fontsize',fsize+2);
    end
    if strcmp(variable,'NOX')
        ylabel('mol/mol','fontsize',fsize+2);
    elseif strcmp(variable,'O3')
        ylabel('ND','fontsize',fsize+2);
    end
    if i == length(presindex)
        lh = legend(namestitle{2:3});
        set(lh,'fontsize',fsize,'box','off')
    end
    
    title([monthstitle{months},', ',num2str(preslev(i)),' hPa, ',num2str(abs(lats)),'{\circ}S'],'fontsize',fsize+4)
    set(gca,'xtick',1:2:20,'xticklabel',2000:2:2020,'fontsize',fsize);
end
%pmtit = mtit([variable,', ',namestitle{5}],'fontsize',fsize+6,'color','k','xoff',0,'yoff',.045);

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/TimeSeries/',variable,'_',monthstitle{months},'1007050hPa_0985S.pdf'];
export_fig(fig,filename);
