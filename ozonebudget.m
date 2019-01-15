%WACCM ozone budget
clear all
clc

reactionrates = 0;

if reactionrates
    variable = 'OddOx_lossrates'; %O3, NO2, NO, NOX, CLO
    allvars = {'OddOx_CLOxBROx_Loss','OddOx_HOx_Loss','OddOx_Loss_Tot',...
        'OddOx_NOx_Loss','OddOx_Ox_Loss','OddOx_Prod_Tot'};

    allvarslegend = {'60-35S','35-60N','Halogens','HOx','NOx','Ox'};

    directory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/reactionrates/'];
    files = dir([directory,variable,'*']);
else
    variable = 'O3';
    directory = ['/Users/kanestone/work/projects/WACCM/netcdffiles/',variable,'/'];
    files = dir([directory,variable,'*']);
end

TAdirectory = '/Users/kanestone/work/projects/WACCM/netcdffiles/Temperature/';
TAfiles = dir([TAdirectory,'TA*']);

HCdirectory = '/Users/kanestone/work/projects/WACCM/netcdffiles/hybrid_ht_conversion/';
HCfiles = dir([HCdirectory,'hybrid*']);

datedir = '/Users/kanestone/work/projects/WACCM/netcdffiles/date/';
datefiles = dir([datedir,'date*']);

for i = 1:length(HCfiles)
    if ~isempty(regexp(HCfiles(i).name,'1990', 'once'))
        HCfiles(i) = [];       
        break
    end
end

for i = 1:length(TAfiles)
    if ~isempty(regexp(TAfiles(i).name,'1990', 'once'))
        TAfiles(i) = [];    
        break
    end
end

for i = 1:length(datefiles)
    if ~isempty(regexp(datefiles(i).name,'1990', 'once'))
        datefiles(i) = []; 
        break
    end
end

%colormaps
cbrew = cbrewer('qual','Set1',10);
cbrew([5,7],:) = [];
cbrewfade = cbrewer('qual','Pastel1',10);
cmap = flipud(cbrewer('div','RdBu',10));

names = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
namestitle = {'CCMI','MAM','VC-MAM','Chem-only','Chem-only-fixedSSTs','Chem-only-noleap'};
%names = {'CCMI','ChemonlyfixedSSTs','Chemonlynoleap'};
type = 'latitude_pressure'; %trend, difference, raw,'latitude_pressure
yearmin = 20000201;
yearmax = 20151201;
remove2002 = 'Southern';

latitudes = [[-90 -80];[80 90]];
szlat = size(latitudes);
percent = 1;
convert_to_ND = 0;
volcanic = 0;

montouse = 1:12;
mons = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
if length(montouse) == 1
    montitle = mons(montouse,:);
else montitle = 'yearly';
end
zonalmean = zeros(size(latitudes,1),12,6,88,16);
zonalmean (zonalmean == 0) = NaN;

for i = 1:length(files)
    %read in O3
    [info.(names{i}), data.(names{i}), attributes.(names{i})] = ...
        Read_in_netcdf([directory,files(i).name]);        
    
    %read in temperature
    [TAinfo.(names{i}), TAdata.(names{i}), TAattributes.(names{i})] = ...
        Read_in_netcdf([TAdirectory,TAfiles(i).name]);
    
    %read in hybrid_ht conversion parameters
    [HCinfo.(names{i}), HCdata.(names{i}), HCattributes.(names{i})] = ...
        Read_in_netcdf([HCdirectory,HCfiles(i).name]);
    
    %read in dates
    [dateinfo.(names{i}), datedata.(names{i}), dateattributes.(names{i})] = ...
        Read_in_netcdf([datedir,datefiles(i).name]);
    
    %finding latitudes
    for p = 1:szlat(1)
        latindex(p).p = find(data.(names{i}).lat >= latitudes(p,1) & ...
            data.(names{i}).lat <= latitudes(p,2));
    end
        
    %calculating to pressure
    Pressure.(names{i}) = zeros(size(HCdata.(names{i}).PS,1),...
        size(HCdata.(names{i}).hyam,1),size(HCdata.(names{i}).PS,2));    
    hyam = permute(repmat(HCdata.(names{i}).hyam,[1,size(HCdata.(names{i}).PS)]),[2,1,3]);
    hybm = permute(repmat(HCdata.(names{i}).hybm,[1,size(HCdata.(names{i}).PS)]),[2,1,3]);
    P0 = permute(repmat(HCdata.(names{i}).P0,[size(HCdata.(names{i}).hyam,1),size(HCdata.(names{i}).PS)]),[2,1,3]);
    PS = double(permute(repmat(HCdata.(names{i}).PS,[1,1,size(HCdata.(names{i}).hyam,1)]),[1,3,2]));
    Pressure.(names{i}) = hyam .* P0 + hybm .* PS;
    
    %Altitude.(names{i}) = -(log(Pressure.(names{i})./PS)*1.38066e-23.*TAdata.(names{i}).T)./((28.97*1.66e-27)*9.80616);        
    
    Altitude.(names{i}) = -(log(Pressure.(names{i})./PS)*8.31.*273.15)./(9.806*.0289);        
    
    %converting to number density
    if convert_to_ND
        for j = 1:length(allvars)        
            convert.(names{i}) = vmr2conc(data.(names{i}).(variable),TAdata.(names{i}).T,Pressure.(names{i}),variable,'conc');
        end
    else 
        convert.(names{i}) = data.(names{i});
    end
    %extracting variables 
    Pressure20002014.(names{i}) = Pressure.(names{i})(:,:,...
            datedata.(names{i}).date >= yearmin & datedata.(names{i}).date < yearmax);
        
    Altitude20002014.(names{i}) = Altitude.(names{i})(:,:,...
            datedata.(names{i}).date >= yearmin & datedata.(names{i}).date < yearmax);
    
    Pressure20002014.(names{i}) = cat(3,Pressure20002014.(names{i}),Pressure20002014.(names{i})(:,:,end-1:end));
    Altitude20002014.(names{i}) = cat(3,Altitude20002014.(names{i}),Altitude20002014.(names{i})(:,:,end-1:end));
        
    for j = 1:length(allvars)
        %extracting dates
        var20002014.(names{i}).(allvars{j}) = convert.(names{i}).(allvars{j})(:,:,...
            datedata.(names{i}).date >= yearmin & datedata.(names{i}).date < yearmax);
        datetemp = datedata.(names{i}).date(datedata.(names{i}).date >= yearmin & ...
            datedata.(names{i}).date < yearmax);
                        
        switch remove2002
            case 'all_lats'
                var20002014.(names{i}).(allvars{j})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
                Pressure20002014.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
                Altitude20002014.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            case 'Southern'        
                 var20002014.(names{i}).(allvars{j})(data.(names{i}).lat < 0,...
                    datetemp >= 20021001 & datetemp < 20031001) = NaN; 
                Altitude20002014.(names{i})(data.(names{i}).lat < 0,...
                    datetemp >= 20021001 & datetemp < 20031001) = NaN; 
        end
        var20002014.(names{i}).(allvars{j}) = cat(3,var20002014.(names{i}).(allvars{j}),var20002014.(names{i}).(allvars{j})(:,:,end-1:end));
    end            
end

%% 
cbrewqual = zeros(10,3);
cbrewqual(3:10,:) = cbrewer('qual','Set1',8);

%cbrewqual(6,1:2) = .6
sollats = [-90 -80; 80 90];
vartouse = 'ChemonlyfixedSSTs';
vartousetitle = 'Chem-only-fixedSSTs';


if strcmp(variable,'O3')
    yinterval = 1;
else
    yinterval = 5;
end

%% 
if reactionrates
for k = 1:4
    if k == 1
        %solpres = [[.0005 .005]; [2.5 5]; [25 50]; [25 50]]*100;
        solpres = [[25 50]; [25 50]]*100;
    elseif k == 2
        solpres = [[2.5 5]; [2.5 5]]*100;
    elseif k == 3
        solpres = [[.25 .5]; [.25 .5]]*100;
    elseif k == 4
        solpres = [[.025 .05]; [.025 .05]]*100;
    end
    for l = 1:length(allvars)
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
            
            solpreslatmean(i,:,:) = nanmean(Pressure20002014.(vartouse)(latindex,:,:),1);

            solpindex(i).t = find(nanmean(solpreslatmean(i,:,:),3) > solpres(i,1) & ...
                nanmean(solpreslatmean(i,:,:),3) < solpres(i,2));
            
            %Calculate fractional Ox loss from catalytic cycle            
            fracloss.(vartouse).(allvars{l}) = var20002014.(vartouse).(allvars{l})./var20002014.(vartouse).(allvars{3});        
            
            solpressure(i,l).t = squeeze(nanmean(fracloss.(vartouse).(allvars{l})(:,solpindex(i).t,:),2));
            
            [yearmeananomaly(k,l,i,:), yearmean(k,l,i,:)] = TCOanomaly(solpressure(i,l).t, ...
                sollats(i,:), data.(vartouse).lat, [2000 2014], [2000 2015], montouse);
            
            sollatmean(k,l,i,:) = squeeze(nanmean(solpressure(i,l).t(latindex,:)));
            
            
        end
    end
end
end

%% Vertical Profiles
pres = [1000:-50:150 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.05:.15 ... 
    .1:-.005:.015 ...
    .01:-.0005:.0015 ...
    .001:-.00005:.00015];

preslog = log(pres);
%prestick = [1000 500 250 100 50 25 10 5 2.5 1 .5 .25 .1 .05 .025 .01 .005 .0025 ...
%    .001 .0005 .00025 .0001];
prestick = [1000 100 10 1 .1 .01 ...
    .001 .0001];
presticklog = log(prestick);

clearvars latindex

latindex(1,:) = find(data.(vartouse).lat >= sollats(1,1) & ...
    data.ChemonlyfixedSSTs.lat <= sollats(1,2));     
latindex(2,:) = find(data.(vartouse).lat >= sollats(2,1) & ...
    data.ChemonlyfixedSSTs.lat <= sollats(2,2));     
Presforinterp(1,:,:) = nanmean(Pressure20002014.(vartouse)(latindex(1,:),:,:),1);
Presforinterp(2,:,:) = nanmean(Pressure20002014.(vartouse)(latindex(2,:),:,:),1);

Presforinterp200608(1,:) = nanmean(Presforinterp(1,:,8*12+1:11*12),3);
Presforinterp200608(2,:) = nanmean(Presforinterp(2,:,8*12+1:11*12),3);
Presforinterp201315(1,:) = nanmean(Presforinterp(1,:,12*12+1:15*12),3);
Presforinterp201315(2,:) = nanmean(Presforinterp(2,:,12*12+1:15*12),3);
Presforinterp200002(1,:) = nanmean(Presforinterp(1,:,0*12+1:3*12),3);
Presforinterp200002(2,:) = nanmean(Presforinterp(2,:,0*12+1:3*12),3);

for l = 1:length(allvars)
    frac = fracloss.(vartouse).(allvars{l});    
    
    chem200608 = squeeze(nanmean(frac(:,:,8*12+1:11*12),3));
    chem201315 = squeeze(nanmean(frac(:,:,12*12+1:15*12),3));
    chem200002 = squeeze(nanmean(frac(:,:,0*12+1:3*12),3));    
        
    
    for i = 1:length(sollats);
        lats = data.(vartouse).lat;
        
        W = cosd(lats);
        W (lats < sollats(i,1) | lats > sollats(i,2)) = [];
        
        for j = 1:size(frac,2)
            
            temp = chem200608(:,j);
            temp2 = chem201315(:,j);
            temp3 = chem200002(:,j); 
            
            temp(lats < sollats(i,1) | lats > sollats(i,2)) = []; 
            temp2(lats < sollats(i,1) | lats > sollats(i,2)) = []; 
            temp3(lats < sollats(i,1) | lats > sollats(i,2)) = []; 
            
            W1 = W;           
            W1 (isnan(temp)) = NaN;
            
            frac200608(l,i,j) = nansum(temp.*W1)./nansum(W1);
            frac201315(l,i,j) = nansum(temp2.*W1)./nansum(W1);
            frac200002(l,i,j) = nansum(temp3.*W1)./nansum(W1);                        
        end
        %interpolate to regular pressure levels
        frac200608interp(l,i,:) = interp1(log(squeeze(Presforinterp200608(i,:)/100)),squeeze(frac200608(l,i,:)),preslog);
        frac201315interp(l,i,:) = interp1(log(squeeze(Presforinterp200608(i,:)/100)),squeeze(frac201315(l,i,:)),preslog);
        frac200002interp(l,i,:) = interp1(log(squeeze(Presforinterp200608(i,:)/100)),squeeze(frac200002(l,i,:)),preslog);
    end    
end

frac200608interp = cat(1,frac200608interp(1:2,:,:),frac200608interp(4:5,:,:));
frac201315interp = cat(1,frac201315interp(1:2,:,:),frac201315interp(4:5,:,:));
frac200002interp = cat(1,frac200002interp(1:2,:,:),frac200002interp(4:5,:,:));

fracdiff1 = frac200002interp-frac200608interp;
fracdiff2 = frac201315interp-frac200608interp;
fracdiff3 = frac201315interp-frac200002interp;

diff = 1;
if diff
    fracinterp = cat(4,fracdiff1,fracdiff2,fracdiff3);
else
    fracinterp = cat(4,frac200608interp,frac201315interp,frac200002interp);
end

%% Plotting vertical profiles
if diff
    years = {'2000 to 2002 - 2008 to 2010','2012 to 2014 - 2008 to 2010','2012 to 2014 - 2000 to 2002'};
else
    years = {'2000 - 2002','2008 - 2010',' 2012 - 2014'};
end
fsize = 14;
figure;
fig = gcf;
set(fig,'color','white','position',[100 100 800 800]);
count = 1;
for j = 1:3
    for i = 1:2;
        subplot(3,2,count);
        if diff            
            plot(zeros(1,length(preslog)),preslog,'--k','LineWidth',1)
            hold on
        end
        ph = plot(squeeze(fracinterp(:,i,:,j)),preslog,'LineWidth',2);           
        set(gca,'YDir','reverse','ytick',fliplr(presticklog),'yticklabel',fliplr(prestick),'fontsize',fsize);
        if diff
            xlim([-.07 .07])  
            set(gca,'xtick',-.06:.02:.06);
        else
            xlim([0 1])
        end
        ylim([presticklog(end) presticklog(1)])
        if i == 1;
            title([years{j},', ',num2str(abs(sollats(i,1))),' - ',num2str(abs(sollats(i,2))),'S'],'fontsize',fsize);
        else
            title([years{j},', ',num2str(abs(sollats(i,1))),' - ',num2str(abs(sollats(i,2))),'N'],'fontsize',fsize);
        end
        
        if mod(i,2) ~= 0
            ylabel('Pressure (hPa)','fontsize',fsize)
        end
        if count == 5 || count == 6;
            xlabel('Fraction loss rate','fontsize',fsize)
        end
        
        count = count+1;
    % ylim([
    end
end

if diff
    pmtit = mtit([vartousetitle,' - ', 'differences in fractional loss rates'],...
             'fontsize',fsize+4,'color','k',...
             'xoff',0,'yoff',.05);
else
    pmtit = mtit([vartousetitle,' - ', 'differences in fractional loss rates'],...
             'fontsize',fsize+4,'color','k',...
             'xoff',0,'yoff',.05);
end

%creating axes for title
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
    'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
vertlegend = {'Halogens','HOx','NOx','Ox'};
lh = legend(ph,vertlegend);        
set(lh,'fontsize',14,'box','off','Orientation','Horizontal','position',[0.45, .02, .01, .01]);

if diff
    filename = strcat('/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/Solarcycleanalysis/',...
    vartouse,'_',variable,'_',montitle,'_frac_diff.pdf');
else
filename = strcat('/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/lineplots/Solarcycleanalysis/',...
    vartouse,'_',variable,'_',montitle,'_frac.pdf');
end
export_fig(fi,filename);


%% Plotting time series
cbrewqual = cbrewer('qual','Set1',6);
plot_frac_time_series = 0;
if plot_frac_time_series

    figure;
    fi = gcf;
    set(fi,'color','white','position',[100 100 1000 800]);

    yearmean2 = cat(2,yearmean(:,1:2,:,:),yearmean(:,4:5,:,:));
    for k = 1:size(yearmean2,1)     
        if k == 1
            %solpres = [[.0005 .005]; [2.5 5]; [25 50]; [25 50]]*100;
            solpres = [[25 50]];
        elseif k == 2
            solpres = [[2.5 5]];
        elseif k == 3
            solpres = [[.25 .5]];
        elseif k == 4
            solpres = [[.025 .05]];
        end


            sp = subplot(2,2,k);
            sp_pos = get(sp,'position');      
            for l = 1:size(yearmean2,2)
                if l == 1
                    ph1(:,1) = plot(squeeze(yearmean2(k,l,1,:))','color','k','LineWidth',2);  
                    hold on
                    ph1(:,2) = plot(squeeze(yearmean2(k,l,2,:))','--','color','k','LineWidth',2);                    
                end
                ph1(:,l+2) = plot(squeeze(yearmean2(k,l,1,:))','color',cbrewqual(l,:),'LineWidth',2);                                                       
                hold on
                ph2(:,l+2) = plot(squeeze(yearmean2(k,l,2,:))','--','color',cbrewqual(l,:),'LineWidth',2);                          

            end

            plot(repmat([1,2,6,7],10,1),repmat(0:.1:.9,4,1)','--k');

            if k == 1
                set(sp,'position',[sp_pos(1) sp_pos(2)-.025 sp_pos(3:4)]);                    
                ylabel('Fractional Ox loss','fontsize',16);
            elseif k == 2
                set(sp,'position',[sp_pos(1) sp_pos(2)-.025 sp_pos(3:4)]);                    
            elseif k == 3            
                set(sp,'position',[sp_pos(1) sp_pos(2)+.025 sp_pos(3:4)]);                    
                ylabel('Fractional Ox loss','fontsize',16);
                xlabel('Year','fontsize',16);
            elseif k == 4
                set(sp,'position',[sp_pos(1) sp_pos(2)+.025 sp_pos(3:4)]);                    
                xlabel('Year','fontsize',16);
            end

            set(gca,'xtick',1:2:16,'xticklabel',2000:2:2050,'fontsize',14);
            set(gca,'ytick',0:.1:1)


            xlim([0 16]);
            title([num2str(solpres(1)),' to ',num2str(solpres(2)),' hPa'],'fontsize',18);
            if k == 4
                %lh = columnlegend(2,solarlinenames,'fontsize',18,'boxoff','location','user');   
            end


    end
    pmtit = mtit([vartousetitle,' - yearly average ', 'fractional loss rates'],...
             'fontsize',22,'color','k',...
             'xoff',0,'yoff',.075);

    %creating axes for title
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
        'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

    lh = legend(ph1,allvarslegend);        
    set(lh,'fontsize',14,'box','off','Orientation','Horizontal','position',[0.45, .05, .01, .01]);

    
    filename = [vartouse,'_',variable,'_',montitle,'_frac.pdf'];
    filepath = ['/Users/kanestone/Dropbox/Work_Share/MITwork/solarcycle/lineplots/',filename]
    export_fig(fi,filepath);
    scriptname = mfilename;
    
    creatlogfile(filename,scriptname,'/Users/kanestone/Dropbox/Work_Share/MITwork/solarcycle/logfiles/');

end
