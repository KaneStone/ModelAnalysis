%Calculates vertical trends in both %/year and DU/km/year

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

%naming
names = {'CCMI','MAM','VCMAM','Chemonly','ChemonlyfixedSSTs','Chemonlynoleap'};
namestitle = {'CCMI','MAM','VC-MAM','Chem-only','Chem-only-fixedSSTs','Chem-only-noleap'};

yearmin = 20000201;
yearmax = 20151001;
remove2002 = 'Southern';

%switches
percent = 1;
convert_to_ND = 1;
latboundaries = [-60 60];
mons = 1:12;

if latboundaries(1) < 0
    ns{1} = 'S';
else ns{1} = 'N';
end
if latboundaries(2) >= 0
    ns{2} = 'N';
else ns{2} = 'S';
end

if length(mons) == 12
    montit = 'Yearly';
elseif length(mons) == 1
    montit = sprintf('%02d',mons);
else montit = sprintf('%02d',mons);
end

%% Begin read in data
for i = 1:length(files)
    
    %read in variable
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
    
    %calculating pressure    
    Pressure.(names{i}) = zeros(size(HCdata.(names{i}).PS,1),...
        size(HCdata.(names{i}).hyam,1),size(HCdata.(names{i}).PS,2));    
    hyam = permute(repmat(HCdata.(names{i}).hyam,[1,size(HCdata.(names{i}).PS)]),[2,1,3]);
    hybm = permute(repmat(HCdata.(names{i}).hybm,[1,size(HCdata.(names{i}).PS)]),[2,1,3]);
    P0 = permute(repmat(HCdata.(names{i}).P0,[size(HCdata.(names{i}).hyam,1),size(HCdata.(names{i}).PS)]),[2,1,3]);
    PS = double(permute(repmat(HCdata.(names{i}).PS,[1,1,size(HCdata.(names{i}).hyam,1)]),[1,3,2]));
    Pressure.(names{i}) = hyam .* P0 + hybm .* PS;
    
    %calculating altitude
    Altitude.(names{i}) = -(log(Pressure.(names{i})./101300.25)*8.31.*273.15)./(9.806*.0289);
        
    %converting to number density
    if convert_to_ND
        convert.(names{i}) = vmr2conc(data.(names{i}).(variable),TAdata.(names{i}).T,Pressure.(names{i}),variable,'conc');
    else 
        convert.(names{i}) = data.(names{i}).(variable);
    end
    
    %extracting date ranges   
    var20002014.(names{i}) = convert.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date < yearmax);    
    datetemp = datedata.(names{i}).date(datedata.(names{i}).date >= yearmin & ...
        datedata.(names{i}).date < yearmax);    
    Pressure20002014.(names{i}) = Pressure.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date < yearmax);    
    Altitude20002014.(names{i}) = Altitude.(names{i})(:,:,...
        datedata.(names{i}).date >= yearmin & datedata.(names{i}).date < yearmax);
    
    %Removing 2002
    switch remove2002
        case 'all_lats'
            var20002014.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN; 
            Pressure20002014.(names{i})(:,datetemp >= 20021001 & datetemp < 20031001) = NaN;             
        case 'Southern'        
             var20002014.(names{i})(data.(names{i}).lat < 0,...
                datetemp >= 20021001 & datetemp < 20031001) = NaN;             
    end
    
end

%% For pressure interpolation

pres = [1000:-50:150 ...
    100:-5:15 ...
    10:-.5:1.5 ...
    1:-.05:.15 ... 
    .1:-.005:.015 ...
    .01:-.0005:.0015 ...
    .001:-.00005:.00015];

preslog = log(pres);

prestick = [1000 500 250 100 50 25 10 5 2.5 1 .5 .25 .1 .05 .025 .01 .005 .0025 ...
    .001 .0005 .00025 .0001];

presticklog = log(prestick);

%% Calculating vertical ozone trends

% Interp onto regular km levels       
alts = 1000:1000:60000;
for k = 1:length(names)
    for i = 1:size(var20002014.(names{k}),1)
        for j = 1:size(var20002014.(names{k}),3)                                        
            datainterpalt(k,i,:,j) = interp1(squeeze(Altitude20002014.(names{k})(i,:,j)),squeeze(var20002014.(names{k})(i,:,j)),alts,'linear');                                        
        end
    end    
end
DU_coeff = 1e5*1.38e-21*1e3*(273.1/10.13);
datainterpaltDU = datainterpalt*DU_coeff;
% Calculate vertical anomalies        
for i = 1:length(names)         
    %input = squeeze(datainterpaltDU(i,:,:,9:12:end));
    input = squeeze(datainterpaltDU(i,:,:,:));
    [yearmeananomaly(i,:,:), yearmean(i,:,:)] = Verticalanomaly(input,...
        alts, latboundaries, data.(names{i}).lat, [2000 2014], [2000 2015], mons,1);
end

%% Calculate vertical percent trends
for k = 1:length(names)
    for i = 1:length(alts)
        [bvert(k,i,:), bvertint(k,i,:,:),~,~] = regress(squeeze(yearmean(k,:,i))',...
            [ones(1,size(yearmean,2)); 1:size(yearmean,2)]',.8);
    end
end
bvert (bvert == 0) = NaN;


%plotting
lwidth = 2;
fsize = 20;
figure;
set(gcf,'color','white','position',[100 100 1000 700]);
ph = plot(bvert([2,3,5],:,2)',1:60,'LineWidth',lwidth);
hold on
plot(zeros(1,2),[-1,61],'--k');
set(gca,'ytick',0:4:30,'yticklabel',0:4:30,'fontsize',fsize);
%set(gca,'xtick',-.1:.01:.1,'xticklabel',-.1:.01:.1);
set(gca,'xtick',-.1:.01:.1,'xticklabel',-.1:.01:.1);
xlabel('DU/km/year','fontsize',fsize+2)
ylabel('Altitude (km)','fontsize',fsize+2)
title([montit,' ',num2str(latboundaries(1)),'-',num2str(latboundaries(2)),'{\circ}N,',' 2000-2014 trends'],...
    'fontsize',fsize+4);
xlim([-.1 .1])
ylim([0 31])

lnames = {'CCMI','MAM','VC-MAM','Chem-only','Chem-only-fixedSSTs','Chem-only-noleap'};

lh = legend(ph,lnames{[2,3,5]},'location','SouthWest');
set(lh,'fontsize',fsize+2,'box','off');

filename = ['/Users/kanestone/Dropbox/Work_Share/MITwork/WACCMO3trends/verticalplots/Trends/',...
    variable,'_',num2str(abs(latboundaries(1))),ns{1},'_',num2str(abs(latboundaries(2))),ns{2},'_',montit,'_verticaltrends.pdf'];

export_fig(filename,'-pdf');
