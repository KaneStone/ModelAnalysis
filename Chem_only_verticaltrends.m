% create vertical profile linear trends in DU/km for chem-only
clear all

Stimeperiod = [2000 2014];
[~,SWOOSH,~] = Read_in_netcdf('/Volumes/MyBook/work/data/SWOOSH/O3/combinedo3q_swoosh-v02.6-198401-201611-latpress-10deg-L31.nc');
SWOOSHyears = repmat(1984:2015,[12,1]);
SWOOSHyears = [SWOOSHyears(:);ones(11,1)*2016];
SWOOSHextract = permute(SWOOSH.combinedo3q(:,:,SWOOSHyears >= Stimeperiod(1) & SWOOSHyears <= Stimeperiod(2)),[2,1,3]);

%% import highCl
highCllevel = [SWOOSH.level;[.9,.8,.7,.6,.5,.4,.3,.2,.1]'];

%% import Chem-only and MAM
SDtimeperiod = [2000 2014];
[SDWaccmData,Latitudes,~,~,Altitudes,ND,Pressure] = ReadInSDWACCM(highCllevel,SDtimeperiod);

alt_chemonly = Altitudes.Chemonlynoleap(:,:,13:end-11);
ND_chemonly = ND.Chemonlynoleap(:,:,13:end-11);
NO2_chemonly = permute(SDWaccmData(2).NO2nointerp,[2,1,3]);
Pressure_chemonly = Pressure.Chemonlynoleap(:,:,13:end-11);
%% use regression to remove SPE and solar cycle


% for i = 1:size(SDWaccmData(2).O3nointerp,1)
%     for j = 1:size(SDWaccmData(2).O3nointerp,2)
%         for k = 1:12
%             if j < 49
%                 NO2 = squeeze(SDWaccmData(2).NO2nointerp(j,i,k:12:end));
%                 solar = squeeze(SDWaccmData(2).solar(k:12:end));
%                 
%                 NO2norm = 2.*(NO2 - min(NO2))./(max(NO2) - min(NO2))-1;
%                 solarnorm = 2.*(solar - min(solar))./(max(solar) - min(solar))-1;
%                 
%                 predictors = [ones(15,1),NO2norm,solarnorm];
%             else
%                 predictors = [ones(15,1),squeeze(SDWaccmData(2).solar(k:12:end))];
%             end
%             [~,~,r(i,j,k,:)] = regress(squeeze(SDWaccmData(2).O3nointerp(i,j,k:12:end))-...
%                 nanmean(squeeze(SDWaccmData(2).O3nointerp(i,j,k:12:end)),3),predictors);
%         end
%     end
% end

%% convert to DU/km

DU_coeff = 1e5*1.38e-21*1e3*(273.1/10.13);

alts = 1000:1000:60000;

for i = 1:size(ND_chemonly,1)
    for j = 1:size(ND_chemonly,3)                                        
        datainterpalt(i,:,j) = interp1(squeeze(alt_chemonly(i,:,j)),squeeze(ND_chemonly(i,:,j)),alts,'linear','extrap');                                        
        datainterpaltNO2(i,:,j) = interp1(squeeze(alt_chemonly(i,:,j)),squeeze(NO2_chemonly(i,:,j)),alts,'linear','extrap');                                        
        pressure(i,:,j) = exp(interp1(squeeze(alt_chemonly(i,:,j)),log(squeeze(Pressure_chemonly(i,:,j))),alts,'linear','extrap'))./100;
    end
end    
datainterpaltDU = datainterpalt(:,:,:).*DU_coeff;
datainterpaltDUNO2 = datainterpaltNO2(:,:,:).*DU_coeff;

%% average over latitudes
lats = [-90 -63];
latind = Latitudes >= lats(1) & Latitudes <= lats(2);
for i = 1:size(datainterpaltDU,2)
    Chem_wa(i,:) = weightedaverage(squeeze(datainterpaltDU(latind,i,:)),Latitudes(latind));
    Chem_wa_no2(i,:) = weightedaverage(squeeze(datainterpaltDUNO2(latind,i,:)),Latitudes(latind));
   
end
Chem_wa2(:,:) = nanmean(squeeze(datainterpaltDU(latind,:,:)),1);

%% put into monthly and take linear trends

for i = 1:size(Chem_wa,1)
    for j = 1:12
        b(i,j,:) = regress(squeeze(Chem_wa(i,j:12:end))',[ones(size(Chem_wa,2)/12,1),[1:size(Chem_wa,2)/12]']);
        b2(i,j,:) = regress(squeeze(Chem_wa2(i,j:12:end))',[ones(size(Chem_wa,2)/12,1),[1:size(Chem_wa,2)/12]']);
    end
end

%% remove solar and NO2 from 

for i = 1:size(datainterpaltDU,1)
    for j = 1:size(datainterpaltDU,2)
        for k = 1:12
            solar = squeeze(SDWaccmData(2).solar(k:12:end));
            solarnorm = 2.*(solar - min(solar))./(max(solar) - min(solar))-1;
            if j > 25
                NO2 = squeeze(datainterpaltDUNO2(i,j,k:12:end));
                                
                NO2norm = 2.*(NO2 - min(NO2))./(max(NO2) - min(NO2))-1;                
                
                predictors = [ones(15,1),NO2norm,solarnorm];
                %predictors = [ones(15,1),NO2norm];%,solarnorm];
            else
                NO2 = detrend(squeeze(datainterpaltDUNO2(i,j,k:12:end)));
                NO2norm = 2.*(NO2 - min(NO2))./(max(NO2) - min(NO2))-1;                
                predictors = [ones(15,1),NO2norm,solarnorm];
                %predictors = [ones(15,1),NO2norm];%,solarnorm];
            end
            [~,~,r(i,j,k,:)] = regress(squeeze(datainterpaltDU(i,j,k:12:end))-...
                nanmean(squeeze(datainterpaltDU(i,j,k:12:end))),predictors);
            
            b_speremoved(i,j,k,:) = regress(squeeze(r(i,j,k,:)),[ones(15,1),[1:15]']);
        end
    end
end

%% take weighted polar average of residual trends

for i = 1:size(b_speremoved,2)
    bwa(i,:) = weightedaverage(squeeze(b_speremoved(latind,i,:,2)),Latitudes(latind));
    pwa(i,:) = weightedaverage(squeeze(pressure(latind,i,:)),Latitudes(latind));
end

%%
fsize = 20;
mon = [8,9,10,11];
prestick = [1000:-100:100,90:-10:10,9:-1:1];
presticklabel = {1000,[],[],[],[],500,[],[],200,100,[],[],[],[],50,[],[],20,10,[],[],[],[],...
    5,[],[],2,1};
logprestick = log(prestick);

createfig('largeportrait','on');

for i = 1:length(mon)
    subplot(2,2,i)
    plot(b(1:60,mon(i),2),log(pwa(:,1)),'LineWidth',3)
    hold on
    plot(bwa(1:60,mon(i)),log(pwa(:,1)),'LineWidth',3)
    set(gca,'YDir','reverse');
    set(gca,'ytick',fliplr(logprestick),'yticklabel',fliplr(presticklabel))
    xlim([-.02 .15])
    ylim([log(.08) log(1000)])
    set(gca,'fontsize',fsize);
    ylabel('Pressure (hPa)','fontsize',fsize+2)
    xlabel('O_3 (DU/km/year)','fontsize',fsize+2)
    title(monthnames(mon(i),0,0));
    TCOtrends = sum(b(1:60,mon(i),2));
    TCOtrends2 = sum(bwa(1:60,mon(i)));
    lh = legend(['Chem-only, TCO trends = ',sprintf('%.2f',TCOtrends),' DU/year'],['solar cycle/SPEs removed, TCO trends = ',sprintf('%.2f',TCOtrends2),' DU/year']);
    set(lh,'fontsize',fsize-4,'box','off');
    
    
end

 annotation('textbox',[0 .905 1 .1],'String','90-63S, 2000-2014 ozone profile trends','FitBoxToText','on','VerticalAlignment','top','HorizontalAlignment','center','fontsize',fsize+4,...
        'EdgeColor','none','fontweight','bold')
    
 filename = '/Users/kanestone/Dropbox (MIT)/Work_Share/MITWork/solarcycleandtrends/LatPres/ProfileTrends_2000-2014.pdf';
 
 export_fig(filename,'-pdf');

