%% CESM Z3 SAM analysis

clear all
directory = '/Volumes/My Book for Mac/work/data/CESM-CCMI/Z3/Z3SAM/';
files = dir([directory,'*.nc*']);

tic;
[datedata,PSdata] = ReadinCESMDateandPS;
toc;

% read in hybrid_height coordinates
[~,hhc,~] = Read_in_netcdf('/Volumes/My Book for Mac/work/data/CESM-CCMI/hybrid_coords/hhc.nc');

for i = 1:length(files);    
    [data(i).attributes, data(i).data,data(i).info] = Read_in_netcdf([directory,files(i).name]);
    data(i).data.Z3 = squeeze(data(i).data.Z3);            

    ap = permute(repmat(hhc.ap(end-6:end).*100000,1,size(data(i).data.Z3,1),size(data(i).data.Z3,2),size(data(i).data.Z3,4)),[2,3,1,4]);
    b = permute(repmat(hhc.b(end-6:end),1,size(data(i).data.Z3,1),size(data(i).data.Z3,2),size(data(i).data.Z3,4)),[2,3,1,4]);
    PS = permute(repmat(PSdata(i).p.PS(:,1:size(data(i).data.Z3,2),:),1,1,1,length(data(i).data.lev)),[1,2,4,3]);        
    Pressure(i).p =  ap + b.*PS;
        
end

pres = 850;
    

%% first put onto regular grid by interpolating pressure
a = 1;
for i = 1:length(data)
    for j = 1:size(data(i).data.Z3,1)
        for k = 1:size(data(i).data.Z3,2)
            %for l = 1:size(data(i).data.Z3,4)   
                a = a+1
                dataregpres1(i).d(j,k,:)  = interp1(log(squeeze(double(Pressure(i).p(j,k,:,:)))),...
                    squeeze(data(i).data.Z3(j,k,:,:)),log(double(pres)),'linear','extrap');
            %end
        end
    end
end
