if write_netcdf
%% write weighted average to netcdf
for i = 1:length(files)
    mode = netcdf.getConstant('NETCDF4');
    mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
    mode = bitor(mode,netcdf.getConstant('CLOBBER'));
    %mode = bitor(mode,netcdf.getConstant('NC_WRITE'));
    ncid = netcdf.create(strcat('/Volumes/My Book for Mac/work/data/CESM-CCMI/Z3/6090mean/',files(i).name),mode);    
    
    dimid = netcdf.defDim(ncid,'p',size(Pressureweighted(i).wa,1));
    dimid2 = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));

    varid_Z3 = netcdf.defVar(ncid,'Z3','double',[dimid dimid2]);
    varid_pressure = netcdf.defVar(ncid,'p','double',[dimid dimid2]);
    varid_time = netcdf.defVar(ncid,'time','double',[dimid2]);
    varid_date = netcdf.defVar(ncid,'date','int',[dimid2]);
    
    netcdf.endDef(ncid)

    %Writing netcdf variables
    netcdf.putVar(ncid,varid_Z3,[0,0],size(Z3weighted(i).wa),Z3weighted(i).wa);
    netcdf.putVar(ncid,varid_pressure,[0,0],size(Pressureweighted(i).wa),Pressureweighted(i).wa);    
    netcdf.putVar(ncid,varid_time,data(i).data.time);       
    netcdf.putVar(ncid,varid_date,data(i).data.date);

    netcdf.reDef(ncid)

    %Writing netcdf Attributes
    %SZA
    netcdf.putAtt(ncid,varid_Z3,'long_name','Geopotential Height (above sea level)')    
    netcdf.putAtt(ncid,varid_Z3,'units','m')
    netcdf.putAtt(ncid,varid_Z3,'mdims',1)
    netcdf.putAtt(ncid,varid_Z3,'cell_methods','time: mean')
    
    %time    
    netcdf.putAtt(ncid,varid_time,'long_name','time')    
    netcdf.putAtt(ncid,varid_time,'bounds','time_bounds')
    netcdf.putAtt(ncid,varid_time,'units','days since 1955-1-1 00:00:00')
    netcdf.putAtt(ncid,varid_time,'calendar','365_day')
    netcdf.putAtt(ncid,varid_time,'axis','T')
    
    %date    
    netcdf.putAtt(ncid,varid_date,'long_name','current date (YYYYMMDD)')    
    
    %pressure   
    netcdf.putAtt(ncid,varid_pressure,'long_name','pressure')    
    netcdf.putAtt(ncid,varid_pressure,'units','Pa')
    
    netcdf.close(ncid);
    %netcdf.putAtt(ncid,varid_SZA,'_FillValue',-9999)    
end

%% ensmean netcdf save
for i = 1:length(ensfiles)
    mode = netcdf.getConstant('NETCDF4');
    mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
    mode = bitor(mode,netcdf.getConstant('CLOBBER'));
    %mode = bitor(mode,netcdf.getConstant('NC_WRITE'));
    ncid = netcdf.create(strcat('/Volumes/My Book for Mac/work/data/CESM-CCMI/Z3/6090mean/ensmean/',ensfiles(i).name),mode);    
    
    dimid = netcdf.defDim(ncid,'p',size(ensPressureweighted(i).wa,1));
    dimid2 = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));

    varid_Z3 = netcdf.defVar(ncid,'Z3','double',[dimid dimid2]);
    varid_pressure = netcdf.defVar(ncid,'p','double',[dimid dimid2]);
    varid_time = netcdf.defVar(ncid,'time','double',[dimid2]);
    varid_date = netcdf.defVar(ncid,'date','int',[dimid2]);
    
    netcdf.endDef(ncid)

    %Writing netcdf variables
    netcdf.putVar(ncid,varid_Z3,[0,0],size(ensZ3weighted(i).wa),ensZ3weighted(i).wa);
    netcdf.putVar(ncid,varid_pressure,[0,0],size(ensPressureweighted(i).wa),ensPressureweighted(i).wa);    
    netcdf.putVar(ncid,varid_time,ensdata(i).data.time);       
    netcdf.putVar(ncid,varid_date,ensdata(i).data.date);

    netcdf.reDef(ncid)

    %Writing netcdf Attributes
    %SZA
    netcdf.putAtt(ncid,varid_Z3,'long_name','Geopotential Height (above sea level)')    
    netcdf.putAtt(ncid,varid_Z3,'units','m')
    netcdf.putAtt(ncid,varid_Z3,'mdims',1)
    netcdf.putAtt(ncid,varid_Z3,'cell_methods','time: mean')
    
    %time    
    netcdf.putAtt(ncid,varid_time,'long_name','time')    
    netcdf.putAtt(ncid,varid_time,'bounds','time_bounds')
    netcdf.putAtt(ncid,varid_time,'units','days since 1955-1-1 00:00:00')
    netcdf.putAtt(ncid,varid_time,'calendar','365_day')
    netcdf.putAtt(ncid,varid_time,'axis','T')
    
    %date    
    netcdf.putAtt(ncid,varid_date,'long_name','current date (YYYYMMDD)')    
    
    %pressure   
    netcdf.putAtt(ncid,varid_pressure,'long_name','pressure')    
    netcdf.putAtt(ncid,varid_pressure,'units','Pa')
    
    netcdf.close(ncid);
    %netcdf.putAtt(ncid,varid_SZA,'_FillValue',-9999)    
end
end
