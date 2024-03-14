function [outvar] = load_wrf_var(filename,varname)

ncid = netcdf.open(filename ,'NOWRITE');     

varid1 =  netcdf.inqVarID(ncid,varname);

%Read variables
outvar = squeeze(squeeze(netcdf.getVar(ncid,varid1,'double')));

%Close inputfiles
netcdf.close(ncid);     