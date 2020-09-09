# Print a summary of surrogate dataset
library(ncdf4)

fnm = 'MERRA/CONUS_Output/interpolated_merra2_for_SARTA_two_slab_2019_MAM_CONUS_NCA_09UTC.nc'

nc1 = nc_open(fnm)
# Print all variable names
print(nc1)

for (j in seq(1,length(nc1$var)) ) {
    print(nc1$var[[j]]$name)
    tmdt = ncvar_get(nc1,nc1$var[[j]]$name)
    # If 4D, summarize by vertical level
    print(dim(tmdt))
    print(summary(as.vector(tmdt)))
    if (length(dim(tmdt)) == 4 ) {
      nlv = dim(tmdt)[3]
      for (k in seq(1,nlv)) {
          print(summary(as.vector(tmdt[,,k,])))    
      }
    }
    
}

nc_close(nc1)






