# Process MERRA for CONUS templates, with cloud slab computations

import numpy
from netCDF4 import Dataset
import os
import pandas
import datetime
import merra_support

dst = datetime.date(2019,6,1)
dfn = datetime.date(2019,9,1)
ddf = dfn - dst
totdy = ddf.days
dysq = numpy.arange(totdy,dtype=numpy.float64)

hrs = 9
sstr = 'JJA'

print(totdy)

dtdr = '/home/jonhobbs/Documents/Research/AIRS/MERRA/Hour%02d' % (hrs)

# AIRS/SARTA levels
rnm = '/home/jonhobbs/Documents/Research/AIRS/MERRA/AIRS_Levels_Quantiles.nc'
f = Dataset(rnm,'r')
airs_sarta_levs = f.variables['level'][:]
f.close()

# Static MERRA
rnm = '/home/jonhobbs/Documents/Research/AIRS/MERRA/MERRA2_101.const_2d_ctm_Nx.00000000.nc4'
f = Dataset(rnm,'r')
lon = f.variables['lon'][:]
lat = f.variables['lon'][:]
frlnd = f.variables['FRLAND'][0,:,:]
f.close()

# CONUS Mask
rnm = '/home/jonhobbs/Documents/Research/AIRS/MERRA/merra2_NCA_mask.nc'
f = Dataset(rnm,'r')
lonmsk = f.variables['lon'][:]
latmsk = f.variables['lat'][:]
conus = f.variables['NCA_mask'][:,:]
f.close()

lnsq = numpy.arange(lonmsk.shape[0])
ltsq = numpy.arange(latmsk.shape[0])

ussmln = numpy.sum(conus,axis=0)
lnsb = lnsq[ussmln > 0]
print(lonmsk[lnsb])

ussmlt = numpy.sum(conus,axis=1)
ltsb = ltsq[ussmlt > 0]
print(latmsk[ltsb])

nlvout = airs_sarta_levs.shape[0]

print(lnsb)
print(ltsb)
print(conus.shape)

lt1 = ltsb[0]
lt2 = numpy.amax(ltsb) + 1
ln1 = lnsb[0]
ln2 = numpy.amax(lnsb) + 1
ussb = conus[lt1:lt2,ln1:ln2]

outdr = '/home/jonhobbs/Documents/Research/AIRS/MERRA/CONUS_Output'
outnm = '%s/interpolated_merra2_cloud_slab_%04d_%s_CONUS_NCA_%02dUTC.nc' % (outdr,dst.year,sstr,hrs)

for t1 in range(0,totdy):
    d2 = dst + datetime.timedelta(days=t1)

    slbfrm = merra_support.merra_find_cloud_slab(dtdr, d2, lnsb, ltsb, outnm, t1)

