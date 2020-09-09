# Match input MERRA dataset to AIRS granules 

from netCDF4 import Dataset
import numpy
import numpy.ma as ma
import datetime
import os
import pandas
import merra_support 

ddir = '/home/jhobbs/AIST16/CONUS_DF/2019/Data'
mkoutput = True

dst = datetime.date(2019,3,1)
dfn = datetime.date(2019,6,1)
ddf = dfn - dst
totdy = ddf.days
dysq = numpy.arange(totdy,dtype=numpy.float64)

hrs = 9
sstr = 'MAM'

merfl = '%s/interpolated_merra2_for_SARTA_two_slab_%04d_%s_CONUS_NCA_%02dUTC.nc' % (ddir,dst.year, sstr, hrs)

f = Dataset(merfl)
lon = f.variables['lon'][:]
lat = f.variables['lat'][:]
tms = f.variables['time'][:]
tmunit = f.variables['time'].units
f.close() 

print(tmunit)

outfl = '%s/CONUS_AIRS_CldFrc_Match_%04d_%s_%02dUTC.nc' % (ddir,dst.year,sstr,hrs)
if mkoutput:
    # Set up output file
    merra_support.setup_airs_cloud(outfl,tms,lat,lon,tmunit)

f = Dataset(outfl)
msgvl = f.variables['AIRS_CldFrac_1'].missing_value
f.close()

print(msgvl)

d0 = datetime.datetime(dst.year,dst.month,dst.day,hrs,0,0)
for t in range(totdy):
    dcr = d0 + datetime.timedelta(days=t)
    print(dcr)
    merra_support.airs_cfrac_match_merra(outfl,t,dcr,lat,lon)
    
