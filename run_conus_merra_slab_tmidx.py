# Process MERRA for CONUS templates, with cloud slab computations
# Use selected time index in diurnal cycle

import numpy
from netCDF4 import Dataset
import os
import pandas
import datetime
import merra_support

dst = datetime.date(2017,6,2)
dfn = datetime.date(2017,9,2)
ddf = dfn - dst
totdy = ddf.days
dysq = numpy.arange(totdy,dtype=numpy.float64)

tmchc = 0
hrs = 0
sstr = 'JJA'

print(totdy)

dtdr = '/Users/jhobbs/Documents/AIST16_UQ/Data/MERRA/Diurnal_Data'

# AIRS/SARTA levels
rnm = '/Users/jhobbs/Documents/AIST16_UQ/Data/AIRS/AIRS_Levels_Quantiles.nc'
f = Dataset(rnm,'r')
airs_sarta_levs = f.variables['level'][:]
f.close()

# Static MERRA
rnm = '/Users/jhobbs/Documents/AIST16_UQ/Data/MERRA/MERRA2_101.const_2d_ctm_Nx.00000000.nc4'
f = Dataset(rnm,'r')
lon = f.variables['lon'][:]
lat = f.variables['lon'][:]
frlnd = f.variables['FRLAND'][0,:,:]
f.close()

# CONUS Mask
rnm = '/Users/jhobbs/Documents/AIST16_UQ/Data/MERRA/merra2_NCA_mask.nc'
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
tmpout = numpy.zeros( (totdy,nlvout,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
qvout = numpy.zeros( (totdy,nlvout,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
rhout = numpy.zeros( (totdy,nlvout,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
clqdout = numpy.zeros( (totdy,nlvout,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
ciceout = numpy.zeros( (totdy,nlvout,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
cfrcout = numpy.zeros( (totdy,nlvout,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
htout = numpy.zeros( (totdy,nlvout,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
psfout = numpy.zeros( (totdy,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
t2mout = numpy.zeros( (totdy,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
tsfcout = numpy.zeros( (totdy,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
zsfmout = numpy.zeros( (totdy,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
cftotout = numpy.zeros( (totdy,ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.
frlndout = numpy.zeros( (ltsb.shape[0],lnsb.shape[0]), dtype=numpy.float32) - 9999.

print(lnsb)
print(ltsb)
print(conus.shape)

lt1 = ltsb[0]
lt2 = numpy.amax(ltsb) + 1
ln1 = lnsb[0]
ln2 = numpy.amax(lnsb) + 1
ussb = conus[lt1:lt2,ln1:ln2]
frlndout[:,:] = frlnd[lt1:lt2,ln1:ln2]

for t1 in range(0,totdy):
    d2 = dst + datetime.timedelta(days=t1)

    tmparr = merra_support.merra_conv_temp_prof(dtdr, d2, lnsb, ltsb, airs_sarta_levs, tmidx=tmchc)
    tmpout[t1,:,:,:] = tmparr[:,:,:]

    qvarr = merra_support.merra_conv_shum_prof(dtdr, d2, lnsb, ltsb, airs_sarta_levs, tmidx=tmchc, vrnm='QV')
    qvout[t1,:,:,:] = qvarr[:,:,:]

    rharr = merra_support.rh_from_qv_prof(qvarr, tmparr, airs_sarta_levs)
    rhout[t1,:,:,:] = rharr[:,:,:]

    #clqdarr = merra_support.merra_conv_shum_prof(dtdr, d2, lnsb, ltsb, airs_sarta_levs, vrnm='QL')
    #clqdout[t1,:,:,:] = clqdarr[:,:,:]

    #cicearr = merra_support.merra_conv_shum_prof(dtdr, d2, lnsb, ltsb, airs_sarta_levs, vrnm='QI')
    #ciceout[t1,:,:,:] = cicearr[:,:,:]

    #cfrcarr = merra_support.merra_conv_cfrac_prof(dtdr, d2, lnsb, ltsb, airs_sarta_levs, vrnm='CLOUD')
    #cfrcout[t1,:,:,:] = cfrcarr[:,:,:]

    #slbfrm = merra_support.merra_find_cloud_slab(dtdr, d2, lnsb, ltsb)

    psfc = merra_support.merra_subset_2d(dtdr, d2, lnsb, ltsb, 'PS', tmidx=tmchc)
    psfout[t1,:,:] = 0.01 * psfc[:,:]
    t2m = merra_support.merra_subset_2d(dtdr, d2, lnsb, ltsb, 'T2M', tmidx=tmchc)
    t2mout[t1,:,:] = t2m[:,:]
    tsfc = merra_support.merra_subset_2d(dtdr, d2, lnsb, ltsb, 'TS', tmidx=tmchc)
    tsfcout[t1,:,:] = t2m[:,:]
    zsfc = merra_support.merra_subset_2d(dtdr, d2, lnsb, ltsb, 'PHIS', tmidx=tmchc, mflstr = 'inst3_3d_asm')
    zsfmout[t1,:,:] = zsfc[:,:] / 9.8
    #cftot = merra_support.merra_subset_2d(dtdr, d2, lnsb, ltsb, 'CLDTOT', mflstr = 'tavg1_2d_rad')
    #cftotout[t1,:,:] = cftot[:,:] 

    htarr = merra_support.merra_conv_heights(dtdr, d2, lnsb, ltsb, airs_sarta_levs, tmparr, tmidx=tmchc)
    htout[t1,:,:,:] = htarr[:,:,:]

# Static fields/masks

# Output
outdr = '/Users/jhobbs/Documents/AIST16_UQ/Data/MERRA/CONUS_Output'
outnm = '%s/interpolated_merra2_cloud_slab_%04d_%s_CONUS_NCA_%02dUTC.nc' % (outdr,dst.year,sstr,hrs)

fout = Dataset(outnm,'w')

dimlv = fout.createDimension('lev',nlvout)
dimlt = fout.createDimension('lat',ltsb.shape[0])
dimln = fout.createDimension('lon',lnsb.shape[0])
dimtm = fout.createDimension('time',size=None)

varlv = fout.createVariable('lev','f8',['lev'], fill_value=-9999)
varlv[:] = airs_sarta_levs
varlv.long_name = 'vertical level'
varlv.units = 'hPa'
varlv.missing_value = -9999

varln = fout.createVariable('lon','f8',['lon'], fill_value=-9999)
varln[:] = lonmsk[lnsb]
varln.long_name = 'longitude'
varln.units = 'degrees_east'
varln.missing_value = -9999

varlt = fout.createVariable('lat','f8',['lat'], fill_value=-9999)
varlt[:] = latmsk[ltsb]
varlt.long_name = 'latitude'
varlt.units = 'degrees_north'
varlt.missing_value = -9999

vartm = fout.createVariable('time','f8',['time'], fill_value=-9999)
vartm[:] = dysq
vartm.long_name = 'time'
vartm.units = 'days since %s 00:00:00' % (dst.strftime('%Y-%m-%d'))
vartm.missing_value = -9999

varmsk = fout.createVariable('NCA_mask','i4',['lat','lon'], fill_value=-99)
varmsk[:] = ussb 
varmsk.long_name = 'National Climate Assessment region mask'
varmsk.missing_value = -99

varpsf = fout.createVariable('spres','f4',['time','lat','lon'], fill_value=-9999)
varpsf[:] = psfout
varpsf.long_name = 'surface air pressure'
varpsf.units = 'hPa'
varpsf.missing_value = -9999

varzsf = fout.createVariable('salti','f4',['time','lat','lon'], fill_value=-9999)
varzsf[:] = zsfmout
varzsf.long_name = 'surface geopotential height'
varzsf.units = 'm'
varzsf.missing_value = -9999

vartsf = fout.createVariable('stemp','f4',['time','lat','lon'], fill_value=-9999)
vartsf[:] = tsfcout
vartsf.long_name = 'surface skin temperature'
vartsf.units = 'K'
vartsf.missing_value = -9999

#varcft = fout.createVariable('cfrac','f4',['time','lat','lon'], fill_value=-9999)
#varcft[:] = cftotout
#varcft.long_name = 'total cloud area fraction'
#varcft.units = 'None'
#varcft.missing_value = -9999

vartp = fout.createVariable('ptemp','f4',['time','lev','lat','lon'], fill_value=-9999)
vartp[:] = tmpout
vartp.long_name = 'air temperature'
vartp.units = 'K'
vartp.missing_value = -9999

#varqv = fout.createVariable('qv','f8',['time','lev','lat','lon'], fill_value=-9999)
#varqv[:] = qvout
#varqv.long_name = 'specific humidity'
#varqv.units = 'kg kg-1'
#varqv.missing_value = -9999

varrh = fout.createVariable('rh','f4',['time','lev','lat','lon'], fill_value=-9999)
varrh[:] = rhout
varrh.long_name = 'relative humidity'
varrh.units = 'unitless'
varrh.missing_value = -9999

#varrh = fout.createVariable('cfrac_profile','f4',['time','lev','lat','lon'], fill_value=-9999)
#varrh[:] = cfrcout
#varrh.long_name = 'cloud fraction for radiation'
#varrh.units = 'unitless'
#varrh.missing_value = -9999

#varcl = fout.createVariable('cloud_liquid','f4',['time','lev','lat','lon'], fill_value=-9999)
#varcl[:] = clqdout
#varcl.long_name = 'mass fraction of cloud liquid water'
#varcl.units = 'kg kg-1'
#varcl.missing_value = -9999

#varci = fout.createVariable('cloud_ice','f4',['time','lev','lat','lon'], fill_value=-9999)
#varci[:] = ciceout
#varci.long_name = 'mass fraction of cloud ice'
#varci.units = 'kg kg-1'
#varci.missing_value = -9999

varht = fout.createVariable('palts','f4',['time','lev','lat','lon'], fill_value=-9999)
varht[:] = htout
varht.long_name = 'geopotential height'
varht.units = 'm'
varht.missing_value = -9999

varlf = fout.createVariable('landfrac','f4',['lat','lon'], fill_value=-9999)
varlf[:] = frlndout
varlf.long_name = 'land fraction'
varlf.units = 'unitless'
varlf.missing_value = -9999

varctp1 = fout.createVariable('ctype1','f4',['time','lat','lon'], fill_value=-9999)
varctp1.long_name = 'cloud slab 1 type'
varctp1.comment = 'Type 101.0: Liquid. Type 201.0: Ice'
varctp1.units = 'None'
varctp1.missing_value = -9999

varctp2 = fout.createVariable('ctype2','f4',['time','lat','lon'], fill_value=-9999)
varctp2.long_name = 'cloud slab 2 type'
varctp2.comment = 'Type 101.0: Liquid. Type 201.0: Ice'
varctp2.units = 'None'
varctp2.missing_value = -9999

varpsz1 = fout.createVariable('cpsize1','f4',['time','lat','lon'], fill_value=-9999)
varpsz1.long_name = 'cloud slab 1 particle size'
varpsz1.units = 'micron'
varpsz1.missing_value = -9999

varpsz2 = fout.createVariable('cpsize2','f4',['time','lat','lon'], fill_value=-9999)
varpsz2.long_name = 'cloud slab 2 particle size'
varpsz2.units = 'micron'
varpsz2.missing_value = -9999

varpbt1 = fout.createVariable('cprbot1','f4',['time','lat','lon'], fill_value=-9999)
varpbt1.long_name = 'cloud slab 1 bottom pressure'
varpbt1.units = 'hPa'
varpbt1.missing_value = -9999

varpbt2 = fout.createVariable('cprbot2','f4',['time','lat','lon'], fill_value=-9999)
varpbt2.long_name = 'cloud slab 2 bottom pressure'
varpbt2.units = 'hPa'
varpbt2.missing_value = -9999

varptp1 = fout.createVariable('cprtop1','f4',['time','lat','lon'], fill_value=-9999)
varptp1.long_name = 'cloud slab 1 top pressure'
varptp1.units = 'hPa'
varptp1.missing_value = -9999

varptp2 = fout.createVariable('cprtop2','f4',['time','lat','lon'], fill_value=-9999)
varptp2.long_name = 'cloud slab 2 top pressure'
varptp2.units = 'hPa'
varptp2.missing_value = -9999

varctt1 = fout.createVariable('cstemp1','f4',['time','lat','lon'], fill_value=-9999)
varctt1.long_name = 'cloud slab 1 top temperature'
varctt1.units = 'K'
varctt1.missing_value = -9999

varctt2 = fout.createVariable('cstemp2','f4',['time','lat','lon'], fill_value=-9999)
varctt2.long_name = 'cloud slab 2 top temperature'
varctt2.units = 'K'
varctt2.missing_value = -9999

varngw1 = fout.createVariable('cngwat1','f4',['time','lat','lon'], fill_value=-9999)
varngw1.long_name = 'cloud slab 1 non gas water'
varngw1.units = 'g m^-2'
varngw1.missing_value = -9999

varngw2 = fout.createVariable('cngwat2','f4',['time','lat','lon'], fill_value=-9999)
varngw2.long_name = 'cloud slab 2 non gas water'
varngw2.units = 'g m^-2'
varngw2.missing_value = -9999

fout.close()

# Construct and save slabs
for t1 in range(0,totdy):
    d2 = dst + datetime.timedelta(days=t1)

    slbfrm = merra_support.merra_find_cloud_slab(dtdr, d2, lnsb, ltsb, outnm, t1)


