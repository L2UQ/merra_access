# Given regional reference granules, find AIRS scan row overlap 

import merra_support
import h5py
from netCDF4 import Dataset
import csv
import pandas
import numpy

# Read config containing ref granule info 
df = pandas.read_csv("2018/CONUS_DJF_Templates_2018_Row.csv", \
                     dtype = {'MaskVal':int, 'Region':str, 'Abbrev':str, \
                              'Hour':int, 'Year':int, 'Season':str, \
                              'CtrLon':float, 'CtrLat':float, 'Month':int, 'Day':int, 'Granule':int, \
                              'MinRow':int, 'ModRow':int, 'MaxRow':int, 'Status':str })
ntmplt = df.shape[0]
flflt = numpy.array([-9999.],dtype=numpy.float32)
sstr = 'DJF'

qdr = '/home/jhobbs/AIST16'
mtdr = '/home/jhobbs/AIST16/CONUS_DF/%04d/Data' % (df['Year'].values[0])
csdr = '/home/jhobbs/AIST16/CONUS_DF/%04d/Outputs' % (df['Year'].values[0])
airdr = mtdr 
outdr = '/home/jhobbs/AIST16/CONUS_DF/%04d/Quantile' % (df['Year'].values[0])

for j in range(ntmplt):
    mtfl = '%s/interpolated_merra2_cloud_slab_%04d_%s_CONUS_NCA_%02dUTC.nc' % (mtdr, \
                    df['Year'].values[j], sstr, df['Hour'].values[j])
    # Deal with year boundary
    if ( (sstr ==  'DJF')  and (df['Month'].values[j] != 12) ):
        yrvl = df['Year'].values[j] + 1
    else:
        yrvl = df['Year'].values[j]

    rwsout = merra_support.airs_granule_overlap(mtfl, yrvl, df['Month'].values[j], df['Day'].values[j], \
                                                df['Granule'].values[j], 'NCA_mask', df['MaskVal'].values[j])
    df['MinRow'].values[j] = rwsout[0]
    df['ModRow'].values[j] = rwsout[1]
    df['MaxRow'].values[j] = rwsout[2]

# Re-write templates  
df.to_csv("2018/CONUS_DJF_Templates_2018_Row.csv",index=False)

