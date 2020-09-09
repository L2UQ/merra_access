# merra_access
Support for accessing MERRA-2 reference data

For the AIST effort, MERRA2 was assembled on the NCCS systems. 
Moving forward, standard GES-DISC subsetting and download will probably be needed. Subsetting within python is documented in this [GES DISC notebook](https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Use%20the%20Web%20Services%20API%20for%20Subsetting%20MERRA-2%20Data)

MERRA Products Needed (instantaneous): Accessible with Python subset interface

* Meteorology, assimilated 3d fields on pressure levels 3-hourly: M2I3NPASM, MERRA-2 inst3_3d_asm_Np  
Access with `merra_conus_gesdisc.py`
    - QV: specific humidity
    - T: temperature
    - H: geopotential height
    - QI: cloud ice
    - QL: cloud liquid water
    - PHIS: surface geopotential height
* Met Diagnostics, 2d assimilated single-level diagnostics 1-hourly: M2I1NXASM: MERRA-2 inst1_2d_asm_Nx  
Access with `merra_conus_gesdisc_sfc.py`
    - TS: skin temperature
    - PS: surface pressure
    - T2M: 2-meter temperature

MERRA time-average radiation diagnostics (times are offset 90 minutes from instantaneous fields): 

* Radiation, time-averaged 3d radiation diagnostics 3-hourly: M2T3NPRAD, MERRA-2 tavg3_3d_rad_Np  
Access with `merra_conus_gesdisc_rad.py`
    - CLDTOT: total cloud area fraction
    - EMIS: surface emissivity
* Cloud, time-averaged 3d cloud diagnostics 3-hourly: M2T3NPCLD, MERRA-2 tavg3_3d_cld_Np  
Access with `merra_conus_gesdisc_rad3d.py`
    - CLOUD: cloud fraction for radiation
    
MERRA constant fields, single download    

* Constant fields, available in `MERRA2_101.const_2d_ctm_Nx.00000000.nc4`
    - FRLAND: land fraction
    
Datasets can be optionally extracted with [simple subset wizard](https://disc.gsfc.nasa.gov/SSW)

