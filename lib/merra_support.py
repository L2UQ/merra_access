# Subroutines for MERRA processing to SARTA inputs

import numpy
from netCDF4 import Dataset
import datetime
import os
from scipy import stats
from numpy import ndarray, ma
import pandas
import math

def sfclvl(psfc, levarr):
    ### Return array with surface level indicator: the lowest vertical level above the surface pressure
    # Assume psfc is 2D (lat, lon)

    nlt = psfc.shape[0]
    nln = psfc.shape[1]

    nz = levarr.shape[0]
    psq = numpy.arange(nz)
    slvs = numpy.zeros((nlt,nln),dtype=numpy.int16)

    for j in range(nlt):
        for i in range(nln):
            psb = psq[levarr <= psfc[j,i]]
            slvs[j,i] = psb[-1]
    return slvs

def sfclvl_rev_met(psfc, levarr, metarr):
    ### Return array with surface level indicator: the lowest vertical level above the surface pressure
    ### Here the levarr is descending in pressure
    ### Also check where metarr is not masked
    # Assume psfc is 2D (lat, lon)
    # levarr:  1D array of pressure levels
    # metarr:  3D met array (masked)

    nlt = psfc.shape[0]
    nln = psfc.shape[1]

    nz = levarr.shape[0]
    psq = numpy.arange(nz)
    slvs = numpy.zeros((nlt,nln),dtype=numpy.int16)
    for j in range(nlt):
        for i in range(nln):
            psb = psq[(levarr <= psfc[j,i]) & numpy.logical_not(metarr[:,j,i].mask)]
            slvs[j,i] = psb[0]
    return slvs

def merra_subset_2d( srchdr, srchdt, lnsq, ltsq, varnm, mflstr = 'inst1_2d_asm'):
    ### Extract a MERRA 2D variable from spatial subset
    flst = os.listdir(srchdr)
    dtstr = srchdt.strftime('%Y%m%d')

    f2d = -1
    f3d = -1
    j = 0

    while ( (j < len(flst)) and ((f2d < 0) or (f3d < 0) ) ):
        if ( (mflstr in flst[j]) and (dtstr in flst[j]) ):
            f2d = j
            mer2d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer2d,'r')
            vararr = f.variables[varnm][0,ltsq,lnsq]
            f.close()
        j = j + 1

    return vararr


def merra_conv_temp_prof( srchdr, srchdt, lnsq, ltsq, lvout, msgval=-9999.):
    ### Convert MERRA temperature profile to SARTA/AIRS pressure levels
    ### srchdr:     Directory with MERRA daily files
    ### srchdt:     Desired date (a python datetime object)
    ### lnsq:       Longitude subset sequence
    ### ltsq:       Latitude subset sequence
    ### lvout:      Array of pressure levels

    flst = os.listdir(srchdr)
    dtstr = srchdt.strftime('%Y%m%d')

    f2d = -1
    f3d = -1
    j = 0

    nlvout = lvout.shape[0]
    lwwt = numpy.zeros(nlvout,)
    hiwt = numpy.zeros(nlvout,)
    lwidx = numpy.zeros(nlvout,dtype=numpy.int32)
    hiidx = numpy.zeros(nlvout,dtype=numpy.int32)

    print(dtstr)
    while ( (j < len(flst)) and ((f2d < 0) or (f3d < 0) ) ):
        if ( ('inst1_2d_asm' in flst[j]) and (dtstr in flst[j]) ):
            f2d = j
            mer2d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer2d,'r')
            t2m = f.variables['T2M'][0,ltsq,lnsq]
            psfc = f.variables['PS'][0,ltsq,lnsq]
            #tskn = f['TS'][0,ltsq,lnsq]
            f.close()
        if ( ('inst3_3d_asm' in flst[j]) and (dtstr in flst[j]) ):
            f3d = j
            mer3d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer3d,'r')
            tmpmer = f.variables['T'][0,:,ltsq,lnsq]
            lvmer = f.variables['lev'][:]
            f.close()
        j = j + 1

    # Develop weights
    nlvmr = lvmer.shape[0]
    for k in range(nlvout):
        lmr = nlvmr - 1
        kfd = -1
        # No need to loop in upper atm
        if ( lvout[k] < lvmer[lmr]):
            kfd = lmr
            lwidx[k] = lmr
            hiidx[k] = lmr + 1
            lwwt[k] = 1.0
            hiwt[k] = 0.0
        # Loop otherwise
        while ( (lmr > 0) and (kfd < 0) ):
            if ((lvout[k] >= lvmer[lmr]) and (lvout[k] < lvmer[lmr-1]) ):
                kfd = lmr
                lwidx[k] = lmr - 1
                hiidx[k] = lmr
                lwwt[k] = (numpy.log(lvout[k]) - numpy.log(lvmer[lmr])) / \
                          (numpy.log(lvmer[lmr-1]) - numpy.log(lvmer[lmr]))
                hiwt[k] = 1.0 - lwwt[k]
            lmr = lmr - 1

    srtsfc = sfclvl(psfc*0.01, lvout)
    mersfc = sfclvl_rev_met(psfc*0.01, lvmer, tmpmer)
    tmpmer[tmpmer > 1.0e10] = msgval

    # Set up output temp profile
    tprfout = numpy.zeros((nlvout,ltsq.shape[0],lnsq.shape[0]),numpy.float32) + msgval

    # Subset levels at upper atmos, linear profile in upper mesosphere
    hsq = numpy.arange(nlvout)
    hsb = hsq[hiidx == nlvmr]
    nupr = hsb.shape[0]

    # Model fit
    nxy = (ltsq.shape[0] * lnsq.shape[0])
    lva = nlvmr
    lvb = nlvmr-3
    tmn = numpy.mean(tmpmer[lva:lvb,:,:],axis=0)
    tmnovr = numpy.mean(tmpmer[lva:lvb,:,:])

    for p in range(nlvmr-1,nlvmr-4,-1):
        prsrp = numpy.tile(numpy.log(lvmer[p]),nxy)
        tflt = tmpmer[p,:,:].flatten()
        if p == (nlvmr-1):
            prsreg =  numpy.zeros((nxy,),prsrp.dtype)
            prsreg[:] = prsrp
            tpreg = numpy.zeros((nxy,),tflt.dtype)
            tpreg[:] = tflt
        else:
            prsreg = numpy.append(prsreg,prsrp.flatten())
            tpreg = numpy.append(tpreg,tflt)
    slp, itcpt, r2, pvl, stderr = stats.linregress(prsreg,tpreg)
    sstr = 'T = %.3f + %.4e logp, R2 = %.4f' % (itcpt, slp, r2)
    print(sstr)

    #tprfout[plvsb,q,p] = itcpt + slp * lvout[plvsb]
    tbsln = itcpt + slp * numpy.log(lvout[hsb])
    # Getting close to thermosphere
    tbsln[0:2] = itcpt + slp * numpy.log(lvout[3]) + slp / 2.0 * (numpy.log(lvout[0:2]) - numpy.log(lvout[3]))

    for k in range(nupr):
        tprfout[k,:,:] = tbsln[k] + (tmn - tmnovr)

    for k in range(hsb[nupr-1]+1,nlvout):
        lcmsk = numpy.zeros((ltsq.shape[0],lnsq.shape[0]),dtype=numpy.float32)
        # Upper atmosphere
        #if (hiidx[k] == nlvmr):
        #    tprfout[k,:,:] = tmpmer[nlvmr-1,:,:]
        # Masking tricks for the rest
        if (hiidx[k] > 0):
            airspsfc = psfc * 0.01
            lcmsk[ (srtsfc >= (k)) & (airspsfc > lvmer[lwidx[k]]) & (mersfc <= lwidx[k]) ] = 1.0
            tprfout[k,:,:] =  lcmsk * (lwwt[k] * tmpmer[lwidx[k],:,:] + hiwt[k] * tmpmer[hiidx[k],:,:]) + msgval * (1.0 - lcmsk)

    # Loop through all locations for sfc behavior
    pref = 1000.0
    kappa = 0.286
    pvlds = numpy.arange(nlvout)
    for q in range(ltsq.shape[0]):
        for p in range(lnsq.shape[0]):
            # Identify levels needed for extrapolation
            airspsfc = psfc[q,p] * 0.01
            plvsb = pvlds[ (pvlds <= (srtsfc[q,p]+1)) &  ( (pvlds > srtsfc[q,p]) | (hiidx <= mersfc[q,p] ) ) ]

            # Average potential temp
            prs2 = numpy.array([ lvmer[mersfc[q,p]+2], lvmer[mersfc[q,p]+1], lvmer[mersfc[q,p]], airspsfc ])
            tmp2 = numpy.array([ tmpmer[mersfc[q,p]+2,q,p], tmpmer[mersfc[q,p]+1,q,p], tmpmer[mersfc[q,p],q,p], t2m[q,p] ])
            prt2 = pref / prs2
            thet2 = tmp2 * numpy.power(prt2,kappa)

            thmn = numpy.mean(thet2)

            #slp, itcpt, r2, pvl, stderr = stats.linregress(prs2,tmp2)
            #tprfout[plvsb,q,p] = itcpt + slp * lvout[plvsb]

            prtinv = lvout[plvsb] / pref
            tprfout[plvsb,q,p] = thmn * numpy.power(prtinv,kappa)

            #if ( (p == 32) and (q == 42) ):
            #    mrlv = 'Plvs Needed to interpolate 42, 32: %d' % (plvsb.shape)
            #    print(mrlv)
            #    print(plvsb)
            #    print(lvout[plvsb])
            #    print(tprfout[plvsb[0]-1,q,p])
            #    print(prs2)
            #    print(thet2)

    #mrlv = 'MERRA Sfc Lvl 42, 32: %d, %.4f, %.4f, %.2f' % (mersfc[42,32],lvmer[mersfc[42,32]],psfc[42,32],tmpmer[mersfc[42,32],42,32])
    #print(mrlv)
    #print(hiidx)

    return tprfout

def merra_conv_shum_prof( srchdr, srchdt, lnsq, ltsq, lvout, vrnm = 'QV', msgval=-9999.):
    ### Convert MERRA specific humidity profile to SARTA/AIRS pressure levels
    ### srchdr:     Directory with MERRA daily files
    ### srchdt:     Desired date (a python datetime object)
    ### lnsq:       Longitude subset sequence
    ### ltsq:       Latitude subset sequence
    ### lvout:      Array of pressure levels
    ### vrnm:       Variable name for 3D variable to evaluate

    flst = os.listdir(srchdr)
    dtstr = srchdt.strftime('%Y%m%d')

    f2d = -1
    f3d = -1
    j = 0

    nlvout = lvout.shape[0]
    lwwt = numpy.zeros(nlvout,)
    hiwt = numpy.zeros(nlvout,)
    lwidx = numpy.zeros(nlvout,dtype=numpy.int32)
    hiidx = numpy.zeros(nlvout,dtype=numpy.int32)

    print(dtstr)
    while ( (j < len(flst)) and ((f2d < 0) or (f3d < 0) ) ):
        if ( ('inst1_2d_asm' in flst[j]) and (dtstr in flst[j]) ):
            f2d = j
            mer2d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer2d,'r')
            psfc = f.variables['PS'][0,ltsq,lnsq]
            f.close()
        if ( ('inst3_3d_asm' in flst[j]) and (dtstr in flst[j]) ):
            f3d = j
            mer3d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer3d,'r')
            qvmer = f.variables[vrnm][0,:,ltsq,lnsq]
            lvmer = f.variables['lev'][:]
            f.close()
        j = j + 1

    # Develop weights
    nlvmr = lvmer.shape[0]
    for k in range(nlvout):
        lmr = nlvmr - 1
        kfd = -1
        # No need to loop in upper atm
        if ( lvout[k] < lvmer[lmr]):
            kfd = lmr
            lwidx[k] = lmr
            hiidx[k] = lmr + 1
            lwwt[k] = 1.0
            hiwt[k] = 0.0
        # Loop otherwise
        while ( (lmr > 0) and (kfd < 0) ):
            if ((lvout[k] >= lvmer[lmr]) and (lvout[k] < lvmer[lmr-1]) ):
                kfd = lmr
                lwidx[k] = lmr - 1
                hiidx[k] = lmr
                lwwt[k] = (numpy.log(lvout[k]) - numpy.log(lvmer[lmr])) / \
                          (numpy.log(lvmer[lmr-1]) - numpy.log(lvmer[lmr]))
                hiwt[k] = 1.0 - lwwt[k]
            lmr = lmr - 1

    srtsfc = sfclvl(psfc*0.01, lvout)
    mersfc = sfclvl_rev_met(psfc*0.01, lvmer, qvmer)
    qvmer[qvmer > 1.0e10] = msgval

    # Set up output qv profile
    qvout = numpy.zeros((nlvout,ltsq.shape[0],lnsq.shape[0]),numpy.float32) + msgval

    for k in range(nlvout):
        lcmsk = numpy.zeros((ltsq.shape[0],lnsq.shape[0]),dtype=numpy.float32)
        # Upper atmosphere
        if (hiidx[k] == nlvmr):
            qvout[k,:,:] = qvmer[nlvmr-1,:,:]
        # Masking tricks for the rest
        elif (hiidx[k] > 0):
            airspsfc = psfc * 0.01
            lcmsk[ (srtsfc >= (k)) & (airspsfc > lvmer[lwidx[k]]) & (mersfc <= lwidx[k]) ] = 1.0
            qvout[k,:,:] =  lcmsk * (lwwt[k] * qvmer[lwidx[k],:,:] + hiwt[k] * qvmer[hiidx[k],:,:]) + msgval * (1.0 - lcmsk)

    # Loop through all locations for sfc behavior
    pvlds = numpy.arange(nlvout)
    for q in range(ltsq.shape[0]):
        for p in range(lnsq.shape[0]):
            # Identify levels needed for extrapolation
            airspsfc = psfc[q,p] * 0.01
            plvsb = pvlds[ (pvlds <= (srtsfc[q,p]+1)) &  ( (pvlds > srtsfc[q,p]) | (hiidx <= mersfc[q,p] ) ) ]

            # Average QV
            prs2 = numpy.array([ lvmer[mersfc[q,p]+2], lvmer[mersfc[q,p]+1], lvmer[mersfc[q,p]] ])
            tmp2 = numpy.array([ qvmer[mersfc[q,p]+2,q,p], qvmer[mersfc[q,p]+1,q,p], qvmer[mersfc[q,p],q,p] ])
            thmn = numpy.mean(tmp2)
            #prs2 = numpy.array([ lvmer[mersfc[q,p]], lvmer[mersfc[q,p]+1] ])
            #tmp2 = numpy.array([ qvmer[mersfc[q,p],q,p], qvmer[mersfc[q,p]+1,q,p] ])
            #slp, itcpt, r2, pvl, stderr = stats.linregress(prs2,tmp2)
            #qvrslt = itcpt + slp * lvout[plvsb]
            qvrslt = numpy.tile(thmn,plvsb.shape[0])
            qvrslt[qvrslt < 0.0] = 0.0
            qvout[plvsb,q,p] = qvrslt

                        #if qvout[k,q,p] < 0.0:
                        #    qvout[k,q,p] = 0.0
                        #if qvout[k,q,p] < 0.0:
                        #    qvout[k,q,p] = 0.0

    return qvout

def rh_from_qv_prof( qvprf, tprf, plvs, msgval=-9999.):
    ### Compute relative humidity profile from specific humidity and temperature
    ### qvprf:      Specific humidity profile
    ### plvs:       Vector of pressure levels (hPa)
    ### msgval:     Missing value


    nlvout = plvs.shape[0]
    ny = tprf.shape[1]
    nx = tprf.shape[2]

    # Set up output RH profile
    rhout = numpy.zeros((nlvout,ny,nx),numpy.float32) + msgval

    for k in range(nlvout):
        # Set up a locmask
        lcmsk = numpy.zeros((ny,nx),dtype=numpy.float32)
        qvtmp = qvprf[k,:,:]
        lcmsk[qvtmp != msgval] = 1.0

        # Calculate sat vap pres in hPa, AMS Glossary
        es = lcmsk * (6.112 * numpy.exp(17.67 * (tprf[k,:,:]-273.15) / (tprf[k,:,:] - 29.65)) )  + msgval * (1.0 - lcmsk)
        mmrs = lcmsk * (0.622 * es / plvs[k]) + msgval * (1.0 - lcmsk)
        rhout[k,:,:] = lcmsk * (qvprf[k,:,:] / mmrs) + msgval * (1.0 - lcmsk)

    return rhout

def merra_conv_heights( srchdr, srchdt, lnsq, ltsq, lvout, tprf, vrnm = 'H', sfcht = 'PHIS', msgval=-9999.):
    ### Convert MERRA geopotential heights to SARTA grid
    ### srchdr:     Directory with MERRA daily files
    ### srchdt:     Desired date (a python datetime object)
    ### lnsq:       Longitude subset sequence
    ### ltsq:       Latitude subset sequence
    ### lvout:      Array of pressure levels
    ### tprf:       Temperature profile on output grid (for thickness calcs)
    ### vrnm:       Variable name for 3D variable to evaluate
    ### sfcht:      Name of surface geopotential variable

    flst = os.listdir(srchdr)
    dtstr = srchdt.strftime('%Y%m%d')

    f2d = -1
    f3d = -1
    j = 0

    nlvout = lvout.shape[0]
    lwwt = numpy.zeros(nlvout,)
    hiwt = numpy.zeros(nlvout,)
    lwidx = numpy.zeros(nlvout,dtype=numpy.int32)
    hiidx = numpy.zeros(nlvout,dtype=numpy.int32)

    print(dtstr)
    while ( (j < len(flst)) and ((f2d < 0) or (f3d < 0) ) ):
        if ( ('inst1_2d_asm' in flst[j]) and (dtstr in flst[j]) ):
            f2d = j
            mer2d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer2d,'r')
            psfc = f.variables['PS'][0,ltsq,lnsq]
            f.close()
        if ( ('inst3_3d_asm' in flst[j]) and (dtstr in flst[j]) ):
            f3d = j
            mer3d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer3d,'r')
            htmer = f.variables[vrnm][0,:,ltsq,lnsq]
            lvmer = f.variables['lev'][:]
            phisfc = f[sfcht][0,ltsq,lnsq] / 9.8
            f.close()
        j = j + 1

    # Develop weights
    nlvmr = lvmer.shape[0]
    for k in range(nlvout):
        lmr = nlvmr - 1
        kfd = -1
        # No need to loop in upper atm
        if ( lvout[k] < lvmer[lmr]):
            kfd = lmr
            lwidx[k] = lmr
            hiidx[k] = lmr + 1
            lwwt[k] = 1.0
            hiwt[k] = 0.0
        # Loop otherwise
        while ( (lmr > 0) and (kfd < 0) ):
            if ((lvout[k] >= lvmer[lmr]) and (lvout[k] < lvmer[lmr-1]) ):
                kfd = lmr
                lwidx[k] = lmr - 1
                hiidx[k] = lmr
                lwwt[k] = (numpy.log(lvout[k]) - numpy.log(lvmer[lmr])) / \
                          (numpy.log(lvmer[lmr-1]) - numpy.log(lvmer[lmr]))
                hiwt[k] = 1.0 - lwwt[k]
            lmr = lmr - 1

    srtsfc = sfclvl(psfc*0.01, lvout)
    mersfc = sfclvl_rev_met(psfc*0.01, lvmer, htmer)
    htmer[htmer > 1.0e10] = msgval

    # Set up output height profile
    htout = numpy.zeros((nlvout,ltsq.shape[0],lnsq.shape[0]),numpy.float32) + msgval

    # Subset levels at upper atmos
    Rd = 287.0
    hsq = numpy.arange(nlvout)
    hsb = hsq[hiidx == nlvmr]
    print(hsb)
    nupr = hsb.shape[0]
    for k in range(hsb[nupr-1],-1,-1):
        tmpmn = (tprf[hsb[k],:,:] + tprf[hsb[k]+1,:,:]) / 2.0
        if (k == (nupr-1)):
            # Work from MERRA
            pupr = lvout[hsb[k]]
            plwr = lvmer[nlvmr-1]
            thk = Rd * tmpmn / 9.8 * (numpy.log(plwr/pupr))
            htout[hsb[k],:,:] = htmer[nlvmr-1,:,:] + thk
        else:
            # Work from prvs
            pupr = lvout[hsb[k]]
            plwr = lvmer[nlvmr-1]
            thk = Rd * tmpmn / 9.8 * (numpy.log(plwr/pupr))
            htout[hsb[k],:,:] = htout[hsb[k]+1,:,:] + thk


    for k in range(hsb[nupr-1]+1,nlvout):
        lcmsk = numpy.zeros((ltsq.shape[0],lnsq.shape[0]),dtype=numpy.float32)
        # Masking tricks for the rest
        if (hiidx[k] > 0):
            airspsfc = psfc * 0.01
            lcmsk[ (srtsfc >= (k)) & (airspsfc > lvmer[lwidx[k]]) & (mersfc <= lwidx[k]) ] = 1.0
            htout[k,:,:] =  lcmsk * (lwwt[k] * htmer[lwidx[k],:,:] + hiwt[k] * htmer[hiidx[k],:,:]) + msgval * (1.0 - lcmsk)

    # Loop through all locations for sfc behavior
    pvlds = numpy.arange(nlvout)
    for q in range(ltsq.shape[0]):
        for p in range(lnsq.shape[0]):
            # Identify levels needed for extrapolation
            airspsfc = psfc[q,p] * 0.01
            plvsb = pvlds[ (pvlds <= (srtsfc[q,p])) &  ( (pvlds > srtsfc[q,p]) | (hiidx <= mersfc[q,p] ) ) ]

            # Average potential temp
            prs2 = numpy.array([ lvmer[mersfc[q,p]+1], lvmer[mersfc[q,p]], airspsfc ])
            tmp2 = numpy.array([ htmer[mersfc[q,p]+1,q,p], htmer[mersfc[q,p],q,p], phisfc[q,p] ])
            slp, itcpt, r2, pvl, stderr = stats.linregress(prs2,tmp2)
            httmp = itcpt + slp * lvout[plvsb]
            httmp[httmp < 0] = 0.0
            htout[plvsb,q,p] = httmp


    return htout

def merra_conv_cfrac_prof( srchdr, srchdt, lnsq, ltsq, lvout, vrnm = 'CLOUD', msgval=-9999.):
    ### Convert MERRA cloud fraction profile to SARTA/AIRS pressure levels
    ### srchdr:     Directory with MERRA daily files
    ### srchdt:     Desired date (a python datetime object)
    ### lnsq:       Longitude subset sequence
    ### ltsq:       Latitude subset sequence
    ### lvout:      Array of pressure levels
    ### vrnm:       Variable name for 3D variable to evaluate

    flst = os.listdir(srchdr)
    dtstr = srchdt.strftime('%Y%m%d')

    f2d = -1
    f3d = -1
    j = 0

    nlvout = lvout.shape[0]
    lwwt = numpy.zeros(nlvout,)
    hiwt = numpy.zeros(nlvout,)
    lwidx = numpy.zeros(nlvout,dtype=numpy.int32)
    hiidx = numpy.zeros(nlvout,dtype=numpy.int32)

    print(dtstr)
    while ( (j < len(flst)) and ((f2d < 0) or (f3d < 0) ) ):
        if ( ('inst1_2d_asm' in flst[j]) and (dtstr in flst[j]) ):
            f2d = j
            mer2d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer2d,'r')
            psfc = f.variables['PS'][0,ltsq,lnsq]
            f.close()
        if ( ('tavg3_3d_rad' in flst[j]) and (dtstr in flst[j]) ):
            f3d = j
            mer3d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer3d,'r')
            cfmer = f.variables[vrnm][0,:,ltsq,lnsq]
            lvmer = f.variables['lev'][:]
            f.close()
        j = j + 1

    # Develop weights
    nlvmr = lvmer.shape[0]
    for k in range(nlvout):
        lmr = nlvmr - 1
        kfd = -1
        # No need to loop in upper atm
        if ( lvout[k] < lvmer[lmr]):
            kfd = lmr
            lwidx[k] = lmr
            hiidx[k] = lmr + 1
            lwwt[k] = 1.0
            hiwt[k] = 0.0
        # Loop otherwise
        while ( (lmr > 0) and (kfd < 0) ):
            if ((lvout[k] >= lvmer[lmr]) and (lvout[k] < lvmer[lmr-1]) ):
                kfd = lmr
                lwidx[k] = lmr - 1
                hiidx[k] = lmr
                lwwt[k] = (numpy.log(lvout[k]) - numpy.log(lvmer[lmr])) / \
                          (numpy.log(lvmer[lmr-1]) - numpy.log(lvmer[lmr]))
                hiwt[k] = 1.0 - lwwt[k]
            lmr = lmr - 1

    srtsfc = sfclvl(psfc*0.01, lvout)
    mersfc = sfclvl_rev_met(psfc*0.01, lvmer, cfmer)
    cfmer[cfmer > 1.0e10] = msgval

    # Set up output qv profile
    cfout = numpy.zeros((nlvout,ltsq.shape[0],lnsq.shape[0]),numpy.float32) + msgval

    for k in range(nlvout):
        lcmsk = numpy.zeros((ltsq.shape[0],lnsq.shape[0]),dtype=numpy.float32)
        # Upper atmosphere
        if (hiidx[k] == nlvmr):
            cfout[k,:,:] = cfmer[nlvmr-1,:,:]
        # Masking tricks for the rest
        elif (hiidx[k] > 0):
            airspsfc = psfc * 0.01
            lcmsk[ (srtsfc >= (k)) & (airspsfc > lvmer[lwidx[k]]) & (mersfc <= lwidx[k]) ] = 1.0
            cfout[k,:,:] =  lcmsk * (lwwt[k] * cfmer[lwidx[k],:,:] + hiwt[k] * cfmer[hiidx[k],:,:]) + msgval * (1.0 - lcmsk)

    # Loop through all locations for sfc behavior
    pvlds = numpy.arange(nlvout)
    for q in range(ltsq.shape[0]):
        for p in range(lnsq.shape[0]):
            # Identify levels needed for extrapolation
            airspsfc = psfc[q,p] * 0.01
            plvsb = pvlds[ (pvlds <= (srtsfc[q,p]+1)) &  ( (pvlds > srtsfc[q,p]) | (hiidx <= mersfc[q,p] ) ) ]

            # Average cfrac
            prs2 = numpy.array([ lvmer[mersfc[q,p]+2], lvmer[mersfc[q,p]+1], lvmer[mersfc[q,p]] ])
            tmp2 = numpy.array([ cfmer[mersfc[q,p]+2,q,p], cfmer[mersfc[q,p]+1,q,p], cfmer[mersfc[q,p],q,p] ])
            thmn = numpy.mean(tmp2)
            cfrslt = numpy.tile(thmn,plvsb.shape[0])
            cfrslt[cfrslt < 0.0] = 0.0
            cfrslt[cfrslt > 1.0] = 1.0
            cfout[plvsb,q,p] = cfrslt

    return cfout

def setup_airs_cloud(flnm, tms, lats, lons, tmunit = 'Seconds since 1993-01-01 00:00:00'):
    # Set up matched AIRS/MERRA cloud file
    # flnm:    Name of output file
    # tms:     Time variable array
    # lats:    Latitude variable array
    # lons:    Longitude variable array

    ntm = tms.shape[0]
    nlat = lats.shape[0]
    nlon = lons.shape[0]

    # Create Output file
    qout = Dataset(flnm,'w')

    dimln = qout.createDimension('lon',nlon)
    dimlt = qout.createDimension('lat',nlat)
    dimtm = qout.createDimension('time',ntm)
    dimtrk = qout.createDimension('AIRSFOV',9)

    if (lons.dtype == numpy.float32):
        lntp = 'f4'
    else:
        lntp = 'f8'
    varlon = qout.createVariable('lon',lntp,['lon'], fill_value = -9999)
    varlon[:] = lons
    varlon.long_name = 'longitude'
    varlon.units='degrees_east'
    varlon.missing_value = -9999

    if (lats.dtype == numpy.float32):
        lttp = 'f4'
    else:
        lttp = 'f8'
    varlat = qout.createVariable('lat',lttp,['lat'], fill_value = -9999)
    varlat[:] = lats
    varlat.long_name = 'latitude'
    varlat.units='degrees_north'
    varlat.missing_value = -9999

    if (tms.dtype == numpy.float32):
        tmtp = 'f4'
    else:
        tmtp = 'f8'
    vartm = qout.createVariable('time',lttp,['time'], fill_value = -9999)
    vartm[:] = tms
    vartm.long_name = 'time'
    vartm.units = tmunit
    vartm.missing_value = -9999

    # Other output variables
    varcfrc1 = qout.createVariable('AIRS_CldFrac_1','f4',['time','lat','lon','AIRSFOV'], fill_value = -9999)
    varcfrc1.long_name = 'AIRS cloud fraction, upper level'
    varcfrc1.units = 'unitless'
    varcfrc1.missing_value = -9999

    varcfrc2 = qout.createVariable('AIRS_CldFrac_2','f4',['time','lat','lon','AIRSFOV'], fill_value = -9999)
    varcfrc2.long_name = 'AIRS cloud fraction, lower level'
    varcfrc2.units = 'unitless'
    varcfrc2.missing_value = -9999

    varcqc1 = qout.createVariable('AIRS_CldFrac_QC_1','i2',['time','lat','lon','AIRSFOV'], fill_value = -99)
    varcqc1.long_name = 'AIRS cloud fraction quality flag, upper level'
    varcqc1.units = 'unitless'
    varcqc1.missing_value = -99

    varcqc2 = qout.createVariable('AIRS_CldFrac_QC_2','i2',['time','lat','lon','AIRSFOV'], fill_value = -99)
    varcqc2.long_name = 'AIRS cloud fraction quality flag, lower level'
    varcqc2.units = 'unitless'
    varcqc2.missing_value = -99

    varncld = qout.createVariable('AIRS_nCld','i2',['time','lat','lon','AIRSFOV'], fill_value = -99)
    varncld.long_name = 'AIRS number of cloud layers'
    varncld.units = 'unitless'
    varncld.missing_value = -99

    qout.close()

    return

def airs_cfrac_match_merra(flnm, tmidx, tmday, lats, lons,  msgvl = -9999, \
                           l2srch = '/archive/AIRSOps/airs/gdaac/v6'):
    # Set up matched AIRS/MERRA cloud file
    # flnm:    Name of output file
    # tms:     Time index in output
    # tmday:   Datetime object with time information
    # lats:    Longitude variable array
    # lons:    Longitude variable array

    # Search AIRS Level 2
    airsdr = '%s/%04d/%02d/%02d/airs2ret' % (l2srch,tmday.year,tmday.month,tmday.day)

    dsclst = []
    asclst = []

    nlat = lats.shape[0]
    nlon = lons.shape[0]

    lonmn = lons[0] - 5.0
    lonmx = lons[nlon-1] + 5.0
    latmn = lats[0] - 5.0
    latmx = lats[nlat-1] + 5.0
    d0 = datetime.datetime(1993,1,1,0,0,0)
    ddif = tmday - d0
    bsdif = ddif.total_seconds()

    # Set up reference frame
    ltrp = numpy.repeat(lats,nlon)
    ltidx = numpy.repeat(numpy.arange(nlat),nlon)
    lnrp = numpy.tile(lons,nlat)
    lnidx = numpy.tile(numpy.arange(nlon),nlat)
    merfrm = pandas.DataFrame({'GridLonIdx': lnidx, 'GridLatIdx': ltidx, \
                               'GridLon': lnrp, 'GridLat': ltrp})

    if (os.path.exists(airsdr)):
        fllst = os.listdir(airsdr)
        #print(fllst)

        for j in range(len(fllst)):
            lncr = len(fllst[j])
            l4 = lncr - 4
            if (fllst[j][l4:lncr] == '.hdf'):
                l2fl = '%s/%s' % (airsdr,fllst[j])
                ncl2 = Dataset(l2fl)
                slrzn = ncl2.variables['solzen'][:,:]
                l2lat = ncl2.variables['Latitude'][:,:]
                l2lon = ncl2.variables['Longitude'][:,:]
                l2tm = ncl2.variables['Time'][:,:]
                ncl2.close()

                # Check lat/lon ranges and asc/dsc
                l2tmdf = numpy.absolute(l2tm - bsdif)
                l2mntm = numpy.min(l2tmdf)

                # Within 4 hours
                if l2mntm < 14400.0:
                   ltflt = l2lat.flatten()
                   lnflt = l2lon.flatten()
                   latsb = ltflt[(ltflt >= latmn) & (ltflt <= latmx)]
                   lonsb = lnflt[(lnflt >= lonmn) & (lnflt <= lonmx)]
                   if ( (latsb.shape[0]  > 0) and (lonsb.shape[0] > 0) ):
                       asclst.append(fllst[j])
                       sstr = '%s %.2f' % (fllst[j], l2mntm)
                       print(sstr)


    # Set up outputs
    cld1arr = numpy.zeros((nlat,nlon,9),dtype=numpy.float32) + msgvl
    cld2arr = numpy.zeros((nlat,nlon,9),dtype=numpy.float32) + msgvl
    cld1qc = numpy.zeros((nlat,nlon,9),dtype=numpy.int16) - 99
    cld2qc = numpy.zeros((nlat,nlon,9),dtype=numpy.int16) - 99
    ncldarr = numpy.zeros((nlat,nlon,9),dtype=numpy.int16) - 99

    #print(asclst)
    tmch = 0
    if (len(asclst) > 0):
        # Start matchups
        for j in range(len(asclst)):
            l2fl = '%s/%s' % (airsdr,asclst[j])
            ncl2 = Dataset(l2fl)
            l2lat = ncl2.variables['Latitude'][:,:]
            l2lon = ncl2.variables['Longitude'][:,:]
            cfrcair = ncl2.variables['CldFrcStd'][:,:,:,:,:]
            cfrcaqc = ncl2.variables['CldFrcStd_QC'][:,:,:,:,:]
            ncldair = ncl2.variables['nCld'][:,:,:,:]
            ncl2.close()

            nairtrk = l2lat.shape[0]
            nairxtk = l2lat.shape[1]

            # Data Frame
            tkidx = numpy.repeat(numpy.arange(nairtrk),nairxtk)
            xtidx = numpy.tile(numpy.arange(nairxtk),nairtrk)
            l2lnflt = l2lon.flatten().astype(numpy.float64)
            l2ltflt = l2lat.flatten().astype(numpy.float64)
            l2frm = pandas.DataFrame({'L2LonIdx': xtidx, 'L2LatIdx': tkidx, \
                                      'L2Lon': l2lnflt, 'L2Lat': l2ltflt})
            l2frm['GridLon'] = numpy.around(l2frm['L2Lon']/0.625) * 0.625
            l2frm['GridLat'] = numpy.around(l2frm['L2Lat']/0.5) * 0.5

            l2mrg = pandas.merge(l2frm,merfrm,on=['GridLon','GridLat'])
            print(l2mrg.shape)
            tmch = tmch + l2mrg.shape[0]

            #if j  == 0:
            #    print(asclst[j])
            #    print(l2mrg[0:15])
            #mrggrp = l2mrg.groupby(['GridLatIdx','GridLonIdx']).count()
            #gtxt = 'Group: %d' % mrggrp.shape[0]

            # Output data if available
            for k in range(l2mrg.shape[0]):
                yidxout = l2mrg['GridLatIdx'].values[k]
                xidxout = l2mrg['GridLonIdx'].values[k]
                yidxl2 = l2mrg['L2LatIdx'].values[k]
                xidxl2 = l2mrg['L2LonIdx'].values[k]
                cld1arr[yidxout,xidxout,:] = cfrcair[yidxl2,xidxl2,:,:,0].flatten().astype(numpy.float32)
                cld2arr[yidxout,xidxout,:] = cfrcair[yidxl2,xidxl2,:,:,1].flatten().astype(numpy.float32)
                cld1qc[yidxout,xidxout,:] = cfrcaqc[yidxl2,xidxl2,:,:,0].flatten().astype(numpy.int16)
                cld2qc[yidxout,xidxout,:] = cfrcaqc[yidxl2,xidxl2,:,:,1].flatten().astype(numpy.int16)
                ncldarr[yidxout,xidxout,:] = ncldair[yidxl2,xidxl2,:,:].flatten().astype(numpy.int16)

                if (cfrcair[yidxl2,xidxl2,1,1,0] < 0.0):
                    print(cfrcair[yidxl2,xidxl2,1,1,0])
                    frcstr = '%d, %d: %.4f: ' % (yidxout,xidxout,cld1arr[yidxout,xidxout,4])
                    print(frcstr)

        c1chk = cld1arr[:,:,4].flatten()
        c1sb = c1chk[c1chk >= 0.0]
        print(c1sb.shape)

    print(tmch)

    # Output
    qout = Dataset(flnm,'r+')

    varcfrc1 = qout.variables['AIRS_CldFrac_1']
    varcfrc1[tmidx,:,:,:] = cld1arr[:,:,:]

    varcfrc2 = qout.variables['AIRS_CldFrac_2']
    varcfrc2[tmidx,:,:,:] = cld2arr[:,:,:]

    varcfqc1 = qout.variables['AIRS_CldFrac_QC_1']
    varcfqc1[tmidx,:,:,:] = cld1qc[:,:,:]

    varcfqc2 = qout.variables['AIRS_CldFrac_QC_2']
    varcfqc2[tmidx,:,:,:] = cld2qc[:,:,:]

    varncld = qout.variables['AIRS_nCld']
    varncld[tmidx,:,:,:] = ncldarr[:,:,:]

    qout.close()

    return

def merra_find_cloud_slab( srchdr, srchdt, lnsq, ltsq, outfl, tidx, msgval=-9999.):
    ### Identify cloud slabs from MERRA cloud content profiles
    ### srchdr:     Directory with MERRA daily files
    ### srchdt:     Desired date (a python datetime object)
    ### lnsq:       Longitude subset sequence
    ### ltsq:       Latitude subset sequence
    ### outfl:      Output file name
    ### tidx:       Output time index

    flst = os.listdir(srchdr)
    dtstr = srchdt.strftime('%Y%m%d')

    mincldvl = 1.0e-7

    f2d = -1
    f3d = -1
    j = 0

    print(dtstr)
    while ( (j < len(flst)) and ((f2d < 0) or (f3d < 0) ) ):
        if ( ('inst1_2d_asm' in flst[j]) and (dtstr in flst[j]) ):
            f2d = j
            mer2d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer2d,'r')
            psfc = f.variables['PS'][0,ltsq,lnsq]
            f.close()
        if ( ('inst3_3d_asm' in flst[j]) and (dtstr in flst[j]) ):
            # Temp, QI, QL in 3D files
            f3d = j
            mer3d = '%s/%s' % (srchdr,flst[j])
            f = Dataset(mer3d,'r')
            qlmer = f.variables['QL'][0,:,ltsq,lnsq]
            qimer = f.variables['QI'][0,:,ltsq,lnsq]
            tmpmer = f.variables['T'][0,:,ltsq,lnsq]
            lvmer = f.variables['lev'][:]
            f.close()
        j = j + 1

    #slbtyp = numpy.array([-99,-99],dtype=numpy.int16)
    # Set up reference frame
    nlon = lnsq.shape[0]
    nlat = ltsq.shape[0]
    #ltidx = numpy.repeat(ltsq,nlon)
    #lnidx = numpy.tile(lnsq,nlat)

    # Liquid, ice indicators
    qimer[qimer > 1.0e10] = 0.0
    qlmer[qlmer > 1.0e10] = 0.0

    qiind = (qimer >= mincldvl)
    #qiind = ((qimer >= mincldvl) & (qimer > qlmer))
    qiind.dtype = numpy.int8
    qlind = (qlmer >= mincldvl)
    #qlind = ((qlmer >= mincldvl) & (qlmer > qimer))
    qlind.dtype = numpy.int8

    qism = numpy.sum(qiind,axis=0)
    qlsm = numpy.sum(qlind,axis=0)

    # Output arrays
    nslb = numpy.zeros(qism.shape,dtype=numpy.int16) + 1
    nslb[(qism == 0) & (qlsm == 0)] = 0

    ctyp1 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval
    ctyp2 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval
    cpsz1 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval
    cpsz2 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval

    cpbt1 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval
    cpbt2 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval
    cptp1 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval
    cptp2 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval

    cttp1 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval
    cttp2 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval
    cngwt1 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval
    cngwt2 = numpy.zeros(qism.shape,dtype=numpy.float32) + msgval

    # Difference booleans
    print('Differencing')
    qidif = numpy.diff(qiind,n=1,axis=0)
    print(qidif.shape)
    print(numpy.amax(qidif))
    print(numpy.amin(qidif))
    qldif = numpy.diff(qlind,n=1,axis=0)

    nlvmer = lvmer.shape[0]
    lvsq = numpy.arange(lvmer.shape[0])
    lvsqdf = numpy.arange(1,lvmer.shape[0])

    for j in range(nlat):
        for i in range(nlon):
            ltspt = j + ltsq[0]
            lnspt = i + lnsq[0]
            slbct = 0
            if (nslb[j,i] > 0):
                # Freezing level
                lvfr = lvsq[(tmpmer[:,j,i] < 273.15) & (tmpmer[:,j,i] > 0.0)]
                if (lvfr.shape[0] == lvsq.shape[0]):
                    frzlvl = psfc[j,i] * 0.01
                    #print('Sfc Frz Level')
                else:
                    dllgp = (numpy.log(lvmer[lvfr[0]-1]) - numpy.log(lvmer[lvfr[0]])) * \
                                    (273.15 - tmpmer[lvfr[0],j,i]) / \
                                    (tmpmer[lvfr[0]-1,j,i] - tmpmer[lvfr[0],j,i])
                    frzlvl = lvmer[lvfr[0]] * numpy.exp(dllgp)
                    #if ( (tmpmer[lvfr[0]-1,j,i] == tmpmer[lvfr[0],j,i]) ):
                    #    str1 = 'Equal temp at freezing level %d' (lvfr)
                    #    print(str1)
                    #if ( (lvmer[lvfr[0]] == 0.0) or (lvmer[lvfr[0]-1] == 0.0)):
                    #    str1 = 'Equal temp at freezing level %d' (lvfr)
                    #    print(str1)
                    #    print(tmpmer[:,j,i])
                nlqd = 0
                nice = 0
                if (qism[j,i] > 1):
                    # Ice slab
                    litp = lvsqdf[qidif[:,j,i] == -1]
                    libt = lvsqdf[qidif[:,j,i] == 1]
                    if (litp.shape[0] != libt.shape[0]):
                        # Surface case
                        libt = numpy.append([0],libt)
                    nice = litp.shape[0]
                    iwp = numpy.zeros(nice,dtype=numpy.float32)
                    pbt = numpy.zeros(nice,dtype=numpy.float32)
                    ptp = numpy.zeros(nice,dtype=numpy.float32)
                    ctt = numpy.zeros(nice,dtype=numpy.float32)
                    pstrt = (psfc[j,i]-1.0) * 0.01
                    for sl1 in range(nice):
                        # Find bottom pres
                        if libt[sl1] > 0:
                            lst1 = libt[sl1] - 1
                            dllgp = (numpy.log(lvmer[lst1]) - numpy.log(lvmer[lst1+1])) * \
                                    (mincldvl - qimer[lst1+1,j,i]) / \
                                    (qimer[lst1,j,i] - qimer[lst1+1,j,i])
                            pchk = lvmer[lst1+1] * numpy.exp(dllgp)
                        else:
                            lst1 = libt[sl1]
                            # Cloud at bottom level
                            lst1 = llbt[sl1]
                            pchk = (psfc[j,i]-1.0) * 0.01
                        if (pchk < pstrt):
                            pstrt = pchk
                        pbt[sl1] = pstrt

                        # Top pres
                        lst2 = litp[sl1] - 1
                        # Find top pres
                        dllgp = (numpy.log(lvmer[lst2]) - numpy.log(lvmer[lst2+1])) * \
                                (mincldvl - qimer[lst2+1,j,i]) / \
                                (qimer[lst2,j,i] - qimer[lst2+1,j,i])
                        pchk = lvmer[lst2+1] * numpy.exp(dllgp)
                        if (pchk > pbt[sl1]):
                            pchk = pbt[sl1] - 10.0
                        ptp[sl1] = pchk
                        pstrt = pchk

                        # Temperature
                        lwwt = (numpy.log(lvmer[lst2+1]) - numpy.log(ptp[sl1])) / \
                               (numpy.log(lvmer[lst2]) - numpy.log(lvmer[lst2+1]))
                        hiwt = 1.0 - lwwt
                        ctt[sl1] = lwwt * tmpmer[lst2,j,i] + hiwt * tmpmer[lst2+1,j,i]

                        # Integrate IWP (kg m^-2)
                        iwp[sl1] = 0.0
                        for slw in range(lst1,lst2):
                            dlprs = (lvmer[slw] - lvmer[slw+1]) * 100.0
                            iwp[sl1] = iwp[sl1] + (qimer[slw,j,i] + qimer[slw+1,j,i]) * dlprs / 9.8
                    icefrm = pandas.DataFrame({'PresBot': pbt, 'PresTop': ptp, 'CTTemp': ctt, \
                                               'WtrPath': iwp})
                    icefrm['CldType'] = 'Ice'
                if (qlsm[j,i] > 1):
                    # Water slabs
                    lltp = lvsqdf[qldif[:,j,i] == -1]
                    llbt = lvsqdf[qldif[:,j,i] == 1]
                    if (lltp.shape[0] != llbt.shape[0]):
                        # Surface case
                        llbt = numpy.append([0],llbt)
                    if ( (j == 21) and (i == 11) ):
                        print(lltp)
                        print(llbt)
                    nlqd = lltp.shape[0]
                    pbt = numpy.zeros(nlqd,dtype=numpy.float32)
                    ptp = numpy.zeros(nlqd,dtype=numpy.float32)
                    ctt = numpy.zeros(nlqd,dtype=numpy.float32)
                    lwp = numpy.zeros(nlqd,dtype=numpy.float32)
                    pstrt = (psfc[j,i]-1.0) * 0.01
                    for sl1 in range(nlqd):
                        if llbt[sl1] > 0:
                            lst1 = llbt[sl1] - 1
                            dllgp = (numpy.log(lvmer[lst1]) - numpy.log(lvmer[lst1+1])) * \
                                    (mincldvl - qlmer[lst1+1,j,i]) / \
                                    (qlmer[lst1,j,i] - qlmer[lst1+1,j,i])
                            pchk = lvmer[lst1+1] * numpy.exp(dllgp)
                        else:
                            # Cloud at bottom level
                            lst1 = llbt[sl1]
                            pchk = (psfc[j,i]-1.0) * 0.01

                        # Find bottom pres
                        if (pchk < pstrt):
                            pstrt = pchk
                        pbt[sl1] = pstrt

                        # Top pres
                        lst2 = lltp[sl1] - 1
                        # Find top pres
                        dllgp = (numpy.log(lvmer[lst2]) - numpy.log(lvmer[lst2+1])) * \
                                (mincldvl - qlmer[lst2+1,j,i]) / \
                                (qlmer[lst2,j,i] - qlmer[lst2+1,j,i])
                        pchk = lvmer[lst2+1] * numpy.exp(dllgp)
                        if (pchk > pbt[sl1]):
                            pchk = pbt[sl1] - 10.0
                        ptp[sl1] = pchk
                        pstrt = pchk

                        # Temperature
                        lwwt = (numpy.log(ptp[sl1]) - numpy.log(lvmer[lst2+1])) / \
                               (numpy.log(lvmer[lst2]) - numpy.log(lvmer[lst2+1]))
                        hiwt = 1.0 - lwwt
                        ctt[sl1] = lwwt * tmpmer[lst2,j,i] + hiwt * tmpmer[lst2+1,j,i]

                        # Integrate LWP (kg m^-2)
                        lwp[sl1] = 0.0
                        for slw in range(lst1,lst2):
                            dlprs = (lvmer[slw] - lvmer[slw+1]) * 100.0
                            lwp[sl1] = lwp[sl1] + (qlmer[slw,j,i] + qlmer[slw+1,j,i]) * dlprs / 9.8
                    wtrfrm = pandas.DataFrame({'PresBot': pbt, 'PresTop': ptp, 'CTTemp': ctt, \
                                               'WtrPath': lwp})
                    wtrfrm['CldType'] = 'Water'
                if ( (nice > 0) and (nlqd >  0)):
                    cldfrm = icefrm.append(wtrfrm,ignore_index=True)
                elif (nice > 0):
                    cldfrm = icefrm
                elif (nlqd > 0):
                    cldfrm = wtrfrm
                else:
                    nstr = 'No slab found for lat %d, lon %d' % (j,i)
                    slbct = 0

                if ( (nice > 0) or (nlqd > 0) ):
                    cldfrm['FrzLvl'] = frzlvl
                    cldfrm['DPCloud'] = cldfrm['PresBot'] - cldfrm['PresTop']
                    cldfrm['DPPhase'] = 0.0
                    cldfrm.loc[ (cldfrm['CldType'] == 'Water') &  \
                                (cldfrm['PresTop'] < cldfrm['FrzLvl']),'DPPhase' ] = cldfrm['PresTop']-cldfrm['FrzLvl']
                    cldfrm.loc[ (cldfrm['CldType'] == 'Ice') &  \
                                (cldfrm['PresBot'] > cldfrm['FrzLvl']),'DPPhase' ] = cldfrm['FrzLvl']-cldfrm['PresBot']
                    cldfrm['AdjWtrPath'] = (cldfrm['DPCloud'] + cldfrm['DPPhase']) * cldfrm['WtrPath'] / cldfrm['DPCloud']
                    cldfrm = cldfrm.sort_values(by=['AdjWtrPath'],ascending=[False])

                    if ( (j == 21) and (i == 11) ):
                        print(psfc[j,i]*0.01)
                        print(cldfrm)
                        print(lwwt)
                        print(qlmer[:,j,i])
                        print(qimer[:,j,i])

                    # Final slab selection
                    if cldfrm.shape[0] >= 2:
                        cldslc = cldfrm[0:2]
                        slbct = 2
                        # Check overlap and re-calc
                        chg2 = False
                        if ((cldslc['PresBot'].values[0] > cldslc['PresBot'].values[1]) and \
                            (cldslc['PresTop'].values[0] < cldslc['PresBot'].values[1])):
                            cldslc['PresBot'].values[1] = cldslc['PresTop'].values[0] - 10.0
                            chg2 = True
                        elif ((cldslc['PresBot'].values[0] < cldslc['PresBot'].values[1]) and \
                              (cldslc['PresTop'].values[1] < cldslc['PresBot'].values[0])):
                            cldslc['PresTop'].values[1] = cldslc['PresBot'].values[0] + 10.0
                            chg2 = True
                        if chg2:
                            if cldslc['CldType'].values[1] == 'Ice':
                                cpth = qimer[:,j,i]
                            elif cldslc['CldType'].values[1] == 'Water':
                                cpth = qlmer[:,j,i]
                            pbt2 = cldslc['PresBot'].values[1]
                            ptp2 = cldslc['PresTop'].values[1]
                            lvrg = lvsq[(lvmer <= pbt2) & (lvmer >= ptp2)]
                            if lvrg.shape[0] == 0:
                                # Remove the slab!
                                cldslc = cldslc[0:1]
                                slbct = 1
                            else:
                                lst1 = lvrg[0]

                                # Temperature
                                lst2 = lvrg[lvrg.shape[0]-1]
                                lwwt = (numpy.log(ptp2) - numpy.log(lvmer[lst2+1])) / \
                                       (numpy.log(lvmer[lst2]) - numpy.log(lvmer[lst2+1]))
                                hiwt = 1.0 - lwwt
                                cldslc['CTTemp'].values[1] = lwwt * tmpmer[lst2,j,i] + hiwt * tmpmer[lst2+1,j,i]
                                if ( (j == 25) and (i == 27) ):
                                    strt = 'Temp adjust, Lst: %d, Pres[lst]: %.4f, PTP: %.4f' % (lst2,lvmer[lst2], ptp2)
                                    print(strt)
                                    print(lwwt)

                                # Integrate LWP (kg m^-2)
                                pth2 = 0.0
                                for slw in range(lst1,lst2):
                                    dlprs = (lvmer[slw] - lvmer[slw+1]) * 100.0
                                    pth2 = pth2 + (cpth[slw] + cpth[slw+1]) * dlprs / 9.8
                                cldslc['WtrPath'].values[1] = pth2
                        # Finally sort vertically
                        cldslc = cldslc.sort_values(by=['PresBot'],ascending=[False])
                    elif cldfrm.shape[0] == 1:
                        cldslc = cldfrm
                        slbct = 1
                    cldslc = cldslc.reset_index(drop=True)
                    if ( (j == 21) and (i == 11) ):
                        print(cldslc)
                        #print(tmpmer[:,j,i])

            # Output arrays
            nslb[j,i] = slbct
            if slbct >= 1:
                if cldslc['CldType'].values[0] == 'Water':
                    ctyp1[j,i] = 101.0
                    cpsz1[j,i] = 20.0
                elif cldslc['CldType'].values[0] == 'Ice':
                    ctyp1[j,i] = 201.0
                    cpsz1[j,i] = 80.0
                cngwt1[j,i] = cldslc['WtrPath'].values[0]
                cpbt1[j,i] = cldslc['PresBot'].values[0]
                cptp1[j,i] = cldslc['PresTop'].values[0]
                cttp1[j,i] = cldslc['CTTemp'].values[0]
            if slbct == 2:
                if cldslc['CldType'].values[1] == 'Water':
                    ctyp2[j,i] = 101.0
                    cpsz2[j,i] = 20.0
                elif cldslc['CldType'].values[1] == 'Ice':
                    ctyp2[j,i] = 201.0
                    cpsz2[j,i] = 80.0
                cngwt2[j,i] = cldslc['WtrPath'].values[1]
                cpbt2[j,i] = cldslc['PresBot'].values[1]
                cptp2[j,i] = cldslc['PresTop'].values[1]
                cttp2[j,i] = cldslc['CTTemp'].values[1]

    # Output
    qout = Dataset(outfl,'r+')

    varctp1 = qout.variables['ctype1']
    varctp1[tidx,:,:] = ctyp1[:,:]

    varctp2 = qout.variables['ctype2']
    varctp2[tidx,:,:] = ctyp2[:,:]

    varpsz1 = qout.variables['cpsize1']
    varpsz1[tidx,:,:] = cpsz1[:,:]

    varpsz2 = qout.variables['cpsize2']
    varpsz2[tidx,:,:] = cpsz2[:,:]

    varpbt1 = qout.variables['cprbot1']
    varpbt1[tidx,:,:] = cpbt1[:,:]

    varpbt2 = qout.variables['cprbot2']
    varpbt2[tidx,:,:] = cpbt2[:,:]

    varptp1 = qout.variables['cprtop1']
    varptp1[tidx,:,:] = cptp1[:,:]

    varptp2 = qout.variables['cprtop2']
    varptp2[tidx,:,:] = cptp2[:,:]

    varctt1 = qout.variables['cstemp1']
    varctt1[tidx,:,:] = cttp1[:,:]

    varctt2 = qout.variables['cstemp2']
    varctt2[tidx,:,:] = cttp2[:,:]

    # Convert ngwat to g m^-2
    cngwt1[cngwt1 > 0.0] = cngwt1[cngwt1 > 0] * 1000.0
    varngw1 = qout.variables['cngwat1']
    varngw1[tidx,:,:] = cngwt1[:,:]

    cngwt2[cngwt2 > 0.0] = cngwt2[cngwt2 > 0] * 1000.0
    varngw2 = qout.variables['cngwat2']
    varngw2[tidx,:,:] = cngwt2[:,:]

    qout.close()

    return

def airs_granule_overlap(mtfl, yrchc, mnchc, dychc, grnchc, mskvr, mskvl, \
                           l2srch = '/archive/AIRSOps/airs/gdaac/v6'):

    # Find range of AIRS scan rows overlapping a template region
    # mtfl:    MERRA file with mask information
    # yrchc:   Year
    # mnchc:   Month
    # dychc:   Day
    # grnchc:  Granule
    # mskvr:   Name of region mask variable
    # mskvl:   Value of region mask for Region Choice

    # Mask, lat, lon
    f = Dataset(mtfl,'r')
    mask = f.variables[mskvr][:,:]
    latmet = f.variables['lat'][:]
    lonmet = f.variables['lon'][:]
    tminf = f.variables['time'][:]
    tmunit = f.variables['time'].units[:]
    f.close()

    mskind = numpy.zeros((mask.shape),dtype=mask.dtype)
    print(mskvl)
    mskind[mask == mskvl] = 1
    lnsq = numpy.arange(lonmet.shape[0])
    ltsq = numpy.arange(latmet.shape[0])

    # Subset a bit
    lnsm = numpy.sum(mskind,axis=0)
    ltsm = numpy.sum(mskind,axis=1)

    lnmn = numpy.amin(lnsq[lnsm > 0])
    lnmx = numpy.amax(lnsq[lnsm > 0]) + 1
    ltmn = numpy.amin(ltsq[ltsm > 0])
    ltmx = numpy.amax(ltsq[ltsm > 0]) + 1

    stridx = 'Lon Range: %d, %d\nLat Range: %d, %d \n' % (lnmn,lnmx,ltmn,ltmx)
    print(stridx)

    nx = lnmx - lnmn
    ny = ltmx - ltmn
    nzout = 101

    lnrp = numpy.tile(lonmet[lnmn:lnmx],ny)
    ltrp = numpy.repeat(latmet[ltmn:ltmx],nx)
    mskblk = mskind[ltmn:ltmx,lnmn:lnmx]
    mskflt = mskblk.flatten()

    # Set up reference frame
    #ltrp = numpy.repeat(lats,nlon)
    #ltidx = numpy.repeat(numpy.arange(nlat),nlon)
    #lnrp = numpy.tile(lons,nlat)
    #lnidx = numpy.tile(numpy.arange(nlon),nlat)

    merfrm = pandas.DataFrame({'GridLon': lnrp, 'GridLat': ltrp, 'MaskInd': mskflt})
    print(merfrm.shape)
    mersb = merfrm[ merfrm['MaskInd'] == 1]
    print(mersb.shape)

    # Find reference granule
    # Search AIRS Level 2
    airsdr = '%s/%04d/%02d/%02d/airs2ret' % (l2srch,yrchc,mnchc,dychc)
    l2fd = -1
    if (os.path.exists(airsdr)):
        fllst = os.listdir(airsdr)

        ldstr = 'AIRS.%04d.%02d.%02d.%03d' % (yrchc, mnchc, dychc, grnchc)
        for j in range(len(fllst)):
            lncr = len(fllst[j])
            l4 = lncr - 4
            if ((fllst[j][l4:lncr] == '.hdf') and (ldstr in fllst[j])):
                l2fl = '%s/%s' % (airsdr,fllst[j])
                ncl2 = Dataset(l2fl)
                l2lat = ncl2.variables['Latitude'][:,:]
                l2lon = ncl2.variables['Longitude'][:,:]
                ncl2.close()
                l2fd = j
                print(l2lat[22,15])
                print(l2lon[22,15])

    rwsq = numpy.arange(45)
    rwrp = numpy.repeat(rwsq,30)
    if l2fd >= 0:
        l2lnflt = l2lon.flatten().astype(numpy.float64)
        l2ltflt = l2lat.flatten().astype(numpy.float64)
        l2frm = pandas.DataFrame({'L2RowIdx': rwrp, \
                                  'L2Lon': l2lnflt, 'L2Lat': l2ltflt})
        l2frm['GridLon'] = numpy.around(l2frm['L2Lon']/0.625) * 0.625
        l2frm['GridLat'] = numpy.around(l2frm['L2Lat']/0.5) * 0.5

        l2mrg = pandas.merge(l2frm,mersb,on=['GridLon','GridLat'])
        #print(l2mrg.shape)
        #print(l2mrg[0:10])

        grprw = l2mrg.groupby('L2RowIdx')
        rwcts = grprw.apply(df_ct_rows,grpvr='L2RowIdx')
        rwcts.reset_index(drop=True,inplace=True)

        rwctsb = rwcts[rwcts['Count'] >= 3]
        mxrw = rwctsb['L2RowIdx'].max()
        mnrw = rwctsb['L2RowIdx'].min()
        mdrw = math.floor( (mxrw + mnrw) / 2.0)

        arout = numpy.array([mnrw, mdrw, mxrw], dtype=numpy.int16)
        arout = arout + 1

    return arout


def df_ct_rows(dtfrm,grpvr):
    # Count of data frame rows
    fct = dtfrm.shape[0]
    xvl = dtfrm[grpvr].values[0]
    dfout = pandas.DataFrame({grpvr: xvl, \
                              'Count': fct}, index=[0])
    return(dfout)
