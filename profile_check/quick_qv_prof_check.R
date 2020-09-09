# Check QV profile

# Check on temperature profs

library(ncdf4)
library(reshape2)
library(ggplot2)
library(plyr)
library(colorspace)

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 12
theme_mat$axis.text.x$size = 11
theme_mat$axis.title.y$size = 12
theme_mat$axis.text.y$size = 11
theme_mat$plot.title$size = 14
theme_mat$plot.title$hjust = 0.5
theme_mat$legend.text$size = 11
theme_mat$legend.title$size = 11

xidx = 25
yidx = 25
tidx = 31

# Pathological cases [x,y,t]
#  25, 25, 31
#  33, 43, 2
#  21, 28, 27

nc1 = nc_open("MERRA/CONUS_Output/interpolated_merra2_for_SARTA_two_slab_2019_MAM_CONUS_NCA_09UTC.nc")
qvprf = ncvar_get(nc1,"qv")
rhprf = ncvar_get(nc1,"rh")
psfc = ncvar_get(nc1,"spres")
lvsrt = as.vector(ncvar_get(nc1,"lev"))
lnsrt = as.vector(ncvar_get(nc1,"lon"))
ltsrt = as.vector(ncvar_get(nc1,"lat"))
nc_close(nc1)

lvsqsrt = seq(51,98)

lnidx = array(rep(1:92,48),c(92,48))
ltidx = array(rep(1:48,each=92),c(92,48))

ttmp = qvprf[,,90,31]
ptmp = psfc[,,31]
lnsb = lnidx[ (ttmp < 0) & !is.na(ttmp)]
ltsb = ltidx[ (ttmp < 0) & !is.na(ttmp)]

# Specific comparison
merfl = sprintf("MERRA/Hour09/MERRA2_400.inst3_3d_asm_Np.201903%02d.SUB.nc",tidx)
nc2 = nc_open(merfl)
lvmer = ncvar_get(nc2,"lev")
lnmer = ncvar_get(nc2,"lon")
ltmer = ncvar_get(nc2,"lat")
qvmer = ncvar_get(nc2,"QV")
nc_close(nc2)

xsq = seq(1,length(lnmer))
ysq = seq(1,length(ltmer))


xsb = xsq[lnmer == lnsrt[xidx]]
ysb = ysq[ltmer == ltsrt[yidx]]
lvsqmer = seq(1,24)

# Plot
srtfrm = data.frame(QV=qvprf[xidx,yidx,lvsqsrt,tidx],Level=lvsrt[lvsqsrt],RH=rhprf[xidx,yidx,lvsqsrt,tidx])
srtfrm$LogRev = 7 - log(srtfrm$Level)
merfrm = data.frame(QV=qvmer[xsb,ysb,lvsqmer],Level=lvmer[lvsqmer])
merfrm$LogRev = 7 - log(merfrm$Level)


pspts = c(1000,700,500,300,200)
pspts2 = c(1000,925,850,700)
lspts = 7 - log(pspts)
lspts2 = 7 - log(pspts2)

gttl = sprintf("QV Profile %.2f, %.2f",lnsrt[xidx],ltsrt[yidx])
gqv = ggplot(srtfrm,aes(x=QV,y=LogRev)) + geom_point() + 
  geom_point(data=merfrm,shape=4) + 
  xlab("Specific humidity") + scale_y_continuous("Pressure",breaks=lspts,labels=pspts) + 
  theme_mat + ggtitle(gttl)

gttl = sprintf("RH Profile %.2f, %.2f",lnsrt[xidx],ltsrt[yidx])
grh = ggplot(srtfrm,aes(x=RH,y=LogRev)) + geom_point() + 
  xlab("Relative humidity") + scale_y_continuous("Pressure",breaks=lspts,labels=pspts) + 
  theme_mat + ggtitle(gttl)

