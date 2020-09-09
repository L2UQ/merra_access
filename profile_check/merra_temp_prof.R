library(ncdf4)
library(colorspace)
library(ggplot2)
library(plyr)
library(reshape2)

theme_mat = theme_bw() 
theme_mat$axis.title.x$size = 14
theme_mat$axis.text.x$size = 14
theme_mat$axis.title.y$size = 14
theme_mat$axis.text.y$size = 14
theme_mat$plot.title$size = 18
theme_mat$plot.title$hjust = 0.5

xidx = 21
yidx = 28
tidx = 27

nc1 = nc_open("MERRA/CONUS_Output/interpolated_merra2_for_SARTA_two_slab_2019_MAM_CONUS_NCA_09UTC.nc")
lvsrt = as.vector(ncvar_get(nc1,"lev"))
lnsrt = as.vector(ncvar_get(nc1,"lon"))
ltsrt = as.vector(ncvar_get(nc1,"lat"))
nc_close(nc1)

# Specific comparison
merfl = sprintf("MERRA/Hour09/MERRA2_400.inst3_3d_asm_Np.201903%02d.SUB.nc",tidx)
nc2 = nc_open(merfl)
lvmer = ncvar_get(nc2,"lev")
lnmer = ncvar_get(nc2,"lon")
ltmer = ncvar_get(nc2,"lat")
tmpmer = ncvar_get(nc2,"T")
nc_close(nc2)

xsq = seq(1,length(lnmer))
ysq = seq(1,length(ltmer))

xsb = xsq[lnmer == lnsrt[xidx]]
ysb = ysq[ltmer == ltsrt[yidx]]

set.seed(352523)

XSpt = seq(xsb-21,xsb+21)
YSpt = seq(ysb-15,ysb+15)
tmrsb = tmpmer[XSpt,YSpt,]

dimnames(tmrsb)[3] = list(lvmer)
dimnames(tmrsb)[2] = list(YSpt)
dimnames(tmrsb)[1] = list(XSpt)
tmpmlt = melt(tmrsb,varnames=c("XSpt","YSpt","Level"),value.name = "Temp")

merrf = expand.grid(XSpt,YSpt)
names(merrf) = c("XSpt","YSpt")
msmp = sort(sample(seq(1,nrow(merrf)),size=210,replace=FALSE  )  )
mersmp = merrf[msmp,]
mersmp$SdSpt = (mersmp$YSpt-1) * 100 + mersmp$XSpt

mermrg = merge(mersmp,tmpmlt)
mermrg$LogRev = 7 - log(mermrg$Level)
mermrg$SdSpt = factor(mermrg$SdSpt)
mermrg = mermrg[order(mermrg$SdSpt,mermrg$Level),]

pspts = c(1000,500,200,100,50,20,10,1,0.5,0.1)
lspts = 7 - log(pspts)
gttl = sprintf("MERRA Temp Profile %.2f, %.2f",lnsrt[xidx],ltsrt[yidx])
gtmp = ggplot(mermrg,aes(x=Temp,y=LogRev,group=SdSpt)) + geom_path(color="#333333") + 
  xlab("Temperature [K]") + scale_y_continuous("Pressure",breaks=lspts,labels=pspts) + 
  theme_mat + ggtitle(gttl)
pfnm = "MERRA_TempProf_CONUS.pdf"
pdf(pfnm,width=10,height=8)
print(gtmp)
dev.off()

