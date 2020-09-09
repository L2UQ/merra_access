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

nc1 = nc_open("AIRS_L2_Ascending_2013.08.05.224.nc")
tprf = ncvar_get(nc1,"TAirSup")
tqc = ncvar_get(nc1,"TAirSup_QC")
lvs = ncvar_get(nc1,"XtraPressureLev")
nc_close(nc1)

set.seed(352523)

dimnames(tprf)[1] = list(lvs)
tmpmlt = melt(tprf,varnames=c("Level","XSpt","YSpt"),value.name = "Temp")

dimnames(tqc)[1] = list(lvs)
qcmlt = melt(tqc,varnames=c("Level","XSpt","YSpt"),value.name = "TempQC")

tmrg = merge(tmpmlt,qcmlt)
tmrg = tmrg[ (tmrg$TempQC < 2) & (!is.na(tmrg$Temp)),]

XSpt = seq(1,30)
YSpt = seq(1,45)
merrf = expand.grid(XSpt,YSpt)
names(merrf) = c("XSpt","YSpt")
msmp = sort(sample(seq(1,nrow(merrf)),size=200,replace=FALSE  )  )
mersmp = merrf[msmp,]
mersmp$SdSpt = (mersmp$YSpt-1) * 100 + mersmp$XSpt

mermrg = merge(mersmp,tmrg)
mermrg$LogRev = 7 - log(mermrg$Level)
mermrg$SdSpt = factor(mermrg$SdSpt)
mermrg = mermrg[order(mermrg$SdSpt,mermrg$Level),]

pspts = c(1000,500,200,100,50,20,10,1,0.5,0.1,0.05,0.01)
lspts = 7 - log(pspts)
gttl = "AIRS Profiles"
gtmp = ggplot(mermrg,aes(x=Temp,y=LogRev,group=SdSpt)) + geom_path(color="#333333") + 
  xlab("Temperature [K]") + scale_y_continuous("Pressure",breaks=lspts,labels=pspts) + 
  theme_mat + ggtitle(gttl)
pfnm = "AIRS_TempProf_MAGIC.pdf"
pdf(pfnm,width=10,height=8)
print(gtmp)
dev.off()

