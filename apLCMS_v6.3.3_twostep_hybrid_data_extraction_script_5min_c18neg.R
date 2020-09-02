library(apLCMS)
library(gtools)
library(xMSanalyzer)

print(packageVersion("apLCMS"))

####


min.run<-c(3,4)
min.pres<-c(0.8,0.5)
mz.tol<-2.5e-6
align.mz.tol<-2.5e-6
align.chr.tol<-10
target.mz.tol.ppm=2.5
target.time.tol=10
cdfloc="/home/ravneet.kaur/Projects/QSTD/c18neg/adjustedIROA/"
known.table.common.pos="/home/ravneet.kaur/Projects/QSTD/scripts/IROA_known.table.c18neg.txt"
apLCMSoutloc="/home/ravneet.kaur/Projects/QSTD/apLCMSnocalib2step_c18neg_adj/"
batch_info="/home/ravneet.kaur/Projects/QSTD/scripts/batchinfo_qstd3_c18neg.txt"


numnodes<-6


known.table.common.pos<-read.table(known.table.common.pos,sep= "\t", header=TRUE)

dir.create(apLCMSoutloc)
#batch_info<-as.character(batch_info)
batch_info<- read.table(batch_info,sep= "\t", header=TRUE)

batch_info$Batch<-as.numeric(factor(batch_info$Batch))

batch_info$Batch<-as.numeric(factor(batch_info$Batch))

known.table.common.pos<-known.table.common.pos[!is.na(known.table.common.pos[,6]),]
known.table.common.pos[,7:12]<-NA
known.table.common.pos<-cbind(known.table.common.pos, known.table.common.pos[,7:12])
file.pattern=".mzXML"
cdf.files<-list.files(cdfloc,file.pattern)
match_tol_ppm=5
baseline.correct.noise.percentile = 0.25
shape.model = "bi-Gaussian"
baseline.correct = NA
peak.estim.method = "moment"
min.bw = NA
max.bw = NA
sd.cut = c(0.125, 15)
sigma.ratio.lim = c(0.33,3)
subs = NA
max.align.mz.diff = 0.01
pre.process = FALSE
recover.mz.range = 10e-6
recover.chr.range = 45
use.observed.range = FALSE
recover.min.count = 1
new_feature_min_count=4
component.eliminate= 0.03
moment.power=2
min.within.batch.num.detect=3
min.within.batch.num.report=3
##min.within.batch.prop=0.001
min.batch.prop=0.1
min.within.batch.prop.detect = 0.1
min.within.batch.prop.report = 0.1
batch.align.mz.tol=align.mz.tol
batch.align.chr.tol=align.chr.tol

for(ii in 1:length(min.pres)){
  setwd(cdfloc)
  aligned.hyb<-two.step.hybrid(folder=cdfloc,
                               n.nodes=numnodes, file.pattern=file.pattern,
                               known.table=known.table.common.pos, sd.cut=sd.cut,sigma.ratio.lim=sigma.ratio.lim, component.eliminate=component.eliminate,
                               moment.power=moment.power, min.pres=min.pres[ii], min.run=min.run[ii], mz.tol=mz.tol, baseline.correct.noise.percentile=baseline.correct.noise.percentile,
                               align.mz.tol=align.mz.tol, align.chr.tol=align.chr.tol, max.align.mz.diff=max.align.mz.diff, recover.mz.range=recover.mz.range,
                               recover.chr.range=recover.chr.range,use.observed.range=use.observed.range, shape.model=shape.model,
                               new.feature.min.count=new_feature_min_count,
                               recover.min.count=recover.min.count,
                               #batchwise parameter
                               info=batch_info,
                               #min.within.batch.num.detect=min.within.batch.num.detect, min.within.batch.num.report=min.within.batch.num.report,
                                min.within.batch.prop.detect=min.within.batch.prop.detect,
                                min.within.batch.prop.report=min.within.batch.prop.report,
                                min.batch.prop=min.batch.prop, batch.align.mz.tol=batch.align.mz.tol, batch.align.chr.tol= batch.align.chr.tol)
  finalfeattable<-aligned.hyb$final.ftrs
  setwd(apLCMSoutloc)
  fname1<-paste("apLCMShybridv6.3.3_min.run",min.run[ii],"min.pres",min.pres[ii],"min.within.batch.num.detect",min.within.batch.num.detect,"mztol",mz.tol,"alignmztol",align.mz.tol,"alignchr",align.chr.tol,".txt",sep="")
  fname2<-paste("apLCMShybridv6.3.3_min.run",min.run[ii],"min.pres",min.pres[ii],"min.within.batch.num.detect",min.within.batch.num.detect,"mztol",mz.tol,"alignmztol",align.mz.tol,"alignchr",align.chr.tol,".rds",sep="")
  write.table(finalfeattable,file=fname1,sep="\t",row.names=FALSE)
  print("Dimension of feature table")
  print(dim(aligned.hyb$final.ftrs))
  saveRDS(aligned.hyb,file=fname2)
}
