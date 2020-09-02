library(apLCMS)

###########Parameters to change##############################
numnodes<-4

#create known table


###apLCMS code
cdfloc<-"F:\\ken\\c18\\"
apLCMSoutloc<-"F:\\ken\\c18\\aplcms\\"
#apLCMSoutloc= paste(cdfloc, "apLCMSoutput/", sep= "")
dir.create(apLCMSoutloc)

#run_list<- read.table("F:/Benzene_QEHF/Complete_sample_list_BEW.txt", sep= "\t", header=TRUE)
#batch_info<- run_list[seq(1,nrow(run_list),2),]


file.pattern = ".cdf"

cdf.files<-list.files(cdfloc,file.pattern)
min.pres<-c(0.5,0.8)   
min.run<- c(4,3)
min.exp<-round(0.1*length(cdf.files),0)
mz.tol<-10e-6
align.mz.tol<-10e-6
align.chr.tol = 45
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
component.eliminate= 0.01
moment.power=2




for(ii in 1:length(min.pres)){
#call cdf.to.ftr for extraction
print(paste("Settings", ii, "of", length(min.pres), sep=" "))

aligned.hyb<-cdf.to.ftr(folder=cdfloc,min.exp=min.exp,min.pres = min.pres[ii], min.run = min.run[ii],
mz.tol = mz.tol, align.mz.tol = align.mz.tol, align.chr.tol = align.chr.tol, n.nodes=numnodes, file.pattern=file.pattern, baseline.correct.noise.percentile=baseline.correct.noise.percentile,
shape.model = shape.model, baseline.correct = baseline.correct, peak.estim.method = peak.estim.method, 
min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, sigma.ratio.lim = sigma.ratio.lim, subs = subs,  
max.align.mz.diff = max.align.mz.diff, pre.process = pre.process, recover.mz.range = recover.mz.range, 
recover.chr.range = recover.chr.range, use.observed.range = use.observed.range, recover.min.count = recover.min.count, component.eliminate=0.01, moment.power=2)

finalfeattable<-aligned.hyb$final.ftrs

setwd(apLCMSoutloc)

#fname<-paste("apLCMShybrid_aligned_feature_table_min.run",min.run,"min.pres",min.pres,"min.exp",min.exp,".txt",sep="")
fname1<-paste("apLCMSuntargetedv",packageVersion("apLCMS"),"_min.run",min.run[ii],"min.pres",min.pres[ii],"min.exp",min.exp,"mztol",mz.tol,"alignmztol",align.mz.tol,"alignchr",align.chr.tol,"recover.mz",recover.mz.range,"recover.chr",recover.chr.range,".txt",sep="")
fname2<-paste("apLCMSuntargetedv",packageVersion("apLCMS"),"_min.run",min.run[ii],"min.pres",min.pres[ii],"min.exp",min.exp,"mztol",mz.tol,"alignmztol",align.mz.tol,"alignchr",align.chr.tol,"recover.mz",recover.mz.range,"recover.chr",recover.chr.range,".Rda",sep="")

write.table(finalfeattable,file=fname1,sep="\t",row.names=FALSE)
save(aligned.hyb,file=fname2)
print("Dimension of feature table")
print(dim(aligned.hyb$final.ftrs))

}















