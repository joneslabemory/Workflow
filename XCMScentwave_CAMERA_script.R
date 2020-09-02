library(xcms)
library(parallel)
library(CAMERA)
###########Parameters to change##############################
numnodes<-4 #detectCores()-2
numnodes_recovery<-0.5*(numnodes)

cdfloc<-"C:\\cdfloc\\"

XCMS.outloc<-"C:\\xcmsOutloc\\"

bw<-2 #recommended Nature Protocol #bw5

mzdiff<-(-0.001)

ppmthresh<-2.5 

snthresh<-3

minfrac.val<-0.1 #percent files in which a non-missing signal is present

mzwid.val=0.015 #recommended Nature Protocol

max=50;

retcor.method="obiwarp";

peakwidth=c(5,20); #recommended Nature Protocol #c(10,60)

prefilter.list=c(3,1000)

polarity="positive" #"negative"

groupval.method="medret"
step.list=c(0.01);
nSlaves=numnodes
minsamp.val=1; sleep.val=0; run.order.file=NA;subs=NA;retcor.family="symmetric";
retcor.plottype="deviation";target.mz.list = NA

dir.create(XCMS.outloc)

setwd(cdfloc) #location of CDF files
dir.create(XCMS.outloc,showWarnings=FALSE) #location of XCMS output
cdf_files=list.files(cdfloc,".cdf|.mzxml|.mzML",ignore.case=TRUE)

	
m<-mzdiff
b<-bw
p<-ppmthresh
t<-snthresh
s<-step.list[1]
numsamp=length(cdf_files)
        
#added on 7/2/2019
snowparam <- SnowParam(workers = numnodes, type = "SOCK")

#1) xcmsSet
fname=paste(XCMS.outloc, "/XCMSstep1centwave","_snthresh", t,"_mzdiff", m,"_pw",peakwidth[1],"_",peakwidth[2],"_bw",b,"_ppm",p,".Rda", sep="")
check_if_exists<-suppressWarnings(try(load(fname),silent=TRUE))

if(is(check_if_exists,"try-error")){

xset=xcmsSet(cdf_files, method="centWave", ppm=p, snthresh=t,mzdiff=m,peakwidth=peakwidth,prefilter=prefilter.list,
					integrate=1, verbose.columns=TRUE,fitgauss=FALSE, BPPARAM = snowparam) #nSlaves=nSlaves)
					
save(xset,file=fname)
}else{
	print(paste("Loading: ",fname),sep="")
	
	load(fname)
}

####2) group.density
fname=paste(XCMS.outloc, "/XCMSstep2centwave","_snthresh", t,"_mzdiff", m,"_pw",peakwidth[1],"_",peakwidth[2],"_bw",b,"_ppm",p,".Rda", sep="")
check_if_exists<-suppressWarnings(try(load(fname),silent=TRUE))

if(is(check_if_exists,"try-error")){

xset <- group(xset, bw=b, minfrac=minfrac.val, minsamp=minsamp.val,  mzwid=mzwid.val, max=max, sleep=sleep.val)

save(xset,file=fname)
}else{

	print(paste("Loading: ",fname),sep="")
	
	load(fname)
}

fname=paste(XCMS.outloc, "/XCMSstep3centwave","_snthresh", t,"_mzdiff", m,"_pw",peakwidth[1],"_",peakwidth[2],"_bw",b,"_ppm",p,".Rda", sep="")
check_if_exists<-suppressWarnings(try(load(fname),silent=TRUE))

if(is(check_if_exists,"try-error")){
	if(retcor.method=="loess"){
					
					#####3) retention time correction using loess; mdevden
                                        xset2 <- retcor(xset, method=retcor.method,family = retcor.family, plottype = retcor.plottype)
	}else{
						if(retcor.method=="obiwarp"){
						###3) retention time correction using obiwarp method
						xset2 <- retcor(xset, method=retcor.method,plottype = retcor.plottype, profStep=s)
						}
						else{
						stop("Please enter obiwarp or loess as retention time correction method.")
						}
	}


	save(xset2,file=fname)
}else{
	print(paste("Loading: ",fname),sep="")
	
	load(fname)
}

#4) group

fname=paste(XCMS.outloc, "/XCMSstep4centwave","_snthresh", t,"_mzdiff", m,"_pw",peakwidth[1],"_",peakwidth[2],"_bw",b,"_ppm",p,".Rda", sep="")
check_if_exists<-suppressWarnings(try(load(fname),silent=TRUE))

if(is(check_if_exists,"try-error")){
xset2 <- group(xset2, bw=b, minfrac=minfrac.val, minsamp=minsamp.val,  mzwid=mzwid.val, max=max, sleep=sleep.val)

save(xset2,file=fname)
}else{
	print(paste("Loading: ",fname),sep="")
	
	load(fname)
}

#5) fillPeaks

fname=paste(XCMS.outloc, "/XCMSstep5centwave","_snthresh", t,"_mzdiff", m,"_pw",peakwidth[1],"_",peakwidth[2],"_bw",b,"_ppm",p,".Rda", sep="")
check_if_exists<-suppressWarnings(try(load(fname),silent=TRUE))

if(is(check_if_exists,"try-error")){
xset3=suppressWarnings(fillPeaks(xset2,method="chrom"))
save(xset3,file=fname)

}else{
	print(paste("Loading: ",fname),sep="")
	
	load(fname)
}

print("Getting sample intensities")
finalfeatmat={}
					
finalfeatmat<- peakTable(object=xset3)

fname=paste(XCMS.outloc, "/XCMS_centwave","_snthresh", t,"_mzdiff", m,"_max",max,"_bw",b,"_ppm",p,".Rda", sep="")
save(xset3,file=fname)

fname=paste(XCMS.outloc, "/XCMS_centwave","_snthresh", t,"_mzdiff", m,"_max",max,"_bw",b,"_ppm",p,".txt", sep="")
cnames=colnames(finalfeatmat)
cnames[1]="mz"
cnames[4]="time"
colnames(finalfeatmat)=c(cnames[1:8],cdf_files)
cnames<-tolower(cnames)
cnames<-gsub(".cdf", "", cnames)

write.table(finalfeatmat,file=fname,sep="\t",row.names=FALSE)

#run CAMERA
xsa <- xsAnnotate(xset3)
#Group after RT value of the xcms grouped peak
xsaF <- groupFWHM(xsa, perfwhm=0.6)
#Verify grouping
xsaC <- groupCorr(xsaF)
#Annotate isotopes, could be done before groupCorr
xsaFI <- findIsotopes(xsaC)
#Annotate adducts
xsaFA <- findAdducts(xsaFI, polarity=polarity)

fname=paste(XCMS.outloc, "/CAMERAresult_centwave","_snthresh", t,"_mzdiff", m,"_max",max,"_bw",b,"_ppm",p,".csv", sep="")
#Get final peaktable and store on harddrive
write.csv(getPeaklist(xsaFA),file=fname)
		
fname=paste(XCMS.outloc, "/CAMERAresult_centwave","_snthresh", t,"_mzdiff", m,"_max",max,"_bw",b,"_ppm",p,".Rda", sep="")		
save(xsaFA,file=fname)
					setwd(XCMS.outloc)
					suppressWarnings(try(unlink("XCMSstep*.Rda",force=TRUE,recursive=TRUE),silent=TRUE))
