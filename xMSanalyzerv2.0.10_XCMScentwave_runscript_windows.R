
#load xMSanalyzer package
library(xMSanalyzer)
library(xMSannotator)
library(RColorBrewer)
data(keggCompMZ)


source("G:\\Medicine\\Pulmonary_ISILON\\Research\\Jones_Lab\\Orbitrap\\Rscripts_data_extraction\\Windows\\xMSanalyzer_2.0.10.R")
data(example_target_list_pos)
data(example_target_list_neg)
###Input parameters###########

#1) cdfloc: The folder where all CDF files to be processed are located. For example "C:/CDF/"
# Note: set cdfloc=NA if the cdf files are already aligned using XCMS and the results 
# exist in XCMS.outloc
cdfloc="/home/kuppal2/Documents/Projects/KenLiu/Germ_Free_Mouse_Mito/neg_c18/"

#2) Note: Feature table at each individual parameter setting (just like XCMS)
XCMSoutloc="/home/kuppal2/Documents/Projects/KenLiu/Germ_Free_Mouse_Mito/neg_c18/xcmsoutput/"

#Note: Merged feature table (like XCMS, but with feature quality summary)
#3) xMSanalyzer.outloc: The folder where xMSanalyzer output will be written. For example "C:/xMSanalyzeroutput/"
xMSanalyzeroutloc="/home/kuppal2/Documents/Projects/KenLiu/Germ_Free_Mouse_Mito/neg_c18/xMSanalyzerv2.0.7xcmscentwavev12/"

#4) Sequence file path
sample_info_file<-NA #"~/Documents/Projects/KenLiu/Rotifer/neg_c18_sample_mapping.txt" #"~/Documents//Projects/faahKO/sequence_file_pos.txt"

#5) reference chemicals; use NA for the example_target_list provided with the package
reference_chemicals_file<-NA #"~/Documents/Projects/xMSanalyzer/valid_chem_mz.txt"

#6) charge type
charge_type="pos"

#7) calibration method
calibration.method="median.adjustment" #"multiplicative.signal.correction"

#########################################END of Input parameters##########################################################

dir.create(XCMSoutloc,showWarnings=FALSE)
dir.create(xMSanalyzeroutloc,showWarnings=FALSE)

########xMSanalyzer usage##################

par(mfrow=c(2,2))
pdf("Rplots.pdf")
	result<-try(
	{

      res.list<-xMSwrapper.XCMS.centWave(cdfloc=cdfloc, XCMS.outloc=XCMSoutloc,xMSanalyzer.outloc=xMSanalyzeroutloc,ppm.list=c(2.5), mz.diff.list=c(-0.00005), sn.thresh.list=c(3), prefilter.list=c(3,1000), bw.val=c(5,10),groupval.method="medret", 
step.list=c(0.1),max=50,minfrac.val=0.1, minsamp.val=1, mzwid.val=0.015, sleep.val=0,retcor.method="obiwarp",retcor.family="symmetric", retcor.plottype="deviation", peakwidth=c(10,60), 
numnodes=4,run.order.file=NA,max.mz.diff=5,max.rt.diff=300, merge.eval.pvalue=0.05,mergecorthresh=0.7,deltamzminmax.tol=100,
num_replicates=3,subs=NA, mz.tolerance.dbmatch=10, adduct.list=NA, samp.filt.thresh=0.70,
feat.filt.thresh=100,cormethod="pearson",mult.test.cor=TRUE,missingvalue=0,ignore.missing=TRUE,
sample_info_file=sample_info_file,refMZ=reference_chemicals_file,refMZ.mz.diff=10,refMZ.time.diff=NA,void.vol.timethresh=30,
replacezeroswithNA=TRUE,scoreweight=30,filepattern=".mzXML",
db_name=c("KEGG","HMDB","LipidMaps"),charge_type=charge_type,
data.norm.pipeline="AC",calibration.method=calibration.method)

	})

dev.off()
setwd(xMSanalyzeroutloc)
save(res.list,file="xMSwrapper_XCMS_centwave.Rda")


##############
