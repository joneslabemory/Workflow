
library(apLCMS)
library(RColorBrewer)
library(xMSanalyzer)
library(xMSannotator)
library(parallel)
data(keggCompMZ)


#please see the manual for description of functions and parameters

source("G:\\Medicine\\Pulmonary_ISILON\\Research\\Jones_Lab\\Orbitrap\\Rscripts_data_extraction\\Windows\\xMSanalyzer_2.0.10.R")


data(example_target_list_pos)
data(example_target_list_neg)
###Input parameters###########
#1) cdfloc: The folder where all CDF files to be processed are located. For example "C:/CDF/"
# Note: set cdfloc=NA if the cdf files are already aligned using apLCMS and the results 
# exist in apLCMS.outloc

cdfloc=NA #"E:\\24_25_18\\c18neg\\"

#Note: Feature table at each individual parameter setting (just like apLCMS)
#2) apLCMS.outloc: The folder where alignment output will be written. For example "C:/CDFoutput/"

apLCMSoutloc="C:\\Users\\kuppal2\\Downloads\\apLCMSout\\" 


#Note: Merged feature table (like apLCMS, but with feature quality summary)
#3) xMSanalyzer.outloc: The folder where xMSanalyzer output will be written. For example "C:/xMSanalyzeroutput/"

xMSanalyzeroutloc="C:\\Users\\kuppal2\\Downloads\\xMSchdwbposC\\" 

#4) Sequence file path; Need for batch-effect evaluation; eg: "C:/Documents/Emory/JonesLab/Projects/pos/sequence_file_pos.txt"
#Column A: Names matching .cdf or .mzXML files
#Column B: Sample ID/name
#Column C: Batch (column should be labeled "Batch")

sample_info_file<-"G:\\Medicine\\Pulmonary_ISILON\\Research\\Jones_Lab\\Orbitrap\\cdf\\ViLinh\\HFQE\\CHDWB_Obese\\aug 2019\\CHDWB_obese_mapping_hilicpos.txt" 
 
#5) reference chemicals; use NA for the example_target_list provided with the package
# eg:"C:/Documents/Emory/JonesLab/Projects/xMSanalyzer/valid_chem_mz.txt"
reference_chemicals_file<-"C:\\Users\\kuppal2\\Downloads\\5minhilicpos_ref_list_xmsanalyzer_calibration.txt"


#6) Ionization mode: use "pos" for positive and AE; use "neg" for negative and C18
charge_type="pos"

#7) Length of chromatography: 300 for 5 min method; 600 for 10 min method
max_retention_time<-300

#8) file pattern: ".cdf" or ".mzXML"
filepattern=".mzXML"

#9) calibration method
calibration.method="median.adjustment" #"multiplicative.signal.correction"

#########################################END of Input parameters##########################################################

dir.create(apLCMSoutloc,showWarnings=FALSE)
dir.create(xMSanalyzeroutloc,showWarnings=FALSE)

numnodes <- detectCores() - 1
numnodes<-4 #round(numnodes*0.6)

########xMSanalyzer usage##################
if(max_retention_time>300){
alignchrtol = 45
}else{
alignchrtol = 30
}

	result<-try(
	{
	par(mfrow=c(2,2))
	pdf("Rplots.pdf")	
		res.list<-xMSwrapper.apLCMS(cdfloc=cdfloc, apLCMS.outloc=apLCMSoutloc, xMSanalyzer.outloc=xMSanalyzeroutloc, 
		
		#apLCMSv6.3 or above arguments
		min.run.list = c(4,3), min.pres.list = c(0.5,0.8), minexp.pct = 0.1, mztol = 2.5e-6, alignmztol = 10e-6,       
		alignchrtol = alignchrtol, numnodes = numnodes, apLCMSmode="untargeted", baseline.correct.noise.percentile = 0.25,                  
		shape.model = "bi-Gaussian", baseline.correct = NA, peak.estim.method = "moment", 
		min.bw = NA, max.bw = NA, sd.cut = c(0.125, 15), sigma.ratio.lim = c(0.33,3), subs = NA, 
		max.align.mz.diff = 0.01, pre.process = FALSE, recover.mz.range = 10e-6,     
		recover.chr.range = alignchrtol, use.observed.range = FALSE, recover.min.count = 1, new_feature_min_count=4,component.eliminate= 0.01,
		moment.power=2, known_table=NA, match_tol_ppm=5,
	
		#xMSanalyzerv2.0.8 arguments
                max.mz.diff = 5, max.rt.diff = 0.5*(max_retention_time), merge.eval.pvalue = 0.05, mergecorthresh = 0.7, deltamzminmax.tol = 100,
		num_replicates = 3, 
		mz.tolerance.dbmatch = 5, adduct.list = c("M+H"), samp.filt.thresh = 0.7,  
		feat.filt.thresh = 75, cormethod = "pearson", mult.test.cor = FALSE,         
		missingvalue = 0, ignore.missing = TRUE, filepattern = filepattern,
		sample_info_file=sample_info_file,refMZ=reference_chemicals_file,refMZ.mz.diff=5,refMZ.time.diff=30,void.vol.timethresh=30,
		replacezeroswithNA=TRUE,charge_type=charge_type,plotEICs="target",rawprofileloc=cdfloc,
		peak.score.thresh=NA,
		reference_sample_index = NA, merge.pairwise = FALSE, summarize.replicates=TRUE,
		summary.method="median",max.number.of.replicates.with.missingvalue=1,summary.na.replacement="zeros",
		db_name=c("KEGG","HMDB","LipidMaps"),qc_label="q3june2014",data.norm.pipeline="AC",calibration.method=calibration.method
		
		#end
		)                
                 
                 
     try(dev.off(),silent=TRUE)

	})

