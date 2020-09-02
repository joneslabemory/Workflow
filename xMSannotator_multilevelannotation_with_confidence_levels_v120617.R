
library(xMSannotator)
library(xmsPANDA)

#Package data files
data(example_data)  #example peak intensity matrix
data(adduct_table)
data(adduct_weights)
data(customIDs)  #example for custom IDs
data(customDB)   #example for custom DB
#data(hmdbAllinf) 
#data(keggotherinf)
#data(t3dbotherinf)


###########Parameters to change##############
#dataA<-example_data #use example data provided with the package
#######################################

peak_intensity_matrix<-"G:\\Medicine\\Pulmonary_ISILON\\Research\\MaHPIC_Metabolomics\\ProcessedData\\Experiment13_5timepoints\\AE\\Experiment13_selected_time_points_feature_table_only_AE_vjan132014_PANDA.txt"

#output location
outloc<-"G:\\Medicine\\Pulmonary_ISILON\\Research\\MaHPIC_Metabolomics\\ProcessedData\\Experiment13_5timepoints\\AE\\xMSannotatortestv1.3.2\\"

max.mz.diff<-10  #mass search tolerance for DB matching in ppm
max.rt.diff<-10 #retention time tolerance between adducts/isotopes
corthresh<-0.7 #correlation threshold between adducts/isotopes
max_isp=5 #maximum number of isotopes to search for
mass_defect_window=0.01 #mass defect window for isotope search

num_nodes<-2   #number of cores to be used; 2 is recommended for desktop computers due to high memory consumption

db_name="HMDB" #other options: HMDB,Custom,KEGG, LipidMaps, T3DB
status=NA #other options: "Detected", NA, "Detected and Quantified", "Expected and Not Quantified"
num_sets<-300 #number of sets into which the total number of database entries should be split into;

ionization_mode<-"pos" #pos or neg, options for ionization mode

#provide list of database IDs (depending upon selected database) for annotating only specific metabolites
customIDs<-NA #c("HMDB15352","HMDB60549","HMDB00159","HMDB00222"); read.csv("/Users/mzmatch_95stdmx_HMDBIDs.csv")

#provide your own custom database to be used for annotation
#set db_name="Custom" if you use this option
#Format: ID, Name, Formula, MonoisotopicMass
customDB<-NA #read.table("/Users/karanuppal/Documents/Emory/JonesLab/Projects/xMSannotator/IROA/IROA_customDB_xMSannotator_plate1.txt",sep="\t",header=TRUE) #custom database; default NA

#number of technical replicates
num_replicates<-3

#########################

if(num_replicates==1){
dataA<-read.table(peak_intensity_matrix, sep="\t", header=TRUE)
}else{

dataA<-read.table(peak_intensity_matrix, sep="\t", header=TRUE)
avg_data_res<-data_preprocess(X=dataA,feature_table_file=NA,parentoutput_dir=outloc,class_labels_file=NA,num_replicates=num_replicates,feat.filt.thresh=NA,summarize.replicates=TRUE,summary.method="median",all.missing.thresh=0.01,group.missing.thresh=NA,log2transform=FALSE,medcenter=FALSE,znormtransform=FALSE,quantile_norm=FALSE,lowess_norm=FALSE,madscaling=FALSE,missing.val=0,samplermindex=NA, rep.max.missing.thresh=0.5,summary.na.replacement="zeros")


dataA<-avg_data_res$data_matrix_afternorm_scaling #avg_data5$data_matrix_prescaling
}

#########################


if(ionization_mode=="neg"){
filter.by=c("M-H")
queryadductlist=c("M-H","M-H2O-H","M+Na-2H","M+Cl","M+FA-H")

}else{
filter.by=c("M+H")
queryadductlist=c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H","M+2Na-H","M+H-H2O","M+H-2H2O")

}


print(paste("Input feature table file is: ",peak_intensity_matrix,sep=""))
print(paste("Ionization mode is: ",ionization_mode,sep=""))
print(paste("Number of replicates is: ",num_replicates,sep=""))
print(paste("Adducts used for annotations: ",sep=""))
print(queryadductlist)


dataA<-unique(dataA)
print(dim(dataA))
print(format(Sys.time(), "%a %b %d %X %Y"))

system.time(annotres<-multilevelannotation(dataA=dataA,max.mz.diff=max.mz.diff,max.rt.diff=max.rt.diff,cormethod="pearson",num_nodes=num_nodes,queryadductlist=queryadductlist,
mode=ionization_mode,outloc=outloc,db_name=db_name, adduct_weights=adduct_weights,num_sets=num_sets,allsteps=TRUE,
corthresh=corthresh,NOPS_check=TRUE,customIDs=customIDs,missing.value=NA,deepsplit=2,networktype="unsigned",
minclustsize=10,module.merge.dissimilarity=0.2,filter.by=filter.by,biofluid.location=NA,origin=NA,status=status,boostIDs=NA,max_isp=max_isp,
customDB=customDB,
HMDBselect="union",mass_defect_window=mass_defect_window,pathwaycheckmode="pm",mass_defect_mode="pos")
)


print(format(Sys.time(), "%a %b %d %X %Y"))
