#Author: Chunyu Ma
#Date: Jan 28, 2019
#This script is used to automatically separate the sequencing files to two indivival files based on columns

library('gtools')

#1) raw sequencing files location.
raw_files_loc="E:\\LatentTB_Collins\\sequencefile\\"

#2) output file location.
out_files_loc="E:\\LatentTB_Collins\\sequencefile\\"

#3) Input file format. Options: 'TXT' or 'CSV'
informat="CSV"

#4) Extracted Columns.
columns=c('File Name','Sample ID','Comment')

#5) Output file format. Options: 'TXT' or 'CSV'
outformat="TXT" 

#6) Output file names. #1 for odd number #2 for even number 
filename=c('LatentTB_mapping_hilicpos.txt','LatentTB_mapping_c18neg.txt')

############################No changes needed below this line###########################
dir.create(out_files_loc,showWarnings=FALSE)
files=mixedsort(list.files(raw_files_loc))
files=paste(raw_files_loc,'\\',mixedsort(list.files(raw_files_loc)),sep='')

for(i in 1:length(files)){
	
	if(informat=='CSV'){
		tmpfile=read.table(files[i],sep=',',header=TRUE,skip=1,stringsAsFactors=FALSE,check.names=FALSE)
		tmpfile=tmpfile[,columns]
		
		if(i==1){
			odd=tmpfile[seq(1,dim(tmpfile)[1],2),]
			even=tmpfile[seq(2,dim(tmpfile)[1],2),]
		}else{
			odd=rbind(odd,tmpfile[seq(1,dim(tmpfile)[1],2),])
			even=rbind(even,tmpfile[seq(2,dim(tmpfile)[1],2),])
		}
	
	}else{
		tmpfile=read.table(files[i],sep='\t',header=TRUE,skip=1,stringsAsFactors=FALSE,check.names=FALSE)
		tmpfile=tmpfile[,columns]
		
		if(i==1){
			odd=tmpfile[seq(1,dim(tmpfile)[1],2),]
			even=tmpfile[seq(2,dim(tmpfile)[1],2),]
		}else{
			odd=rbind(odd,tmpfile[seq(1,dim(tmpfile)[1],2),])
			even=rbind(even,tmpfile[seq(2,dim(tmpfile)[1],2),])
		}
		
	}
	

}

colnames(odd)[grep('Comment',colnames(odd))]='Batch'
colnames(even)[grep('Comment',colnames(even))]='Batch'

for(i in 1:2){
			
	if(i==1){
		out=odd
	}else{
		out=even
	}
	
	if(outformat=='CSV'){
				
		write.table(out,paste(out_files_loc,filename[i],sep='\\'),sep=',',row.names=FALSE)
	
	}else{
				
		write.table(out,paste(out_files_loc,filename[i],sep='\\'),sep='\t',row.names=FALSE)

	}
}