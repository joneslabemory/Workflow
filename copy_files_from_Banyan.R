

#sequence file with column A as the instrument file name
#floc=read.csv("F:\\QSTDs\\HFQE\\NIST_QSTD_CHEAR\\seq_data_nist_chear_qstd.csv")

BanyanPath="/home/data/Bassig_NCI_HFQE/"

isilon_path_rawfiles="G:\\Medicine\\Pulmonary_ISILON\\Research\\Jones_Lab\\_Orbitrap_RAW\\HFQE_Raw\\"



dir.create(local_path)

num_batches=2

#files<-paste(BanyanPath,"/",floc[,1],".raw",sep="")

library(ssh)

username_hostipaddress<-"kuppal2@170.140.231.74"

password_val<-"19inchnails"

session <- ssh_connect(username_hostipaddress,passwd=password_val)


res<-lapply(1:num_batches,function(j)
{

dirname_1=paste(BanyanPath,"batch",j,sep="")

scp_download(session, dirname_1, to = local_path, verbose = TRUE)
}
)