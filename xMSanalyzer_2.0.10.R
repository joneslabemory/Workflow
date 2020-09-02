find.peak.regions<-function (x, nups = 1, ndowns = nups, zero = "0", peakpat = NULL, 
    minpeakheight = -Inf, minpeakdistance = 1, threshold = 0, 
    npeaks = 0, sortstr = FALSE) 
{

    #pracma
    stopifnot(is.vector(x, mode = "numeric") || length(is.na(x)) == 
        0)
    if (!zero %in% c("0", "+", "-")) 
        stop("Argument 'zero' can only be '0', '+', or '-'.")
    xc <- paste(as.character(sign(diff(x))), collapse = "")
    xc <- gsub("1", "+", gsub("-1", "-", xc))
    if (zero != "0") 
        xc <- gsub("0", zero, xc)
    if (is.null(peakpat)) {
        peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
    }
    rc <- gregexpr(peakpat, xc)[[1]]
    if (rc[1] < 0) 
        return(NULL)
    x1 <- rc
    x2 <- rc + attr(rc, "match.length")
    attributes(x1) <- NULL
    attributes(x2) <- NULL
    n <- length(x1)
    xv <- xp <- numeric(n)
    for (i in 1:n) {
        xp[i] <- which.max(x[x1[i]:x2[i]]) + x1[i] - 1
        xv[i] <- x[xp[i]]
    }
    inds <- which(xv >= minpeakheight & xv - pmax(x[x1], x[x2]) >= 
        threshold)
    X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])
    if (minpeakdistance < 1) 
        warning("Handling 'minpeakdistance < 1' is logically not possible.")
    if (sortstr || minpeakdistance > 1) {
        sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
        X <- X[sl, , drop = FALSE]
    }
    if (length(X) == 0) 
        return(c())
    if (minpeakdistance > 1) {
        no_peaks <- nrow(X)
        badpeaks <- rep(FALSE, no_peaks)
        for (i in 1:no_peaks) {
            ipos <- X[i, 2]
            if (!badpeaks[i]) {
                dpos <- abs(ipos - X[, 2])
                badpeaks <- badpeaks | (dpos > 0 & dpos < minpeakdistance)
            }
        }
        X <- X[!badpeaks, , drop = FALSE]
    }
    if (npeaks > 0 && npeaks < nrow(X)) {
        X <- X[1:npeaks, , drop = FALSE]
    }
    return(X)
}

call.supervised.peak.detection<-function(filenum,filenames,targetfeature_list,ppm_error,mzXML.loc,plot.peak=FALSE,snr.threshold=1,basepeak.height.threshold=1000,quality.score.threshold=1,void.time=20,
step=0.001,basepeak.select=FALSE,iqr.prop.nonmissing.thresh=0.1,min.peak.length.rt=3,deltatime=30,profstep=0.1,sample_names=NA,min.pw=10,max.pw=60,peak.method="custom"){

setwd(mzXML.loc)
Sys.sleep(1)

tempobject_fname<-paste(filenames[filenum],"temp.Rda",sep="")

if(is.na(sample_names)==FALSE){
cur_filename=gsub(filenames[filenum],pattern=".mzXML|.cdf|.CDF",replacement="")

sample_names[,1]=gsub(sample_names[,1],pattern=".mzXML|.cdf|.CDF",replacement="")

sample_name=sample_names[which(sample_names[,1]==cur_filename),2]

}

#save(filenum,filenames,targetfeature_list,ppm_error,mzXML.loc,plot.peak,snr.threshold,basepeak.height.threshold,quality.score.threshold,void.time,step,basepeak.select,iqr.prop.nonmissing.thresh,min.peak.length.rt,file=tempobject_fname)

#rm(list=ls())

#load(tempobject_fname)

xraw_fname<-paste(filenames[filenum],"xraw.Rda",sep="")

check_xrawfiles<-list.files(".","xraw.Rda$")

m1<-grep(xraw_fname,check_xrawfiles)

if(length(m1)>0){
	check_fload<-try(load(xraw_fname),silent=TRUE)

	if(is(check_fload,"try-error")){

		xraw<-xcmsRaw(filenames[filenum],includeMSn=FALSE,profstep=profstep)
	#	save(xraw,file=xraw_fname)
	}
}else{
	xraw<-xcmsRaw(filenames[filenum],includeMSn=FALSE,profstep=profstep)
	#save(xraw,file=xraw_fname)
}

mzrange<-xraw@mzrange
filename1<-filenames[filenum]


targetfeature_list<-as.data.frame(targetfeature_list)

if(ncol(targetfeature_list)==1){
	targetmz_list<-targetfeature_list[which(targetfeature_list>=mzrange[1] & targetfeature_list<=mzrange[2]),1]
	targettime_list=rep(NA,length(targetmz_list))

}else{
	if(ncol(targetfeature_list)>=2){
	targetmz_list<-targetfeature_list[which(targetfeature_list>=mzrange[1] & targetfeature_list<=mzrange[2]),1]
	targettime_list<-targetfeature_list[which(targetfeature_list>=mzrange[1] & targetfeature_list<=mzrange[2]),2]
	}
}

if(peak.method=="centwave"){
			
					
			print("Using the CentWave algorithm for peak detection")
}else{
			print("Using time and intensity patterns for peak detection")
}

data_list<-lapply(1:length(targetmz_list),function(t1){

targetmz<-targetmz_list[t1]
targettime<-targettime_list[t1]


if(is.na(targetmz)==FALSE){
peaks<-supervised.peak.detection(raw_fname=filenames[filenum],targetmz=targetmz,targettime=targettime,ppm_error=ppm_error,xcms.raw=xraw,plot.peak=plot.peak,snr.threshold=snr.threshold,basepeak.height.threshold=basepeak.height.threshold,
quality.score.threshold=quality.score.threshold,void.time=void.time,step=step,basepeak.select=basepeak.select,iqr.prop.nonmissing.thresh=iqr.prop.nonmissing.thresh,min.peak.length.rt=min.peak.length.rt,deltatime=deltatime,profstep=profstep,
sample_name=sample_name,min.pw=min.pw,max.pw=max.pw,peak.method=peak.method)

rm(xraw)


if(length(peaks)<1){

	peaks<-c(targetmz,rep(0,14))
	peaks<-t(peaks)
	peaks<-as.data.frame(peaks)
	
}else{

if(length(peaks)>0){
	if(nrow(peaks)<1){

		peaks<-c(targetmz,rep(0,14))
		peaks<-t(peaks)
		peaks<-as.data.frame(peaks)
		
	}else{
	
		#peaks<-peaks[order(peaks$basepeak.height,decreasing=TRUE),]
	}
	}
}
colnames(peaks)<-c("mz","maxo","rt","min.rt","max.rt","rtmin.1","rtmax.1","into","intb","sn","peak_area","sn.ratio","basepeak.height","basepeak.rt","peak.quality.score")

peaks<-peaks[,-c(6:7)]

peaknum=t1
peaks<-cbind(filenum,peaknum,peaks)

peaks<-as.data.frame(peaks)



#print(mem_used())

tempobject_fname<-paste(filenames[filenum],"temp.Rda",sep="")
#save(,file=tempobject_fname)

return(peaks)
}
})

xraw_fname<-paste(filenames[filenum],"_",filenum,"datalist.Rda",sep="")
save(data_list,file=xraw_fname)
#print(mem_used())
return(data_list)
#rm(list=ls())

}


find.peak.regions2<-function (x, thresh = 0) 
{
    pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 
        0) + 2
    if (!missing(thresh)) {
        if (sign(thresh) < 0) 
            thresh <- -thresh
        pks[x[pks - 1] - coredata(x[pks]) > thresh]
    }
    else pks
}


#supervised.peak.detection(raw_fname="F005_161012_001.mzXML",targetmz=175.1189,targettime=NA,ppm_error=5)
#source("C:\\Users\\kuppal2\\Documents\\xMSanalyzer\\MSMS_process\\MSMS_process\\xMSanalyzer_2.0.8.25.R")
#pdf("Test_methionine.pdf")
#supervised.peak.detection(raw_fname="F005_161012_001.mzXML",targetmz=150.0580,targettime=NA,ppm_error=5,xcms.raw=xraw,plot.peak=TRUE)
#dev.off()

supervised.peak.detection<-function(xcms.raw=NA,raw_fname,targetmz,targettime=NA,ppm_error=5,plot.peak=FALSE,snr.threshold=1,basepeak.height.threshold=1000,
quality.score.threshold=1,void.time=20,step=0.001,basepeak.select=FALSE,iqr.prop.nonmissing.thresh=0.9,min.peak.length.rt=3,deltatime=30,profstep=0.1,sample_name=NA,min.pw=10,max.pw=60,peak.method="centwave")
{
min_mz<-targetmz-(ppm_error*(10^-6)*targetmz)
max_mz<-targetmz+(ppm_error*(10^-6)*targetmz)

if(is.na(iqr.prop.nonmissing.thresh)==TRUE){
	iqr.prop.nonmissing.thresh=0
}


###Below code for extracting peak data####

#print(mem_used())
check_class<-class(xcms.raw)

if(check_class=="logical"){
xraw<-xcmsRaw(raw_fname,includeMSn=FALSE,profstep=0.1)
}else{
xraw<-xcms.raw
}

intvec<-getTIC(targetmz=targetmz,targettime=targettime,deltappm=ppm_error,deltatime=NA,profstep=profstep,object=xraw)
     

peak.method=tolower(peak.method)
rm(xraw)
#print(mem_used())
#if(length(intvec@eic[[1]][[1]])>0){
if(length(intvec)>0){



	xvec<-intvec[,1]
	yvec<-intvec[,2]

if(max(yvec)>basepeak.height.threshold){

	#plot(xvec,yvec)
	fnametemp=paste("peaks",targetmz,".Rda",sep="")
	#save(xvec,intvec,yvec,file=fnametemp)
	
	# define function that returns the SSE
calcSSE <- function(x,xvec,yvec){
  
  loessMod <- try(loess(yvec ~ xvec, span=x), silent=T)
  res <- try(loessMod$residuals, silent=T)
  
 
  if(class(res)!="try-error"){
    if((sum(res, na.rm=T) > 0)){
      sse <- sum(res^2)  
    }
  }else{
    sse <- 99999
  }
  return(sse)
}


	

	fnametemp=paste("peaks",targetmz,raw_fname,".Rda",sep="")
	
	#print(fnametemp)
	
	check_if_exists<-0 #suppressWarnings(try(load(fnametemp),silent=TRUE))

	#if(is(check_if_exists,"try-error"))
	
	if(check_if_exists==0)
	{
	
	
	if(peak.method=="centwave"){
		
		peaks<-matrix(1,nrow=3,ncol=5)
	
		
	}else{
	
	peaks<-find.peak.regions(yvec,nups=1,ndowns=1,threshold=(1),minpeakheight=(10000),minpeakdistance=(1))
		
	}
	
	#save(intvec,yvec,min.pw,max.pw,peaks,file=fnametemp)
	if(length(peaks)>0)
	{
		
		if(nrow(peaks)>1)
		{
		
		
			if(peak.method=="centwave"){
			
					
			#print("Using the CentWave algorithm for peak detection")
			
			peaks <- peaksWithCentWave(intvec[,2], intvec[,1],peakwidth=c(min.pw,max.pw))
			
			peaks1<-as.data.frame(peaks)
			
			}else{
			
					#print("Using time and intensity patterns for peak detection")
			
					peaks1<-peaks[order(peaks[,3]),]
				
				
					peaks1<-cbind(seq(1,nrow(peaks1)),peaks1)
					
					
					colnames(peaks1)<-c("ID","baseheight","basetime","start","end")
					
					peaks1<-as.data.frame(peaks1)
					median_base_height=1000
					for(p1 in 2:nrow(peaks1)){
					
							if((peaks1$start[p1]-peaks1$end[p1-1])<2 & peaks1$baseheight[p1]>=median_base_height){
								
									peaks1$ID[p1]=peaks1$ID[p1-1]
									
									start_range<-min(peaks1$start[which(peaks1$ID==peaks1$ID[p1-1])])
									end_range<-max(peaks1$end[which(peaks1$ID==peaks1$ID[p1-1])])
									median_base_height<-sd(peaks1$baseheight[which(peaks1$ID==peaks1$ID[p1])]) #max(peaks1$baseheight[which(peaks1$ID==peaks1$ID[p1-1])])-min(peaks1$baseheight[which(peaks1$ID==peaks1$ID[p1-1])]) #quantile(yvec,0.1) #[start_range:end_range],0.1,na.rm=TRUE) #==peaks1$ID[p1-1])])
									
							}else{
							
							if((peaks1$start[p1]-peaks1$end[p1-1])<2 & peaks1$baseheight[p1]<=peaks1$baseheight[p1-1]){
								
									peaks1$ID[p1]=peaks1$ID[p1-1]
									
									start_range<-min(peaks1$start[which(peaks1$ID==peaks1$ID[p1-1])])
									end_range<-max(peaks1$end[which(peaks1$ID==peaks1$ID[p1-1])])
									median_base_height<-sd(peaks1$baseheight[which(peaks1$ID==peaks1$ID[p1])]) #max(peaks1$baseheight[which(peaks1$ID==peaks1$ID[p1-1])])-min(peaks1$baseheight[which(peaks1$ID==peaks1$ID[p1-1])]) #quantile(yvec,0.1) #[start_range:end_range],0.1,na.rm=TRUE) #==peaks1$ID[p1-1])])
									
							}
								
							}
					}
					
					
					peaks1<-as.data.frame(peaks1)
					
					#save(intvec,yvec,xvec,min.pw,max.pw,peaks,peaks1,targetmz,file=fnametemp)
				
						peaks_clust<-peaks1 %>%
								  dplyr::group_by(ID) %>%
								  dplyr::summarise(snr=(max(baseheight,na.rm=TRUE)[1]-mean(baseheight,na.rm=TRUE))/sd(baseheight,na.rm=TRUE), stdev=sd(baseheight),start = dplyr::first(start), end = dplyr::last(end),basert=basetime[which(baseheight==max(baseheight,na.rm=TRUE))],peakwidth=abs(start-end),baseheight=max(baseheight))
					
					peaks_clust<-as.data.frame(peaks_clust)
					
					peaks_clust<-cbind(peaks_clust,xvec[peaks_clust$start],xvec[peaks_clust$end])
					
					
					
					
					if(FALSE){
					
				
					peaks1 <- peaksWithCentWave(intvec[,2], intvec[,1],peakwidth=c(min.pw,max.pw))
					
					peaks1<-as.data.frame(peaks1)
					
						#peaks1 <- peaksWithCentWave(intvec[,2], intvec[,1],peakwidth=c(min.pw,max.pw))
						
			
						#peaks1 <- peaksWithCentWave(intvec[,2], intvec[,1],peakwidth=c((min(c(min.pw,peakwidth_range[1]))),max(c(peakwidth_range[2],max.pw))))
					
					peaks<-as.data.frame(peaks1)
					}
					
					#peaks2<-cbind(peaks_clust,xvec[peaks_clust$basert],xvec[peaks_clust$start],xvec[peaks_clust$end])
				
					peaks<-cbind(xvec[peaks_clust$basert],peaks_clust[,9],peaks_clust[,10],peaks_clust[,8],peaks_clust[,8],peaks_clust[,8],peaks_clust$snr)
					
					colnames(peaks)<-c("rt","rtmin","rtmax","maxo","maxo","maxo","sn")
					peaks1<-peaks
			}
			
			
	
	#rm(yvec1)
	
	
	
	
	if(length(peaks1)>0){
		
		if(nrow(peaks1)>0){
		
		if(nrow(peaks)>1){
	
			peaks<-peaks[order(peaks[,1],decreasing=TRUE),]
		}
		
		#baseheight, basert,start,end,seq,xvecsstart,xvecend,xvecbase
		#peaks<-cbind(peaks,seq(1,nrow(peaks)),xvec[peaks[,3]],xvec[peaks[,4]],xvec[peaks[,2]])
		
		peaks<-cbind(peaks[,c(6,1,2:3)],seq(1,nrow(peaks)),peaks[,c(2:3)],peaks[,c(1)],peaks[,c(4,5,7)])
		
		

		peak_groupnum=1
		peak_groups={}
		peak_groups=c(peak_groups,peak_groupnum)


		peaks<-as.data.frame(peaks)
		peaks<-peaks[which(peaks[,1]>basepeak.height.threshold),]

		
		#save(intvec,yvec,xvec,min.pw,max.pw,peaks,peaks1,targetmz,file=fnametemp)
		
		
		
		area_vec<-lapply(1:nrow(peaks),function(clustnum){
		
		lb=peaks[which(peaks[,5]==clustnum),6]
		ub=peaks[which(peaks[,5]==clustnum),7]
		i <- xvec >= lb & xvec <= ub

		peak_area=0
		zig_zag_sum=0
		window=3
		zig_zag_index=0
		mcq_index=0
		Gaussian_similarity=0
		sharpness=0
		peak_significance=0
		TPASR=0
		peak.quality.score=0
		sn.ratio<-0
		basepeak.height<-0
		basepeak.rt<-0
		
		
		
	
		
		#if(length(i)>0){
		if(length(which(i==TRUE))>0){
		
		
		
		yvectemp=replace(yvec[i],which(yvec[i]==0),NA)
		
		sum_vec<-summary(xvec[i])
		
		low_limit_time<-lb
		
		high_limit_time<-ub #sum_vec[5] #+((sum_vec[5]-sum_vec[2]))
		
		iqr_index<-which(xvec[i]>=low_limit_time & xvec[i]<=high_limit_time)
		
		prop_nonmissing<-length(which(yvec[i][iqr_index]>0))/length(yvec[i][iqr_index])
		
		length_peak_time<-diff(range(xvec[i]))
		
		if(is.na(prop_nonmissing)==TRUE){
		
			prop_nonmissing=iqr.prop.nonmissing.thresh
		}
		
		if(prop_nonmissing>=iqr.prop.nonmissing.thresh){
		
		
		
		if(length_peak_time>min.peak.length.rt & length_peak_time<=max(xvec,na.rm=TRUE)*1){

		
		yvect<-yvectemp[which(yvectemp<quantile(yvectemp,0.75,na.rm=TRUE))]
		
		if(length(yvect)<2){
				sn.ratio<-1 #(max(yvectemp,na.rm=TRUE)[1]-mean(yvect,na.rm=TRUE))
		}else{
			sn.ratio<-(max(yvectemp,na.rm=TRUE)[1]-mean(yvect,na.rm=TRUE))/(sd(yvect,na.rm=TRUE))
		}
		sn.ratio<-round(sn.ratio,2)
		
	
		
		fnametemp=paste("area",targetmz,"_",clustnum,".Rda",sep="")
		#save(xvec,yvec,lb,ub,i,sn.ratio,file=fnametemp)

		peak_area<-MESS::auc(xvec[i],yvec[i], type = 'linear')
		#peak_area<-MESS::auc(xvec[i],yvec[i], type = 'spline')
		peak_area<-round(peak_area,1)

		basepeak.height<-max(yvec[i],na.rm=TRUE)[1]
		
		if(basepeak.height==0){
			basepeak.rt<-0
		}else{
			
			basepeak.rt<-xvec[i][which(yvec[i]==basepeak.height)][1]
		
		}


		eic=yvec[i]



		
		
		if(sn.ratio>snr.threshold){

		#1. zig_zag_index
		zig_zag_sum=0

	
		
		
		epi=max(eic)-((eic[1]+eic[2]+eic[length(eic)-1]+eic[length(eic)])/4)

		 local_zig_zag<-lapply(2:(length(eic)-1),function(ii){
			  local_zig_zag=(2*eic[ii]-eic[ii-1]-eic[ii+1])^2
			  
		 })
		 local_zig_zag<-unlist(local_zig_zag)
		 zig_zag_sum=sum(local_zig_zag)
		zig_zag_index=zig_zag_sum/(epi^2*length(eic))
			
		 #2. mcq
			window=3
			eic_scale=eic/sqrt(sum(eic^2))
			
			pad_num = floor(window/2)
			eic_padding=c(rep(eic[1],pad_num),eic,rep(eic[length(eic)],window-pad_num))
			eic_window=eic

			
			eic_window<-lapply(1:length(eic),function(ii){
					mean(eic_padding[ii:(ii+window)])
			})
			eic_window<-unlist(eic_window)
			eic_window_mean  = mean(eic_window)
			eic_window_sd   = sd(eic_window)
			eic_window_standarization = (eic_window-eic_window_mean)/eic_window_sd
			eic_window_standarization_scale=eic_window_standarization/sqrt(sum(eic_window_standarization^2))  
			mcq_index=sum(eic_scale*eic_window_standarization_scale)
			
		#3.gaussian similarity
			Gaussian_similarity=0
			x<-1:length(eic)
			y<-eic
			p=try(polyfit(x,log(y),2),silent=TRUE)
			if(is(p,"try-error")){
				Gaussian_similarity=0
			}else{
					fit_sigma=suppressWarnings(sqrt(-1/(2*p[1])))
					fit_mu=p[2]*fit_sigma^2
					fit_a=exp(p[3]+fit_mu^2/(2*fit_sigma^2))
					
					Gaussian_Peak=fit_a*exp((c(1:length(eic))-fit_mu)^2/(-2*fit_sigma^2))
					Gaussian_Peak_standard=(Gaussian_Peak-mean(Gaussian_Peak))/sd(Gaussian_Peak)
					Gaussian_Peak_Scale=Gaussian_Peak_standard/sqrt(sum(Gaussian_Peak_standard^2))
					eic_standard=(eic-mean(eic))/sd(eic)
					eic_scale=eic_standard/sqrt(sum(eic_standard^2))
					Gaussian_similarity=sum(Gaussian_Peak_Scale*eic_scale)
					if(is.na(Gaussian_similarity)==TRUE){
					Gaussian_similarity=0
					}else{
					Gaussian_similarity=round(Gaussian_similarity,2)
					}
			}
		#4. 
		   sharpness=0

			apex_index=which(eic==max(eic,na.rm=TRUE))[1]
			sharpness<-lapply(1:length(eic),function(ii){	
			  if(ii < apex_index){
			    return((eic[ii+1]-eic[ii])/max(1,eic[ii]))
			  }else{
			    if(ii<length(eic)){
			     return((eic[ii]-eic[ii+1])/max(1,eic[ii+1]))
			    }
			  }
			  
			})
		      sharpness<-unlist(sharpness)
			sharpness<-sum(sharpness)
			
		#5. peak_significance
			peak_significance=0
			 apex_index=which(eic==max(eic,na.rm=TRUE))[1]
			
			if(apex_index!=1 && apex_index!=length(eic) && length(eic)>4){
			  sum1=mean(eic[c(1,2,length(eic)-1,length(eic))],na.rm=TRUE)
			  sum2=mean(c(eic[apex_index-1],max(eic)[1],eic[apex_index+1]),na.rm=TRUE)
			  peak_significance =sum2/(max(sum1,1)[1])
			}

			#6. TPASR       
			peak_level=max(eic,na.rm=TRUE)[1]
			base_level=0
			
			if(apex_index!=1 && apex_index!=length(eic) && length(eic)>4){
			  base_level=(sum(eic[c(1,2,length(eic)-1,length(eic))]))/4
			  peak_level=(eic[apex_index-1]+ max(eic) +eic[apex_index+1])/3
			}
			
			trangular_area=0.5*length(eic)*(peak_level-min(eic[1],eic[length(eic)]))
			#peak_area_1=sum(eic)-min(eic[1],eic[length(eic)])[1]*length(eic)
			TPASR=abs(trangular_area-peak_area)/trangular_area
			
			
			if(is.na(sn.ratio)==TRUE){
				sn.ratio=0
			}
			#sharpness
			peak.quality.score=(sn.ratio+3*Gaussian_similarity+peak_significance+mcq_index)-(TPASR+3*zig_zag_index)
		
			peak_area<-round(peak_area,1)
			basepeak.height<-round(basepeak.height,1)
			peak.quality.score<-round(peak.quality.score,1)
			sn.ratio<-round(sn.ratio,2)
		}	
				
		}
	  }
		if(plot.peak==TRUE){
		
				if(peak.quality.score>quality.score.threshold){
				
				
					if(is.na(basepeak.rt)==TRUE){
					
						
						basepeak.rt=mean(c(lb,ub))
					}
					
					if(is.na(deltatime)==TRUE){
						xvec<-intvec[,1] # [which(intvec[,2]>0),1]
							yvec<-intvec[,2] #[which(intvec[,2]>0),2]

							i_peak=i
							i=1:nrow(intvec)
							
							plot(xvec[i],yvec[i],xlab="Time",ylab="Intensity",cex.axis=0.7,
							main=paste("mz: ",targetmz," time range: (",round(lb,1),",",round(ub,1),") \n file: ",
							raw_fname,"(sample: ",sample_name,")\n peak area: ",peak_area," base peak height: ",
							basepeak.height,"\nSNR: ",sn.ratio," peak.quality.score:",peak.quality.score,sep=""),cex.main=0.7,type="p",cex=0.1)
							
							
							polygon(c(lb,xvec[i_peak],ub),c(0,yvec[i_peak],0),col="brown",cex=0.1)
					}else{
					
						
						
						if(is.na(deltatime)==FALSE && is.na(targettime)==FALSE){
							if(abs(targettime-basepeak.rt)<deltatime){
						
								xvec<-intvec[,1] # [which(intvec[,2]>0),1]
								yvec<-intvec[,2] #[which(intvec[,2]>0),2]
								
						
								i_peak=i
								i=1:nrow(intvec)
								plot(xvec[i],yvec[i],xlab="Time",ylab="Intensity",cex.axis=0.7,
								main=paste("mz: ",targetmz," time range: (",round(lb,1),",",round(ub,1),") \n file: ",
								raw_fname,"(sample: ",sample_name,")\n peak area: ",peak_area," base peak height: ",
								basepeak.height,"\nSNR: ",sn.ratio," peak.quality.score:",peak.quality.score,sep=""),cex.main=0.7,type="p",cex=0.1)
								
				
								
								polygon(c(lb,xvec[i_peak],ub),c(0,yvec[i_peak],0),col="brown",cex=0.1)
							}
						}else{
							xvec<-intvec[,1] # [which(intvec[,2]>0),1]
								yvec<-intvec[,2] #[which(intvec[,2]>0),2]
								
								
								i_peak=i
								i=1:nrow(intvec)
								plot(xvec[i],yvec[i],xlab="Time",ylab="Intensity",cex.axis=0.7,
								main=paste("mz: ",targetmz," time range: (",round(lb,1),",",round(ub,1),") \n file: ",
								raw_fname,"(sample: ",sample_name,")\n peak area: ",peak_area," intb: ",peaks[which(peaks[,5]==clustnum)[1],9][1],"\n base peak height: ",
								basepeak.height,"\nSNR: ",sn.ratio," peak.quality.score:",peak.quality.score,sep=""),cex.main=0.7,type="p",cex=0.1)
								
					
								
								polygon(c(lb,xvec[i_peak],ub),c(0,yvec[i_peak],0),col="brown",cex=0.1)
							
						}
					}
					
					}

				}
	}		
			
			
		return(list(peak_area=peak_area,sn.ratio=sn.ratio,basepeak.height=basepeak.height,basepeak.rt=basepeak.rt,peak.quality.score=peak.quality.score))


		})
		
		#save(intvec,yvec,min.pw,max.pw,peaks,peaks1,area_vec,file=fnametemp)
	
		
	
		#print(mem_used())
		
		if(length(peaks)>0){
		if(nrow(peaks)>0){

		if(length(area_vec)>0){
		
		
		
		peak_info<-ldply(area_vec,as.data.frame)
		

		peaks<-cbind(rep(targetmz,nrow(peaks)),peaks[,-c(5,8)],peak_info)
	
		
		
		 peaks<-as.data.frame(peaks)
		 
		 cnames1<-colnames(peaks)
		
		cnames1[1]<-"mz"
		cnames1[4]<-"min.rt"
		cnames1[5]<-"max.rt"
		
		colnames(peaks)<-cnames1
		
		# save(peaks,file="peak1.rda")

		rm(peak_info)
		 
		if(length(peaks)>0){
		 peaks<-peaks[which(peaks$sn.ratio>snr.threshold & peaks$basepeak.rt>void.time & peaks$peak.quality.score>quality.score.threshold & peaks$peak_area>0),]

		 }
		 
		if(length(peaks)>0){
			peaks$basepeak.height<-round(peaks$basepeak.height,1)
			peaks$min.rt<-round(peaks$min.rt,2)
			peaks$max.rt<-round(peaks$max.rt,2)
			peaks$basepeak.rt<-round(peaks$basepeak.rt,2)
			
			if(basepeak.select==TRUE){
			
				peaks<-peaks[1,]
			}
		  }
		}
		}
		
		#rm(area_vec)
		
		}else{
			peaks<-{}
		}
		#print(mem_used())
		return(peaks)
		}
	} 
	
		}
	}
	
}

}



}

}



match_msms_data<-function(q,query_mz_list,mz.thresh,time.thresh,msms_data,msms.file.location,msnraw.filename,DatabaseSearchName="PubChem",DatabaseSearchRelativeMassDeviation=5, FragmentPeakMatchAbsoluteMassDeviation=0.001,
FragmentPeakMatchRelativeMassDeviation=5,
plot.relative.intensity=TRUE,localcsv_file_path=NA,num_nodes=4,max_metfrag_hits=25,outloc=outloc,fileindex=NA,min.intensity.threshold=1000){

    query_mz=query_mz_list[q,1] #round(query_mz_list[q,1],4)
    
    setwd(outloc)

    log_fname=msnraw.filename
    
    res_fname<-paste("tmpres",query_mz,"_",query_mz_list[q,4],"_",fileindex,".Rda",sep="")
	
    #res_fname<-paste("tmpres",query_mz,"_",fileindex,".Rda",sep="")

    check_fname<-try(load(res_fname),silent=TRUE)
 
   
    
    plot_res_list<-new("list")
    if(is(check_fname,"try-error")){
    
    mzmin=query_mz-(10^(-6)*query_mz*mz.thresh)
    mzmax=query_mz+(10^(-6)*query_mz*mz.thresh)

    resloc=msms.file.location
    
    msms_data_mztime<-msms_data[,c("mz","rt")]
    colnames(msms_data_mztime)<-c("mz","time")
    

    
    matching_parent_scans<-getVenn(dataB=msms_data_mztime,name_b="msms", dataA=query_mz_list[q,c(1:2)],
    name_a="querydata",mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=NA, xMSanalyzer.outloc=outloc,use.unique.mz=FALSE,plotvenn=FALSE,nearest.time.match=FALSE)
    

    
    if(nrow(matching_parent_scans$common)>0)
    {
    
        
        matching_msms_data<-as.data.frame(msms_data[matching_parent_scans$common$index.B,])
        
        fname<-paste("matching_msms_data",query_mz,".Rda",sep="")
        
        ms1_data<-matching_msms_data[which(matching_msms_data$msLevel==1),]
        
        if(nrow(ms1_data)>0)
        {
	
	
	
	ms1_data<-ms1_data[order(ms1_data$rt),]
	
	yvec1<-ms1_data$intensity
	xvec<-ms1_data$rt
	
	save(ms1_data,file="ms1_data.Rda")
	

		ms1_data_clust<-cbind(ms1_data,seq(1,nrow(ms1_data)))
		
	
	
	if(length(ms1_data_clust)>0){
		
		if(nrow(ms1_data_clust)>=1){
		
		
		ms1_data_clust<-as.data.frame(ms1_data_clust)
		colnames(ms1_data_clust)<-c(colnames(ms1_data),"clusterlabels")
		ms1_data_clust<-as.data.frame(ms1_data_clust)
		
	    
	   setwd(outloc)
	   clusterlevels=seq(1,nrow(ms1_data_clust))
	
            plot_res_list<-lapply(1:nrow(ms1_data_clust),function(snum_val)
            {
                mainlab1<-paste("Peak clusters based on time and intensity for mz:",round(query_mz,5),"",sep=" ")
                
              
                
                select_int_scan<-ms1_data_clust[snum_val,]
                
                ind_peakid<-which(msms_data$MSnParentPeakID%in%select_int_scan$peakID)
                
                max_int_scan<-select_int_scan
                
		
                rda_fname<-paste(query_mz,"_",ind_peakid[1],log_fname,"clust",snum_val,"matching.Rda",sep="")
              # save(matching_ms1_data,file=rda_fname)
                
                rda_fname<-paste(query_mz,"_",ind_peakid[1],log_fname,"clust",snum_val,"ms1clust.Rda",sep="")
             # save(ms1_data_clust,file=rda_fname)
                
                
               	# save(msms_data,file="msms_data.Rda")
                
                
                
                if(length(ind_peakid)>0){
                    
                    msms_max_int_scan<-msms_data[ind_peakid,]
                    
                    msms_max_int_scan<-msms_max_int_scan[order(msms_max_int_scan$intensity,decreasing=TRUE),]
                    
                    rda_fname<-paste(query_mz,"_",ind_peakid[1],log_fname,"clust",snum_val,".Rda",sep="")
                    # save(msms_max_int_scan,file=rda_fname)
                    if(length(which(duplicated(round(msms_max_int_scan$mz,1))==TRUE))>0){
                        msms_max_int_scan<-msms_max_int_scan[-which(duplicated(round(msms_max_int_scan$mz,1))==TRUE),]
                    }
                    
                    if(is.na(min.intensity.threshold)==FALSE){
		    
			msms_max_int_scan<-msms_max_int_scan[which(msms_max_int_scan$intensity>min.intensity.threshold),]
			
		    }
		    
		     if(nrow(msms_max_int_scan)>0){
		     
                    data_m<-msms_max_int_scan[,c("mz","intensity")]
		    
		    
                    
                    scan_num<-paste("parentPeakID",msms_max_int_scan$MSnParentPeakID[1],sep="")
                    
                    #pdffname<-paste(query_mz,"_",scan_num,log_fname,"spectra.pdf",sep="")
                    
                    
                      setwd(outloc)
                        
                   
                        data_m<-data_m[order(data_m[,2],decreasing=TRUE),]
                        data_m[,2]<-100*(data_m[,2]/max(data_m[,2]))
                        data_m<-data_m[which(data_m[,2]>1),]
                        name<-rep("U",length(data_m[,1]))
                        filename<-rep("U1",length(data_m[,1]))
                        display<-rep("FALSE",length(data_m[,1]))
                        identity<-rep("NA",length(data_m[,1]))
                        
                        data_m2<-cbind(data_m[,1],data_m[,2])
                        colnames(data_m2)<-c("mz","intensity")
                        
                        setwd(outloc)
                        dir.create("Metlin")
                        setwd("Metlin")
                        fname<-paste(query_mz,"_",scan_num,log_fname,".csv",sep="")
                        
                        write.table(data_m2,file=fname,sep=",",row.names=FALSE)
                        
                        setwd(outloc)
                        
                        dir.create("MetFrag")
                        setwd("MetFrag")
                        fname<-paste(query_mz,"_",scan_num,log_fname,".txt",sep="")
                        
                        write.table(data_m2,file=fname,sep="\t",row.names=FALSE)
                        
                        setwd(outloc)
                        data_m2<-cbind(name,filename,data_m[,1],data_m[,2],display,identity)
                        colnames(data_m2)<-c("compound","filename","mz","intensity","display","identity")
                        
                        sample<-rep("U",length(data_m[,1]))
                        instrument<-rep("Fusion",length(data_m[,1]))
                        category1<-rep("U1",length(data_m[,1]))
                        category2<-rep("U2",length(data_m[,1]))
                        category3<-rep("U3",length(data_m[,1]))
                        rt1d<-rep("U4",length(data_m[,1]))
                        rt2d<-rep("U5",length(data_m[,1]))
                        formula<-rep("U6",length(data_m[,1]))
                        comment<-rep("NA",length(data_m[,1]))
                        
                        
                        meta<-cbind(name,filename,sample,instrument,category1,category2,category3,rt1d,rt2d,formula,comment)
                        
                        colnames(meta)<-c("compound","filename","sample","instrument","category.1","category.2","category.3","rt.1d","rt.2d","formula","comment")
                        
                        fname<-paste("xcms/",query_mz,"_",scan_num,log_fname,"meta1.csv",sep="")
                        
                        data_m2<-as.data.frame(data_m2)
                        meta<-as.data.frame(meta)
                        
                        #write.table(data_m2,file=fname,sep=";",row.names=FALSE)
                        
                        fname<-paste("xcms/",query_mz,"_",scan_num,log_fname,"meta2.csv",sep="")
                        #write.table(meta,file=fname,sep=";",row.names=FALSE)
                        
                        setwd(outloc)
                        
                        
                        dir.create("mzCloud")
                        
                        setwd("mzCloud")
                        
                        fname<-paste(query_mz,"_",scan_num,log_fname,".msp",sep="")
                        WriteMspFile(spectra=data_m2,metadata=meta,filename=fname)
                        
                        setwd(outloc)
                        
                        ####code for searching in MetFrag
                         fname<-paste(query_mz,"_",scan_num,log_fname,"debugmetfrag.Rda",sep="")
			#save(data_m,DatabaseSearchName,DatabaseSearchRelativeMassDeviation,FragmentPeakMatchAbsoluteMassDeviation,FragmentPeakMatchRelativeMassDeviation,query_mz_list,localcsv_file_path,file=fname)
                       # print(head(data_m))
			
                        DatabaseName<-DatabaseSearchName
                        
                        settingsObject<-list()
                        #
                        # set database parameters to select candidates
                        #
		
                        settingsObject[["DatabaseSearchRelativeMassDeviation"]]<-DatabaseSearchRelativeMassDeviation
                        settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]]<-FragmentPeakMatchAbsoluteMassDeviation
                        settingsObject[["FragmentPeakMatchRelativeMassDeviation"]]<-FragmentPeakMatchRelativeMassDeviation
                        settingsObject[["MetFragDatabaseType"]]<-DatabaseName
                        
                        if(dim(query_mz_list)[2]>2){
                            settingsObject[["NeutralPrecursorMass"]]<-as.numeric(as.character(query_mz_list[q,3]))
                            settingsObject[["NeutralPrecursorMolecularFormula"]]<-as.character(query_mz_list[q,4]) #"C23H45N1O4"
                        }
                        
                        if(DatabaseName=="LocalCSV"){
                            settingsObject[["MetFragDatabaseType"]]<-"LocalCSV"
                            settingsObject[["LocalDatabasePath"]]<-localcsv_file_path
                        }
                        
                        
                        data_m[,1]<-round(data_m[,1],4)
                        
                        settingsObject[["PeakList"]]<-cbind(data_m[,1],data_m[,2])
			options(digits=10)
			scored.candidates<-try(run.metfrag(settingsObject),silent=TRUE)
                       
                        rm(settingsObject)
			
			#print(head(scored.candidates))
                            
                            scan_num<-msms_max_int_scan$MSnParentPeakID[1]
                            maintext=paste("MS/MS spectra for selected precursor m/z: ",round(select_int_scan$mz[1],5)," time: ",round(select_int_scan$rt[1],2)," (cluster: ",clusterlevels[snum_val]," peakID: ",scan_num,")",sep="")
                           #  print("here")
                            if(plot.relative.intensity==TRUE){
                                msms_max_int_scan$intensity<-100*(msms_max_int_scan$intensity/max(msms_max_int_scan$intensity,na.rm=TRUE))
                                ms1_data_clust$intensity<-100*(ms1_data_clust$intensity)/max(ms1_data_clust$intensity,na.rm=TRUE)
                                ylab1="Relative Intensity (%)"
                                yinc=3.5
                            }else{
                                ylab="Intensity"
                                maxint<-max(msms_max_int_scan$intensity,na.rm=TRUE)
                                minint<-min(msms_max_int_scan$intensity,na.rm=TRUE)
                                if(abs(maxint-minint)>10000){
                                    yinc=median(msms_max_int_scan$intensity)
                                }else{
                                    yinc=min(100,abs(maxint-minint))
                                }
                            }
			    
			 #   print(head(msms_max_int_scan))
			    
                            mytheme <- gridExtra::ttheme_default(
                            core = list(fg_params=list(cex = 0.5,hjust=0,x=0)),
                            colhead = list(fg_params=list(cex = 0.55,hjust=0,x=0)),
                            rowhead = list(fg_params=list(cex = 0.5,hjust=0,x=0)))
                           # save(ms1_data,ms1_data_clust,msms_max_int_scan,mainlab1,clusterlabels,ylab1,maintext,msnraw.filename,file="debug_2.Rda")
                            #maintext<-c(maintext,"\n",top_metfrag_hits_string,sep="")
			    
			    
                            x=msms_max_int_scan$mz
                            y=msms_max_int_scan$intensity
                            x=round(x,4)
                            y=round(y,2)
                            tmp1<-cbind(x,y)
                            tmp1<-as.data.frame(tmp1)
                         #   save(tmp1,file="tmp1.Rda")
			    mzmin=query_mz-(10^(-6)*query_mz*DatabaseSearchRelativeMassDeviation)
		            mzmax=query_mz+(10^(-6)*query_mz*DatabaseSearchRelativeMassDeviation)
			    setwd(msms.file.location)
			    xraw<-xcmsRaw(msnraw.filename,includeMSn=TRUE)
			    intvec<-getEIC(xraw,mzrange=c(mzmin,mzmax),step=0.001) #$intensity
			    
			    intvec=intvec@eic[[1]][[1]]
			    
			    rm(xraw)
			#print(head(intvec))
			intvec=as.data.frame(intvec)
			mainlab<-paste("EIC for mz:",round(query_mz,5),"(file:",msnraw.filename,")",sep=" ")
			if(plot.relative.intensity==TRUE){			
				intvec$intensity<-100*(intvec$intensity)/max(intvec$intensity,na.rm=TRUE)
				ylab1="Relative Intensity(%)"
			}else{
				ylab1="Intensity"
			}
			intvec=cbind(intvec$rt,intvec$intensity)
			intvec=as.data.frame(intvec)
			colnames(intvec)<-c("rt","intensity")
		      
	               	plots_list=ggplot(data=intvec,aes(x=rt, y=intensity)) + labs(x="time",y=ylab1) + ggtitle(mainlab) +geom_area(aes()) + theme(axis.text=element_text(size=7),
                                                        axis.title=element_text(size=7), plot.title=element_text(size=7))
                            p1=ggplot(data=ms1_data_clust,aes(x=rt, y=intensity)) + labs(x="time",y=ylab1) + ggtitle(mainlab1) +geom_area(aes(fill=clusterlabels,group=clusterlabels)) + geom_text(data=ms1_data_clust,aes(x=rt, y=intensity,label=clusterlabels),size=2) + theme(axis.text=element_text(size=7),
                            axis.title=element_text(size=7), plot.title=element_text(size=7))
                            
                            p2=ggplot(data=tmp1,aes(x,y)) + labs(x="m/z",y=ylab1)  + ggtitle(maintext) + geom_segment(x=x,xend=x,y=0,yend=y,data=tmp1) + geom_text(data=tmp1,aes(x=x,y=y,label=x),size=2) + theme(axis.text=element_text(size=7),
                            axis.title=element_text(size=7), plot.title=element_text(size=7))
 
                            myt <-{} #gridExtra::tableGrob(top_metfrag_hits_string, theme = mytheme)
                            
                            tmp2<-new("list")
                           top_metfrag_hits_string<-{}
			 if(is(scored.candidates,"try-error")){
                                tmp2=NULL
                        }else{

			
                        scored.candidates<-as.data.frame(scored.candidates)

				
                       
                         fname<-paste(query_mz,"_",scan_num,log_fname,"scoredmetfrag.Rda",sep="")
                        #save(scored.candidates,file=fname)

			     if(nrow(scored.candidates)>0){
			     
			    

				scored.candidates$InChIKey<-as.character(scored.candidates$InChIKey)
                           # data(hmdbAllinf)
                            scored.candidates<-scored.candidates[order(scored.candidates$Score,decreasing=TRUE),]
                           # hmdbAllinf$InChiKey=gsub(hmdbAllinf$InChiKey,pattern="InChIKey=",replacement="")
                           # hmdbAllinf$InChiKey<-as.character(hmdbAllinf$InChiKey)
                           # hmdbAllinf<-hmdbAllinf[,c("HMDBID","Name","InChiKey","KEGGID","BioFluidLocation","CellularLocation","Origin","BioFunctions","HMDBStatus")]
                            m1<-scored.candidates #merge(scored.candidates,hmdbAllinf,by.x="InChIKey",by.y="InChiKey")
                            rm(scored.candidates)
			    m1<-m1[order(m1$Score,decreasing=TRUE),]
                     
			       top_scores<-m1[which(m1$Score>0),]
                            dup_inchi_check<-which(duplicated(top_scores$InChIKey)==TRUE)
                            if(length(dup_inchi_check)>0){
                                top_scores<-top_scores[-dup_inchi_check,]
                            }
                            top_scores<-top_scores[order(top_scores$Score,decreasing=TRUE),]
                            if(nrow(top_scores)>max_metfrag_hits){
                                top_scores<-top_scores[1:max_metfrag_hits,]
                            }
                            top_metfrag_hits_string<-{}
                            for(mrownum in 1:nrow(top_scores)){
                                expl_peaks<-gsub(top_scores$ExplPeaks[mrownum],pattern="_[0-9.]*",replacement="")
                                top_metfrag_hits_string<-rbind(top_metfrag_hits_string,cbind(as.character(top_scores$Identifier[mrownum]),as.character(top_scores$Name[mrownum]),round(top_scores$Score[mrownum],2),round(top_scores$FragmenterScore[mrownum],2),expl_peaks))
                            }
                            top_metfrag_hits_string<-as.data.frame(top_metfrag_hits_string)
                            colnames(top_metfrag_hits_string)<-c("ID", "Name", "Score (0 to1)","FragmenterScore","Explained Peaks")
                            rownames(top_metfrag_hits_string)<-NULL
                            setwd(outloc)
                            setwd("MetFrag")
                            fname<-paste(query_mz,"_",scan_num,log_fname,"candidates.csv",sep="")
                            write.csv(m1,fname,row.names=FALSE)
                            #fname<-paste(query_mz,"_",scan_num,log_fname,"candidatesHMDB.csv",sep="")
                            #write.csv(m1,fname,row.names=FALSE)
 

           		     myt <- gridExtra::tableGrob(top_metfrag_hits_string, theme = mytheme)
 
                           top_metfrag_hits_string<-cbind(top_metfrag_hits_string,query_mz,scan_num,log_fname)
                            top_metfrag_hits_string<-as.data.frame(top_metfrag_hits_string)
                            colnames(top_metfrag_hits_string)<-c("ID", "Name", "Score (0 to1)","FragmenterScore","Explained Peaks","query_mz","ScanPeakID","FileName")
                            rownames(top_metfrag_hits_string)<-NULL
			    

			    }                            
			}                            
                            tmp2<-list(p1=p1,p2=p2,myt=myt,top_metfrag_hits_string=top_metfrag_hits_string,p3=plots_list)
                            fname<-paste(query_mz,"tmp2.Rda",sep="")
                            #save(tmp2,file=fname)
			  
			  rm(p1)
			  rm(p2)
			  rm(myt)
			  rm(top_metfrag_hits)
			  rm(plots_list)
			  rm(intvec)
			  rm(data_m2)
			  rm(matching_ms1_data)
			  rm(tmp1)
			  rm(m1)
			  rm(meta)
			  
			  setwd(outloc)

			   return(tmp2)
                                
                             
                            
                        }
                        
                    
                    #dev.off()
                    
                    setwd(resloc)
                    
                }
                
                
            }) #end lapply
	    
	    rm(matching_msms_data)
            rm(msms_data)
	    rm(ms1_data_clust)
            #plot_res[[q]]<-plot_res_list
	    
	#    fname<-paste("tmpres",query_mz,"_",log_fname,".Rda",sep="")
	   # fname<-paste("tmpres",query_mz,"_",query_mz_list[q,4],"_",log_fname,".Rda",sep="")
	   
	    fname<-paste("tmpres",query_mz,"_",query_mz_list[q,4],"_",fileindex,".Rda",sep="")
	    save(plot_res_list,file=fname)
            return(plot_res_list)
	    
	    
	    }
	    }
        }
        
    }
    
    }else{
        #print("m/z not found")
	load(res_fname)
	return(plot_res_list)
    }
    
}

get.EIC.target<-function(MS_data_directory,input_mz_list_fname, ppm_mz_thresh=10, sec_time_thresh=60, scans_per_second=1.3,
method.process="xcms",
rawfile.conversion.mode="centroid",
plot.relative.intensity=TRUE,ms.inputfile.format=".raw",num_nodes=4,outloc=NA){

setwd(MS_data_directory)
conversion_mode<-rawfile.conversion.mode
resloc<-MS_data_directory
setwd(resloc)
query_mz_list<-read.table(input_mz_list_fname,sep="\t",header=TRUE)

mz.thresh=ppm_mz_thresh
time.thresh=sec_time_thresh

        dir.create(outloc)
        ms.inputfile.format<-tolower(ms.inputfile.format)
        if(ms.inputfile.format==".raw"){
                        dirfiles<-list.files(".",".raw$")
                        dirfiles<-paste(resloc,"",dirfiles,sep="")
                        dirfiles<-gsub(dirfiles,pattern="\\\\",replacement="/")
                        write.table(dirfiles,file="Rawfiles.txt",sep="\t",row.names=FALSE,quote = FALSE,col.names=FALSE)
                        if(conversion_mode=="centroid"){
                        system2("C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.9098\\msconvert.exe", args=c("-f  Rawfiles.txt","--64","--zlib","--mzXML","--filter", "\"peakPicking true 1-2\"","--filter", "\"msLevel 1-2\""))
                        }else{
                                if(conversion_mode=="profile"){
                                system2("C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.9098\\msconvert.exe", args=c("-f  Rawfiles.txt","--64","--zlib","--mzXML","--filter", "\"msLevel 1-\""))
                                }
                        }
        }
        dirfiles<-list.files(".",".mzXML$")
        dir.create(outloc)
        setwd(outloc)
        cl <- parallel::makeCluster(num_nodes)
           clusterExport(cl, "getVenn")
          clusterExport(cl, "WriteMspFile")
         clusterExport(cl, "run.metfrag")
          clusterExport(cl, "rawEIC")
          clusterExport(cl, "find_peaks")
          clusterExport(cl, "find.Overlapping.mzs")
	      clusterExport(cl, "find.Overlapping")
	       clusterExport(cl, "find.overlapping.single")
          clusterExport(cl, "vennCounts")
           clusterEvalQ(cl, library(xMSannotator))
         clusterEvalQ(cl, library(xcms))
          clusterEvalQ(cl, library(ggplot2))
         clusterEvalQ(cl, library(gridExtra))
         process_raw_xcms<-function(file.num,dirfiles,resloc){
                        setwd(resloc)
                        msnraw.filename=dirfiles[file.num]
                        xraw<-xcmsRaw(msnraw.filename)
                        return(xraw)
         }
         generate_eic_plot<-function(q,query_mz_list,mz.thresh,xraw,dirfiles){
                query_mz=query_mz_list[q,1]
                        mzmin=query_mz-(10^(-6)*query_mz*mz.thresh)
                        mzmax=query_mz+(10^(-6)*query_mz*mz.thresh)
                
			
			plots_list<-lapply(1:length(dirfiles),function(j){
			
			msnraw.filename=dirfiles[j]
			intvec<-rawEIC(xraw[[j]] ,mzrange=c(mzmin,mzmax)) #$intensity
                        mainlab1<-paste("EIC for mz:",round(query_mz,5),"Name:",query_mz_list[q,2],"\n(file:",msnraw.filename,")",sep=" ")
                        #intvec$intensity<-log10(intvec$intensity+1)

			if(plot.relative.intensity==TRUE){			
				intvec$intensity<-100*(intvec$intensity)/max(intvec$intensity)
				ylab1="Relative Intensity(%)"
			}else{
				ylab1="Intensity"
			}
                      
			ms1_data_clust=cbind(intvec$scan,intvec$intensity)
			ms1_data_clust=as.data.frame(ms1_data_clust)
			colnames(ms1_data_clust)<-c("rt","intensity")
		        #plots_list=qplot(intvec$scan, intvec$intensity,xlab="Scan",ylab=ylab,main=mainlab1,geom=c("line"))
	               	plots_list=ggplot(data=ms1_data_clust,aes(x=rt, y=intensity)) + labs(x="time",y=ylab1) + ggtitle(mainlab1) +geom_area(aes()) + theme(axis.text=element_text(size=8),
                                                        axis.title=element_text(size=8), plot.title=element_text(size=8))


			 return(plots_list)
			})
			return(plots_list)
        }
         xraw<-parLapply(cl,1:length(dirfiles),process_raw_xcms,dirfiles=dirfiles,resloc=resloc)
        plot_res<-parLapply(cl,1:dim(query_mz_list)[1],generate_eic_plot,query_mz_list,mz.thresh,xraw,dirfiles=dirfiles)
        #do.call(grid.arrange,p)
        save(plot_res,file="plot_res_ms1.Rda")
        stopCluster(cl)
	pdffname<-"EIC_results.pdf"
                                                        pdf(pdffname,onefile=TRUE)
                                                        par(mfrow=c(3,1))

        for(j in 1:length(plot_res)){

                        if(length(plot_res[[j]])>0){
                                        #try(do.call(grid.arrange,plot_res[[j]]$p1),silent=TRUE)
        					
	                		#do.call(grid.arrange,plot_res[[j]])
					lapply(1:length(plot_res[[j]]),function(k){
						plot(plot_res[[j]][[k]])
					#	do.call(grid.arrange,plot_res[[j]][[k]])
					})
			}
        }
        dev.off()
}



find_peaks <- function (x, m = 1){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1)
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
     pks <- unlist(pks)
     pks
}

findPeaks1 <-
function(x, thresh=0) {
  pks <- which(diff(sign(diff(x, na.pad=FALSE)),na.pad=FALSE) < 0) + 1
  if( !missing(thresh) ) {
    if(sign(thresh) < 0)
      thresh <- -thresh
    pks[x[pks-1]-coredata(x[pks]) > thresh]
  } else pks
}


findValleys1 <-
function(x, thresh=0) {

  
  pks <- c(which(diff(sign(diff(x, na.pad=FALSE)),na.pad=FALSE) > 0) + 1,which(x==min(x)))
  pks<-unique(pks)
  if( !missing(thresh) ) {
    if(sign(thresh) > 0)
      thresh <- -thresh
    pks[x[pks-1]-coredata(x[pks]) < thresh]
  } else pks
}



#rawfile.conversion.mode: centroid, profile, NA
#rawfile.conversion.mode: centroid, profile, NA
raw.file.convert<-function(MSn_data_directory,rawfile.conversion.mode="centroid",proteowizardversion="3.0.9098",numnodes=2,overwrite=TRUE)
{

		setwd(MSn_data_directory)


		conversion_mode<-rawfile.conversion.mode


		
		
		dirfiles<-list.files(".",".raw$")


		
		if(overwrite==FALSE){
			l2<-list.files(".",".mzXML$")

			lcom<-which(gsub(dirfiles,pattern=".raw|.RAW",replacement="")%in%gsub(l2,pattern=".mzXML",replacement=""))

		
			if(length(lcom)>0){
				print("The following files are already converted:")
				print(dirfiles[lcom])
				dirfiles<-dirfiles[-lcom]
			}
			
		}
		dirfiles<-paste(MSn_data_directory,"\\",dirfiles,sep="")
		
		dirfiles<-gsub(dirfiles,pattern="\\\\",replacement="/")
		
		max_num_file_per_list<-2
		
		list_num<-seq(1,length(dirfiles),max_num_file_per_list)

		list_index<-new("list")
		
		list_index<-lapply(1:length(list_num),function(k){
				
					return(seq(list_num[k],(max_num_file_per_list+list_num[k]-1)))
		})
		
		
		convert_files<-function(list_num,conversion_mode,dirfiles,list_index,proteowizardversion)
		{
					
					temp_dirfiles<-dirfiles[list_index[[list_num]]]
					
					temp_fname<-paste("Rawfiles",list_num,".txt",sep="")
					write.table(temp_dirfiles,file=temp_fname,sep="\t",row.names=FALSE,quote = FALSE,col.names=FALSE)
					
					
					
					if(conversion_mode=="centroid"){
					
					system2(paste("C:\\Program Files\\ProteoWizard\\ProteoWizard ",proteowizardversion,"\\msconvert.exe",sep=""), args=c(paste("-f  ",temp_fname,sep=""),"--64","--zlib","--mzXML","--filter", "\"peakPicking true 1-2\"","--filter", "\"msLevel 1-2\""))
					
					#"C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.9098\\msconvert.exe"
					}else{
						if(conversion_mode=="profile"){
						#system2(paste("C:\\Program Files\\ProteoWizard\\ProteoWizard ",proteowizardversion,"\\msconvert.exe",sep=""), args=c("-f  Rawfiles.txt","--64","--zlib","--mzXML","--filter", "\"msLevel 1-\""))
						system2(paste("C:\\Program Files\\ProteoWizard\\ProteoWizard ",proteowizardversion,"\\msconvert.exe",sep=""), args=c(paste("-f  ",temp_fname,sep=""),"--64","--zlib","--mzXML","--filter", "\"msLevel 1-\""))
						
						}
					}
		}
		cl<-parallel::makeCluster(numnodes)
		clusterEvalQ(cl,"convert_files")
		res<-parLapply(cl,1:length(list_num),convert_files,conversion_mode,dirfiles,list_index,proteowizardversion)
		
		stopCluster(cl)

}

MSMS.process<-function(MSMS_data_directory,input_mz_list_fname, ppm_mz_thresh=10, sec_time_thresh=60, scans_per_second=1.3, 
method.process="xcms",
rawfile.conversion.mode="centroid",
DatabaseSearchName="PubChem",DatabaseSearchRelativeMassDeviation=5, FragmentPeakMatchAbsoluteMassDeviation=0.001,
FragmentPeakMatchRelativeMassDeviation=5,
plot.relative.intensity=TRUE,localcsv_file_path=NA,msms.inputfile.format=".mzXML",num_nodes=4,max_metfrag_hits=25,xcms.method="centWave",xcms.ppm=2.5,
xcms.snthresh=1,outloc=NA,xcms.peakwidth=c(5,20),xcms.mzdiff=0.001,min.intensity.threshold=1000,xcms.prefilter=c(3,100),proteowizard_version="ProteoWizard 3.0.18264.d4d63fbdd",process.raw.files=TRUE){

setwd(MSMS_data_directory)


conversion_mode<-rawfile.conversion.mode

resloc<-MSMS_data_directory
setwd(resloc)

query_mz_list<-read.table(input_mz_list_fname,sep="\t",header=TRUE)

method.process<-tolower(method.process)

#query_mz_list<-query_mz_list[1:3,]
mz.thresh=ppm_mz_thresh
time.thresh=sec_time_thresh
if(method.process=="xcms"){

	
	msms.inputfile.format<-tolower(msms.inputfile.format)
	
	if(msms.inputfile.format==".raw"){
	
			dirfiles<-list.files(".",".raw$")
			
			print(dirfiles)
			
			dirfiles<-paste(resloc,"",dirfiles,sep="")
			
			dirfiles<-gsub(dirfiles,pattern="\\\\",replacement="/")

			
			write.table(dirfiles,file="Rawfiles.txt",sep="\t",row.names=FALSE,quote = FALSE,col.names=FALSE)
			
			sys_cmd<-paste("C:\\Program Files\\ProteoWizard\\",proteowizard_version,"\\msconvert.exe", sep="")
			
			if(conversion_mode=="centroid"){
			
			system2(sys_cmd, args=c("-f  Rawfiles.txt","--64","--zlib","--mzXML","--filter", "\"peakPicking true 1-2\"","--filter", "\"msLevel 1-2\""))
			
			
			}else{
				if(conversion_mode=="profile"){
				system2(sys_cmd, args=c("-f  Rawfiles.txt","--64","--zlib","--mzXML","--filter", "\"msLevel 1-\""))
				}
			}

	}
	
	dirfiles<-list.files(".",".mzXML$") 
	#round(0.5*parallel::detectCores())
	snowparam <- SnowParam(workers =num_nodes, type = "SOCK")
	
	if(xcms.method=="MS1"){
	
	
	xs <- xcmsSet(dirfiles, method = "MS1",BPPARAM = snowparam)
	}else{
		if(xcms.method=="matchedFilter"){
			
			xs <- xcmsSet(dirfiles, method = "matchedFilter",step=0.01,mzdiff=xcms.mzdiff,snthresh=xcms.snthresh,max=50,BPPARAM = snowparam)
			
		}else{
			#xs <- xcmsSet(msnraw.filename, method = xcms.method, ppm=xcms.ppm, peakwidth=xcms.peakwidth)
			
			xs<-xcmsSet(dirfiles, method="centWave", ppm=xcms.ppm, snthresh=xcms.snthresh,mzdiff=xcms.mzdiff,peakwidth=xcms.peakwidth,prefilter=xcms.prefilter,
					integrate=1, verbose.columns=TRUE,fitgauss=FALSE, BPPARAM = snowparam)
		}
	}
	
	xfrag <- xcmsFragments(xs,snthresh=xcms.snthresh)

	#xfrag@peaks<-rbind(xfrag@peaks,xfrag2@peaks)

 
        plot_res<-lapply(1:length(dirfiles),function(file.num)
       {
					setwd(resloc)
					msnraw.filename=dirfiles[file.num]
					
					print(msnraw.filename)
					
					
					#xraw<-xcmsRaw(msnraw.filename,includeMSn=TRUE)

					#peaks <- findPeaks(xraw, method=xcms.method) #method = "MS1",
					

					
					msms_data<-xfrag@peaks[which(xfrag@peaks[,7]==file.num),]
					
					msms_data<-as.data.frame(msms_data)
					
					
					msms_data$mz<-round(msms_data$mz,4)
					msms_data$rt<-round(msms_data$rt,1)
					
					msms_data<-as.data.frame(msms_data)
					msms_data_mztime<-msms_data[,c("mz","rt")]
					dir.create(outloc)
					setwd(outloc)
					fname<-paste("scans",xcms.method,"ppm",ppm_mz_thresh,"pw",xcms.peakwidth[1],"_",xcms.peakwidth[2],"sn",xcms.snthresh,msnraw.filename,".csv",sep="")
					write.csv(msms_data,file=fname,row.names=FALSE)
					
					colnames(msms_data_mztime)<-c("mz","time")
					
					log_fname=msnraw.filename
				
					
					
					cl <- parallel::makeCluster(num_nodes)
					   clusterExport(cl, "getVenn")
					  clusterExport(cl, "WriteMspFile")  
					 clusterExport(cl, "run.metfrag")
					  clusterExport(cl, "rawEIC")
					  clusterExport(cl, "find_peaks")
					  clusterExport(cl, "find.Overlapping.mzs")
					   clusterExport(cl, "find.overlapping.single")
					   clusterExport(cl, "find.Overlapping")
					  clusterExport(cl, "vennCounts")
					   clusterEvalQ(cl, library(xMSannotator))
					   clusterExport(cl, "find.peak.regions")
					 clusterEvalQ(cl, library(xcms))
					  clusterEvalQ(cl, library(ggplot2))
					 clusterEvalQ(cl, library(gridExtra))
					 

					 
				if(FALSE)
				{
				#for(qindval in 1:dim(query_mz_list)[1])
					for(qindval in 1:2)
					{

						
						print(qindval)
						print(query_mz_list[qindval,])
						plot_res<-match_msms_data(qindval,query_mz_list=query_mz_list,mz.thresh=mz.thresh,time.thresh=time.thresh,msms_data=msms_data,msms.file.location=resloc,msnraw.filename=msnraw.filename,DatabaseSearchName=DatabaseSearchName,DatabaseSearchRelativeMassDeviation=DatabaseSearchRelativeMassDeviation, FragmentPeakMatchAbsoluteMassDeviation=FragmentPeakMatchAbsoluteMassDeviation,
						FragmentPeakMatchRelativeMassDeviation=FragmentPeakMatchRelativeMassDeviation,
						plot.relative.intensity=plot.relative.intensity,localcsv_file_path=localcsv_file_path,num_nodes=num_nodes,max_metfrag_hits=max_metfrag_hits,outloc=outloc,fileindex=file.num,min.intensity.threshold=min.intensity.threshold)
					}
				}

					   
				#if(FALSE)
				{
					plot_res<-parLapply(cl,1:dim(query_mz_list)[1],match_msms_data,query_mz_list=query_mz_list,mz.thresh=mz.thresh,time.thresh=time.thresh,msms_data=msms_data,msms.file.location=resloc,msnraw.filename=msnraw.filename,DatabaseSearchName=DatabaseSearchName,DatabaseSearchRelativeMassDeviation=DatabaseSearchRelativeMassDeviation, FragmentPeakMatchAbsoluteMassDeviation=FragmentPeakMatchAbsoluteMassDeviation,
					FragmentPeakMatchRelativeMassDeviation=FragmentPeakMatchRelativeMassDeviation,
					plot.relative.intensity=plot.relative.intensity,localcsv_file_path=localcsv_file_path,num_nodes=num_nodes,max_metfrag_hits=max_metfrag_hits,outloc=outloc,fileindex=file.num,min.intensity.threshold=min.intensity.threshold)
					    #end parLapply
				}	
				stopCluster(cl)
				plot_res<-{}
			
					    
	    
	    })
	    


    setwd(outloc)
	
	
	tmpres_files<-list.files(".","tmpres")
	
	pdffname<-"MSMS_results.pdf"
							pdf(pdffname,onefile=TRUE)
							par(mfrow=c(4,1))
							
							tt3 <- ttheme_minimal(
							  core=list(bg_params = list(fill = blues9[1:4], col=NA),
								    fg_params=list(fontface=3)),
							  colhead=list(fg_params=list(col="navyblue", fontface=4L)),
							  rowhead=list(fg_params=list(col="orange", fontface=3L)))

                            tmpmat1<-{}
	for(find1 in 1:length(tmpres_files))
	{
	
			check_fname<-try(load(tmpres_files[find1]),silent=TRUE)
			if(is(check_fname,"try-error")){
				next
			}else{
					plot_res<-plot_res_list
					for(j in 1:length(plot_res)){
					
						
						
							if(length(plot_res[[j]])>0){
						
							
						    check_metfrag_hits<-try(plot_res[[j]]$top_metfrag_hits_string,silent=TRUE)
						    
						    if(is(check_metfrag_hits,"try-error")){
							check_metfrag_hits<-{}
						    }else{
						    
								 try(grid.arrange(plot_res[[j]]$p3,plot_res[[j]]$p2,plot_res[[j]]$myt,layout_matrix = rbind(c(1),c(2),c(3),c(4))),silent=TRUE)
							    if(length(plot_res[[j]]$top_metfrag_hits_string)>0){
								
								
								if(is.na(plot_res[[j]]$top_metfrag_hits_string[1,1])==FALSE){
								
								tmpmat1<-rbind(tmpmat1,plot_res[[j]]$top_metfrag_hits_string)
								    
								# try(grid.arrange(plot_res[[j]]$p3,plot_res[[j]]$p1+theme(legend.position='hidden'),plot_res[[j]]$p2,plot_res[[j]]$myt,layout_matrix = rbind(c(1),c(2),c(3),c(4))),silent=TRUE)
									
							    
								
								
								}
							    }
						    
						
							}
						}
					}
			}
	}
	dev.off()
	
	
    if(length(tmpmat1)>0){
        if(nrow(tmpmat1)>0){
            write.table(tmpmat1,file="MetFragresults.txt",sep="\t",row.names=FALSE)
	    
	    
	    number_explained_peaks<-apply(tmpmat1,1,function(x){num_peaks<-length(unlist(strsplit(as.character(x[4]),";")))})
	    
	    tmpmat1<-cbind(tmpmat1,number_explained_peaks)
	    
	    write.table(tmpmat1,file="MetFragresults.txt",sep="\t",row.names=FALSE)
	    
	    
        }
    }
    
   #unlink("*.Rda")

    
}


if(method.process=="deconmsn"){


resloc<-MSMS_data_directory
#setwd(resloc)

dir.create(outloc)
setwd(outloc)
	dirfiles<-list.files(resloc,".raw$")
	print(dirfiles)

	#dirfiles<-paste(resloc,"",dirfiles,sep="")
	pdf("MSn_spectrum_deconmsn.pdf")
	
	if(process.raw.files==TRUE){
		for(filename in dirfiles){
		 
			#F1: start at scan 1; -L400: last scan #400; -C1: charge 1; -B50: minimum m/z; -T2000: maximum m/z; -I1: minimum number of ions in a scan
			#my $cmd=q("C:\\Program Files (x86)\\DeconMSn\\DeconMSn.exe" -XDTA -C1 -F1 -L400 -B50 -T2000 -I1 ).$MSMS_data_directory.$file;
			
			
			system2("C:\\Program Files (x86)\\DeconMSn\\DeconMSn.exe", args=c("-XDTA","-C1","-B1","-I1",paste(MSMS_data_directory,filename,sep="")))
			
			#system2("C:\\Program Files (x86)\\DeconMSn\\DeconMSn.exe", args=c("-XDTA","-B1","-I0",paste(MSMS_data_directory,filename,sep="")))
			
			
		}
		
	}



match_check=c(rep(0,length(query_mz_list)))


group_numbers<-c() #c("MSMS_ChearRef_hilic_5G") #c(group_numbers1,group_numbers2)


fileheader1=""; #"MSMS_ChearRef_hilic_5G_"; #"kl_20160212_"; 
fileheader2="_DeconMSn_log";


setwd(outloc)
log_fnames<-list.files(resloc,"log.txt$") #,full.names=TRUE)

#pdf("MetFrag_results.pdf")
top_metfrag_hits_string_all<-{}
for(q in 1:dim(query_mz_list)[1])
{
	query_mz_dirname<-query_mz_list[q,1]
	mzdir<-paste("mz",query_mz_dirname,sep="")
	
	
	query_mz<-query_mz_list[q,1]
	
	print("query")
	print(query_mz)

	for(gnum in 1:length(log_fnames))
	{

	
	log_fname<-log_fnames[gnum] 
	
	print("Processing file:")
	setwd(resloc)
	print(log_fname)
	decon_msn_log<-read.table(log_fname,header=TRUE,skip=8)
	
	setwd(outloc)
	decon_msn_log<-decon_msn_log[order(decon_msn_log$Mono_Mz),]
	
	#select scans that match the query mz within a defined threshold
	select_rows<-which(10^6*(abs(decon_msn_log$Parent_Mz-query_mz)/query_mz)<mz.thresh)

	select_rows2<-which(10^6*(abs(decon_msn_log$Mono_Mz-query_mz)/query_mz)<mz.thresh)

	select_rows<-c(select_rows,select_rows2)
	select_rows<-unique(select_rows)
	
	if(length(select_rows)>0){

	print("Matching scans based on mono m/z")
	print(head(select_rows))
	decon_msn_log_mz<-decon_msn_log[select_rows,]
	

	max_int<-max(decon_msn_log_mz$Parent_Intensity)[1]
	
	scan_num<-decon_msn_log_mz[which(decon_msn_log_mz$Parent_Intensity==max(decon_msn_log_mz$Parent_Intensity)[1])[1],1]
	
	int_summary<-summary(decon_msn_log_mz$Parent_Intensity)
	
	int_thresh<-int_summary[5]
	
	scan_nums<-decon_msn_log_mz[which(decon_msn_log_mz$Parent_Intensity>=max_int),1]
	
	print(scan_nums)
	
	for(scan_num in scan_nums)
	{

	
	cur_logfile<-gsub(log_fname,pattern="_DeconMSn_log.txt",replacement="")
	#cur_logfile<-gsub(cur_logfile,pattern="_",replacement="")
	if(scan_num<100){
		regpat=paste(cur_logfile,"(.)*(00)(",scan_num,")(.)*(dta)",sep="")
	}else{
		if(scan_num<1000){
			regpat=paste(cur_logfile,"(.)*(0)(",scan_num,")(.)*(dta)",sep="")
		}
	}
	
	setwd(resloc)
	#dir.create("deconmsn")

	mzfiles<-list.files(".",pattern=regpat) #glob2rx("*sca*USD*")) #(".",".dta")

	
	print("Using file:")
	print(mzfiles)
	
	setwd(outloc)
	dir.create("deconmsn")
	
	for(m in 1:length(mzfiles))
	{
	
		
		{
			setwd(resloc)
			data_m<-read.table(mzfiles[m])
			setwd(outloc)
			
			ppmb=(mz.thresh)*(query_mz/1000000)
		
			
			if(abs(data_m[1,1]-query_mz)<=ppmb)
			{	
				print(paste("match: ",query_mz," : ",mzfiles[m],sep=""))
			
				match_check[q]<-1
				data_m<-data_m[order(data_m[,2],decreasing=TRUE),]
				data_m[,2]<-100*(data_m[,2]/max(data_m[,2]))
				data_m<-data_m[which(data_m[,2]>1),]
				name<-rep("U",length(data_m[,1]))
				filename<-rep("U1",length(data_m[,1]))
				display<-rep("FALSE",length(data_m[,1]))
				identity<-rep("NA",length(data_m[,1]))
				
				data_m2<-cbind(data_m[,1],data_m[,2])
				colnames(data_m2)<-c("mz","intensity")
				
				print(head(data_m2))
				print("here1")
				msms_max_int_scan<-data_m2
				
				msms_max_int_scan<-as.data.frame(msms_max_int_scan)
				
				fname<-paste("deconmsn/",round(query_mz,2),"_",scan_num,cur_logfile,"Metlin.csv",sep="")
				
				write.table(data_m2,file=fname,sep=",",row.names=FALSE)
				
				fname<-paste("deconmsn/",round(query_mz,2),"_",scan_num,cur_logfile,"MetFrag.csv",sep="")
				
				write.table(data_m2,file=fname,sep="\t",row.names=FALSE)
				print("here2")
				data_m2<-cbind(name,filename,data_m[,1],data_m[,2],display,identity)
				colnames(data_m2)<-c("compound","filename","mz","intensity","display","identity")
				
				sample<-rep("U",length(data_m[,1]))
				instrument<-rep("Fusion",length(data_m[,1]))
				category1<-rep("U1",length(data_m[,1]))
				category2<-rep("U2",length(data_m[,1]))
				category3<-rep("U3",length(data_m[,1]))
				rt1d<-rep("U4",length(data_m[,1]))
				rt2d<-rep("U5",length(data_m[,1]))
				formula<-rep("U6",length(data_m[,1]))
				comment<-rep("NA",length(data_m[,1]))
				
				
				meta<-cbind(name,filename,sample,instrument,category1,category2,category3,rt1d,rt2d,formula,comment)
				
				colnames(meta)<-c("compound","filename","sample","instrument","category.1","category.2","category.3","rt.1d","rt.2d","formula","comment")
				
				fname<-paste("deconmsn/",round(query_mz,2),"_",scan_num,gnum,"meta1.csv",sep="")
				
				data_m2<-as.data.frame(data_m2)
				meta<-as.data.frame(meta)
				
				write.table(data_m2,file=fname,sep=";",row.names=FALSE)
				
				fname<-paste("deconmsn/",round(query_mz,2),"_",scan_num,gnum,"meta2.csv",sep="")
				write.table(meta,file=fname,sep=";",row.names=FALSE)
				
				fname<-paste("deconmsn/",round(query_mz,2),"_",scan_num,gnum,".msp",sep="")
				WriteMspFile(spectra=data_m2,metadata=meta,filename=fname)
				
				
				####code for searching in MetFrag
				
				
					                 DatabaseName<-DatabaseSearchName
                        
                        settingsObject<-list()
                        #
                        # set database parameters to select candidates
                        #
		
                        settingsObject[["DatabaseSearchRelativeMassDeviation"]]<-DatabaseSearchRelativeMassDeviation
                        settingsObject[["FragmentPeakMatchAbsoluteMassDeviation"]]<-FragmentPeakMatchAbsoluteMassDeviation
                        settingsObject[["FragmentPeakMatchRelativeMassDeviation"]]<-FragmentPeakMatchRelativeMassDeviation
                        settingsObject[["MetFragDatabaseType"]]<-DatabaseName
                        
                        if(dim(query_mz_list)[2]>2){
                            settingsObject[["NeutralPrecursorMass"]]<-as.numeric(as.character(query_mz_list[q,3]))
                            settingsObject[["NeutralPrecursorMolecularFormula"]]<-as.character(query_mz_list[q,4]) #"C23H45N1O4"
                        }
                        
                        if(DatabaseName=="LocalCSV"){
                            settingsObject[["MetFragDatabaseType"]]<-"LocalCSV"
                            settingsObject[["LocalDatabasePath"]]<-localcsv_file_path
                        }
                        
                     
			
                        data_m[,1]<-round(data_m[,1],4)
                        
                        settingsObject[["PeakList"]]<-cbind(data_m[,1],data_m[,2])
			options(digits=10)
			scored.candidates<-try(run.metfrag(settingsObject),silent=TRUE)
                     
                        rm(settingsObject)
			msms_max_int_scan<-as.data.frame(msms_max_int_scan)
                            
                      #maintext=paste("Precursor query m/z: ",round(query_mz,5)," time: ",round(query_mz_list[q,2],1)," \n Scan #: ",scan_num," \n File: ",log_fname,sep="")
				maintext=paste("MS/MS spectra for selected precursor m/z: ",round(query_mz,5)," peakID/scan: ",scan_num,sep="")
                            
                            if(plot.relative.intensity==TRUE){
                                msms_max_int_scan$intensity<-100*(msms_max_int_scan$intensity/max(msms_max_int_scan$intensity,na.rm=TRUE))
                                
                                ylab1="Relative Intensity (%)"
                                yinc=3.5
                            }else{
                                ylab="Intensity"
                                maxint<-max(msms_max_int_scan$intensity,na.rm=TRUE)
                                minint<-min(msms_max_int_scan$intensity,na.rm=TRUE)
                                if(abs(maxint-minint)>10000){
                                    yinc=median(msms_max_int_scan$intensity)
                                }else{
                                    yinc=min(100,abs(maxint-minint))
                                }
                            }
                            mytheme <- gridExtra::ttheme_default(
                            core = list(fg_params=list(cex = 0.5,hjust=0,x=0)),
                            colhead = list(fg_params=list(cex = 0.55,hjust=0,x=0)),
                            rowhead = list(fg_params=list(cex = 0.5,hjust=0,x=0)))
                         
			
                            x=msms_max_int_scan$mz
                            y=msms_max_int_scan$intensity
                            x=round(x,4)
                            y=round(y,2)
                            tmp1<-cbind(x,y)
                            tmp1<-as.data.frame(tmp1)
                   
			   if(nrow(tmp1)>10){
			   
				tmp1<-tmp1[order(msms_max_int_scan$intensity,decreasing=TRUE)[1:10],]
				
				minyval<-min(tmp1[,2])
				
				x=tmp1[,1]
				y=tmp1[,2]
			   }else{
				minyval=0
			   }
		
			 setwd(resloc)
			    mzmin=query_mz-(10^(-6)*query_mz*DatabaseSearchRelativeMassDeviation)
		            mzmax=query_mz+(10^(-6)*query_mz*DatabaseSearchRelativeMassDeviation)
			  msnraw.filename<-gsub(log_fname,pattern="_DeconMSn_log.txt",replacement=".mzXML")
			 xraw<-xcmsRaw(msnraw.filename,includeMSn=TRUE)
			 intvec<-getEIC(xraw,mzrange=c(mzmin,mzmax),step=0.001) #$intensity
			    intvec=intvec@eic[[1]][[1]]
			 #   intvec=cbind(decon_msn_log_mz[,1],decon_msn_log_mz$Parent_Intensity) #intvec@eic[[1]][[1]]
			    colnames(intvec)<-c("rt","intensity")
			    intvec<-as.data.frame(intvec)
			    
			    rm(xraw)
			#print(head(intvec))
	
			mainlab<-paste("EIC for mz:",round(query_mz,5),"(file:",msnraw.filename,")",sep=" ")
			if(plot.relative.intensity==TRUE){			
				intvec$intensity<-100*(intvec$intensity)/max(intvec$intensity,na.rm=TRUE)
				ylab1="Relative Intensity(%)"
			}else{
				ylab1="Intensity"
			}
			intvec=cbind(intvec$rt,intvec$intensity)
			intvec=as.data.frame(intvec)
			colnames(intvec)<-c("rt","intensity")
		        
			
			 setwd(outloc)
			  
			#  rawfile <- openMSfile(raw_fname,backend='pwiz') 
			 # scan_peak_list <- mzR::peaks(rawfile) 
  
		#	  p1=ggplot(data=ms1_data_clust,aes(x=rt, y=intensity)) + labs(x="time",y=ylab1) + ggtitle(mainlab1) +geom_area(aes(fill=clusterlabels,group=clusterlabels)) + geom_text(data=ms1_data_clust,aes(x=rt, y=intensity,label=clusterlabels),size=2) + theme(axis.text=element_text(size=7),
                     #       axis.title=element_text(size=7), plot.title=element_text(size=7))
			    
	               	plots_list=ggplot(data=intvec,aes(x=rt, y=intensity)) + labs(x="scan",y=ylab1) + ggtitle(mainlab) +geom_area(aes()) + theme(axis.text=element_text(size=7),
                                                        axis.title=element_text(size=7), plot.title=element_text(size=7))
                         #   p1=ggplot(data=msms_max_int_scan,aes(x=rt, y=intensity)) + labs(x="scan",y=ylab1) + ggtitle(mainlab1) +geom_area(aes(fill=clusterlabels,group=clusterlabels)) + geom_text(data=ms1_data_clust,aes(x=rt, y=intensity,label=clusterlabels),size=2) + theme(axis.text=element_text(size=7),
                      #      axis.title=element_text(size=7), plot.title=element_text(size=7))
                            
			    mainlab2<-paste("MS/MS EIC for precursor mz:",round(query_mz,5),"(file:",msnraw.filename,")",sep=" ")
			    
			      intvec2=cbind(decon_msn_log_mz[,1],decon_msn_log_mz$Parent_Intensity)
			      colnames(intvec2)<-c("rt","intensity")
			    intvec2<-as.data.frame(intvec2)
			     	p1=ggplot(data=intvec2,aes(x=rt, y=intensity)) + labs(x="scan",y=ylab1) + ggtitle(mainlab2) +geom_area(aes()) + theme(axis.text=element_text(size=7),
                                                       axis.title=element_text(size=7), plot.title=element_text(size=7))
			    
                            p2=ggplot(data=tmp1,aes(x,y)) + labs(x="m/z",y=ylab1)  + ggtitle(maintext) + geom_segment(x=x,xend=x,y=10,yend=y,data=tmp1) + geom_text(data=tmp1,aes(x=x,y=y,label=x),size=2) + theme(axis.text=element_text(size=7),
                            axis.title=element_text(size=7), plot.title=element_text(size=7))
 
                            myt <-{} #gridExtra::tableGrob(top_metfrag_hits_string, theme = mytheme)
                            
                            tmp2<-new("list")
                           top_metfrag_hits_string<-{}
			 if(is(scored.candidates,"try-error")){
                                tmp2=NULL
                        }else{
                        scored.candidates<-as.data.frame(scored.candidates)

				
                        scored.candidates$InChIKey<-as.character(scored.candidates$InChIKey)
                         fname<-paste(query_mz,"_",scan_num,log_fname,"scoredmetfrag.Rda",sep="")
                        #save(scored.candidates,file=fname)

			     if(nrow(scored.candidates)>0){


                           # data(hmdbAllinf)
                            scored.candidates<-scored.candidates[order(scored.candidates$Score,decreasing=TRUE),]
                           # hmdbAllinf$InChiKey=gsub(hmdbAllinf$InChiKey,pattern="InChIKey=",replacement="")
                           # hmdbAllinf$InChiKey<-as.character(hmdbAllinf$InChiKey)
                           # hmdbAllinf<-hmdbAllinf[,c("HMDBID","Name","InChiKey","KEGGID","BioFluidLocation","CellularLocation","Origin","BioFunctions","HMDBStatus")]
                            m1<-scored.candidates #merge(scored.candidates,hmdbAllinf,by.x="InChIKey",by.y="InChiKey")
                            rm(scored.candidates)
			    m1<-m1[order(m1$Score,decreasing=TRUE),]
                     
			       top_scores<-m1[which(m1$Score>0),]
                            dup_inchi_check<-which(duplicated(top_scores$InChIKey)==TRUE)
                            if(length(dup_inchi_check)>0){
                                top_scores<-top_scores[-dup_inchi_check,]
                            }
                            top_scores<-top_scores[order(top_scores$Score,decreasing=TRUE),]
                            if(nrow(top_scores)>max_metfrag_hits){
                                top_scores<-top_scores[1:max_metfrag_hits,]
                            }
                            top_metfrag_hits_string<-{}
                            for(mrownum in 1:nrow(top_scores)){
                                expl_peaks<-gsub(top_scores$ExplPeaks[mrownum],pattern="_[0-9.]*",replacement="")
                                top_metfrag_hits_string<-rbind(top_metfrag_hits_string,cbind(as.character(top_scores$Identifier[mrownum]),as.character(top_scores$Name[mrownum]),round(top_scores$Score[mrownum],2),round(top_scores$FragmenterScore[mrownum],2),expl_peaks))
                            }
                            top_metfrag_hits_string<-as.data.frame(top_metfrag_hits_string)
                            colnames(top_metfrag_hits_string)<-c("ID", "Name", "Score (0 to1)","FragmenterScore","Explained Peaks")
                            rownames(top_metfrag_hits_string)<-NULL
                            setwd(outloc)
                           # setwd("MetFrag")
                            fname<-paste(query_mz,"_",scan_num,gnum,"candidates.csv",sep="")
                            write.csv(m1,fname,row.names=FALSE)
                            #fname<-paste(query_mz,"_",scan_num,log_fname,"candidatesHMDB.csv",sep="")
                            #write.csv(m1,fname,row.names=FALSE)
 
				print("here3")
				print(top_metfrag_hits_string)

           		     myt <- gridExtra::tableGrob(top_metfrag_hits_string, theme = mytheme)
 
                           top_metfrag_hits_string<-cbind(top_metfrag_hits_string,query_mz,scan_num,log_fname)
                            top_metfrag_hits_string<-as.data.frame(top_metfrag_hits_string)
                            colnames(top_metfrag_hits_string)<-c("ID", "Name", "Score (0 to1)","FragmenterScore","Explained Peaks","query_mz","ScanPeakID","FileName")
                            rownames(top_metfrag_hits_string)<-NULL
			    
				top_metfrag_hits_string_all<-rbind(top_metfrag_hits_string_all,top_metfrag_hits_string)
			    }                            
			}                            
                            tmp2<-list(p1=p1,p2=p2,myt=myt,top_metfrag_hits_string=top_metfrag_hits_string,p3=plots_list)
                            fname<-paste(query_mz,"tmp2.Rda",sep="")
                            save(tmp2,file=fname)
			      #   try(grid.arrange(plots_list,p2,myt),silent=TRUE)
			      #movehere
                          # grid.arrange(plots_list,p2,myt)
			  
			  if(length(tmp2$myt)>0){
			  grid.arrange(tmp2$p3,tmp2$p2,tmp2$myt,layout_matrix=rbind(c(1),c(2),c(3),c(4)))
			  }else{
				 grid.arrange(tmp2$p3,tmp2$p2,layout_matrix=rbind(c(1),c(2),c(3),c(4)))
			  }
			  rm(p1)
			  rm(p2)
			  rm(myt)
			  rm(top_metfrag_hits)
			  rm(plots_list)
			  rm(intvec)
			  rm(data_m2)
			 
			  rm(tmp1)
			  rm(m1)
			  rm(meta)
			  
			  setwd(outloc)

			  # return(tmp2)
                                
                             
				
				
				###################
			}
			
		}
	}
	
	}
	}

	}
	#return(msms_max_int_scan)
}

print("Match found for")
print(query_mz_list[(which(match_check==1)),])
dev.off()
setwd(outloc)
write.csv(top_metfrag_hits_string_all,file="MetFrag_results.csv",row.names=FALSE)

#try(unlink(".","*.Rda"),silent=TRUE)
}

}

getturnpoints<-function (x, calc.proba = TRUE) 
{
    data <- deparse(substitute(x))
    if (is.null(ncol(x)) == FALSE) 
        stop("Only one series can be treated at a time")
    x <- as.vector(x)
    n <- length(x)
    diffs <- c(x[1] - 1, x[1:(n - 1)]) != x
    uniques <- x[diffs]
    n2 <- length(uniques)
    poss <- (1:n)[diffs]
    exaequos <- c(poss[2:n2], n + 1) - poss - 1
    if (n2 < 3) {
        warning("Less than 3 unique values, no calculation!")
        nturns <- NA
        firstispeak <- FALSE
        peaks <- rep(FALSE, n2)
        pits <- rep(FALSE, n2)
        tppos <- NA
        proba <- NA
        info <- NA
    }
    else {
        m <- n2 - 2
        ex <- matrix(uniques[1:m + rep(3:1, rep(m, 3)) - 1], 
            m)
        peaks <- c(FALSE, apply(ex, 1, max, na.rm = TRUE) == 
            ex[, 2], FALSE)
        pits <- c(FALSE, apply(ex, 1, min, na.rm = TRUE) == ex[, 
            2], FALSE)
        tpts <- peaks | pits
        if (sum(tpts) == 0) {
            nturns <- 0
            firstispeak <- FALSE
            peaks <- rep(FALSE, n2)
            pits <- rep(FALSE, n2)
            tppos <- NA
            proba <- NA
            info <- NA
        }
        else {
            tppos <- (poss + exaequos)[tpts]
            tptspos <- (1:n2)[tpts]
            firstispeak <- tptspos[1] == (1:n2)[peaks][1]
            nturns <- length(tptspos)
            if (isTRUE(calc.proba)) {
                if (nturns < 2) {
                  inter <- n2 + 1
                  posinter1 <- tptspos[1]
                }
                else {
                  inter <- c(tptspos[2:nturns], n2) - c(1, tptspos[1:(nturns - 
                    1)]) + 1
                  posinter1 <- tptspos - c(1, tptspos[1:(nturns - 
                    1)])
                }
                posinter2 <- inter - posinter1
                posinter <- pmax(posinter1, posinter2)
                proba <- 2/(inter * gamma(posinter) * gamma(inter - 
                  posinter + 1))
                info <- -log(proba, base = 2)
            }
            else {
                proba = NULL
                info <- NULL
            }
        }
    }
    res <- list(data = data, n = n, points = uniques, pos = (poss + 
        exaequos), exaequos = exaequos, nturns = nturns, firstispeak = firstispeak, 
        peaks = peaks, pits = pits, tppos = tppos, proba = proba, 
        info = info)
    class(res) <- "turnpoints"
    res
}

feature_evaluation<-function(Xmat=NA,featuretable,rawdata,tool=NA,outloc,rawformat='mzXML',mz_tolerance=10,time_tolerance=45, nu_cores=2){
  
  #load the required packages
  library(mzR)
  #library(msdata)
  library(MSnbase)
  library(parallel)
  
  # Calculate the number of cores
  #nu_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- parallel::makeCluster(ceiling(nu_cores*0.6))
  
  ######################### define function ################################## 
  
  ###
  # This is the main function
  scan_feature <- function(i,scan_peak_list,timelist,featlist,mz_tolerance=10,time_tolerance=15,outloc){

    library(ggplot2)
    library(pracma)
    
    ######################### define function ################################## 
    
    ###
    #This function is used to extract the EIC by binning the data point with
    #the binning tolerance around the mass_value 
    eic_extraction_binning<-function(scan_peak_list,timelist,mz,tolerance=10){
      count=0
      totalscan=length(scan_peak_list)
      eic_intensity_list=as.data.frame(matrix(NA,nrow=totalscan,ncol=3))
      colnames(eic_intensity_list)=c("scan","time","intensity")
      for(ii in 1:totalscan){
        peaklist=as.data.frame(scan_peak_list[[ii]])
	
	#######added by Karan#######
	colnames(peaklist)<-c("mz","intensity")
	
	
        find_index=(peaklist$mz>mz-(tolerance*10^-6*mz) & peaklist$mz<mz+(tolerance*10^-6*mz))
        if(sum(find_index)>0){
          eic_intensity_list[ii,"scan"] <- ii
          eic_intensity_list[ii,"time"] <- timelist[ii]
          if(sum(peaklist[find_index,"intensity"])>0){
            eic_intensity_list[ii,"intensity"] <- sum(peaklist[find_index,"intensity"])
            count=count+1
          }else{
            eic_intensity_list[ii,"intensity"] <- NA
          }
        }else{
          eic_intensity_list[ii,"scan"] <- ii
          eic_intensity_list[ii,"time"] <- timelist[ii]
          eic_intensity_list[ii,"intensity"] <- NA
        }
      }
      return(list(eic_intensity_list,count))
    }
    ###
    
    ###
    #This function is used to generate the EIC plot
    plot_eic<-function(plotdata,feature,tolerance,count){
      
      pos=seq(0,max(plotdata$time),by=30)
      if(pos[length(pos)]<max(plotdata$time)){
        pos=c(pos,pos[length(pos)]+30)
      }
      
      plotdata[plotdata$time>=as.numeric(feature[2]-tolerance) & plotdata$time<=as.numeric(feature[2]+tolerance),"col"]<-"red"
      plotdata[!(plotdata$time>=as.numeric(feature[2]-tolerance) & plotdata$time<=as.numeric(feature[2]+tolerance)),"col"]<-"black"
      
      feature=suppressWarnings(as.numeric(as.character(feature)))
      
      if(dim(plotdata[!is.na(plotdata$intensity) & plotdata$col=="red",])[1]>=5){
        subtitle=paste(
          paste("Retention Time(s):",round(feature[2],1)),
	  paste("Number of data points:",feature[3]),
          paste(paste("MCQ index(whole EIC):",round(feature[4],4)),"; ",paste("Zigzag index(whole EIC):",round(feature[5],4)),sep=""),
          paste(paste("MCQ index(Peak area):",round(feature[6],4)),"; ",paste("Zigzag index(Peak area):",round(feature[7],4)),sep=""),
          paste(paste("Sharpness(Peak area):",round(feature[8],4)),"; ",paste("Gaussian Similarity(Peak area):",round(feature[9],4)),sep=""),
          paste(paste("Peak Significance Level(Peak area):",round(feature[10],4)),"; ",paste("TPASR(Peak area):",round(feature[11],4)),"; ",paste("SNR(Peak area):",round(feature[12],4)),sep=""),
          sep="\n")
      }else{
        subtitle=paste("Retention Time(s):",round(feature[2],1))
      }
      
      if(count>0){
      print(ggplot() +
              geom_point(data=plotdata,aes(x=time,y=intensity,color=col),na.rm=TRUE) +
              geom_vline(xintercept=as.numeric(feature[2]-tolerance),color="red",linetype = 2) +
              geom_vline(xintercept=as.numeric(feature[2]+tolerance),color="red",linetype = 2) +
              scale_y_continuous(name="Intensity") +
              scale_x_continuous(name="Retention time(s)",breaks=pos)+
              labs(subtitle=subtitle, 
                   title=paste("EIC plot for ",round(feature[1],5),sep="")
              ) +
              theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                    axis.line.y = element_line(size = 0.5, colour = "black"),
                    axis.line = element_line(size=1, colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.title=element_text(size = 20),
                    plot.subtitle=element_text(size = 10),
                    text=element_text(size = 16),
                    plot.margin=margin(10,5,10,1),
                    axis.text.x=element_text(colour="black", size = 10),
                    axis.text.y=element_text(colour="black", size = 10))+
              scale_colour_manual(values =c("black"="black","red"="red"))+
              guides(color="none"))
      }else{
	print(ggplot() + labs(subtitle=subtitle, 
                   title=paste("EIC plot for ",round(feature[1],5),sep="")
              ) +
              theme(axis.line.x = element_line(size = 0.5, colour = "black"),
                    axis.line.y = element_line(size = 0.5, colour = "black"),
                    axis.line = element_line(size=1, colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    plot.title=element_text(size = 20),
                    plot.subtitle=element_text(size = 10),
                    text=element_text(size = 16),
                    plot.margin=margin(10,5,10,1),
                    axis.text.x=element_text(colour="black", size = 10),
                    axis.text.y=element_text(colour="black", size = 10))+
              scale_colour_manual(values =c("black"="black","red"="red"))+
              guides(color="none"))
      }
    }
    ###
    
    ###
    #This function is used to calculate the MCQ index
    mcq_index_cal <- function(eic_intensity_list,window=3,min_time,max_time){
      datalist=eic_intensity_list[!is.na(eic_intensity_list$intensity),]
      
      datalist=datalist[datalist$time>=min_time & datalist$time<=max_time,]
      
      eic=datalist$intensity
      
      if(length(eic)>=5){
        eic_scale=eic/sqrt(sum(eic^2))
        
        pad_num = floor(window/2)
        eic_padding=c(rep(eic[1],pad_num),eic,rep(eic[length(eic)],window-pad_num))
        eic_window=eic
        for(ii in 1:length(eic)){
          eic_window[ii]=mean(eic_padding[ii:(ii+window)])
        }
        
        eic_window_mean  = mean(eic_window)
        eic_window_sd   = sd(eic_window)
        eic_window_standarization = (eic_window-eic_window_mean)/eic_window_sd
        eic_window_standarization_scale=eic_window_standarization/sqrt(sum(eic_window_standarization^2))
        
        mcq_index=sum(eic_scale*eic_window_standarization_scale)
      }else{
        mcq_index=NA
      }
      
      return(mcq_index)
    }
    ###
    
    ###
    #This function is used to calculate the ZigZag index 
    zig_zag_index_cal <- function(eic_intensity_list,min_time,max_time){
      zig_zag_sum=0
      datalist=eic_intensity_list[!is.na(eic_intensity_list$intensity),]
      
      datalist=datalist[datalist$time>=min_time & datalist$time<=max_time,]
      
      eic=datalist$intensity
      
      if(length(eic)>=5){
        epi=max(eic)-((eic[1]+eic[2]+eic[length(eic)-1]+eic[length(eic)])/4)
        for(ii in 2:(length(eic)-1)){
          local_zig_zag=(2*eic[ii]-eic[ii-1]-eic[ii+1])^2
          zig_zag_sum=zig_zag_sum+local_zig_zag
        }
        
        zig_zag_index=zig_zag_sum/(epi^2*length(eic))
      }else{
        zig_zag_index=NA
      }
      
      return(zig_zag_index)
    }
    ###
    
    ###
    #This function is used to calculate the sharpness
    sharpness_cal <- function(eic_intensity_list,time,tolerance){
      sharpness=0
      datalist=eic_intensity_list[!is.na(eic_intensity_list$intensity),]
      
      datalist=datalist[datalist$time>=time-tolerance & datalist$time<=time+tolerance,]
      
      eic=datalist$intensity
      
      if(length(eic)>=5){
        apex_index=which(eic==max(eic))[1]
        
        for(ii in 1:length(eic)){
          if(ii < apex_index){
            sharpness=sharpness+(eic[ii+1]-eic[ii])/max(0.0001,eic[ii])
          }else{
            if(ii<length(eic)){
              sharpness=sharpness+(eic[ii]-eic[ii+1])/max(0.0001,eic[ii+1])
            }
          }
        }
      }else{
        sharpness=NA
      }
      
      
      return(sharpness)
      
    }
    ###
    
    ###
    #This function is used to calculate the sharpness
    Gaussian_similarity_cal <- function(eic_intensity_list,time,tolerance){
      Gaussian_similarity=0
      datalist=eic_intensity_list[!is.na(eic_intensity_list$intensity),]
      
      datalist=datalist[datalist$time>=time-tolerance & datalist$time<=time+tolerance,]
      
      eic=datalist$intensity
      
      if(length(eic)>=5){
        #### GUSSIAN CURVE FITTING
        x<-1:length(eic)
        y<-eic
        p=polyfit(x,log(y),2)
        fit_sigma=suppressWarnings(sqrt(-1/(2*p[1])))
        fit_mu=p[2]*fit_sigma^2
        fit_a=exp(p[3]+fit_mu^2/(2*fit_sigma^2))
        
        Gaussian_Peak=fit_a*exp((c(1:length(eic))-fit_mu)^2/(-2*fit_sigma^2))
        Gaussian_Peak_standard=(Gaussian_Peak-mean(Gaussian_Peak))/sd(Gaussian_Peak)
        Gaussian_Peak_Scale=Gaussian_Peak_standard/sqrt(sum(Gaussian_Peak_standard^2))
        eic_standard=(eic-mean(eic))/sd(eic)
        eic_scale=eic_standard/sqrt(sum(eic_standard^2))
        Gaussian_similarity=sum(Gaussian_Peak_Scale*eic_scale)
        
      }else{
        Gaussian_similarity=NA
      }
      
      return(Gaussian_similarity)                           
      
    }
    ###
    
    ###
    #This function is used to calculate the peak significance level
    significance_cal <- function(eic_intensity_list,time,tolerance){
      significance=0
      datalist=eic_intensity_list[!is.na(eic_intensity_list$intensity),]
      
      datalist=datalist[datalist$time>=time-tolerance & datalist$time<=time+tolerance,]
      
      eic=datalist$intensity
      
      if(length(eic)>=5){
        apex_index=which(eic==max(eic))[1]
        
        if(apex_index!=1 && apex_index!=length(eic) && length(eic)>4){
          sum1=sum(eic[c(1,2,length(eic)-1,length(eic))])
          sum2=eic[apex_index-1]+ max(eic) +eic[apex_index+1]
          significance =sum2*4/(max(sum1,0.01)*3)
        }
      }else{
        significance=NA
      }
      
      return(significance)
    }
    ###
    
    ###
    #This function is used to calculate the Triangle Peak Area Similarity Ratio (TPASR) 
    TPASR_cal <- function(eic_intensity_list,time,tolerance){
      
      datalist=eic_intensity_list[!is.na(eic_intensity_list$intensity),]
      
      datalist=datalist[datalist$time>=time-tolerance & datalist$time<=time+tolerance,]
      
      eic=datalist$intensity
      
      if(length(eic)>=5){
        apex_index=which(eic==max(eic))[1]
        
        peak_level=max(eic)
        base_level=0
        
        if(apex_index!=1 && apex_index!=length(eic) && length(eic)>4){
          base_level=(sum(eic[c(1,2,length(eic)-1,length(eic))]))/4
          peak_level=(eic[apex_index-1]+ max(eic) +eic[apex_index+1])/3
        }
        
        trangular_area=0.5*length(eic)*(peak_level-min(eic[1],eic[length(eic)]))
        peak_area=sum(eic)-min(eic[1],eic[length(eic)])*length(eic)
        TPASR=abs(trangular_area-peak_area)/trangular_area
      }else{
        TPASR=NA
      }
      
      return(TPASR)
    }
    ###
    
    ###
    #This function is used to calculate the Signal-to-noise ratio (SNR)
    SNR_cal <- function(eic_intensity_list,time,tolerance){
      
      datalist=eic_intensity_list[!is.na(eic_intensity_list$intensity),]
      
      noiselist=datalist[datalist$time<time-tolerance | datalist$time>time+tolerance,]
      
      signallist=datalist[datalist$time>=time-tolerance & datalist$time<=time+tolerance,]
      
      noise_intensity=noiselist$intensity
      
      signal_intensity=signallist$intensity
      
      if(length(signal_intensity)>=5 & length(noise_intensity)>=3){
        SNR=(max(signal_intensity)-mean(noise_intensity))/sd(noise_intensity)
      }else{
        SNR=NA
      }
      
      return(SNR)
      
    }
    ###
    
    #########################################################################
    
    #print(i)
    target_mz=featlist[i,1]
    target_time=featlist[i,2]
    plotname=paste(outloc,"/PDF_plots/",round(target_mz,5),"_",round(target_time,1),".pdf",sep="")
    result <-eic_extraction_binning(scan_peak_list=scan_peak_list,timelist=timelist,mz=target_mz,tolerance=mz_tolerance)
    eic_intensity_list <- result[[1]]
    count <- result[[2]]
    print(count)
    if(count>5){
      mcq_index_whole_area <- mcq_index_cal(eic_intensity_list,3,0,max(eic_intensity_list$time))
      zig_zag_index_whole_area <- zig_zag_index_cal(eic_intensity_list,0,max(eic_intensity_list$time))
      mcq_index_specific_area <- mcq_index_cal(eic_intensity_list,3,target_time-time_tolerance,target_time+time_tolerance)
      zig_zag_index_specific_area <- zig_zag_index_cal(eic_intensity_list,target_time-time_tolerance,target_time+time_tolerance)
      sharpness_specific_area <- sharpness_cal(eic_intensity_list,target_time,time_tolerance)
      Gaussian_similarity_specific_area <- Gaussian_similarity_cal(eic_intensity_list,target_time,time_tolerance)
      significance_specific_area <- significance_cal(eic_intensity_list,target_time,time_tolerance)
      TPASR_specific_area <- TPASR_cal(eic_intensity_list,target_time,time_tolerance)
      SNR_specific_area <- SNR_cal(eic_intensity_list,target_time,time_tolerance)
      #pdf(plotname)
      #plot_eic(eic_intensity_list,featlist[i,],time_tolerance)
     # dev.off()
    }else{
      mcq_index_whole_area<-NA
      zig_zag_index_whole_area<-NA
      mcq_index_specific_area<-NA
      zig_zag_index_specific_area<-NA
      sharpness_specific_area<-NA
      Gaussian_similarity_specific_area<-NA
      significance_specific_area<-NA
      TPASR_specific_area<-NA
      SNR_specific_area<-NA
    }
    return(list("eic_intensity_list"=eic_intensity_list,"featlist"=featlist[i,],time_tolerance=time_tolerance,qmetric=c(count,mcq_index_whole_area,zig_zag_index_whole_area,mcq_index_specific_area,zig_zag_index_specific_area,sharpness_specific_area,Gaussian_similarity_specific_area,significance_specific_area,TPASR_specific_area,SNR_specific_area)))
  }
  ###
  
  
  #########################################################################
  
  #read in the feature table
  if(!is.na(featuretable)){
    featable=read.table(featuretable,sep="\t",header=TRUE)
  }else{
    if(!any(is.na(Xmat))){
      featable=Xmat
    }else{
      stop("No feature table was provided. Please check your parameters again.")
    }
  }
  
  #read in the rawdata
  if(length(grep('mzXML',rawformat,ignore.case=FALSE))!=0){
    rawfile <- openMSfile(rawdata,backend='pwiz') 
  }else{
    if(length(grep('cdf',rawformat,ignore.case=FALSE))!=0){
      rawfile <- openMSfile(rawdata,backend='netCDF') 
    }else{
      stop("Only mzXML or CDF format can be accepted. Please set rawformat to 'cdf' or 'mzXML'.")
    }
  }
  
  #get the feature list
  if(length(grep('apLCMS',tool,ignore.case=TRUE))!=0 || length(grep('xcms',tool,ignore.case=TRUE))!=0){
    featlist=featable[,c("mz","time")]
  }else{
    featlist=featable[,c(1,2)] # assume the first column is mz and second is time
  }
  
  scan_acquisition_time_arry <- mzR::header(rawfile)$retentionTime #The acquisition time of each scan
  scan_peak_list <- mzR::peaks(rawfile) #get all data ponts for each scan
  totalscan <- length(scan_peak_list) # total scan number 
  
  dir.create(paste(outloc,"PDF_plots",sep="/"))
  
  ### scan all features in parallel way
  start_time <- Sys.time()
  
    
#######Edited by Karan#######
  out_all<-parLapply(cl, as.character(1:nrow(featlist)),scan_feature,scan_peak_list=scan_peak_list,featlist=featlist,timelist=scan_acquisition_time_arry,mz_tolerance=mz_tolerance,time_tolerance=time_tolerance,outloc=outloc)

 save(out_all,file="out_all.Rda")
out_all<-ldply(out_all,rbind)

#######End of Karan's edits#######

  colnames(out_all)<-c("count","mcq_index_whole_area","zig_zag_index_whole_area","mcq_index_specific_area","zig_zag_index_specific_area","sharpness_specific_area","Gaussian_similarity_specific_area","significance_specific_area","TPASR_specific_area","SNR_specific_area")
  featlist<-cbind(featlist,out_all)
  end_time <- Sys.time()
  
  diff=end_time - start_time
  print(diff)
  
  stopCluster(cl)
  
  write.table(featlist,paste(outloc,"/EIC_evaluation",rawdata,".csv",sep=""),sep=",",row.names=FALSE)
  
  
}


check_element<-function(curformula,elementname){
    
    curformula<-as.character(curformula)
    strsplit_var<-strsplit(curformula,split="")
    
    strsplit_elem<-strsplit(elementname,split="")
    elem_len<-length(strsplit_elem[[1]])
    
    g3<-gregexpr(curformula,pattern=paste(elementname,"[0-9]*",sep=""))
    
    numelement<-0
    
    if(length(g3)>0){
        if(elem_len>1){
            
            regexp_len3<-attr(g3[[1]],"match.length")[1]
            if(regexp_len3>2){
                regexp_len3<-attr(g3[[1]],"match.length")[1]-elem_len+1
            }else{
                regexp_len3<-attr(g3[[1]],"match.length")[1]-elem_len
            }
        }else{
            regexp_len3<-attr(g3[[1]],"match.length")[1]-elem_len
        }
        
        if(regexp_len3>=0){
            
            if(regexp_len3==0){
                numelement<-1
                
            }else{
                numelement<-paste(strsplit_var[[1]][(g3[[1]][1]+elem_len):(g3[[1]][1]+regexp_len3)],collapse="")
                
                numelement<-as.numeric(numelement)
                
            }
        }else{
            numelement<-0
            
        }
        
    }
    return(numelement)
}



check_golden_rules<-function(curformula,NOPS_check=FALSE){
    
    
    numnitrogens<-check_element(curformula,"N")
    numcarbons<-check_element(curformula,"C")
    
    numoxygens<-check_element(curformula,"O")
    
    numhydrogens<-check_element(curformula,"H")
    
    numphos<-check_element(curformula,"P")
    numsulphur<-check_element(curformula,"S")
    
    
    
    if(numcarbons<1){
        
        bool_check<-0
    }else{
        
        max_hydrogens<-(2*numcarbons)+numnitrogens+2
        
        nitrogen_to_carbon_ratio<-numnitrogens/numcarbons
        
        oxygen_to_carbon_ratio<-numoxygens/numcarbons
        
        phosphorus_to_carbon_ratio<-numphos/numcarbons
        
        sulphur_to_carbon_ratio<-numsulphur/numcarbons
        
        hydrogens_to_carbon_ratio<-numhydrogens/numcarbons
        
        bool_check<-1
        
        
        
        if(hydrogens_to_carbon_ratio<0.1 | hydrogens_to_carbon_ratio>6){
            
            bool_check=0
        }
        
        if(nitrogen_to_carbon_ratio>4){
            
            bool_check=0
        }
        
        if(oxygen_to_carbon_ratio>3){
            
            bool_check=0
        }
        
        if(phosphorus_to_carbon_ratio>2){
            
            bool_check=0
        }
        
        if(sulphur_to_carbon_ratio>3){
            
            bool_check=0
        }
        
        if(NOPS_check==TRUE){
            #NOPS>1
            if(numnitrogens>1 & numoxygens>1 & numphos>1 & numsulphur>1){
                if(numnitrogens>10 | numoxygens>20 | numphos>4 | numsulphur>3){
                    
                    bool_check<-0
                }
                
            }
            
            
            #NOP>3
            if(numnitrogens>3 & numoxygens>3 & numphos>3){
                if(numnitrogens>11 | numoxygens>22 | numphos>6){
                    
                    bool_check<-0
                }
                
            }
            
            #OPS>1
            if(numoxygens>1 & numphos>1 & numsulphur>1){
                if(numoxygens>14 | numphos>3 | numsulphur>3){
                    
                    bool_check<-0
                }
                
            }
            
            #PSN>1
            if(numnitrogens>1 & numphos>1 & numsulphur>1){
                if(numnitrogens>4 | numphos>3 | numsulphur>3){
                    
                    bool_check<-0
                }
                
            }
            
            #NOS>6
            if(numnitrogens>6 & numoxygens>6 & numsulphur>6){
                if(numnitrogens>19 | numoxygens>14 | numsulphur>8){
                    
                    bool_check<-0
                }
                
            }
            
        }
        
        
    }
    res<-cbind(curformula,bool_check)
    res<-as.data.frame(res)
    
    
    return(res)
    
}

simpleAnnotation2<-function(dataA,max.mz.diff=10,max_diff_rt=10,num_nodes=2,queryadductlist=c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H"),
gradienttype="Acetonitrile",mode="pos",outloc,db_name="KEGG", NOPS_check=TRUE)
{
    
    
    if(db_name=="KEGG"){
        
        
        data(keggCompMZ)
        chemCompMZ<-keggCompMZ
        rm(keggCompMZ)
    }else{
        if(db_name=="HMDB"){
            data(hmdbCompMZ)
            chemCompMZ<-hmdbCompMZ
	    rm(hmdbCompMZ)
        }else{
            if(db_name=="T3DB"){
                data(t3dbCompMZ)
                chemCompMZ<-t3dbCompMZ
                rm(t3dbCompMZ)
            }else{
                
                if(db_name=="LipidMaps"){
                    data(lipidmapsCompMZ)
                    chemCompMZ<-lipidmapsCompMZ
                    rm(lipidmapsCompMZ)
                }
            }
            
        }
    }
    chemCompMZ$mz<-round(chemCompMZ$mz,4)
    dataA$mz<-round(dataA$mz,4)
 
    data(adduct_table)
    allowWGCNAThreads(nThreads=num_nodes)
    #rm(adduct_table)
    #data(adduct_table)
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="2M\\+2H",replacement="M+H")
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="3M\\+3H",replacement="M+H")
    
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="2M\\+2Na",replacement="M+Na")
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="3M\\+3Na",replacement="M+Na")
    
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="2M\\+2K",replacement="M+K")
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="3M\\+3K",replacement="M+K")
    
    adduct_table<-unique(adduct_table)
    
    suppressWarnings(dir.create(outloc))
    
    if(queryadductlist=="all" & mode=="pos"){
        
        adduct_names<-adduct_table$Adduct[(adduct_table$Type=="S" & adduct_table$Mode=="positive") | (adduct_table$Type==gradienttype & adduct_table$Mode=="positive")]
        
        adduct_table<-adduct_table[which(adduct_table$Adduct%in%adduct_names),]
        
    }else{
        if(queryadductlist=="all" & mode=="neg"){
            
            adduct_names<-adduct_table$Adduct[(adduct_table$Type=="S" & adduct_table$Mode=="negative") | (adduct_table$Type==gradienttype & adduct_table$Mode=="negative")]
            adduct_table<-adduct_table[which(adduct_table$Adduct%in%adduct_names),]
        }else{
            
            adduct_names<-adduct_table$Adduct[which(adduct_table$Adduct%in%queryadductlist)]
            
            adduct_table<-adduct_table[which(adduct_table$Adduct%in%adduct_names),]
            
        }
    }
    
    adduct_names<-unique(adduct_names)
    
    
    
    
    chemCompMZ<-chemCompMZ[which(chemCompMZ$Adduct%in%adduct_names),]
    
    print("dim of Db")
    print(dim(chemCompMZ))
    
    
    cl<-parallel::makeCluster(num_nodes)
    
    clusterEvalQ(cl, library(XML))
    clusterEvalQ(cl, library(R2HTML))
    clusterEvalQ(cl, library(RCurl))
    clusterEvalQ(cl, library(SSOAP))
    clusterEvalQ(cl, library(limma))
    
    clusterEvalQ(cl, library(plyr))
    
    clusterEvalQ(cl, "processWSDL")
    clusterEvalQ(cl, library(png))
    clusterExport(cl, "Annotationbychemical_IDschild")
    
    clusterExport(cl, "find.Overlapping.mzs")
     clusterExport(cl, "find.overlapping.single")
        clusterExport(cl, "find.Overlapping")
    
	   clusterExport(cl, "getVenn")
       
       
       s1<-seq(1,length(adduct_names))
       print("Mapping m/z to metabolites:")
       l2<-parLapply(cl,s1,Annotationbychemical_IDschild,dataA=dataA, queryadductlist=c(adduct_names),adduct_type=c("S",gradienttype),
       adduct_table=adduct_table,max.mz.diff=max.mz.diff,outloc=outloc,chemCompMZ=chemCompMZ,otherdbs=FALSE,otherinfo=FALSE)
       
       stopCluster(cl)
       
       
       levelB_res<-{};
       for(j in 1:length(l2)){
           if(length(l2[[j]])>1){
               levelB_res<-rbind(levelB_res,l2[[j]])
           }
       }
       
       
       
       rm(l2)
       
       #m2$mz<-as.numeric(as.character(m2$mz))
       #levelB_res[which(levelB_res$mz>175.117 & levelB_res$mz<175.119),1:8]
       
       levelB_res$mz<-as.numeric(as.character(levelB_res$mz))
       
       levelB_res$time<-as.numeric(as.character(levelB_res$time))
       
       levelB_res<-as.data.frame(levelB_res)
       
       
       uniq_formula<-as.character(unique(levelB_res$Formula))
       
       bad_formula<-which(is.na(uniq_formula)==TRUE)
       if(length(bad_formula)>0){
           uniq_formula<-uniq_formula[-c(bad_formula)]
       }
       
       cl<-parallel::makeCluster(num_nodes)
       
       
       clusterExport(cl, "check_golden_rules")
       clusterExport(cl, "check_element")
       #clusterExport(cl, "uniq_formula")
       #clusterExport(cl, "NOPS_check")
       
       levelB_res_check<-parLapply(cl,1:length(uniq_formula),function(j,uniq_formula,NOPS_check){
           
           curformula<-as.character(uniq_formula[j])
           return(check_golden_rules(curformula,NOPS_check=NOPS_check))
           
       },uniq_formula=uniq_formula,NOPS_check=NOPS_check)
       stopCluster(cl)
       
       #save(levelB_res_check,file="xMSannotator_levelB_check.Rda")
       levelB_res_check2<-ldply(levelB_res_check,rbind)
       
       levelB_res_check3<-levelB_res_check2[which(levelB_res_check2[,2]==1),]
       
       
       levelB_res<-levelB_res[which(levelB_res$Formula%in%as.character(levelB_res_check3[,1])),]
       
       water_adducts<-c("M+H-H2O","M+H-2H2O","M-H2O-H")
       
       water_adduct_ind<-which(levelB_res$Adduct%in%water_adducts)
       
       cl<-parallel::makeCluster(num_nodes)
       
       
       clusterExport(cl, "check_element")
       
       
       
       if(length(water_adduct_ind)>0){
           levelB_res2<-levelB_res[c(water_adduct_ind),]
           
           levelB_res<-levelB_res[-c(water_adduct_ind),]
           
           sind1<-seq(1:dim(levelB_res2)[1])
           
           levelB_res_check3<-parLapply(cl,sind1,function(j){
               
               adduct<-as.character(levelB_res2$Adduct[j])
               curformula<-as.character(levelB_res2$Formula[j])
               
               numoxygens<-check_element(curformula,"O")
               
               if(numoxygens>0){
                   bool_check<-1
               }else{
                   bool_check<-0
               }
               
               res<-cbind(curformula,bool_check)
               res<-as.data.frame(res)
               return(res)
               
               
           })
           
           levelB_res_check4<-ldply(levelB_res_check3,rbind)
           
           valid_form<-{}
           
           if(length(which(levelB_res_check4[,2]==1))>0){
               levelB_res_check4<-levelB_res_check4[which(levelB_res_check4[,2]==1),]
               
               
               valid_form<-which(levelB_res2$Formula%in%as.character(levelB_res_check4[,1]))
           }
           if(length(valid_form)>0){
               levelB_res2<-levelB_res2[valid_form,]
               levelB_res<-rbind(levelB_res,levelB_res2)
           }
           
       }
       multiresmat<-levelB_res
       rm(levelB_res)
       dupmz<-multiresmat$mz[which(duplicated(multiresmat$mz)==TRUE)]
       
       
       MatchCategory<-rep("Multiple",dim(multiresmat)[1])
       
       MatchCategory[-which(multiresmat$mz%in%dupmz)]<-"Unique"
       
       levelB_res<-cbind(MatchCategory,multiresmat)
       
       #save(levelB_res,file="xMSannotator_levelB.Rda")
       return(levelB_res)
}




Annotationbychemical_IDschild<-function(adduct_index=NA,dataA,queryadductlist=c("M+H"),adduct_type=c("S","Acetonitrile"),adduct_table,max.mz.diff=10,outloc, otherdbs=FALSE,otherinfo=FALSE,chemCompMZ){
    
    #load("~/Documents/Emory/JonesLab/Projects/xMSannotator/keggCompMZ.Rda")
    dataA<-as.data.frame(dataA)
    adduct_names<-as.character(adduct_table[,1])
    adductlist<-adduct_table[,4]
    mult_charge<-adduct_table[,3]
    num_mol<-adduct_table[,2]
    names(adductlist)<-as.character(adduct_names)
    names(mult_charge)<-as.character(adduct_names)
    names(num_mol)<-as.character(adduct_names)
    alladducts<-adduct_names
    
    
    keggCompMZ<-chemCompMZ
    
    rm(chemCompMZ)
    
    if(is.na(adduct_index)==FALSE){
        
        queryadductlist=queryadductlist[adduct_index]
    }
    
    #load("~/Documents/Emory/JonesLab/Projects/xMSannotator/keggCompMZ.Rda")
    
    alladducts<-adduct_names
    #print(queryadductlist)
    #print(alladducts)
    if(queryadductlist[1]=="positive")
    {
        queryadductlist<-adduct_names[which(adduct_table[,5]=="positive")]
        
    }else{
        if(queryadductlist[1]=="negative")
        {
            
            queryadductlist<-adduct_names[which(adduct_table[,5]=="negative")]
            
        }else{
            if(queryadductlist[1]=="all"){
                
                
                queryadductlist<-alladducts
                
                
            }else{
                if(length(which(queryadductlist%in%alladducts==FALSE))>0){
                    
                    errormsg<-paste("Adduct should be one of:",sep="")
                    for(i in alladducts){errormsg<-paste(errormsg, i,sep=" ; ")}
                    stop(errormsg, "\n\nUsage: feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"M+H\", \"M+Na\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"positive\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"negative\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"all\"), xMSannotator.outloc, numnodes=1)"
                    )
                }
                
            }
        }
    }
    adduct_table<-as.data.frame(adduct_table)
    
    adduct_table<-adduct_table[which(adduct_table$Type%in%adduct_type),] #=="S" | adduct_table$Type=="Acetonitrile",]
    
    #adduct_table<-adduct_table[which(adduct_table$Mode%in%adduct_mode),] #=="positive",]
    
    adduct_table<-adduct_table[which(adduct_table$Adduct%in%queryadductlist),] #adduct_table$Adduct=="M+H" | adduct_table$Adduct=="M+Na",]
    
    
    suppressWarnings(dir.create(outloc))
    
    setwd(outloc)
    #keggres<-KEGG.annotation(dataA=mz_search_list,queryadductlist = c("positive"),xMSannotator.outloc)
    
    
    #cur_fname<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/MaHPIC/Exp2/c18/apLCMS_with_xMSanalyzer_merged_data/apLCMS_feature_list_at_p1_U_p2cor0.7_CV100.txt"
    #dataA<-read.table(cur_fname,sep="\t",header=TRUE)
    
    mz_search_list_1<-as.data.frame(keggCompMZ[which(keggCompMZ$Adduct%in%adduct_table$Adduct),c(1,7)])
    
    
    
    mz_search_list_1<-apply(mz_search_list_1,2,as.numeric)
    
    gcur<-getVenn(dataA=dataA,name_a="Experimental",name_b="DB",dataB=mz_search_list_1,mz.thresh=max.mz.diff,time.thresh=NA,
    xMSanalyzer.outloc=outloc,alignment.tool=NA)
    
    #print("here 1")
    #print(length(gcur$common))
    #print("here 2")
    
    if(length(gcur$common)>0){
        
        mcur<-merge(keggCompMZ[which(keggCompMZ$Adduct%in%adduct_table$Adduct),],gcur$common,by.x="mz",by.y="mz.data.B")
        
        #print(mcur[1,1:4])
        #print(dataA[1,])
        
        
        mcur_2<-merge(mcur,dataA,by.x="mz.data.A",by.y="mz")
        
        mcur_3<-mcur_2[order(mcur_2$Name),]
        
        mcur_3<-mcur_3[,-c(9,10)]
        
        cnames<-colnames(mcur_3)
        cnames[1]<-"mz"
        cnames[2]<-"theoretical.mz"
        colnames(mcur_3)<-as.character(cnames)
        
        mcur_4<-as.data.frame(mcur_3)
        
        suppressWarnings(rm(keggCompMZ))
        rm(dataA)
        rm(mcur_2)
        rm(mcur_3)
        rm(mz_search_list_1)
        
        
        if(dim(mcur)[1]>1){
            
            h1<-table(mcur_4$mz) #Adduct
            
            if(length(h1)>0){
                #u1<-c(u1,which(h1<=1))
                u1<-which(h1<=1)
            }
            match_status<-rep("Multiple",dim(h1)[1])
            
            
            uniq_kegg_matches<-names(u1)
            
            
            
            match_status[u1]<-"Unique"
            
            
            
            match_status_mat<-cbind(rownames(h1),match_status)
            
            
            
            colnames(match_status_mat)<-c("mz","MatchCategory")
            match_status_mat<-as.data.frame(match_status_mat)
            
            mcur_5<-merge(match_status_mat,mcur_4,by="mz")
            
            
            rm(mcur_4)
            
            #h1<-table(mcur_5$chemical_ID,mcur_5$mz)
            
            #s2<-apply(h1,2,sum)
            
            
            mcur_5<-as.data.frame(mcur_5)
            
            #mcur_5$mz<-as.numeric(mcur_5$mz)
            
            mcur_5<-mcur_5[order(mcur_5$mz),]
            
        }else{
            MatchCategory<-"Unique"
            cnames1<-colnames(mcur_4)
            cnames1<-c(cnames1[1],"MatchCategory",cnames[-c(1)])
            mcur_5<-cbind(mcur_4[1,1],MatchCategory,mcur_4[,-c(1)])
            mcur_5<-as.data.frame(mcur_5)
            colnames(mcur_5)<-as.character(cnames1)
            
            
            #mcur_5$mz<-as.numeric(mcur_5$mz)
            
            mcur_5<-mcur_5[order(mcur_5$mz),]
            
        }
        if(otherinfo==TRUE){
            info_mat<-sapply(1:dim(mcur_5)[1],function(j){
                
                
                b1<-keggGet(paste("cpd:",mcur_5[j,1],sep=""))
                brite_inf<-paste(b1[[1]]$BRITE,collapse=";")
                path_inf<-paste(b1[[1]]$PATHWAYS,collapse=";")
                otherdb_inf<-paste(b1[[1]]$DBLINKS,collapse=";")
                r1<-c(as.character(mcur_5[j,1]),as.character(brite_inf),as.character(path_inf),as.character(otherdb_inf))
                
                return(r1)
            })
            
            
            info_mat_1<-as.data.frame(t(info_mat))
            colnames(info_mat_1)<-c("chemical_ID","BriteCategory","Pathways","ExternalLinks")
            
            
            mcur_6<-merge(info_mat_1,mcur_5,by="chemical_ID")
            
            mcur_7<-unique(mcur_6)
            
            rm(mcur_6)
            
            if(otherdbs==TRUE){
                info_mat_2<-sapply(1:dim(mcur_7)[1],function(j){
                    
                    b1<-keggLink(paste("cpd:",mcur_7[j,1],"+-e",sep=""))
                    hmdbID<-"-"
                    lipidmapsID<-"-"
                    
                    link_text<-b1[,2]
                    
                    t2<-gregexpr(pattern="hmdb:",perl=FALSE,text=link_text)
                    
                    if(length(t2)>1){
                        g_ind<-which(t2==1)
                        
                        if(length(g_ind)>0){
                            if(length(g_ind)>1){
                                for(g in g_ind){
                                    t3=t2[[g]]
                                    
                                    hmdbID<-c(hmdbID,gsub(b1[g,2],pattern="hmdb:",replacement=""))
                                }
                                if(length(g_ind)>1){hmdbID<-paste(hmdbID,collapse=";")}
                            }else{
                                
                                hmdbID<-gsub(b1[g_ind,2],pattern="hmdb:",replacement="")
                            }
                        }
                    }
                    
                    
                    t2<-gregexpr(pattern="lipidmaps:",perl=FALSE,text=link_text)
                    
                    if(length(t2)>1){
                        g_ind<-which(t2==1)
                        
                        if(length(g_ind)>0){
                            
                            if(length(g_ind)>1){
                                for(g in g_ind){
                                    t3=t2[[g]]
                                    
                                    lipidmapsID<-c(lipidmapsID,gsub(b1[g,2],pattern="lipidmaps:",replacement=""))
                                    
                                    
                                }
                                lipidmapsID<-paste(lipidmapsID,collapse=";")
                                
                            }else{lipidmapsID<-gsub(b1[g_ind,2],pattern="lipidmaps:",replacement="")}
                            
                        }
                        
                    }
                    
                    
                    return(list(keggid=as.character(mcur_7[j,1]),hmdb=hmdbID,lipidmaps=lipidmapsID))
                    c1<-c(as.character(mcur_7[j,1]),hmdbID,lipidmapsID)
                    c1<-as.data.frame(c1)
                    return(c1)
                })
                
                info_mat_3<-{}
                #for(i in 1:dim(info_mat_2)[1]){
                
                cdata<-rbind(info_mat_2[1,],info_mat_2[2,])
                cdata<-rbind(cdata,info_mat_2[3,])
                cdata<-as.data.frame(cdata)
                info_mat_3<-rbind(info_mat_3,cdata)
                
                
                #}
                
                #info_mat_3<-as.data.frame(t(info_mat_2))
                info_mat_3<-t(info_mat_3)
                colnames(info_mat_3)<-c("chemical_ID","HMDBID","LIPIDMAPS")
                
                mcur_7<-as.data.frame(mcur_7)
                
                mcur_8<-cbind(info_mat_3,mcur_7) #,by="chemical_ID")
                mcur_8<-unique(mcur_8)
                rownames(mcur_8)<-NULL
                return(mcur_8)
            }else{
                mcur_7<-as.data.frame(mcur_7)
                
                rownames(mcur_7)<-NULL
                return(mcur_7)
            }
        }else{
            mcur_5<-unique(mcur_5)
            return(mcur_5)
            
        }
        
    }else{return("no match found.")}
    #}else{return("no match found.")}
}



getSumreplicateschild<-function(curdata,alignment.tool,numreplicates,rep.num.max.missing.thresh,method="mean",missing.val)
{
    #curdata<-t(curdata)
    #write.table(curdata,file="test.txt",sep="\t",row.names=FALSE)
    numfeats=dim(curdata)[1]
    numsamp=dim(curdata)[2]
    # if(FALSE){
    resvec_1<-lapply(1:numfeats,function(r)
    {
        newrow={}
        finalmat={}
        #for(samp in seq(1,(numsamp),numreplicates))
        {
            # i=samp
            #j=i+numreplicates-1
            
            curdata_int=curdata[r,]
            
            #if(is.na(missing.val)==FALSE){
            #           check_zeros=which(curdata_int==missing.val)
            #           }else{
            check_zeros=which(is.na(curdata_int)==TRUE)
            #           	}
            na_thresh=rep.num.max.missing.thresh #round(rep.num.max.missing.thresh*numreplicates)
            
            
            if(length(check_zeros)>na_thresh)
            {
                meanval<-missing.val
            }
            else
            {
                #temporarily replace the missing intensities, set to 0 in apLCMS,
                #with mean intensity value of the corresponding replicates (with non-zero values)
                #curdata_int[check_zeros]=mean(t(curdata_int[-c(check_zeros)]))
                if(length(check_zeros)>0)
                {
                    if(method=="mean"){
                        meanval<-mean(t(curdata_int[-check_zeros]),na.rm=TRUE)
                    }else{
                        meanval<-median(t(curdata_int[-check_zeros]),na.rm=TRUE)
                    }
                }
                else
                {
                    if(method=="mean"){
                        meanval<-mean(t(curdata_int),na.rm=TRUE)
                    }else{
                        meanval<-median(t(curdata_int),na.rm=TRUE)
                    }
                }
                
            }
            newrow<-cbind(newrow,meanval)
        }
        
        
        finalmat<-rbind(finalmat, newrow)
        return(finalmat)
    })
    
    #colnames(final_set)<-colnames_data
    #rownames(final_set)=NULL
    return(resvec_1)
    
    
    
}


getSumreplicates<-function(curdata,alignment.tool,numreplicates,numcluster,rep.num.max.missing.thresh,summary.method="mean",summary.na.replacement="zeros",missing.val=0)
{
    mean_replicate_difference<-{}
    sd_range_duplicate_pairs<-{}
		  #print(alignment.tool)
          if(alignment.tool=="apLCMS")
          {
              col_end=2
          }
          else
          {
              if(alignment.tool=="XCMS")
              {
                  col_end=2
              }
              else
              {
                  stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
          }
          
          curdata_mz_rt_info=curdata[,c(1:col_end)]
          curdata=curdata[,-c(1:col_end)]
          
          
          
          cl<-parallel::makeCluster(numcluster)
          numfeats=dim(curdata)[1]
          numsamp=dim(curdata)[2]
          
          clusterEvalQ(cl, "getSumreplicateschild")
          sub_samp_list<-list()
          
          sampcount=1
          for(samp in seq(1,(numsamp),numreplicates))
          {
              i=samp
              j=i+numreplicates-1
              if(dim(curdata[,c(i:j)])[1]>0){
                  sub_samp_list[[sampcount]]=curdata[,c(i:j)]
              }
              sampcount=sampcount+1
          }
          
          avg.res<-parSapply(cl,sub_samp_list,getSumreplicateschild,alignment.tool=alignment.tool,numreplicates=numreplicates,rep.num.max.missing.thresh=rep.num.max.missing.thresh,method=summary.method,missing.val=missing.val)
          #avg.res<-getAvgreplicateschild(sub_samp_list[[1]],alignment.tool,numreplicates)
          #print("done")
          
          
          stopCluster(cl)
          
          
          
          final_set<-as.data.frame(avg.res)
          colnames_data<-colnames(curdata)
          colnames_data<-colnames_data[seq(1,(numsamp),numreplicates)]
          colnames(final_set)<-colnames_data
          rownames(final_set)=NULL
          #final_set<-cbind(curdata_mz_rt_info,final_set)
          
          final_set<-apply(final_set,2,as.numeric)
          #	write.table(final_set,file="final_Set.txt",sep="\t",row.names=FALSE)
          
          if(summary.na.replacement=="zeros"){
              
              final_set<-replace(final_set,which(is.na(final_set)==TRUE),0)
          }else{
              if(summary.na.replacement=="halfsamplemin"){
                  
                  
                  final_set<-apply(final_set,2,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
              }else{
                  
                  if(summary.na.replacement=="halfdatamin"){
                      
                      
                      min_val<-min(final_set,na.rm=TRUE)*0.5
                      final_set<-replace(final_set,which(is.na(final_set)==TRUE),min_val)
                      
                      #data_m<-apply(data_m,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
                  }else{
                      if(summary.na.replacement=="halffeaturemin"){
                          
                          
                          final_set<-apply(final_set,1,function(x){naind<-which(is.na(x)==TRUE); if(length(naind)>0){x[naind]<-min(x,na.rm=TRUE)/2}; return(x)})
                          final_set<-t(final_set)
                      }
                  }
              }
              
              
          }
          
          return(final_set)
}



data_summarize<-function(Xmat=NA,Ymat=NA,feature_table_file,parentoutput_dir,class_labels_file,num_replicates=3,summarize.replicates=TRUE,summary.method="mean",missing.val=0,rep.num.max.missing.thresh=1,
summary.na.replacement="zeros",fileheader="RAW"){

    options(warn=-1)
    
    #read file; First row is column headers
    if(is.na(Xmat==TRUE)){
        data_matrix<-read.table(feature_table_file,sep="\t",header=TRUE)
    }else{
        data_matrix<-Xmat
        #rm(Xmat)
    }
    #print("signal filter threshold ")
    #print(group.missing.thresh)
    
    print("missing val is")
    print(missing.val)
    
   
    
   
    
    #use only unique records
    data_matrix<-unique(data_matrix)
    
    if(is.na(missing.val)==FALSE){
        
        print("Replacing missing values with NAs.")
        data_matrix<-replace(as.matrix(data_matrix),which(data_matrix==missing.val),NA)
        print(head(data_matrix[1:3,1:4]))
    }
    
    # print(data_matrix[1:10,1:5])
    
    #print("dim of original data matrix")
    #print(dim(data_matrix))
    data_matrix_orig<-data_matrix
    
    
    snames<-colnames(data_matrix)
    
    
    
    
    
    dir.create(parentoutput_dir,showWarnings=FALSE)
    #parentoutput_dir<-paste(parentoutput_dir,"/Stage1/",sep="")
    
    dir.create(parentoutput_dir,showWarnings=FALSE)
    fheader="transformed_log2fc_threshold_"
    setwd(parentoutput_dir)
    
    data_m<-as.matrix(data_matrix[,-c(1:2)])
    
    if(is.na(Xmat)==FALSE){
        
        #   write.table(Xmat,file="organized_featuretableA.txt",sep="\t",row.names=TRUE)
        
        
    }
    
    if(is.na(Ymat)==FALSE){
        # write.table(Ymat,file="organized_classlabelsA.txt",sep="\t",row.names=FALSE)
        
    }
    
    #Step 2) Average replicates
    if(summarize.replicates==TRUE)
    {
        if(num_replicates>1)
        {
            
            data_m<-getSumreplicates(data_matrix,alignment.tool="apLCMS",numreplicates=num_replicates,numcluster=10,rep.num.max.missing.thresh=rep.num.max.missing.thresh,summary.method=summary.method,summary.na.replacement, missing.val=missing.val)
            
            
            data_m<-replace(data_m,which(is.na(data_m)==TRUE),missing.val)
            
            if(summary.method=="mean"){
                print("Replicate averaging done")
                filename<-paste(fileheader,"_mzcalibrated_untargeted_averaged_featuretable.txt",sep="")
            }else{
                if(summary.method=="median"){
                    print("Replicate median summarization done")
                    filename<-paste(fileheader,"_mzcalibrated_untargeted_mediansummarized_featuretable.txt",sep="")
                }
                
            }
            
            data_m_prenorm<-cbind(data_matrix[,c(1:2)],data_m)
            
            write.table(data_m_prenorm, file=filename,sep="\t",row.names=FALSE)
            
            data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
            #num_samps_group[[1]]=(1/num_replicates)*num_samps_group[[1]]
            #num_samps_group[[2]]=(1/num_replicates)*num_samps_group[[2]]
        }
    }
    
    
    
    data_matrix<-cbind(data_matrix[,c(1:2)],data_m)
    
    return(data_matrix)
}


apLCMS.target.EIC<-function(aligned,finalfeatmat,refMZ.mz.diff=10,refMZ.time.diff=NA,target.mz.list = NA,subs=NA,
colors=NA,apLCMS.outloc,runval,presval,minexp,cdfloc,cvvec=NA){

	if(is.na(aligned)==TRUE){

		setwd(apLCMS.outloc)
		fname<-paste(apLCMS.outloc, "/apLCMS_aligned", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.Rda", sep="")


		load(fname)
	}
	if(is.na(aligned)==TRUE){
		stop("No aligned object found while generating EICs")
	}

	finalfeatmat<-as.data.frame(finalfeatmat[,c(1:2)])
	if(is.na(target.mz.list)==TRUE)
	{

		 Name<-paste("mz",seq(1,dim(finalfeatmat)[1]),sep="")
                        stddata<-cbind(finalfeatmat[,c(1:2)],Name)
			stddata<-as.data.frame(stddata)
			
stddata$mz<-as.numeric(as.character(stddata$mz))
stddata$time<-as.numeric(as.character(stddata$time))
			print(head(stddata))
			stddata1<-stddata[,c(1:2)]
				stddata1<-as.data.frame(stddata1)
				stddata1$mz<-as.numeric(as.character(stddata1$mz))				
				#print(stddata1$time[1]+1)
				stddata1$time<-as.numeric(as.character(stddata1$time))
				#stddata1<-apply(stddata1,2,as.numeric)
				stddata1<-as.data.frame(stddata1)
				print(head(stddata1))
				print(print(stddata1$time[1]+1))
				stddata1<-as.data.frame(stddata1)
				print(head(finalfeatmat[1:30,]))	
		
	}else{

		stddata<-target.mz.list

		stddata<-as.data.frame(stddata)
                        print(head(stddata))
			stddata1<-stddata[,c(1:2)]
				stddata1<-as.data.frame(stddata1)
				stddata1$mz<-as.numeric(as.character(stddata1$mz))				
				#print(stddata1$time[1]+1)
			refMZ.time.diff=NA
	}


	{
			
			#eic_fname<-paste(apLCMS.outloc,"/EICrun",runval,"pres",presval,"target.pdf",sep="")
			#pdf(eic_fname)
				
	#		overlapres5ppm<-getVenn(dataA=finalfeatmat, name_a=paste("Expdata",sep=""), dataB=stddata1, name_b="Target", mz.thresh = refMZ.mz.diff, time.thresh=refMZ.time.diff,alignment.tool="apLCMS",xMSanalyzer.outloc=apLCMS.outloc,plotvenn=FALSE)

			
			overlapres5ppm<-find.Overlapping.mzs(dataA=finalfeatmat, dataB=stddata1, mz.thresh=refMZ.mz.diff, time.thresh=refMZ.time.diff, alignment.tool=NA)
		
				print("Common")	
				print(length(unique(overlapres5ppm$index.A)))
				if(length(unique(overlapres5ppm$index.A))>0)
				{
					finalfeatmat<-as.data.frame(finalfeatmat)

					stddata<-as.data.frame(stddata)
					overlap_res<-overlapres5ppm
                                        overlap_res<-as.data.frame(overlap_res)
                                        #print(overlap_res)
					dup_mz_ind<-which(duplicated(overlapres5ppm$index.A)==TRUE)
					if(length(dup_mz_ind)>0){
						overlap_res<-overlap_res[-c(dup_mz_ind),]
					}
					finalfeatmat<-as.data.frame(finalfeatmat)
					time.list=finalfeatmat$time[c((overlap_res$index.A))]
					mz.list=finalfeatmat$mz[c((overlap_res$index.A))]
					
									
					if(is.na(cvvec)==FALSE){
						cvvec=cvvec[c((overlap_res$index.A))]
					}
					chem.names<-stddata$Name[c((overlap_res$index.B))]
					print(dim(aligned$aligned.ftrs))
					print("cdfloc")
					print(cdfloc)
					#print("mzlist")	
					#print(mz.list)	

				#c1<-custom.EIC.plot(aligned, rows = c((overlap_res$index.A)), colors = colors, transform = "none", 
					#	subset = subs, mz.list=mz.list, time.list=time.list,chem.names=chem.names,min.run=runval,min.pres=presval, max.spline.time.points = 1000,cdfloc=cdfloc) 
				overlap_ind<-unique(overlapres5ppm$index.A)			
				min_rt=min(finalfeatmat$time,na.rm=TRUE)
				max_rt=max(finalfeatmat$time,na.rm=TRUE)	
	apLCMS.EIC.plot(aligned, rows = c((overlap_ind)), colors = colors, transform = "none",
    subset =NA, mz.list=mz.list, time.list=time.list,minrt=min_rt, maxrt=max_rt,chem.names=chem.names,min.run=runval,min.pres=presval,
    max.spline.time.points = 1000,rawprofileloc=cdfloc,cvvec=cvvec)
				}
			}
	

}

call.getEIC<-function(filenum,filenames,targetmz_list,ppm_error,mzXML.loc,sample_names=NA,deltatime=10,profstep=0.1){

setwd(mzXML.loc)
#xraw<-xcmsRaw(filenames[filenum],includeMSn=FALSE,profstep=profstep)

xraw <- xcmsRaw(filenames[filenum],profmethod="bin",profstep=profstep)

cur_filename=gsub(filenames[filenum],pattern=".mzXML|.cdf|.CDF",replacement="")

sample_names[,1]=gsub(sample_names[,1],pattern=".mzXML|.cdf|.CDF",replacement="")

sample_name=sample_names[which(sample_names[,1]==cur_filename),2]

if(ncol(targetmz_list)>=2){
	data_list<-lapply(1:nrow(targetmz_list),function(t1){

	targetmz<-targetmz_list[t1,1]
	targettime<-targetmz_list[t1,2]

	peak_mat<-get.rawEIC(object=xraw,targetmz=targetmz,targettime=targettime,deltappm=ppm_error,deltatime=deltatime)

	base_peak_height<-max(peak_mat[,2],na.rm=TRUE)
	
	if(is.na(sample_names)==FALSE){
		main_text=paste("File: ",filenames[filenum],"(Sample:",sample_name,")\n mz:",targetmz, " time: ",targettime,"\n Base peak height: ",base_peak_height,sep="")
	}else{
		main_text=paste("File: ",filenames[filenum],"\n mz:",targetmz, " time: ",targettime,"\n Base peak height: ",base_peak_height,sep="")
	}
	
	plot(peak_mat,xlab="Time",ylab="Intensity",cex.axis=0.7,main=main_text,cex.main=0.7,type="l")
		
						#lines(peak_mat[,1],peak_mat[,2])
						#polygon(c(lb,peak_mat[,1],ub),c(0,peak_mat[,2],0),col="brown")
					#polygon(peak_mat[,1],peak_mat[,2],col="brown")

	})
	
}else{
	data_list<-lapply(1:length(targetmz_list),function(t1){

	targetmz<-targetmz_list[t1,1]


	peak_mat<-get.rawEIC(object=xraw,targetmz=targetmz,targettime=NA,deltappm=ppm_error)

	
	
	if(is.na(sample_names)==FALSE){
		main_text=paste("File: ",filenames[filenum],"(Sample:",sample_name,")\n mz:",sep="")
	}else{
		main_text=paste("File: ",filenames[filenum],"\n mz:",targetmz,sep="")
	}

	plot(peak_mat,xlab="Time",ylab="Intensity",cex.axis=0.7,main=main_text,cex.main=0.7,type="l")
		
						#lines(peak_mat[,1],peak_mat[,2])
						
						#polygon(peak_mat[,1],peak_mat[,2],col="brown")

	})

}
return(data_list)
}


get.binnedEIC <- function(object,targetmz=NA,targettime=NA,ppm_error=5,deltatime=300,step=0.001){
   
	amu_error=10^(-6)*(ppm_error*targetmz)
	mzmin=targetmz-amu_error
	mzmax=targetmz+amu_error
	
		
	if(is.na(targettime)==TRUE){
		timemin=targettime-deltatime
		timemax=targettime+deltatime	
		res<-getEIC(object,mzrange=c(mzmin,mzmax),step=step) #rawEIC(object,targetmz=range(object@env$mz))
		
		intvec=res@eic$xcmsRaw[[1]][,2]

		timevec=res@eic$xcmsRaw[[1]][,1]
		peak_mat<-cbind(timevec,intvec)
	}else{

		

		res<-getEIC(object,mzrange=c(mzmin,mzmax),step=step) #rawEIC(object,mzrange=c(mzmin,mzmax))
		intvec=res@eic$xcmsRaw[[1]][,2]

		timevec=res@eic$xcmsRaw[[1]][,1]
		
		timemax1<-max(timevec,na.rm=TRUE)-10
		
		timemax<-min(c(timemax1,timemax),na.rm=TRUE)
		peak_mat<-cbind(timevec,intvec)

		#print(head(peak_mat))
		#intvec=intvec[which(timevec>timemin & timevec<timemax)]
		#timevec<-timevec[which(timevec>timemin & timevec<timemax)]
		if(is.na(targettime)==FALSE){
		peak_mat<-peak_mat[which(timevec>timemin & timevec<timemax),]
		}
	}
     return(peak_mat)
}

get.rawEIC <- function(file=NA,rtcor=NULL,targetmz=NA,targettime=NA,deltappm=5,deltatime=10,profstep=0.1,object=NA) {
     
      
        if(is.na(file)==FALSE){
	object <- xcmsRaw(file,profmethod="bin",profstep=profstep)
        }
	
	deltaamu=10^(-6)*(deltappm*targetmz)
	
	if(is.na(targetmz)==TRUE){
		
		intvec<-rawEIC(object,mzrange=range(object@env$mz))$intensity
		peak_mat<-cbind(object@scantime,intvec)
	}else{

		mzmin=targetmz-deltaamu
		mzmax=targetmz+deltaamu

		if(is.na(targettime)==FALSE){
		
		if(is.na(deltatime)==TRUE){
			deltatime=max(object@scantime,na.rm=TRUE)
		}
		timemin=targettime-deltatime
		timemax=targettime+deltatime
		
		intvec<-rawEIC(object,mzrange=c(mzmin,mzmax))$intensity
			peak_mat<-cbind(object@scantime,intvec)
			
		#intvec<-rawEIC(object,mzrange=c(mzmin,mzmax),rtrange=c(timemin,timemax))$intensity
		
		peak_mat<-peak_mat[which(peak_mat[,1]>timemin & peak_mat[,1]<timemax),]
		}else{
			
			intvec<-rawEIC(object,mzrange=c(mzmin,mzmax))$intensity
			peak_mat<-cbind(object@scantime,intvec)
		}
	
	}
	
	
	colnames(peak_mat)<-c("time","intensity")
     return(peak_mat)
}


getTIC <- function(file=NA,rtcor=NULL,targetmz=NA,targettime=NA,deltappm=5,deltatime=10,profstep=0.1,object=NA) {
     
      
        if(is.na(file)==FALSE){
	object <- xcmsRaw(file,profmethod="bin",profstep=profstep)
        }
	
	deltaamu=10^(-6)*(deltappm*targetmz)
	
	if(is.na(targetmz)==TRUE){
		
		intvec<-rawEIC(object,mzrange=range(object@env$mz))$intensity
		peak_mat<-cbind(object@scantime,intvec)
	}else{

		mzmin=targetmz-deltaamu
		mzmax=targetmz+deltaamu

		if(is.na(targettime)==FALSE){
		
		if(is.na(deltatime)==TRUE){
			deltatime=max(object@scantime,na.rm=TRUE)
		}
		timemin=targettime-deltatime
		timemax=targettime+deltatime
		
		intvec<-rawEIC(object,mzrange=c(mzmin,mzmax))$intensity
			peak_mat<-cbind(object@scantime,intvec)
			
		#intvec<-rawEIC(object,mzrange=c(mzmin,mzmax),rtrange=c(timemin,timemax))$intensity
		
		peak_mat<-peak_mat[which(peak_mat[,1]>timemin & peak_mat[,1]<timemax),]
		}else{
			
			intvec<-rawEIC(object,mzrange=c(mzmin,mzmax))$intensity
			peak_mat<-cbind(object@scantime,intvec)
		}
	
	}
	
	
	colnames(peak_mat)<-c("time","intensity")
     return(peak_mat)
}

getEICmultiple<- function(targetindex,files=NA,rtcor=NULL,targetmz_list,deltappm=5,deltatime=10,profstep=0.1,xcms.raw=NA,sample_names) {

	if(ncol(targetmz_list)==2){
	

	targetmz<-targetmz_list[targetindex,1]
	targettime<-targetmz_list[targetindex,2]

	}else{
	
		targetmz=targetmz_list[targetindex]
		targettime<-NA
	}
	check_class<-class(xcms.raw)

	

	TIC<-lapply(1:length(files),function(i){
	
		
		if(check_class=="logical"){
			xraw<-xcmsRaw(files[i],includeMSn=FALSE,profstep=profstep)
		}else{
			xraw<-xcms.raw[[i]]
		}
		return(getTIC(file=NA,rtcor=rtcor,targetmz=targetmz,targettime=targettime,deltappm=deltappm,deltatime=deltatime,profstep=profstep,object=xraw))
     
        })
	
	if(is.na(targettime)==TRUE){
	
		 mainlab=paste("Extracted Ion Chromatograms\n mz: ",targetmz,sep="")
	}else{
		 mainlab=paste("Extracted Ion Chromatograms\n mz: ",targetmz,"; time: ",targettime,sep="")
	}
	N=length(TIC)
     cols <- rainbow(N)
      lty = 1:N
      pch = 1:N
      xlim = range(sapply(TIC, function(x) range(x[,1])))
      ylim = range(sapply(TIC, function(x) range(x[,2])))
      plot(0,0, type="n", xlim = xlim, ylim = ylim, main =mainlab, xlab = "Retention Time (s)", ylab = "Intensity")
      sample_name_vec={}
      for (i in 1:N) {

	if(is.na(sample_names)==FALSE){
		cur_filename=gsub(files[i],pattern=".mzXML|.cdf|.CDF",replacement="")

		sample_names[,1]=gsub(sample_names[,1],pattern=".mzXML|.cdf|.CDF",replacement="")

		sample_name_vec=c(sample_name_vec,as.character(sample_names[which(sample_names[,1]==cur_filename),2]))
		}else{
			sample_name_vec=files
		}

      tic <- TIC[[i]]
      points(tic[,1], tic[,2], col = cols[i], pch = pch[i], type="l")
      
      }
      print(sample_name_vec)
      legend("topright",paste(basename(sample_name_vec)), col = cols, lty = lty, pch = pch,cex=0.5)
}

getTICplots <- function(cdfloc, filepattern=".cdf|.mzxml|mXML",mzlist=NA,rtlist=NA,deltappm=5,deltatime=10,subs=NA,profstep=0.1) {
  

	
    file_list <- list.files(cdfloc, pattern = filepattern,full.names = TRUE,ignore.case=TRUE)
  
    file_list1 <- list.files(cdfloc, pattern = filepattern,full.names = FALSE,ignore.case=TRUE) 

if(is.na(subs)==FALSE){

	file_list<-file_list[subs]
	file_list1<-file_list1[subs]
}
  N <- length(file_list)
  TIC <- vector("list",N)
  

  pdfname="TICsingle.pdf"

        pdf(pdfname)
	par(mfrow=c(3,3))
  #f1<-lapply(1:N,function(i){
for(i in 1:N){  
    cat(file_list[i],"\n")
      

      rtcor <- NULL
      
      TIC[[i]] <- getTIC(file_list[i],rtcor=rtcor,mzlist,rtlist,deltappm,deltatime,profstep)
      mainlab1<-paste("TIC for ",file_list1[i],sep=" ")
	if(length(TIC[[i]][,1])>0){
      plot( TIC[[i]][,1], TIC[[i]][,2],xlab="Retention Time",ylab="TIC",main=mainlab1)
 	}
 }
	dev.off()
   pdfname="TICall.pdf"
  pdf(pdfname,w=16,h=10)
      cols <- rainbow(N)
      lty = 1:N
      pch = 1:N
      xlim = range(sapply(TIC, function(x) range(x[,1])))
      ylim = range(sapply(TIC, function(x) range(x[,2])))
      plot(0,0, type="n", xlim = xlim, ylim = ylim, main = "Total Ion Chromatograms", xlab = "Retention Time", ylab = "TIC")
      for (i in 1:N) {

      tic <- TIC[[i]]
      points(tic[,1], tic[,2], col = cols[i], pch = pch[i], type="l")
      }
      legend("topright",paste(basename(file_list)), col = cols, lty = lty, pch = pch)
  dev.off()
  
  invisible(TIC)
}

##
##  overlay TIC from all files in current folder or from xcmsSet, create pdf
##
getTICmultiple <- function(cdfloc, filepattern=".cdf|.mzxml|.mxml",mzlist=NA,rtlist=NA,deltappm=5,deltatime=10,file_list=NA,profstep=0.1) {
  
  if(is.na(file_list)==TRUE){
    file_list <- list.files(cdfloc, pattern = filepattern,full.names = TRUE,ignore.case=TRUE)
  }

  N <- length(file_list)
  TIC <- vector("list",N)

  for (i in 1:N) {
      cat(file_list[i],"\n")
      

      rtcor <- NULL
      TIC[[i]] <- getTIC(files[i],rtcor=rtcor,mzlist,rtlist,deltappm,deltatime,profstep)
  }

  pdf(pdfname,w=16,h=10)
      cols <- rainbow(N)
      lty = 1:N
      pch = 1:N
      xlim = range(sapply(TIC, function(x) range(x[,1])))
      ylim = range(sapply(TIC, function(x) range(x[,2])))
      plot(0, 0, type="n", xlim = xlim, ylim = ylim, main = "Total Ion Chromatograms", xlab = "Retention Time", ylab = "TIC")
      for (i in 1:N) {
      tic <- TIC[[i]]
      points(tic[,1], tic[,2], col = cols[i], pch = pch[i], type="l")
      }
      legend("topright",paste(basename(files)), col = cols, lty = lty, pch = pch)
  dev.off()
  
  invisible(TIC)
}


##################################################################################
#Function Name: apLCMS.align 
#Description: Call apLCMS function, cdf.to.ftr, at different parameter settings
##cdfloc: The folder where all CDF files to be processed are located. For example "C:/CDF/"
##apLCMS.outloc: The folder where alignment output will be written. For example "C:/CDFoutput/"
##min.run.list: List of values for min.run parameter, eg: c(3,6) would run the cdf.to.ftr function at min.run=3 and min.run=6
##min.pres.list: List of values min.pres, eg: c(0.3,0.8) would run the cdf.to.ftr function at min.run=3 and min.run=6
##minexp: If a feature is to be included in the final feature table, it must be present in at least this number of samples, eg: 2
##subs: If not all the CDF files in the folder are to be processed, the user can define a subset using this parameter. For example, subs=15:30, or subs=c(2,4,6,8)
##run.order.file: Name of a tab-delimited file that includes sample names sorted by the order in which they were run (sample names must match the CDF file names)


apLCMS.align<-function(cdfloc, apLCMS.outloc,min.run.list=c(3,4), min.pres.list=c(0.8,0.5), minexp=2, mztol=2.5e-6, alignmztol=10e-6, alignchrtol=45,
numnodes=2, run.order.file=NA,subs=NA,filepattern=".cdf",apLCMSmode="untargeted",known_table,match_tol_ppm=5,refMZ.mz.diff=10,
refMZ.time.diff=NA,target.mz.list = NA,plotEICs="target",peak.score.thresh=0,num_replicates=3,baseline.correct.noise.percentile = 0.25, 

    shape.model = "bi-Gaussian", baseline.correct = NA, peak.estim.method = "moment", 

    min.bw = NA, max.bw = NA, sd.cut = c(0.125, 15), sigma.ratio.lim = c(0.33, 

        3), max.align.mz.diff = 0.01, pre.process = FALSE, recover.mz.range = NA, 

    recover.chr.range = NA, use.observed.range = TRUE, recover.min.count = 3, new_feature_min_count=4, reference_sample_index=NA,
    cvthresh=100,xMSanalyzer.outloc=NA,component.eliminate = 0.01, moment.power = 2)
{
        setwd(cdfloc)
	  cdf.files=list.files(cdfloc,filepattern,ignore.case=TRUE)
	  cdf.files=tolower(cdf.files)
	  dir.create(apLCMS.outloc,showWarnings=FALSE)
	aligned_data_list=new("list")
	filenames_list=new("list")
	peak_score_list=new("list")
	
	pcount=1
	
	
	
	
        if(is.na(subs[1])==FALSE)
        {
                numsamp=length(subs)
		    cdf.files=cdf.files[subs]
		    
        }
        else
        {
                
                numsamp=length(cdf.files)
        }
        
        if(length(min.run.list)!=length(min.pres.list)){
        	
        	stop("Vectors min.run.list and min.pres.list should be of the same length. eg: min.run.list=c(3,3) and min.pres.list=c(0.3,0.8)")
        }
        for(r in 1:length(min.run.list))
        {
		    runval=min.run.list[r]
		    p=r
              
                        features<-new("list")
                        
                        presval=min.pres.list[p]

			#dev.off()
			
			
			
			print(minexp)
			
			reference_sample_index<-NA

			sys_name<-Sys.info()['sysname']
			
			setwd(apLCMS.outloc)
			
			fileload<-{}
			
			  fname_check<-paste("apLCMS",apLCMSmode,"_aligned", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.Rda", sep="")
			  if(length(which(fname_check%in%dir()))>0){
			  
					setwd(apLCMS.outloc)
						print("loading file")
						print(fname_check)
					fileload<-try(load(fname_check),silent=TRUE)
			  }
			  
				print(fileload)
				print(length(fileload))
				
			   if(is(fileload,"try-error") | length(fileload)<1)
			  {
                        							if(apLCMSmode=="untargeted"){
										if(sys_name=="Windows"){
										
										if(is.na(reference_sample_index)==FALSE){
										
										par(mfrow=c(2,2))
			fname<-paste("Rplots",runval,presval,".pdf",sep="")
			pdf(fname)
											aligned<-cdf.to.ftr(cdfloc,subs=subs,min.exp=minexp,min.run=runval,min.pres=presval,mz.tol=mztol, align.mz.tol=alignmztol, align.chr.tol=alignchrtol,
		n.nodes=numnodes,file.pattern=filepattern,  baseline.correct.noise.percentile = baseline.correct.noise.percentile, 

    shape.model = shape.model, baseline.correct = baseline.correct, peak.estim.method = peak.estim.method, 

    min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, sigma.ratio.lim = sigma.ratio.lim, 

    max.align.mz.diff = max.align.mz.diff, pre.process = pre.process, recover.mz.range = recover.mz.range, 

    recover.chr.range = recover.chr.range, use.observed.range = use.observed.range, recover.min.count = recover.min.count,component.eliminate = component.eliminate, moment.power = moment.power,intensity.weighted=intensity.weighted)
    
	dev.off()
					reference_sample<-aligned$final.ftrs
					medianInt<-apply(reference_sample[,-c(1:4)],1,median)
					reference_sample<-cbind(reference_sample[,c(1:2)],medianInt)
					colnames(reference_sample)<-c("mz","time","medianInt")
					
										
										}else{
										
											reference_sample<-NA
										}
										
										par(mfrow=c(2,2))
			fname<-paste("Rplots",runval,presval,".pdf",sep="")
			pdf(fname)
		aligned<-cdf.to.ftr(cdfloc,subs=subs,min.exp=minexp,min.run=runval,min.pres=presval,mz.tol=mztol, align.mz.tol=alignmztol, align.chr.tol=alignchrtol,
		n.nodes=numnodes,file.pattern=filepattern,  baseline.correct.noise.percentile = baseline.correct.noise.percentile, 

    shape.model = shape.model, baseline.correct = baseline.correct, peak.estim.method = peak.estim.method, 

    min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, sigma.ratio.lim = sigma.ratio.lim, 

    max.align.mz.diff = max.align.mz.diff, pre.process = pre.process, recover.mz.range = recover.mz.range, 

    recover.chr.range = recover.chr.range, use.observed.range = use.observed.range, recover.min.count = recover.min.count,component.eliminate = component.eliminate, moment.power = moment.power)
    dev.off()
										
							}else{
							
							
							if(is.na(reference_sample_index)==FALSE){
										
										par(mfrow=c(2,2))
			fname<-paste("Rplots",runval,presval,".pdf",sep="")
			pdf(fname)
											aligned<-cdf.to.ftr(cdfloc,subs=subs,min.exp=minexp,min.run=runval,min.pres=presval,mz.tol=mztol, align.mz.tol=alignmztol, align.chr.tol=alignchrtol,
		n.nodes=numnodes,file.pattern=filepattern,  baseline.correct.noise.percentile = baseline.correct.noise.percentile, 

    shape.model = shape.model, baseline.correct = baseline.correct, peak.estim.method = peak.estim.method, 

    min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, sigma.ratio.lim = sigma.ratio.lim, 

    max.align.mz.diff = max.align.mz.diff, pre.process = pre.process, recover.mz.range = recover.mz.range, 

    recover.chr.range = recover.chr.range, use.observed.range = use.observed.range, recover.min.count = recover.min.count,component.eliminate = component.eliminate, moment.power = moment.power,intensity.weighted=intensity.weighted)
    
    dev.off()
					reference_sample<-aligned$final.ftrs
					medianInt<-apply(reference_sample[,-c(1:4)],1,median)
					reference_sample<-cbind(reference_sample[,c(1:2)],medianInt)
					colnames(reference_sample)<-c("mz","time","medianInt")
					
										
										}else{
										
											reference_sample<-NA
										}
										
										par(mfrow=c(2,2))
			fname<-paste("Rplots",runval,presval,".pdf",sep="")
			pdf(fname)
		aligned<-cdf.to.ftr(cdfloc,subs=subs,min.exp=minexp,min.run=runval,min.pres=presval,mz.tol=mztol, align.mz.tol=alignmztol, align.chr.tol=alignchrtol,
		n.nodes=numnodes,file.pattern=filepattern,  baseline.correct.noise.percentile = baseline.correct.noise.percentile, 

    shape.model = shape.model, baseline.correct = baseline.correct, peak.estim.method = peak.estim.method, 

    min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, sigma.ratio.lim = sigma.ratio.lim, 

    max.align.mz.diff = max.align.mz.diff, pre.process = pre.process, recover.mz.range = recover.mz.range, 

    recover.chr.range = recover.chr.range, use.observed.range = use.observed.range, recover.min.count = recover.min.count,component.eliminate = component.eliminate, moment.power = moment.power,intensity.weighted=intensity.weighted)
							dev.off()
							
							
						}
							
							}else{
								if(apLCMSmode=="hybrid"){
									
									par(mfrow=c(2,2))
			fname<-paste("Rplots",runval,presval,".pdf",sep="")
			pdf(fname)
aligned<-semi.sup(folder=cdfloc, known.table = known_table, match.tol.ppm = match_tol_ppm, subs=subs,min.exp=minexp,
min.run=runval,min.pres=presval,mz.tol=mztol, align.mz.tol=alignmztol, align.chr.tol=alignchrtol, n.nodes=numnodes,
file.pattern=filepattern,baseline.correct.noise.percentile = baseline.correct.noise.percentile, 

    shape.model = shape.model, baseline.correct = baseline.correct, peak.estim.method = peak.estim.method, 

    min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, sigma.ratio.lim = sigma.ratio.lim, 

    max.align.mz.diff = max.align.mz.diff, pre.process = pre.process, recover.mz.range = recover.mz.range, 

    recover.chr.range = recover.chr.range, use.observed.range = use.observed.range, recover.min.count = recover.min.count,
    component.eliminate = component.eliminate, moment.power = moment.power,new.feature.min.count=new_feature_min_count,intensity.weighted=intensity.weighted,BIC.factor=BIC.factor)

dev.off()
									}
								
								}
				}
								print(names(aligned))
			finalfeatmat=aligned$final.ftrs
		
			fname<-paste(apLCMS.outloc, "/apLCMS",apLCMSmode,"_aligned", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.Rda", sep="")
            #save(aligned,file=fname)
			finalfeatmat<-as.data.frame(finalfeatmat)
					
			fname<-paste(apLCMS.outloc, "/apLCMS",apLCMSmode,"_aligned", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.txt", sep="")
			write.table(finalfeatmat,fname,sep="\t",row.names=F)			
						
	

			rows_ind<-seq(1,dim(finalfeatmat)[1])	
				
			
				#pdf(eic_fname)
				time.list=finalfeatmat$time
				mz.list=finalfeatmat$mz
				chem.names<-paste("mz",seq(1,length(rows_ind)),sep="")
					
						
    
   
	
	
	setwd(apLCMS.outloc)
	
	
	
	if(is.na(peak.score.thresh)==FALSE)
	{
		
		print("getting peakscore")
		print(runval)
		print(presval)
		
		    peak_score_vec=rep(-1,nrow(finalfeatmat))
		    feat.eval.result=evaluate.Features(finalfeatmat, numreplicates=num_replicates,min.samp.percent=0.6,alignment.tool="apLCMS",impute.bool=impute.bool,
						peak_scores=peak_score_vec,numnodes=numnodes)
				
			cvvec=feat.eval.result$median
			
		  fname_check<-paste("apLCMS",apLCMSmode,"_peakscore", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.Rda", sep="")
				  if(length(which(fname_check%in%dir()))>0){
				  
						print("loading file")
						print(fname_check)
						load(fname_check)
				  }else{
	 

						PeakScore<-apLCMS.get.peakscore(aligned, rows = rows_ind, colors = NA, transform = "none",
					    subset = NA, mz.list=NA,time.list=NA,minrt=NA, maxrt=NA, min.run=runval,min.pres=presval, max.spline.time.points = 1000,
					    chem.names=NA,rawprofileloc=cdfloc,cvvec=cvvec,plotEIC=FALSE,numnodes=numnodes,cvthresh=cvthresh)
					    
	  
	    

				}
	    fname<-paste(apLCMS.outloc, "/apLCMS",apLCMSmode,"_peakscore", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.Rda", sep="")
	    save(PeakScore,file=fname)
	    print("done with peakscores")
	  
				  if(is(PeakScore,"try-error")){
				  
						#PeakScore<-rep(1,dim(finalfeatmat)[1])
						stop("Error computing peak scores")
				  }

   }else{
	   PeakScore<-rep(1,dim(finalfeatmat)[1])
	    fname<-paste(apLCMS.outloc, "/apLCMS",apLCMSmode,"_peakscore", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.Rda", sep="")
	    save(PeakScore,file=fname)
	    
   }
			

					                         
                        
	
			
            
				num_samples<-dim(finalfeatmat)[2]-4	


					

				
				
				
                        cnames<-colnames(finalfeatmat[,-c(1:4)])
                        cnames<-tolower(cnames)
                        cnames<-gsub(".cdf", "", cnames)

                        if(is.na(run.order.file)==FALSE)
                        {
                                fileorder=read.table(run.order.file, header=FALSE)
                                fileorder=apply(fileorder,1,tolower)
				cnames=tolower(cnames)
                                ordlist=sapply(1:length(fileorder),function(i){which(cnames==fileorder[i])})
                                ordlist=unlist(ordlist)
                                ordlist=ordlist+4
                                finalfeatmat=finalfeatmat[,c(1:4,ordlist)]
                        }

				fname<-paste(apLCMS.outloc, "/apLCMS",apLCMSmode,"_aligned", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.txt", sep="")
				write.table(finalfeatmat,fname,sep="\t",row.names=F)

				if(is.na(xMSanalyzer.outloc)==FALSE){
					setwd(xMSanalyzer.outloc)
				}
				
				
			getTICplots(cdfloc, filepattern=filepattern,mzlist=NA,rtlist=NA,deltaamu=1,deltatime=10)
		
			setwd(apLCMS.outloc)
		if(is.na(subs)==TRUE){
                                        
					rand_set<-{}
				
					uniq_samp_ind<-seq(1,num_samples,num_replicates)
					num_EIC<-min(30,num_samples)
					num_EIC<-(num_EIC)*(1/num_replicates)
					samp_ind<-sample(uniq_samp_ind,size=num_EIC)
					samp_ind<-samp_ind[order(samp_ind)]

					
					samp_ind1<-samp_ind+1
					samp_ind2<-{}        
					if(num_replicates>2){samp_ind2<-samp_ind+2}

					samp_ind<-c(samp_ind,samp_ind1,samp_ind2)
					samp_ind<-samp_ind[order(samp_ind)]

					rand_set<-c(rand_set,samp_ind)

					rand_set<-rand_set[order(rand_set)]
			}else{
				rand_set<-subs
			}
				
				print("Plotting EICs")
				
				

   
	     
				if(plotEICs=="target"){
					
					if(is.na(target.mz.list)==FALSE){
						eic_fname<-paste(apLCMS.outloc,"/EICrun",runval,"pres",presval,"prefiltering_target.pdf",sep="")

						pdf(eic_fname)				
	
						#try(apLCMS.target.EIC(aligned,finalfeatmat,refMZ.mz.diff,refMZ.time.diff,target.mz.list,subs=rand_set,
						#colors,apLCMS.outloc,runval,presval,minexp,cdfloc),silent=TRUE)
						
	try(apLCMS.target.EIC(aligned=aligned,finalfeatmat=aligned$aligned.ftrs[,c(1:2)],refMZ.mz.diff=refMZ.mz.diff,refMZ.time.diff=refMZ.time.diff,
	     target.mz.list = target.mz.list,subs=rand_set,colors=NA,apLCMS.outloc,runval=runval,presval=presval,
	     minexp,cdfloc=cdfloc),silent=TRUE)
						
						dev.off()
                                 	}
				}
					#getTICplots(cdfloc, filepattern=filepattern,mzlist=NA,rtlist=NA,deltaamu=1,deltatime=10)
				
				if(plotEICs=="all"){
				
				
				eic_fname<-paste(apLCMS.outloc,"/EICrun",runval,"pres",presval,"prefiltering_all.pdf",sep="")
		
				rows_ind<-seq(1,dim(finalfeatmat)[1])	
				#rows_ind<-rows_ind[1:20]
			
				pdf(eic_fname)
				time.list=finalfeatmat$time
					mz.list=finalfeatmat$mz
					chem.names<-paste("mz",seq(1,length(rows_ind)),sep="")
					#c1<-custom.EIC.plot(aligned, rows = rows_ind, colors = NA, transform = "none", 
					#	subset = rand_set, mz.list=mz.list, time.list=time.list,chem.names=chem.names,
					#	min.run=runval,min.pres=presval, max.spline.time.points = 1000,cdfloc=cdfloc) 
						
						c1<-apLCMS.EIC.plot(aligned, rows = c((overlap_res$index.A)), colors = colors, transform = "none",
    subset = subs, mz.list=mz.list, time.list=time.list,minrt=NA, maxrt=NA,chem.names=chem.names,min.run=runval,min.pres=presval,
    max.spline.time.points = 1000,rawprofileloc=cdfloc)
				dev.off()

					finalfeatmat<-finalfeatmat[which(c1>peak.score.thresh),]			

				}
			

			filenames_list[[pcount]]<-paste("apLCMS",apLCMSmode,"_aligned", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.txt", sep="")
			peak_score_list[[pcount]]<-PeakScore
			aligned_data_list[[pcount]]<-finalfeatmat
			pcount=pcount+1
                
		}
	

	#rm(apLCMSres)
	return(list(aligned_data_list=aligned_data_list,peak_score_list=peak_score_list,filenames_list=filenames_list))
}


#################################################################################
##Function: XCMS.align
##Description: Runs XCMS wrapper function, cdf.to.ftr, at different parameter settings
##cdfloc: The folder where all CDF files to be processed are located. For example "C:/CDF/"
##XCMS.outloc: The folder where alignment output will be written. For example "C:/CDFoutput/"
##step.list: list containing values for the step size
##mz.diff.list: list containing values for the minimum difference for features with
#	  retention time overlap
##sn.thresh.list: list containing values for signal to noise ratio cutoff variable
##max.list: list containing values for maxnimum number of peaks per EIC variable
##bw.val: bandwidth value
##minfrac.val:  minimum fraction of samples necessary in at least one of the sample
#         groups for it to be a valid group
##minsamp.val: minimum number of samples necessary in at least one of the sample
#          groups for it to be a valid group
##mzwid.val: width of overlapping m/z slices to use for creating peak density chromatograms
#         and grouping peaks across samples
##max.val: maximum number of groups to identify in a single m/z slice
##sleep.val: seconds to pause between plotting successive steps of the
#          peak grouping algorithm. peaks are plotted as points showing
#          relative intensity. identified groups are flanked by dotted
#          vertical lines.
#
##subs: If not all the CDF files in the folder are to be processed, the user can define
# a subset using this parameter. For example, subs=15:30, or subs=c(2,4,6,8)
##run.order.file: Name of a tab-delimited file that includes sample names sorted by the
# order in which they were run (sample names must match the CDF file names)
###################################################################################
XCMS.align.matchedFilter<-function(cdfloc, XCMS.outloc,step.list=c(0.001), mz.diff.list=c(0.1), sn.thresh.list=c(3), max=50, bw.val=c(10),
minfrac.val=0.5, minsamp.val=2, mzwid.val=0.25, sleep.val=0, run.order.file=NA,subs=NA, retcor.family="symmetric", retcor.plottype="mdevden",groupval.method="medret",target.mz.list = NA,xMSanalyzer.outloc=NA,nSlaves=2)
{

        setwd(cdfloc)
        dir.create(XCMS.outloc,showWarnings=FALSE)
        cdf_files=list.files(cdfloc,".cdf|.mzxml|mXML",ignore.case=TRUE)
	aligned_data_list=new("list")
	pcount=1
	
        if(is.na(subs[1])==FALSE)
        {
                cdf_files=cdf_files[subs]
                numsamp=length(subs)


        }
        else
        {

                numsamp=length(cdf_files)
        }

        for(t in sn.thresh.list)
        {
                for(s in step.list)
                {
                        for(m in mz.diff.list)
                        {
                            #for(maxl in max.list)
                                {
					 snow <- SnowParam(workers = nSlaves, type = "SOCK")

					xset=xcmsSet(cdf_files, step=s,snthresh=t,mzdiff=m,max=max,BPPARAM = snow)
					

                                        xset<-group(xset)
					
                                        xset2 <- retcor(xset, family = retcor.family, plottype = retcor.plottype)
					
					for(b in bw.val){
                                        ###  Group peaks together across samples, set bandwitdh, change important m/z parameters here
                                        ###  Syntax: group(object, bw = 30, minfrac = 0.5, minsamp= 1,  mzwid = 0.25, max = 5, sleep = 0)
                                        xset2 <- group.density(xset2, bw=b, minfrac=minfrac.val, minsamp=minsamp.val,  mzwid=mzwid.val, max=max, sleep=sleep.val)
                                        #xset3=fillPeaks(xset2)
       					xset3=suppressWarnings(fillPeaks(xset2))                                 
					
					print("Getting sample intensities")
                                        finalfeatmat={}
					
					if (dim(xset3@groups)[1] > 0) {
					     groupmat <- groups(xset3)
					     
					     
					     #group_intmat<-groupval(xset3, method = groupval.method, intensity = "into") 

					    group_intmat<-groupval(xset3,groupval.method,"into")
					     finalfeatmat<-cbind(groupmat,group_intmat)
					     finalfeatmat<-as.data.frame(finalfeatmat,row.names = NULL)
					     
					  } else{
						if (length(xset3@sampnames) == 1)
						{
							finalfeatmat  <- xset3@peaks
						}else 
						{
							stop ('First argument must be a xcmsSet with group information or contain only one sample.')
						 }
					}


                                        fname=paste(XCMS.outloc, "/XCMS_matchedFilter","_thresh", t,"_step", s, "_mzdiff", m,"_max",max,"_bw",b,".txt", sep="")
                                        cnames=colnames(finalfeatmat)
                                        cnames[1]="mz"
                                        cnames[4]="time"
                                        colnames(finalfeatmat)=c(cnames[1:8],cdf_files)
                                        cnames<-tolower(cnames)
                                        cnames<-gsub(".cdf", "", cnames)

                                        if(is.na(run.order.file)==FALSE)
                                        {
                                                fileorder=read.table(run.order.file, header=FALSE)
                                                fileorder=apply(fileorder,1,tolower)
                                                ordlist=sapply(1:length(fileorder),function(i){which(cnames==fileorder[i])})
                                                ordlist=unlist(ordlist)
                                                ordlist=ordlist+8
                                                
                                                finalfeatmat=finalfeatmat[,c(1:8,ordlist)]
                                        }

                                        write.table(finalfeatmat,file=fname,sep="\t",row.names=FALSE)
					aligned_data_list[[pcount]]<-finalfeatmat
					pcount=pcount+1

	if(is.na(target.mz.list)==FALSE){

	fname=paste(XCMS.outloc, "/EICmatchedFilter","_thresh", t,"_step", s, "_mzdiff", m,"_max",max,"_bw",b,".pdf", sep="")
                                        pdf(fname)
                                if(is.na(target.mz.list[1,1])==FALSE){
stddata<-target.mz.list

#print(head(stddata))
}else{
        Name<-paste("mz",seq(1,dim(finalfeatmat)[1]),sep="")
        stddata<-cbind(finalfeatmat[,c(1)],Name)

        }
overlapres5ppm<-getVenn(dataA=finalfeatmat, name_a=paste("Expdata",sep=""), dataB=stddata, name_b="Target", mz.thresh = 10, time.thresh=NA,alignment.tool=NA,
xMSanalyzer.outloc=XCMS.outloc,plotvenn=FALSE)
                                if(length(unique(overlapres5ppm$common$index.A))>0)
                                {

					print(cdfloc)
                                        num_samples<-dim(finalfeatmat)[2]-8
                                        min_samp<-6
                                        if(num_samples<min_samp){
                                                min_samp<-num_samples
                                        }
                                        rand_sample_set<-sample(size=num_samples,x=num_samples,replace=FALSE)
                                        rand_set<-c(1:min_samp,(num_samples-2):num_samples)
                                        com_ind<-which(rand_sample_set%in%rand_set)
                                        if(length(com_ind)>0){
                                        rand_sample_set<-rand_sample_set[-c(com_ind)]
                                        rand_sample_set<-rand_sample_set[1:min_samp]
                                        rand_set<-c(rand_set,rand_sample_set)
                                        }else{
                                                rand_set<-c(rand_set,rand_sample_set[1:min_samp])
                                                rand_set<-rand_set[order(rand_set)]
                                        }
                                        rand_set<-na.omit(rand_set)
                                        print(rand_set)
                                        #EIC.plot(aligned,rows=c(unique(overlapres5ppm$common$index.A)),min.run=runval,min.pres=presval)
                                        #apLCMS.EIC.plot(aligned, rows = c(unique(overlapres5ppm$common$index.A)), colors = NA, transform = "none",
                                        #subset = rand_set, minrt=NA, maxrt=NA, min.run=runval,min.pres=presval, max.spline.time.points = 1000)
                                        overlap_res<-overlapres5ppm$common
                                        overlap_res<-as.data.frame(overlap_res)
                                        dup_mz_ind<-which(duplicated(overlapres5ppm$common$index.A)==TRUE)
                                        if(length(dup_mz_ind)>0){
                                                overlap_res<-overlap_res[-c(dup_mz_ind),]
                                        }
                                        finalfeatmat<-as.data.frame(finalfeatmat)
                                        time.list=finalfeatmat$time[c((overlap_res$index.A))]
                                        mz.list=finalfeatmat$mz[c((overlap_res$index.A))]
                                        chem.names<-stddata$Name[c((overlap_res$index.B))]
	
					eicraw <- getEIC(xset3, groupidx = c(overlap_res$index.A), rt = "raw",step=0.001)
					for(i in 1:length(overlap_res$index.A)){
						plot(eicraw, xset3, groupidx = i)
					}
				}
				dev.off()
				}
	
				}
                                }
                        }
                }
        }
		if(is.na(xMSanalyzer.outloc)==FALSE){
				setwd(xMSanalyzer.outloc)
				}
				
				
			getTICplots(cdfloc, filepattern=".cdf|.mzxml|mXML",mzlist=NA,rtlist=NA,deltaamu=1,deltatime=10)
		
			setwd(XCMS.outloc)
	return(aligned_data_list)
}


XCMS.align.centWave<-function(cdfloc, XCMS.outloc,ppm.list=c(2.5), mz.diff.list=c(-0.00005), sn.thresh.list=c(3), prefilter.list=c(3,1000), bw.val=c(5),groupval.method="medret", 
step.list=c(0.1),max=50,minfrac.val=0.5, minsamp.val=1, mzwid.val=0.015, sleep.val=0, run.order.file=NA,subs=NA, retcor.method="obiwarp",retcor.family="symmetric", 
retcor.plottype="deviation", peakwidth=c(10,60),nSlaves=2,target.mz.list = NA,xMSanalyzer.outloc=NA)
{

        setwd(cdfloc)
	dir.create(XCMS.outloc,showWarnings=FALSE)
        cdf_files=list.files(cdfloc,".cdf|.mzxml|mXML",ignore.case=TRUE)
	aligned_data_list=new("list")
	pcount=1
	
	
	
        if(is.na(subs[1])==FALSE)
        {
                cdf_files=cdf_files[subs]
                numsamp=length(subs)


        }
        else
        {

                numsamp=length(cdf_files)
        }

        for(t in sn.thresh.list)
        {
                for(p in ppm.list)
                {
                        for(m in mz.diff.list)
                        {
                                #for(maxl in max.list)
                                for(s in step.list)
				{
					 snow <- SnowParam(workers = nSlaves, type = "SOCK")
					 
					 fname=paste(XCMS.outloc, "/XCMSstep1centwave","_snthresh", t,"_step", s, "_mzdiff", m,"_max",max,"_ppm",p,".Rda", sep="")
					check_if_exists<-suppressWarnings(try(load(fname),silent=TRUE))

					if(is(check_if_exists,"try-error")){
	
					xset=xcmsSet(cdf_files, method="centWave", ppm=p, snthresh=t,mzdiff=m,peakwidth=peakwidth,prefilter=prefilter.list,
					integrate=1, verbose.columns=TRUE,fitgauss=FALSE, BPPARAM = snow)
					save(xset,file=fname)
					
					}else{
						print(paste("Loading: ",fname),sep="")
	
						load(fname)
					}
					

					
					for(b in bw.val){
					
					fname=paste(XCMS.outloc, "/XCMSstep2centwave","_snthresh", t,"_mzdiff", m,"_pw",peakwidth[1],"_",peakwidth[2],"_bw",b,"_ppm",p,".Rda", sep="")
					check_if_exists<-suppressWarnings(try(load(fname),silent=TRUE))

					if(is(check_if_exists,"try-error")){

					####group.density
                                        xset<-group(xset, bw=b, minfrac=minfrac.val, minsamp=minsamp.val,  mzwid=mzwid.val, max=max, sleep=sleep.val)
					save(xset,file=fname)
					}else{
							print(paste("Loading: ",fname),sep="")
	
							load(fname)
						
					}
					
					
                                        fname=paste(XCMS.outloc, "/XCMSstep3centwave","_snthresh", t,"_mzdiff", m,"_pw",peakwidth[1],"_",peakwidth[2],"_bw",b,"_ppm",p,".Rda", sep="")
					check_if_exists<-suppressWarnings(try(load(fname),silent=TRUE))

					if(is(check_if_exists,"try-error")){
						if(retcor.method=="loess"){
					
					
							#####retention time correction using loess; mdevden
							xset2 <- retcor(xset, method=retcor.method,family = retcor.family, plottype = retcor.plottype)
						}else{
							if(retcor.method=="obiwarp"){
							###retention time correction using obiwarp method
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

				
					fname=paste(XCMS.outloc, "/XCMSstep4centwave","_snthresh", t,"_mzdiff", m,"_pw",peakwidth[1],"_",peakwidth[2],"_bw",b,"_ppm",p,".Rda", sep="")
					check_if_exists<-suppressWarnings(try(load(fname),silent=TRUE))

					if(is(check_if_exists,"try-error")){
							
							###  Group peaks together across samples, set bandwitdh, change important m/z parameters here
							###  Syntax: group(object, bw = 30, minfrac = 0.5, minsamp= 1,  mzwid = 0.25, max = 5, sleep = 0)
						       
							xset2 <- group(xset2, bw=b, minfrac=minfrac.val, minsamp=minsamp.val,  mzwid=mzwid.val, max=max, sleep=sleep.val)
							save(xset2,file=fname)
					}else{
						print(paste("Loading: ",fname),sep="")
						
						load(fname)
					}
					
										
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

					fname=paste(XCMS.outloc, "/XCMS_centwave","_snthresh", t,"_step", s, "_mzdiff", m,"_max",max,"_bw",b,"_ppm",p,".Rda", sep="")
					save(xset3,file=fname)

                                        fname=paste(XCMS.outloc, "/XCMS_centwave","_snthresh", t,"_step", s, "_mzdiff", m,"_max",max,"_bw",b,"_ppm",p,".txt", sep="")
                                        cnames=colnames(finalfeatmat)
                                        cnames[1]="mz"
                                        cnames[4]="time"
                                        colnames(finalfeatmat)=c(cnames[1:8],cdf_files)
                                        cnames<-tolower(cnames)
                                        cnames<-gsub(".cdf", "", cnames)

                                        if(is.na(run.order.file)==FALSE)
                                        {
                                                fileorder=read.table(run.order.file, header=FALSE)
                                                fileorder=apply(fileorder,1,tolower)
                                                ordlist=sapply(1:length(fileorder),function(i){which(cnames==fileorder[i])})
                                                ordlist=unlist(ordlist)
                                                ordlist=ordlist+8
                                                finalfeatmat=finalfeatmat[,c(1:8,ordlist)]
                                        }

                                        write.table(finalfeatmat,file=fname,sep="\t",row.names=FALSE)
					aligned_data_list[[pcount]]<-finalfeatmat
					
					
					
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

fname=paste("CAMERAresult_centwave","_snthresh", t,"_mzdiff", m,"_max",max,"_bw",b,"_ppm",p,".csv", sep="")
#Get final peaktable and store on harddrive
write.csv(getPeaklist(xsaFA),file=fname)
		

					
					pcount=pcount+1
					if(is.na(target.mz.list)==FALSE){

        fname=paste(XCMS.outloc, "/EICcentwave","_thresh", t,"_step", s, "_mzdiff", m,"_max",max,"_bw",b,".pdf", sep="")
                                        pdf(fname)
                                if(is.na(target.mz.list[1,1])==FALSE){
				stddata<-target.mz.list

				#print(head(stddata))
				}else{
				Name<-paste("mz",seq(1,dim(finalfeatmat)[1]),sep="")
				stddata<-cbind(finalfeatmat[,c(1)],Name)

				}
				overlapres5ppm<-getVenn(dataA=finalfeatmat, name_a=paste("Expdata",sep=""), dataB=stddata, name_b="Target", mz.thresh = p, time.thresh=NA,alignment.tool=NA,
				xMSanalyzer.outloc=XCMS.outloc,plotvenn=FALSE)
				
                                if(length(unique(overlapres5ppm$common$index.A))>0)
                                {
                                        setwd(cdfloc)
                                        num_samples<-dim(finalfeatmat)[2]-8
                                        min_samp<-6
                                        if(num_samples<min_samp){
                                                min_samp<-num_samples
                                        }
                                        rand_sample_set<-sample(size=num_samples,x=num_samples,replace=FALSE)
                                        rand_set<-c(1:min_samp,(num_samples-2):num_samples)
                                        com_ind<-which(rand_sample_set%in%rand_set)
                                        if(length(com_ind)>0){
                                        rand_sample_set<-rand_sample_set[-c(com_ind)]
                                        rand_sample_set<-rand_sample_set[1:min_samp]
                                        rand_set<-c(rand_set,rand_sample_set)
                                        }else{
                                                rand_set<-c(rand_set,rand_sample_set[1:min_samp])
                                                rand_set<-rand_set[order(rand_set)]
                                        }
                                        rand_set<-na.omit(rand_set)
                                        print(rand_set)
                                        overlap_res<-overlapres5ppm$common
                                        overlap_res<-as.data.frame(overlap_res)
                                        dup_mz_ind<-which(duplicated(overlapres5ppm$common$index.A)==TRUE)
                                        if(length(dup_mz_ind)>0){
                                                overlap_res<-overlap_res[-c(dup_mz_ind),]
                                        }
                                        finalfeatmat<-as.data.frame(finalfeatmat)
					time.list=finalfeatmat$time[c((overlap_res$index.A))]
                                        mz.list=finalfeatmat$mz[c((overlap_res$index.A))]
                                        chem.names<-stddata$Name[c((overlap_res$index.B))]

                                        eicraw <- getEIC(xset3, groupidx = c(overlap_res$index.A), rt = "corrected",step=0.001)
                                        for(i in 1:length(overlap_res$index.A)){
                                                plot(eicraw, xset3, groupidx = i)
                                        }
                                }
                                dev.off()
                                }

					


					}
                                }
                        }
                }
        }
	if(is.na(xMSanalyzer.outloc)==FALSE){
				setwd(xMSanalyzer.outloc)
				}
				
				
			getTICplots(cdfloc, filepattern=".cdf|.mzxml|mXML",mzlist=NA,rtlist=NA,deltaamu=1,deltatime=10)
		
			setwd(XCMS.outloc)
	return(aligned_data_list)
}


#################################################################
#Function: evaluate.Samples
#Description: Evaluate sample consistency based on Pearson Correlation
#curdata: feature alignment output matrix from apLCMS or XCMS with
# 	  intensities
#numreplicates: number of replicates per sample
#alignment.tool: name of the feature alignment tool eg: "apLCMS" or "XCMS"
##################################################################
evaluate.Samples<-function(curdata, numreplicates, alignment.tool, cormethod="pearson",missingvalue=0,ignore.missing=TRUE,replace.bad.replicates=TRUE)
{

	if(is.na(alignment.tool)==FALSE){
        if(alignment.tool=="apLCMS")
        {
              col_end=4
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    col_end=8
              }
              else
              {
                  stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
                  #col_end=NA
              }
        }
       }else
              {
                  #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
                  col_end=NA
              }
       
	rnames<-colnames(curdata)
    if(is.na(col_end)==FALSE){
        rnames<-rnames[-c(1:col_end)]
    }
        rnames=rnames[seq(1,length(rnames),numreplicates)]
        finalmat={}
 
         if(is.na(col_end)==FALSE){
        curdata_int=curdata[,-c(1:col_end)]
         }else{
             curdata_int=curdata
             
         }
        numsamp=dim(curdata_int)[2]
       
       curdata<-apply(curdata,2,as.numeric)
       
        for(samp in seq(1,numsamp,numreplicates))
        {
                samp_rep_last=samp+numreplicates-1
                subdata=curdata_int[,samp:samp_rep_last]
                
		if(ignore.missing==TRUE){
		if(is.na(missingvalue)==FALSE){
		
			subdata<-replace(as.matrix(subdata),which(subdata==missingvalue),NA)
			
		}
		
		}
		rmat=cor(subdata, method=cormethod,use="pairwise.complete.obs")
		
                rmat_upper=rmat[upper.tri(rmat)]
		
		good_reppairs<-which(rmat_upper>0.7)
		check_bad_rep<-length(good_reppairs)
		
		rmat2<-rmat
		diag(rmat2)<-0
		bad_reppair<-apply(rmat2,2,max)
		num_bad_reps<-length(which(bad_reppair<0.7))
		
		if(numreplicates>2){
		if(num_bad_reps==1){
			
				bad_rep<-which(bad_reppair==min(bad_reppair))
				
				subdata[,c(bad_rep)]<-apply(subdata[,-c(bad_rep)],1,function(x){mean(x,na.rm=TRUE)})
			
			}
		
		}
		
		if(replace.bad.replicates==TRUE){
		curdata_int[,samp:samp_rep_last]<-subdata
		rmat=cor(subdata, method=cormethod,use="pairwise.complete.obs")
		
                rmat_upper=rmat[upper.tri(rmat)]
		}
		if(numreplicates==2)
		{
			finalmat<-rbind(finalmat, mean(rmat_upper,na.rm=TRUE))
		}
		else
		{
			if(numreplicates>2)
			{
				rmat_vec=c(rmat_upper,mean(rmat_upper,na.rm=TRUE))
				finalmat<-rbind(finalmat,rmat_vec)
			}
			
		}
		
	}
	if(numreplicates==2)
		{
			colnames(finalmat)<-c(paste(cormethod,"Correlation",sep=""))
			cnames<-"PearsonCorrelation"
		}
		else
		{
			if(numreplicates>2)
			{
	
				cnames={}
				for(repnum in seq(1,numreplicates-1,1))
				{
					for(r1 in seq(repnum+1,numreplicates,1))
					{
						cnames<-c(cnames,paste("rep",repnum,"vs","rep",r1,sep=""))
					}
				}	
				cnames<-c(cnames,paste("mean","Correlation",sep=""))
			}
		}
		if(replace.bad.replicates==TRUE){
		
		if(is.na(col_end)==FALSE){
		curdata<-cbind(curdata[,c(1:col_end)],curdata_int)
		}else{
		curdata<-curdata_int
		}
		curdata<-replace(as.matrix(curdata),which(is.na(curdata)==TRUE),0)
		}
		
	colnames(finalmat)<-cnames
        rownames(finalmat)<-rnames
        return(list("cor.matrix"=finalmat,"feature.table"=curdata))
}

#################################################################
#Function: evaluate.Features
#Description: Evaluate feature consistency based on PID or CV
#curdata: feature alignment output matrix from apLCMS or XCMS with
#         intensities
#numreplicates: number of replicates per sample
#alignment.tool: name of the feature alignment tool eg: "apLCMS" or "XCMS"
##################################################################
evaluate.Features<-function(curdata,numreplicates,min.samp.percent=0.6,alignment.tool,  impute.bool,missingvalue=0,peak_scores=NA,numnodes=2)
{
	
        if(numreplicates==2)
        {
		print("**calculating percent intensity difference**")
                eval.feat.results<-getPID(curdata, alignment.tool,missingvalue, numreplicates,peak_scores=peak_scores)
        }
        else
        {
                  if(numreplicates>2)
                  {
			  print("**calculating CV**")
                          eval.feat.results<-getCVreplicates(curdata, alignment.tool, min.samp.percent=min.samp.percent, numreplicates=numreplicates,impute.bool= impute.bool,missingvalue=missingvalue,peak_scores=peak_scores,numnodes=numnodes)
                  }else{
		  
			#	eval.feat.results<-matrix(1,nrow=dim(curdata)[1],ncol=9)
			#colnames(eval.feat.results)<-c("min", "first_quartile", "median", "mean", "third_quartile", "max", "numgoodsamples","Qscore","peak_scores")
	
			print("****")
			eval.feat.results<-getPID(curdata, alignment.tool,missingvalue, numreplicates,peak_scores=peak_scores)
	
			eval.feat.results<-as.data.frame(eval.feat.results)
   
		  }
        
        }
	return(eval.feat.results)
}

#################################################################
#Function: getPID
#Description: Evaluate feature consistency based on PID
#curdata: feature alignment output matrix from apLCMS or XCMS with
#         intensities
#alignment.tool: name of the feature alignment tool eg: "apLCMS" or "XCMS"
##################################################################
getPID<-function(curdata, alignment.tool,missingvalue, numreplicates,peak_scores=NA)
{
        mean_replicate_difference<-{}
        sd_range_duplicate_pairs<-{}
        if(alignment.tool=="apLCMS")
        {
              col_end=4
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    col_end=8
              }
              else
              {
		     col_end=2
		     print("**Using the first two columns as mz and retention time for PID calculation**")
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }
        curdata_mz_rt_info=curdata[,c(1:col_end)]
        curdata=curdata[,-c(1:col_end)]
        numfeats=dim(curdata)[1]
        numsamp=dim(curdata)[2]
        
        maxint<-apply(curdata,1,function(x){max(x,na.rm=TRUE)})
        maxint<-log(maxint,10)
        
       # maxint<-(maxint)/(max(maxint,na.rm=TRUE))
                
        rnames<-colnames(curdata)
        rnames<-gsub(".cdf", "", rnames, ignore.case=TRUE)
        quantcolnames=c("min", "first_quartile", "median", "mean", "third_quartile", "max")
	
	if(numreplicates>1){
        resvec_1<-lapply(1:numfeats,function(r)
        {
                newrow={}
                finalmat={}
		no_value=0
		goodsamps<-0
                for(samp in seq(1,(numsamp),2))
                {
                        i=samp
                        j=i+1
                        int1=curdata[r,i]
                        int2=curdata[r,j]
			
			curdata_int=curdata[r,c(i:j)]
						if(is.na(missingvalue)==TRUE){
							
							check_zeros=which(is.na(curdata_int)==TRUE)
                        
						}else{
							check_zeros=which(curdata_int==missingvalue)
								
						}
                        
			if(length(check_zeros)>0)
                        {
                                                 
                                replicate_diff<-NA
								no_value<-no_value+1
                        }
                        else
                        {
                                #calculate PID
                                replicate_diff<-100*(abs(int1-int2)/mean(c(int1,int2)))
								goodsamps<-goodsamps+1
                        }
                        newrow<-cbind(newrow,replicate_diff)
                }

                #get indices of the PIDs that are NA
                na_ind=which(is.na(newrow)==TRUE)

                #get quantile summary of the percent intensity difference (PID) vector
                #using only the non-NA values
                if(length(na_ind)>0)
                {
                        sumrow=summary(as.vector(newrow[-na_ind]))
                }
                else
                {
                        sumrow=summary(as.vector(newrow))
                }

                #if quantile values are set  to NA
                if(length(sumrow)<6)
                {
                        for(i in 1:6)
                        {
                                sumrow[i]=200
                        }
                }
                names(sumrow)=quantcolnames

		finalmat<-rbind(finalmat, c(unlist(sumrow),goodsamps))
                return(finalmat)
        })
	
        final_set={}
        for(i in 1:length(resvec_1)){
                if(length(resvec_1[[i]])>1)
                {
                        final_set<-rbind(final_set,resvec_1[[i]])
                }
        }

        final_set<-as.data.frame(final_set)
        rownames(final_set)=NULL
	    final_set<-apply(final_set,2,as.numeric)
		final_set<-as.data.frame(final_set)
		
		colnames(final_set)<-c("min", "first_quartile", "median", "mean", "third_quartile", "max", "numgoodsamples")
		numsamp<-numsamp/2
	mz_min_max<-cbind(curdata_mz_rt_info[,1],curdata_mz_rt_info[,1])
	if(alignment.tool=="apLCMS"){
		mz_min_max<-cbind(curdata_mz_rt_info[,3],curdata_mz_rt_info[,4])
	}else{
		if(alignment.tool=="XCMS"){
			mz_min_max<-cbind(curdata_mz_rt_info[,2],curdata_mz_rt_info[,3])
		}
	}
	mz_min_max<-as.data.frame(mz_min_max)
	
	deltappm_res<-apply(mz_min_max,1,get_deltappm)
	delta_cv_range<-as.numeric(final_set$max)-as.numeric(final_set$min)+0.1
	
	#Qscore<-100*((final_set$numgoodsamples)/(delta_cv_range*as.numeric(final_set$median)*numsamp*(deltappm_res+0.1)))
	
	if(is.na(missingvalue)==TRUE){
							num_tot_nonzeros<-apply(curdata,1,function(x){length(which(is.na(x)==FALSE))})
                        
                        
						}else{
							
							num_tot_nonzeros<-apply(curdata,1,function(x){length(which(x>missingvalue))})
                        	
						}
	
	
				
							

				

	#part of xMSanalyzer_v2.0.7
	termA<-final_set$median+0.01+1
        termB<-(final_set$numgoodsamples)/(num_tot_nonzeros*(1/numreplicates))

	percent_samp_zeros<-100-(100*(num_tot_nonzeros/dim(curdata)[2]))

	termB<-(final_set$numgoodsamples*numreplicates)/(percent_samp_zeros+1)
	termC<-replace(peak_scores,which(is.na(peak_scores)==TRUE),0.1)
	termC<-replace(termC,which(termC=="Inf"),1)
	termC<-replace(termC,which(termC=="-Inf"),1)

	
	Qscore<-100*(termC)*(termB)*(1/termA)
	final_set<-cbind(final_set,Qscore,peak_scores)
	}else{
	
				if(is.na(missingvalue)==TRUE){
							num_tot_nonzeros<-apply(curdata,1,function(x){length(which(is.na(x)==FALSE))})
                        
                        
						}else{
							
							num_tot_nonzeros<-apply(curdata,1,function(x){length(which(x>missingvalue))})
                        	
						}
						
						
									
												

									

						#part of xMSanalyzer_v2.0.7
						termA<-1
						termB<-1

						percent_samp_zeros<-100-(100*(num_tot_nonzeros/dim(curdata)[2]))

						termB<-(num_tot_nonzeros)/(percent_samp_zeros+1)
						termC<-1

						
						Qscore<-100*(termC)*(termB)*(1/termA)
						
						final_set<-matrix(1,nrow=dim(curdata)[1],ncol=9)
						colnames(final_set)<-c("min", "first_quartile", "median", "mean", "third_quartile", "max", "numgoodsamples","Qscore","peak_scores")
						
						final_set<-as.data.frame(final_set)
						
						final_set$Qscore<-Qscore
						
						print(head(final_set))
	}
	
	final_set<-as.data.frame(final_set)
        return(final_set)
}

 get_deltappm<-function(x){
	if(is.na(x[1])==FALSE){
	return(10^6*(abs(x[1]-x[2])/x[1]))
	}else{
	return(0)
	}
}

#################################################################
#Function: getCVreplicates
#Description: Evaluate feature consistency based on coefficient of
# 	   Variation
#curdata: feature alignment output matrix from apLCMS or XCMS with
#         intensities
#numreplicates: number of replicates per sample
#alignment.tool: name of the feature alignment tool eg: "apLCMS" or "XCMS"
##################################################################
getCVreplicates<-function(curdata,alignment.tool,numreplicates, min.samp.percent=0.6,impute.bool=TRUE,missingvalue=0,peak_scores=NA,numnodes=2)
{
        mean_replicate_difference<-{}
        sd_range_duplicate_pairs<-{}
	
        if(alignment.tool=="apLCMS")
        {
              col_end=4
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    col_end=8
              }
              else
              {
			  col_end=2
		     print("**Using the first two columns as mz and retention time for PID calculation**")
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }
        
        curdata_mz_rt_info=curdata[,c(1:col_end)]
        curdata=curdata[,-c(1:col_end)]
	
	#curdata<-curdata[1:10,]
        numfeats=dim(curdata)[1]
        numsamp=dim(curdata)[2]
        rnames<-colnames(curdata)
        rnames<-gsub(".cdf", "", rnames, ignore.case=TRUE)
	quantcolnames=c("min", "first_quartile", "median", "mean", "third_quartile", "max","sampleCount")
	if(impute.bool==TRUE)
	{
		min.samp.percent=min.samp.percent
	}
	else
	{
		min.samp.percent=min.samp.percent
	}

        
                newrow={}
                finalmat={}
		
		
		cl<-parallel::makeCluster(numnodes)
			
		
		
		clusterExport(cl, "getCVreplicates.child") 
		#clusterExport(cl, "numsamp")
		#clusterExport(cl, "numreplicates")
		#clusterExport(cl, "min.samp.percent")
		#clusterExport(cl, "impute.bool")
		#clusterExport(cl, "alignment.tool")
		if(numreplicates>1){
		
		cv.res<-parApply(cl,curdata,1,getCVreplicates.child,numsamp=numsamp,numreplicates=numreplicates,min.samp.percent=min.samp.percent,impute.bool=impute.bool,
		missingvalue=missingvalue)
		
		dim(cv.res)=dim(matrix(nrow=7,ncol=numfeats))
	#cv.res<-t(cv.res)
			
		
	stopCluster(cl)
        final_set<-as.data.frame(cv.res)
        rownames(final_set)=NULL
	final_set<-apply(final_set,2,as.numeric)
	final_set<-as.data.frame(t(final_set))
	#print(final_set)
	numsamp=numsamp/numreplicates
	colnames(final_set)<-c("min", "first_quartile", "median", "mean", "third_quartile", "max","numgoodsamples")
	
	#deltappm_res<-apply(curdata_mz_rt_info,1,get_deltappm)
	#delta_cv_range<-as.numeric(final_set$max)-as.numeric(final_set$min)+0.1
	
	#Qscore<-100*((final_set$numgoodsamples)/(delta_cv_range*numsamp*(deltappm_res+0.1)))
	
	
	#Qscore<-100*(final_set$numgoodsamples/(final_set$median*numsamp))

	 mz_min_max<-cbind(curdata_mz_rt_info[,1],curdata_mz_rt_info[,1])	
	if(alignment.tool=="apLCMS"){
		mz_min_max<-cbind(curdata_mz_rt_info[,3],curdata_mz_rt_info[,4])
	}else{
		if(alignment.tool=="XCMS"){
			mz_min_max<-cbind(curdata_mz_rt_info[,2],curdata_mz_rt_info[,3])
		}
	}
	mz_min_max<-as.data.frame(mz_min_max)
	
	deltappm_res<-apply(mz_min_max,1,get_deltappm)
	delta_cv_range<-1 #as.numeric(final_set$max)-as.numeric(final_set$min)+0.1
	
	
	#Qscore<-100*(final_set$numgoodsamples)/((((as.numeric(final_set$median)+1)*numsamp*((10*(deltappm_res+0.1))))+1))
	
	if(is.na(missingvalue)==TRUE){
							num_tot_nonzeros<-apply(curdata,1,function(x){length(which(is.na(x)==FALSE))})
                        
                        
						}else{
							
							num_tot_nonzeros<-apply(curdata,1,function(x){length(which(x>missingvalue))})
                        	
						}
	

	
	#part of xMSanalyzer_v2.0.7
	termA<-final_set$median+0.01+1
        #termB<-(final_set$numgoodsamples)/(num_tot_nonzeros*(1/numreplicates))
	
	percent_samp_zeros<-100-(100*(num_tot_nonzeros/dim(curdata)[2]))

        termB<-(final_set$numgoodsamples*numreplicates)/(percent_samp_zeros+1)

	termC<-replace(peak_scores,which(is.na(peak_scores)==TRUE),0.1)
	termC<-replace(termC,which(termC=="Inf"),1)
	termC<-replace(termC,which(termC=="-Inf"),1)
	#print(termC)
	
	Qscore<-100*(termC)*(termB)*(1/termA)
	
	
	final_set<-cbind(final_set,Qscore,peak_scores)
	}else{
	
	
		final_set<-matrix(1,nrow=dim(curdata)[1],ncol=9)
			colnames(final_set)<-c("min", "first_quartile", "median", "mean", "third_quartile", "max", "numgoodsamples","Qscore","peak_scores")
	}
	
	final_set<-as.data.frame(final_set)
	
	#quantcolnames
        return(final_set)
}

#min.samp.percent: if greater than equal to min.samp.percent CV is set to NA
getCVreplicates.child<-function(curdata,numsamp,numreplicates,min.samp.percent,impute.bool,missingvalue)
{
	
			newrow={}
			#numsamp=length(curdata)
			goodsamps<-0
			for(samp in seq(1,numsamp,numreplicates))
			{
				i=samp
				j=i+numreplicates-1
				
				curdata_int=as.numeric(curdata[c(i:j)])
				
				if(is.na(missingvalue)==TRUE){
							
					check_zeros=which(is.na(curdata_int)==TRUE)
					check_values=which(is.na(curdata_int)==FALSE)
                        
						}else{
							check_zeros=which(curdata_int==missingvalue)
							check_values=which(curdata_int!=missingvalue)	
						}
				#check_zeros=which(curdata_int==0)
				
				na_thresh=round(min.samp.percent*numreplicates)
				
				if(length(check_values)<2){
				
					cvval<-NA
				}else{
				
				print(na_thresh)
				print(length(check_zeros))
				print(curdata_int)
					if(length(check_zeros)>=na_thresh)
					{
						
							cvval<-NA
						
							#newrow<-cbind(newrow,cvval)
							
					}
					else
					{
						#temporarily replace the missing intensities, set to 0 in apLCMS,
						#with mean intensity value of the corresponding replicates (with non-zero values)
						if(length(check_zeros)>0){
						
						if(impute.bool==TRUE){
						
							if(numreplicates==3){
							
							 replicate_diff<-100*(diff(curdata_int[-c(check_zeros)],na.rm=TRUE))/mean(t(curdata_int[-c(check_zeros)]),na.rm=TRUE)
								if(replicate_diff<20){
									curdata_int[check_zeros]=mean(t(curdata_int[-c(check_zeros)]),na.rm=TRUE)
								}else{
							
									curdata_int[check_zeros]=NA #median(t(curdata_int[-c(check_zeros)]))
								} 
							 
								#curdata_int[check_zeros]=mean(t(curdata_int[-c(check_zeros)]),na.rm=TRUE)
							 
							}else{
									curdata_int[check_zeros]=NA
							
							}
							
							
							sdval<-sd(curdata_int,na.rm=TRUE)
							meanval<-mean(t(curdata_int),na.rm=TRUE)
						}else{
							sdval<-sd(curdata_int[-c(check_zeros)],na.rm=TRUE)
							meanval<-mean(t(curdata_int[-c(check_zeros)]),na.rm=TRUE)
						
						}
						
						
						}else{
						
								sdval<-sd(curdata_int,na.rm=TRUE)
								meanval<-mean(t(curdata_int),na.rm=TRUE)
						}
						
						cvval<-100*(sdval/meanval)
						newrow<-cbind(newrow,cvval)
						
						goodsamps<-goodsamps+1
					}
					}
			}
			if(length(newrow)>0)
			{
				na_ind=which(is.na(newrow)==TRUE)
			
				if(length(na_ind)>0)
				{
					sumrow=summary(as.vector(newrow[-na_ind]),na.rm=TRUE)
				}
				else
				{
					sumrow=summary(as.vector(newrow),na.rm=TRUE)
				}
			}
			else{sumrow<-{}}
			if(length(sumrow)<6)
			{
				for(i in 1:6)
				{
				 sumrow[i]=200
				}
			}
			finalmat<-{}

			finalmat<-rbind(finalmat, c(unlist(sumrow),goodsamps))
	
			
			return(finalmat)

}
#################################################################
#Function: Function that merges results from different parameter settings
#Description: Evaluate feature consistency based on PID or CV
#dataA: feature alignment output matrix from apLCMS or XCMS with
#         intensities at parameter settings P1
#dataB: feature alignment output matrix from apLCMS or XCMS with
#         intensities at parameter settings P2
##max.mz.diff: +/- mz tolerance in ppm for feature matching
##max.rt.diff: retention time tolerance for feature matching
##tstatistic: Threshold for defining signifcance level of the paired t-test
#numreplicates: number of replicates per sample
#alignment.tool: name of the feature alignment tool eg: "apLCMS" or "XCMS"

ttest.mz<-function(x, y)
				                                                        {
				                                                                
				                                                                
				                                                               
				                                                                        ttest_res=t.test(t(x), as.matrix(y),paired=T)
													ttest_res=abs(ttest_res$p.value)
				                                                                
				                                                                return(ttest_res)
				                                                        }



#################################################################
#Function: Function that merges results from different parameter settings
#Description: Evaluate feature consistency based on PID or CV
#dataA: feature alignment output matrix from apLCMS or XCMS with
#         intensities at parameter settings P1
#dataB: feature alignment output matrix from apLCMS or XCMS with
#         intensities at parameter settings P2
##max.mz.diff: +/- mz tolerance in ppm for feature matching
##max.rt.diff: retention time tolerance for feature matching
##tstatistic: Threshold for defining signifcance level of the paired t-test
#numreplicates: number of replicates per sample
#alignment.tool: name of the feature alignment tool eg: "apLCMS" or "XCMS"
##################################################################
merge.Results<-function(dataA, dataB, feat.eval.A,feat.eval.B,max.mz.diff=15, max.rt.diff=300,merge.eval.pvalue=0.05,alignment.tool="apLCMS", numnodes=1, mult.test.cor=FALSE,
mergecorthresh=0.7,missingvalue=0)
{
	#dataA<-unique(dataA)
	#dataB<-unique(dataB)

	#feat.eval.A<-unique(feat.eval.A)
	#feat.eval.B<-unique(feat.eval.B)
	
	dataA<-dataA[order(dataA$mz),]
	dataB<-dataB[order(dataB$mz),]
	
	
	feat.eval.A<-feat.eval.A[order(feat.eval.A$mz),]
	
	feat.eval.B<-feat.eval.B[order(feat.eval.B$mz),]
	

	
	  if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
	      
	      dataA<-cbind(dataA,feat.eval.A[,c(sample.col.start:(dim(feat.eval.A)[2]))])
	      dataB<-cbind(dataB,feat.eval.B[,c(sample.col.start:(dim(feat.eval.B)[2]))])
	curdata<-rbind(dataA,dataB)
	curdata<-unique(curdata)
	curdata<-curdata[order(curdata$mz),]
	curdata<-as.data.frame(curdata)
	
			
	      npeaks<-apply(curdata[,sample.col.start:(dim(curdata)[2]-8)],1,countpeaks)
	      curdatatemp<-cbind(curdata[,c(1:4)],npeaks)
	      curdata<-cbind(curdatatemp,curdata[,sample.col.start:dim(curdata)[2]])
	      sample.col.start=6
	     
	      curdata<-as.data.frame(curdata)
	     
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
		    
				    cnames<-colnames(dataA)
                    cnames[1]="mz"

                    cnames[4]="time"
                    colnames(dataA)=cnames
                    colnames(dataB)=cnames
		    
		       	dataA<-cbind(dataA,feat.eval.A[,c(sample.col.start:(dim(feat.eval.A)[2]))])
				dataB<-cbind(dataB,feat.eval.B[,c(sample.col.start:(dim(feat.eval.B)[2]))])
				curdata<-rbind(dataA,dataB)
				curdata<-unique(curdata)	
				curdata<-curdata[order(curdata$mz),]
				curdata<-as.data.frame(curdata)
			
				#peak score 1	
				curdata[,dim(curdata)[2]]<-1
				
			#curdata<-curdata[-which(curdata$min==curdata$max),]
              }
              else
              {
                    stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }
	

	sub_data_a<-new("list")
	sub_data_b<-new("list")
	lindex<-1
	
	min_mz<-min(curdata[,1])
	max_mz<-max(curdata[,1])

	mz_group=10

	mzdefect<-1*((curdata$mz-floor(curdata$mz)))

	#curdata<-cbind(mzdefect,curdata)
	
	curdata<-as.data.frame(curdata)
	curdata<-unique(curdata)
	max_rt_thresh<-max(curdata$time,na.rm=TRUE)
	
	end_ind<-dim(curdata)[2]-9	
	
	
	
		#lindex<-length(sub_data_a)
	diff_mz_num=1
	 #Step 1 Group features by m/z
        mz_groups<-lapply(1:dim(curdata)[1],function(j){
                                commat={}
                                diffmz=new("list")
                                ppmb=(500)*(curdata$mz[j]/1000000)
                                getbind_same<-which(abs(curdata$mz-curdata$mz[j])<=ppmb & abs(curdata$time-curdata$time[j])<=max.rt.diff)
                                diffmz[[diff_mz_num]]=getbind_same #dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })

	options(warn=-1)
	
	  	  
	  del_list<-{}
	
	  for(m in 1:length(mz_groups)){
	  
	   if((m%in%del_list)==FALSE){
		for(n in (m+1):length(mz_groups)){
		
			if(n>length(mz_groups)){
				break;
			}
			com1<-intersect(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
			if(length(com1)>0){
			
				mz_groups[[m]][[1]]<-c(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
				del_list<-c(del_list,n)
			}
		}
		
		mz_groups[[m]][[1]]<-unique(mz_groups[[m]][[1]])
		
		}
	  }
	  
	  
	
	  if(length(del_list)>0){
		 mz_groups<-mz_groups[-del_list]
		}
	 
	  
	mz_groups<-unique(mz_groups)
	
	sub_data_a<-lapply(1:length(mz_groups),function(j){
                                commat={}
                               # diffmz=new("list")
                                getbind_same=mz_groups[[j]][[1]]
                                diffmz=curdata[getbind_same,]
                                
                                return(diffmz)
          })
	rm(mz_groups)
	lindex<-length(sub_data_a)
	
	sub_data_a<-unique(sub_data_a)
	
	library(parallel)
	
	if(lindex>2){
	
	print("start making cluster")
	print(numnodes)
			cl<-parallel::makeCluster(numnodes)
	print("end make cluster")
		
			clusterEvalQ(cl, "merge.Results.child.ttest")
			clusterEvalQ(cl, "compare_intensities_ttest")
			
			merge.res<-parLapply(cl,sub_data_a,merge.Results.child.ttest,max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,
			alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
		
		if(FALSE){
				merge.res<-new("list")
				
				for(lsub in 1:length(sub_data_a)){
					
					tempres<-merge.Results.child.ttest(dataA=sub_data_a[[lsub]],max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
					
					merge.res[[lsub]]<-tempres
					
				}
		}	
	
		stopCluster(cl)
		
		
		end_ind<-dim(curdata)[2]-9
		l1<-lapply(1:length(merge.res),function(x){nrow(merge.res[[x]])})
		l2<-unlist(l1)
		
		if(length(which(l2<1))>0){
		merge.res<-merge.res[-which(l2<1)]	
		}
		
		l1<-lapply(1:length(merge.res),function(x){if(is.na(merge.res[[x]][1])[1]==FALSE){return(as.data.frame(merge.res[[x]]))}})
		final.res<-ldply(l1,rbind)

		if(FALSE){
		l2<-unlist(l1)
		
		if(length(which(l2==TRUE))>0){
		merge.res<-merge.res[-which(l2==TRUE)]	
		}
		
		final.res={}
		
		
		save(merge.res,file="merge.res.Rda")
		final.res<-ldply(merge.res,rbind)
		#final.res={}
		
		if(FALSE){
		for(s in 1:length(merge.res))
		{	
			tempd<-as.data.frame(merge.res[[s]])
			final.res<-rbind(final.res,tempd)
			
		}
		}
		
		}
		
		
		final.res<-unique(final.res)
	
	}else{
		
			final.res<-merge.Results.child.ttest(curdata,max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,
			alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
	}
	
	final.res<-as.data.frame(final.res)
	
	final.res<-unique(final.res)
	
	
	final.res<-final.res[order(final.res$mz),]
	
	
	
	#if(FALSE)
	{
	curdata<-final.res
	curdata<-as.data.frame(curdata)


	if(ncol(curdata[,c(sample.col.start:end_ind)])>500){
	
		set.seed(555)
		samp_ind<-sample(x=seq(sample.col.start:end_ind),size=500)
		system.time(global_cor<-try(WGCNA::cor(t(curdata[,c(samp_ind)]),nThreads=numnodes,use = 'p'),silent=TRUE))

	}else{

		system.time(global_cor<-try(WGCNA::cor(t(curdata[,c(sample.col.start:end_ind)]),nThreads=numnodes,use = 'p'),silent=TRUE))

	}
	if(is(global_cor,"try-error")){

		non_cor_data<-curdata

	}else{
			cl<-parallel::makeCluster(numnodes)
	
		
			
			clusterEvalQ(cl, "count_overlapping_feats")
			
			#function
			count_overlapping_feats<-function(j,mzvec,rtvec,max.mz.diff,max.rt.diff,mergecorthresh,corvec){
			
				diff_mz<-10^6*abs(mzvec[j]-mzvec)/mzvec[j]; 
				diff_rt<-abs(rtvec[j]-rtvec);
				return(length(which(diff_mz<max.mz.diff & diff_rt<max.rt.diff & corvec[j,]>=mergecorthresh)))
			}
			
		#count_cor<-lapply(1:dim(curdata)[1],function(j){diff_mz<-10^6*abs(curdata[j,1]-curdata[,1])/curdata[j,1]; diff_rt<-abs(curdata[j,2]-curdata[,2]); length(which(diff_mz<max.mz.diff & diff_rt<max.rt.diff & global_cor[j,]>=mergecorthresh))})
		
		count_cor<-parLapply(cl,1:dim(curdata)[1],count_overlapping_feats,mzvec=curdata[,1],rtvec=curdata[,2],max.mz.diff,max.rt.diff,mergecorthresh,corvec=global_cor)
		
		stopCluster(cl)
		count_cor<-unlist(count_cor)
		
		
		if(length(which(count_cor<=1))>0){
			non_cor_data<-curdata[which(count_cor<=1),]
			curdata<-curdata[-which(count_cor<=1),]
			#global_cor<-global_cor[-which(count_cor<=1),-which(count_cor<=1)]
			

		}
	
	}
	
	

	if(dim(curdata)[1]>0){
      
		final.res2<-merge.Results.child.ttest(curdata,max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
		
		final.res2<-unique(final.res2)
	    
		final.res<-rbind(final.res2,non_cor_data)
		final.res<-as.data.frame(final.res)
		final.res<-final.res[order(final.res$mz),]
		final.res<-unique(final.res)	
		
	
	}
    
    curdata<-final.res
   
   	rm(curdata) 
  
    
    
    
    
    
	
	}
	
	if(is.na(missingvalue)==FALSE){
		
			final.res<-replace(as.matrix(final.res),which(is.na(final.res)==TRUE),missingvalue)
			
		}

	final.res<-as.data.frame(final.res)
	
	print("dim of final res")
    
	final.res<-unique(final.res)
	print(dim(final.res))
	
	options(warn=0)
	return(final.res)
	
	
	
}



merge.Resultsv2.0.8.32H<-function(dataA, dataB, feat.eval.A,feat.eval.B,max.mz.diff=15, max.rt.diff=300,merge.eval.pvalue=0.05,alignment.tool="apLCMS", numnodes=1, mult.test.cor=FALSE,
mergecorthresh=0.7,missingvalue=0)
{
	#dataA<-unique(dataA)
	#dataB<-unique(dataB)

	#feat.eval.A<-unique(feat.eval.A)
	#feat.eval.B<-unique(feat.eval.B)
	
	dataA<-dataA[order(dataA$mz),]
	dataB<-dataB[order(dataB$mz),]
	
	
	feat.eval.A<-feat.eval.A[order(feat.eval.A$mz),]
	
	feat.eval.B<-feat.eval.B[order(feat.eval.B$mz),]
	
	

	  if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
	      
	      dataA<-cbind(dataA,feat.eval.A[,c(sample.col.start:(dim(feat.eval.A)[2]))])
	      dataB<-cbind(dataB,feat.eval.B[,c(sample.col.start:(dim(feat.eval.B)[2]))])
	#curdata<-rbind(dataA,dataB)
	#curdata<-unique(curdata)
	#curdata<-curdata[order(curdata$mz),]
	#curdata<-as.data.frame(curdata)
	
			
	  #    npeaks<-apply(curdata[,sample.col.start:(dim(curdata)[2]-8)],1,countpeaks)
	    #  curdatatemp<-cbind(curdata[,c(1:4)],npeaks)
	     # curdata<-cbind(curdatatemp,curdata[,sample.col.start:dim(curdata)[2]])
	      #sample.col.start=6
	     
	      #curdata<-as.data.frame(curdata)
	     
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
		    
		   cnames<-colnames(dataA)
                    cnames[1]="mz"

                    cnames[4]="time"
                    colnames(dataA)=cnames
                    colnames(dataB)=cnames
		    
		         	dataA<-cbind(dataA,feat.eval.A[,c(sample.col.start:(dim(feat.eval.A)[2]))])
				dataB<-cbind(dataB,feat.eval.B[,c(sample.col.start:(dim(feat.eval.B)[2]))])
			#	curdata<-rbind(dataA,dataB)
			
			#	curdata<-curdata[order(curdata$mz),]
			#	curdata<-as.data.frame(curdata)
				
			#	curdata[,dim(curdata)[2]]<-1
				
				#curdata<-curdata[-which(curdata$min==curdata$max),]
              }
              else
              {
			#  stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		  
			    sample.col.start=3
		    
			    cnames<-colnames(dataA)
			    cnames[1]="mz"

			    cnames[2]="time"
			    colnames(dataA)=cnames
			    colnames(dataB)=cnames
		    
			    dataA<-cbind(dataA,feat.eval.A[,c(sample.col.start:(dim(feat.eval.A)[2]))])
			   dataB<-cbind(dataB,feat.eval.B[,c(sample.col.start:(dim(feat.eval.B)[2]))])
			   
			 #  curdata<-rbind(dataA,dataB)
			
				#curdata<-curdata[order(curdata$mz),]
				#curdata<-as.data.frame(curdata)
              }
        }
	

	sub_data_a<-new("list")
	sub_data_b<-new("list")
	lindex<-1
	
	
	mz_group=10


	diff_mz_num=1
	
	s1=Sys.time()
	print(s1)
	

	options(warn=-1)
	print("s2")

	#max.mz.diff=15, max.rt.diff=300
	#common_mz_matrix<-find.Overlapping.mzs(dataA=dataA, dataB=dataB, mz.thresh=max.mz.diff, time.thresh=max.rt.diff, alignment.tool=alignment.tool,nearest.time.match=FALSE)

if(FALSE){
	temp_data<-rbind(dataA[,c("mz","time")],dataB[,c("mz","time")])
	
	overlap_res<-getVenn(dataA=dataA[,c("mz","time")],name_a="A", dataB=dataB[,c("mz","time")],name_b="B",mz.thresh=max.mz.diff, time.thresh=max.rt.diff, alignment.tool=alignment.tool, xMSanalyzer.outloc=NA,use.unique.mz=FALSE,plotvenn=FALSE,nearest.time.match=FALSE, numnodes=numnodes)
	
	
	s2=Sys.time()
	
	print((s2-s1))
	

	print("s3")		 
	s3=Sys.time()
	
	print((s3-s2))
	
	cl<-parallel::makeCluster(numnodes)
	
	unique_indexA<-(overlap_res$common$index.A) #unique(overlap_res$common$index.A)
	
	curdata<-rbind(dataA[overlap_res$common$index.A,],dataB[overlap_res$common$index.B,])
	curdata<-as.data.frame(curdata)
	curdata<-unique(curdata)
	min_mz<-min(curdata[,1])
	max_mz<-max(curdata[,1])

	noncordata<-rbind(dataA[-overlap_res$common$index.A,],dataB[-overlap_res$common$index.B,])
	noncordata<-as.data.frame(noncordata)
	noncordata<-unique(noncordata)
	
	print("dim of redundant data")
	print(dim(curdata))
	
	print("dim of non-redundant data")
	print(dim(noncordata))
}

	curdata<-rbind(dataA,dataB)
	curdata<-as.data.frame(curdata)
	curdata<-unique(curdata)
	
	curdata<-curdata[order(curdata$mz),]
	
	min_mz<-min(curdata[,1])
	max_mz<-max(curdata[,1])

	if(alignment.tool=="apLCMS"){
	
	   npeaks<-apply(curdata[,sample.col.start:(dim(curdata)[2]-8)],1,countpeaks)
	     curdatatemp<-cbind(curdata[,c(1:4)],npeaks)
	      curdata<-cbind(curdatatemp,curdata[,sample.col.start:dim(curdata)[2]])
	    
	     

	
	print("done")
		  sample.col.start=6
		  curdata<-as.data.frame(curdata)
		 
	}
	
	#calculate mass defect
	md1=round(abs(curdata$mz-round(curdata$mz,0)),2)
	
	#get row number
	rowindex=seq(1,nrow(curdata))

	
	curdata_temp<-cbind(rowindex,md1,curdata[,c("mz","time")])
	
	curdata_temp<-as.data.frame(curdata_temp)
	
	mz_groups<-aggregate(curdata_temp,by=list(roundmz=curdata_temp$md1),FUN=function(x){return(x)})
	
	rm(curdata_temp)
	
	
	
	sub_data_a<-lapply(1:length(mz_groups[[2]]),function(j){
                                commat={}
                               # diffmz=new("list")
                                getbind_same=mz_groups[[2]][[j]] #[[1]]
                                diffmz=curdata[getbind_same,]
                                
                                return(diffmz)
          })
	
	lindex=length(sub_data_a)

end_ind<-dim(curdata)[2]-9
	
	if(lindex>2){
			cl<-parallel::makeCluster(numnodes)
	
		
			clusterEvalQ(cl, "merge.Results.child.ttest")
			clusterEvalQ(cl, "compare_intensities_ttest")
			
			merge.res<-parLapply(cl,sub_data_a,merge.Results.child.ttest,max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,
			alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)

	
				stopCluster(cl)
				
				
			save(merge.res,file="mergeres.Rda")
				
				
			end_ind<-dim(curdata)[2]-9
			l1<-lapply(1:length(merge.res),function(x){nrow(merge.res[[x]])})
			l2<-unlist(l1)
			
			if(length(which(l2<1))>0){
			merge.res<-merge.res[-which(l2<1)]	
			}
			
			l1<-lapply(1:length(merge.res),function(x){is.na(merge.res[[x]][1])[1]})
			l2<-unlist(l1)
			
			if(length(which(l2==TRUE))>0){
			merge.res<-merge.res[-which(l2==TRUE)]	
			}
			
			final.res={}
			
			
			
			final.res<-do.call(rbind,merge.res) #ldply(merge.res,rbind)
			
			
			
			final.res<-unique(final.res)
	
	}else{
		
			final.res<-merge.Results.child.ttest(curdata,max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,
			alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
		}
	
	final.res<-as.data.frame(final.res)
	#final.res<-na.omit(final.res)
	
	
	print("s5")
	s5=Sys.time()
	
	print((s5-s4))
	
	
	final.res<-unique(final.res)
	
	final.res<-final.res[order(final.res$mz),]
	
	
	print(dim(final.res))
	
#	save(curdata,file="G:\\Medicine\\Pulmonary_ISILON\\Research\\Jones_Lab\\MoTrPAC\\Pilot1_M293_Sep2017\\xMSanalyzervdebug32C\\curdata.Rda")
	
	
#	save(final.res,file="G:\\Medicine\\Pulmonary_ISILON\\Research\\Jones_Lab\\MoTrPAC\\Pilot1_M293_Sep2017\\xMSanalyzervdebug32C\\final.res2.Rda")
	
	#if(FALSE)
	{
		curdata<-final.res
		curdata<-as.data.frame(curdata)
		
		save(curdata,file="curdata.Rda")
		
		#if(FALSE)
		{
		
		#system.time(global_cor<-try(WGCNA::cor(t(curdata[,c(sample.col.start:end_ind)]),nThreads=numnodes,use = 'p'),silent=TRUE))
		
		global_cor<-WGCNA::cor(t(curdata[,c(sample.col.start:end_ind)]),nThreads=numnodes,use = "p")

		save(global_cor,file="global_cor.Rda")
		
		if(is(global_cor,"try-error")){

			non_cor_data<-curdata

		}else{

		#count_cor<-lapply(1:dim(curdata)[1],function(j){length(which(global_cor[j,]>=0.9))})
		
		#count_cor<-lapply(1:dim(curdata)[1],function(j){diff_rt<-abs(curdata[j,2]-curdata[,2]); length(which(global_cor[j,]>=0.9 & diff_rt<10))})
		
		if(FALSE){
			cl<-parallel::makeCluster(numnodes)
			
		
		
	
		clusterExport(cl,curdata)
		
		clusterExport(cl,max.rt.diff)
		clusterExport(cl,max.mz.diff)
		clusterExport(cl,sample.col.start)
		clusterExport(cl,end_ind)
		clusterExport(cl,global_cor)
		clusterExport(cl,mergecorthresh)
		}
		
		count_cor<-lapply(1:dim(curdata)[1],function(j){
			diff_mz<-10^6*abs(curdata[j,1]-curdata[,1])/curdata[j,1];
			diff_rt<-abs(curdata[j,2]-curdata[,2]);
			same_feat_ind<-which(diff_mz<max.mz.diff & diff_rt<max.rt.diff);
			same_feat_ind<-setdiff(same_feat_ind,j)
			
			if(length(same_feat_ind)>0){
			
				com_data<-curdata[same_feat_ind,c(sample.col.start:end_ind)]
				
				
					diff_mat<-sweep(com_data,1,t(curdata[j,c(sample.col.start:end_ind)]))
					diff_mat<-abs(diff_mat)
				
				if(length(diff_mat)>0){
					if(nrow(diff_mat)>=1){
						check_zeros<-apply(diff_mat,1,function(x){len_zero<-length(which(x==0));if(len_zero>=1){return(len_zero)}else{return(0)}})
					}else{
						check_zeros<-0
					}
				}else{
					check_zeros<-0
				}
				
				#more than one value same
				if(check_zeros>0){
					
					same_val_index<-same_feat_ind[which(check_zeros>0)]
					
					same_feat_ind1<-intersect(same_feat_ind,same_val_index)

					if(length(same_feat_ind1)>0){
					
					check_zeros<-length(which(global_cor[j,same_feat_ind1]>=mergecorthresh))  #max(c(length(which(global_cor[j,same_feat_ind]>=mergecorthresh))-1,check_zeros),na.rm=TRUE)
				
				
						if(check_zeros>0){
							dup_index<-c(j,same_feat_ind1)
							return(dup_index)
						}else{
							return((-1))
						}
					
					}else{
						return((-1))
					}
					
				}else{
					if(length(same_feat_ind)>0){
					
					check_zeros<-length(which(global_cor[j,same_feat_ind]>=mergecorthresh))  #max(c(length(which(global_cor[j,same_feat_ind]>=mergecorthresh))-1,check_zeros),na.rm=TRUE)
				
				
						if(check_zeros>0){
							dup_index<-c(j,same_feat_ind)
							return(dup_index)
						}else{
							return((-1))
						}
					
					}else{
						return((-1))
					}
				}
			}else{
				return((-1))
			}
		})
		
		print("here")
		print("length of count_cor")
		print(length(count_cor))
		save(count_cor,file="count_cor.Rda")
		

		#stopCluster(cl)
		#count_cor<-lapply(1:dim(curdata)[1],function(j){diff_mz<-10^6*abs(curdata[j,1]-curdata[,1])/curdata[j,1]; diff_rt<-abs(curdata[j,2]-curdata[,2]); length(which(diff_mz<max.mz.diff & diff_rt<max.rt.diff & global_cor[j,]>=mergecorthresh))})
		
			
		
		count_cor<-unlist(count_cor)
	
		
			if(length(which(count_cor<=0))>0){
				non_cor_data<-curdata[which(count_cor<=0),]
				curdata<-curdata[-which(count_cor<=0),]
				#global_cor<-global_cor[-which(count_cor<=1),-which(count_cor<=1)]
				
			}else{
				non_cor_data<-{}
			}

		}
		}
		
		

		
	print("s6")
	s6=Sys.time()
	
	print((s6-s5))
	if(dim(curdata)[1]>0){
	
			
			print("dim of curdata 2nd round check")
		       print(dim(curdata))
			final.res2<-merge.Results.child.ttest(curdata,max.mz.diff=max.mz.diff, max.rt.diff=max.rt.diff, merge.eval.pvalue=merge.eval.pvalue,alignment.tool=alignment.tool, mult.test.cor=mult.test.cor,mergecorthresh=mergecorthresh,col.rm.index=NA,missingvalue=missingvalue)
			#final.res2<-na.omit(final.res2)
			final.res2<-unique(final.res2)
		    
			final.res<-rbind(final.res2,non_cor_data)
			final.res<-as.data.frame(final.res)
			final.res<-final.res[order(final.res$mz),]
			final.res<-unique(final.res)	

	}
	print("s7")
	s7=Sys.time()
	
	print((s7-s6))
	curdata<-final.res
    
   	rm(curdata) 

	
	}
	
	if(is.na(missingvalue)==FALSE){
		
			final.res<-replace(as.matrix(final.res),which(is.na(final.res)==TRUE),missingvalue)
			
		}

	final.res<-as.data.frame(final.res)
	
	print("dim of final res")
    
    final.res<-unique(final.res)
	print(dim(final.res))
	
	options(warn=0)
	return(final.res)
	#return(sub_data_a)
}


		

  compare_intensities_ttestvind<-function(other_feats,y,merge.eval.pvalue, mergecorthresh=0.6){
		
											ttest_sum=apply(other_feats,1,function(x)
				                                                        {
				                                                        	
												diff_rt<-abs(x[1]-y[1])
												
												x<-x[-c(1)]
												y<-y[-c(1)]
												x<-as.matrix(x)
												y<-as.matrix(y)
				                                
				                                
												yind<-which(is.na(y)==TRUE)
												xind<-which(is.na(x)==TRUE)
												
												naind<-c(yind,xind)
												
												if(dim(x)[1]>dim(y)[1])
												{
													x<-t(as.numeric(x))
												}
												
												if(dim(y)[1]>dim(x)[1])
												{
													y<-t(as.numeric(y))
												}
												
												if(length(naind)>0){
												 x<-x[-naind]
												 y<-y[-naind]
				                                                                 
												 }
												
												 if(length(x)>1 & length(y)>1)
												 {
													num_com<-length(which(abs(x-y)<0.5))
													
													if(num_com<=length(y))
													{
														
												
													mean_x<-mean((x+1),na.rm=TRUE)
													mean_y<-mean((y+1),na.rm=TRUE)
													
											
													log2fc<-abs(log2(mean_x/mean_y))
													
			
													
													if(log2fc>0.1)
													{
													
													
													ttest_res=try(t.test(x,y,paired=TRUE,na.omit=TRUE),silent=TRUE)
													
													#print(ttest_res)
													if (is(ttest_res, "try-error"))
													{
														ttest_pval<-0
														ttest_pval2<-0

													}else
													{
														
														ttest_pval=ttest_res$p.value
														
														if(is.na(ttest_pval)==TRUE){
															print("here 1")
															ttest_pval=1
														}
														
														ttest_res2=try(t.test(x,y,paired=FALSE,na.omit=TRUE),silent=TRUE)
														ttest_pval2=ttest_res2$p.value
														if(is.na(ttest_pval2)==TRUE){
															print("here 2")
															ttest_pval2=1
														}
														
													}
													
													
													
													
													#if(ttest_pval<merge.eval.pvalue)
													
													
													#r6<-kruskal.test(list(x,y))
													
													#r10<-var.test(x,y)
													
													ttest_pval<-max(ttest_pval,ttest_pval2,na.rm=TRUE)
													
													if(length(x)>2 & length(y)>2)
													{
																cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="spearman",use="pairwise.complete.obs"),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																	cortest_pval=cortest_res$p.value
																	cortest_est=cortest_res$estimate
																	
																	if(cortest_est>=mergecorthresh){
																		
																		ttest_pval<-1
																		cortest_est<-1
																		cortest_pval<-0
																		
																	}
																	
																	if(cortest_est<mergecorthresh){
																		cortest_pval=1
																		
																					cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="pearson",use="pairwise.complete.obs"),silent=TRUE)
																				if (is(cortest_res, "try-error")){
																					
																				cortest_pval=1
																				cortest_est<-0
																				}else{
																					
																					cortest_pval=cortest_res$p.value
																					cortest_est=cortest_res$estimate
																					
																					if(cortest_est>=mergecorthresh){
																		
																						ttest_pval<-1
																						cortest_est<-1
																						cortest_pval<-0
																		
																					}
																					
																					if(cortest_est<mergecorthresh){
																						cortest_pval=1
																						cortest_est<-0
																					}
																					
																				}
																	}
																}
													
													
													
													if(diff_rt>30)
													{

													if(cortest_est>=mergecorthresh & ttest_pval>=merge.eval.pvalue)
													{
													
															print("here 4")
															print(cortest_est)
															print(ttest_pval)
															
															ttest_pval=1
													
													}else{
														ttest_pval=0
														}
													}else{
														
														if(cortest_est>=mergecorthresh | ttest_pval>=merge.eval.pvalue)
														{
													
																print("here 5")
																print(cortest_est)
																print(ttest_pval)
																ttest_pval=1
														}
														
														
														
														if(diff_rt<1)
														{
															print("here 6")
															print(diff_rt)
															ttest_pval=1
															
														}
															
														
														}
													
													}
													else{
														
														ttest_pval<-ttest_pval
													}
													#ttest_pval<-max(c(ttest_pval,r6$p.value,r10$p.value),na.rm=TRUE)
													
													
													}
													else
													{
													print("here 7")
														ttest_pval=1
													}
													
													
													}
													else
													{
													print("here 8")
													ttest_pval=1
													}											
													
													
				                                  }
				                                 else
				                                 {
				                                 	
				                                 	print("here 9")
														ttest_pval=1
													
												}
												
				                                                                return(ttest_pval)
				                                                        })
return(ttest_sum)
}

merge.Results.child.ttestv32C<-function(dataA, max.mz.diff=15, max.rt.diff=300, merge.eval.pvalue=0.05,alignment.tool="apLCMS", mult.test.cor=FALSE,mergecorthresh=0.7,col.rm.index=NA,missingvalue=0)
{	
        diff_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}
        cnames=colnames(dataA)
        dataA<-as.data.frame(dataA)
        if(nrow(dataA)<2){
        	
        	return(dataA)
        }
        
        if(is.na(col.rm.index)==FALSE){
       	 dataA<-as.data.frame(dataA[,-c(col.rm.index)]) 
		}else{
			 dataA<-as.data.frame(dataA)
			}
			
		#dataA<-as.data.frame(dataA,2,as.numeric)
         if(alignment.tool=="apLCMS")
        {
              sample.col.start=6
	      timeindex=2
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    cnames[1]="mz"

                    cnames[4]="time"
                    colnames(dataA)=cnames
                    
		    timeindex=4
		    
		    dataA[,dim(dataA)[2]]<-1
              }
              else
              {
                    stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }

		if(is.na(missingvalue)==FALSE){
		
			dataA<-replace(as.matrix(dataA),which(dataA==missingvalue),NA)
			
		}
		
		dataA<-as.data.frame(dataA)
  
compare_intensities_ttest<-function(other_feats,y,merge.eval.pvalue, mergecorthresh=0.6){
		
											ttest_sum=apply(other_feats,1,function(x)
				                                                        {
				                                                        	
												diff_rt<-abs(x[1]-y[1])
												
												x<-x[-c(1)]
												y<-y[-c(1)]
												x<-as.matrix(x)
												y<-as.matrix(y)
				                                
												if(diff_rt<2){
													
													return(1)
												}
												yind<-which(is.na(y)==TRUE)
												xind<-which(is.na(x)==TRUE)
												
												naind<-c(yind,xind)
												
												if(dim(x)[1]>dim(y)[1])
												{
													x<-t(as.numeric(x))
												}
												
												if(dim(y)[1]>dim(x)[1])
												{
													y<-t(as.numeric(y))
												}
												
												if(length(naind)>0){
												 x<-x[-naind]
												 y<-y[-naind]
				                                                                 
												 }
												
												 if(length(x)>1 & length(y)>1)
												 {
													num_com<-length(which(abs(x-y)<0.5))
													
													if(num_com<1)
													{
														
													
													
													mean_x<-mean((x+1),na.rm=TRUE)
													mean_y<-mean((y+1),na.rm=TRUE)
													
											
													log2fc<-abs(log2(mean_x/mean_y))
													
			
													
													if(log2fc>0.1)
													{
													
													
													ttest_res=try(t.test(x,y,paired=TRUE,na.omit=TRUE),silent=TRUE)
													
												
													if (is(ttest_res, "try-error"))
													{
														ttest_pval<-0
														ttest_pval2<-0

													}else
													{
														
														ttest_pval=ttest_res$p.value
														
														if(is.na(ttest_pval)==TRUE){
													
															ttest_pval=1
														}
														
														ttest_res2=try(t.test(x,y,paired=FALSE,na.omit=TRUE),silent=TRUE)
														ttest_pval2=ttest_res2$p.value
														if(is.na(ttest_pval2)==TRUE){
														
															ttest_pval2=1
														}
														
													}
													
													
													
													
													#if(ttest_pval<merge.eval.pvalue)
													
													
													#r6<-kruskal.test(list(x,y))
													
													#r10<-var.test(x,y)
													
													ttest_pval<-max(ttest_pval,ttest_pval2,na.rm=TRUE)
													
													if(length(x)>2 & length(y)>2)
													{
																cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="spearman",use="pairwise.complete.obs"),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																	cortest_pval=cortest_res$p.value
																	cortest_est=cortest_res$estimate
																	
																	if(cortest_est>=mergecorthresh){
																		
																		ttest_pval<-1
																		cortest_est<-1
																		cortest_pval<-0
																		
																	}
																	
																	if(cortest_est<mergecorthresh){
																		cortest_pval=1
																		
																				cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="pearson",use="pairwise.complete.obs"),silent=TRUE)
																				if (is(cortest_res, "try-error")){
																					
																				cortest_pval=1
																				cortest_est<-0
																				}else{
																					
																					cortest_pval=cortest_res$p.value
																					cortest_est=cortest_res$estimate
																					
																					if(cortest_est>=mergecorthresh){
																		
																							ttest_pval<-1
																							cortest_est<-1
																							cortest_pval<-0
																		
																					}
																					
																					if(cortest_est<mergecorthresh){
																						cortest_pval=1
																						cortest_est<-0
																					}
																					
																				}
																	}
																}
													
													
													
															if(diff_rt>30)
															{

																if(cortest_est>=mergecorthresh && ttest_pval>=merge.eval.pvalue)
																{
																
																
																		ttest_pval=1
																
																}else{
																		ttest_pval=0
																}
															}else{
																
																if(cortest_est>=mergecorthresh && ttest_pval>=merge.eval.pvalue)
																{
															
																		
																		
																		ttest_pval=1
																}
																else{
																	if(diff_rt<5 && cortest_est>=mergecorthresh)
																	{
																		
																		ttest_pval=1
																		
																	}
																	
																}
																
															}
													
													}
													else{
														
														ttest_pval<-ttest_pval
													}
													#ttest_pval<-max(c(ttest_pval,r6$p.value,r10$p.value),na.rm=TRUE)
													
													
													}
													else
													{
													
														ttest_pval=1
													}
													
													
													}
													else
													{
													
														ttest_pval=1
													}											
													
													
				                                  }
												 else
												 {
				                                 	
				                                 	
														ttest_pval=1
													
												}
												
				                                                                return(ttest_pval)
				                                                        })
							return(ttest_sum)
						}

#return(dataA)
  #Step 1 Group features by m/z
        mz_groups<-lapply(1:dim(dataA)[1],function(j){
                                commat={}
                                diffmz=new("list")
                                ppmb=(max.mz.diff)*(dataA$mz[j]/1000000)
                                getbind_same<-which(abs(dataA$mz-dataA$mz[j])<=ppmb)
                                diffmz[[diff_mz_num]]=getbind_same #dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
	})
	  
	#return(mz_groups)  
	
	del_list<-{}
	for(m in 1:length(mz_groups)){
	  
	   if((m%in%del_list)==FALSE){
		for(n in (m+1):length(mz_groups)){
		
			if(n>length(mz_groups)){
				break;
			}
			com1<-intersect(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
			if(length(com1)>0){
			
				mz_groups[[m]][[1]]<-c(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
				del_list<-c(del_list,n)
			}
		}
		
		mz_groups[[m]][[1]]<-unique(mz_groups[[m]][[1]])
		
		}
	  }
	  if(length(del_list)>0){
		 mz_groups<-mz_groups[-del_list]
		}
	 
	  
	mz_groups<-unique(mz_groups)
	
	diff_mz_num=1
	
	
	mz_groups<-lapply(1:length(mz_groups),function(j){
                                commat={}
                                diffmz=new("list")
                                getbind_same=mz_groups[[j]][[1]]
                                diffmz[[diff_mz_num]]=dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })
          
    fname<-paste("mz_groups",max.mz.diff,".Rda",sep="")
     
  # save(mz_groups,file=fname)

	if(mult.test.cor==TRUE){
		mz_group_size<-lapply(1:length(mz_groups),function(j){
		num_rows<-dim(mz_groups[[j]][[1]])
		n=num_rows[[1]][[1]]
		num_comp=(n*(n-1))/2
		})

	
		num_comparisons<-sum(unlist(mz_group_size))
	}else{
		num_comparisons=1
	}

	
	
        #Step 2 Sub-group features from Step 1 by Retention time
        #find the features with RT values within the defined range as compared to the query feature
         diffmat={} #length(mz_groups)
       # for(j in 1:length(mz_groups))
       diffmatres<-lapply(1:length(mz_groups),function(j)
        {
	
		
        	
    			
               temp_diff={}
                tempdata=mz_groups[[j]][[1]]
		
			
		#fname<-paste("mz_groups",round(tempdata$mz,4),tempdata$time,"_",j,".Rda",sep="")
     
             #  save(tempdata,file=fname)
		
                if(dim(tempdata)[1]>1)
                {
			
				
                                tempdata<-tempdata[order(tempdata$Qscore,decreasing=TRUE),]
				
				 tempdata<-unique(tempdata)
				rem_ind<-{}
				
				rem_ind<-which(rem_ind==TRUE)
								if(length(rem_ind)>0){
								
									
									 tempdata<-tempdata[-rem_ind,]
								}
				allrownames<-rownames(tempdata)
					
				dup_list<-{}
				temp_diff={}
				dup_ind<-{}
				curdupnames<-{}
				diffmat<-{}
                                for(d in 1:dim(tempdata)[1])
                                {
                                		
										if(d%in%dup_ind==FALSE & d%in%curdupnames==FALSE)
										{
										rowname<-rownames(tempdata[d,])                                      
										cur_feat=tempdata[d,]
										 other_feats=tempdata
									
										getbind_rtsame<-which(abs(other_feats$time-cur_feat$time)<=max.rt.diff)
					                                        commat={}
										 same_feat_ind={}
										 
											if(length(getbind_rtsame)>1)
											{
												other_feats<-other_feats[getbind_rtsame,]
		
												y=cur_feat[c(timeindex,sample.col.start:(dim(tempdata)[2]-9))]
						
												diff_rt<-abs(other_feats[,timeindex]-y[timeindex])

												ttest_sum<-compare_intensities_ttest(other_feats[,c(timeindex,sample.col.start:(dim(tempdata)[2]-9))],y,merge.eval.pvalue,mergecorthresh)
											
												same_feat_ind=which(ttest_sum>(merge.eval.pvalue/num_comparisons))
												same_feat_ind=c(same_feat_ind,which(is.na(ttest_sum)))
												dup_list<-c(dup_list,names(same_feat_ind))
												curnames<-names(same_feat_ind)
												curdupnames<-which(allrownames%in%curnames)
												
												curnames<-curnames[-which(curnames==rowname)]
												
												
													
												if(length(which(curdupnames%in%dup_ind==TRUE))<1)
													{
																dup_ind<-c(dup_ind,curdupnames)
														
																if(length(same_feat_ind)>1)
																{
																				same_ind<-as.numeric(same_feat_ind)						
																	other_feats<-as.data.frame(other_feats)
													
																	commat=other_feats[c(same_ind),]
																       maxint<-apply(commat,1,function(x){max(x,na.rm=TRUE)})
																	
																	maxint<-log((maxint+1),10)
																	maxint<-maxint+0.001
																	maxint<-maxint/max(maxint,na.rm=TRUE)
																	
																	commat_qscore<-as.numeric(commat$Qscore)*maxint #cv_range/as.numeric(commat$numgoodsamples)*as.numeric(deltappm_res)
																					#inverse of the actual qscore
																	best_level_index=which(as.numeric(commat_qscore)==max(as.numeric(commat_qscore),na.rm=TRUE))
																	if(length(best_level_index)>0)
																	{
																		best_level_index=best_level_index[1]
																	}
																	best_level_data=commat[best_level_index,]
																}else
																{
																	best_level_index=d
																	best_level_data=cur_feat
									
																}
																
																		
																		
																diffmat=rbind(diffmat, best_level_data)
																#diffmat<-as.data.frame(best_level_data)

																temp_diff=rbind(temp_diff,best_level_data)
														
																dup_ind<-c(dup_ind,best_level_index)
														
													} else{
                                											if(d%in%dup_ind==FALSE){
				                                               					 diffmat<-rbind(diffmat,cur_feat)
				                                               					 #diffmat<-as.data.frame(cur_feat)

				                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
																				  dup_ind<-c(dup_ind,d)
				                                              					  }
                                				}
												
					
                                }
                                else{
                                			if(d%in%dup_ind==FALSE){
				                                               					 diffmat<-rbind(diffmat,cur_feat)
				                                               					 #diffmat<-as.data.frame(cur_feat)

				                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
																				  dup_ind<-c(dup_ind,d)
				                                              					  }
                                	}
                               }
                               
							}#for loop
                	
                }
                else{
								cur_feat<-as.matrix(mz_groups[[j]][[1]])
                	
                        
                                               					 diffmat<-rbind(diffmat,cur_feat)
                                               					 # diffmat<-as.data.frame(cur_feat)
                                              					 temp_diff=rbind(temp_diff,tempdata)
                        	}
                        
                        
                        
                
                diffmat<-unique(diffmat)
                 best_level_data={}

		return(diffmat)
        }
        )
       
        diffmat<-do.call(rbind,diffmatres)
        diffmat=unique(diffmat)
        return(diffmat)
        
}



merge.Results.child.ttest<-function(dataA, max.mz.diff=15, max.rt.diff=300, merge.eval.pvalue=0.05,alignment.tool="apLCMS", mult.test.cor=FALSE,mergecorthresh=0.7,col.rm.index=NA,missingvalue=0)
{	
        diff_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}
        cnames=colnames(dataA)
        dataA<-as.data.frame(dataA)
        if(nrow(dataA)<2){
        	
        	return(dataA)
        }
        
        if(is.na(col.rm.index)==FALSE){
       	 dataA<-as.data.frame(dataA[,-c(col.rm.index)]) 
		}else{
			 dataA<-as.data.frame(dataA)
			}
			
		#dataA<-as.data.frame(dataA,2,as.numeric)
         if(alignment.tool=="apLCMS")
        {
              sample.col.start=6
	      timeindex=2
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    cnames[1]="mz"

                    cnames[4]="time"
                    colnames(dataA)=cnames
                    
		    timeindex=4
		    
		    dataA[,dim(dataA)[2]]<-1
              }
              else
              {
                    stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }

		if(is.na(missingvalue)==FALSE){
		
			dataA<-replace(as.matrix(dataA),which(dataA==missingvalue),NA)
			
		}
		
		dataA<-as.data.frame(dataA)
  
  compare_intensities_ttest<-function(other_feats,y,merge.eval.pvalue, mergecorthresh=0.6){
		
											ttest_sum=apply(other_feats,1,function(x)
				                                                        {
				                                                        	
												diff_rt<-abs(x[1]-y[1])
												
												x<-x[-c(1)]
												y<-y[-c(1)]
												x<-as.matrix(x)
												y<-as.matrix(y)
				                                
												if(diff_rt<2){
													
													return(1)
												}
												yind<-which(is.na(y)==TRUE)
												xind<-which(is.na(x)==TRUE)
												
												naind<-c(yind,xind)
												
												if(dim(x)[1]>dim(y)[1])
												{
													x<-t(as.numeric(x))
												}
												
												if(dim(y)[1]>dim(x)[1])
												{
													y<-t(as.numeric(y))
												}
												
												if(length(naind)>0){
												 x<-x[-naind]
												 y<-y[-naind]
				                                                                 
												 }
												
												 if(length(x)>1 & length(y)>1)
												 {
													num_com<-length(which(abs(x-y)<0.5))
													
													if(num_com<1)
													{
														
													
													
													mean_x<-mean((x+1),na.rm=TRUE)
													mean_y<-mean((y+1),na.rm=TRUE)
													
											
													log2fc<-abs(log2(mean_x/mean_y))
													
			
													
													if(log2fc>0.1)
													{
													
													
													ttest_res=try(t.test(x,y,paired=TRUE,na.omit=TRUE),silent=TRUE)
													
												
													if (is(ttest_res, "try-error"))
													{
														ttest_pval<-0
														ttest_pval2<-0

													}else
													{
														
														ttest_pval=ttest_res$p.value
														
														if(is.na(ttest_pval)==TRUE){
													
															ttest_pval=1
														}
														
														ttest_res2=try(t.test(x,y,paired=FALSE,na.omit=TRUE),silent=TRUE)
														ttest_pval2=ttest_res2$p.value
														if(is.na(ttest_pval2)==TRUE){
														
															ttest_pval2=1
														}
														
													}
													
													
													
													
													#if(ttest_pval<merge.eval.pvalue)
													
													
													#r6<-kruskal.test(list(x,y))
													
													#r10<-var.test(x,y)
													
													ttest_pval<-max(ttest_pval,ttest_pval2,na.rm=TRUE)
													
													if(length(x)>2 & length(y)>2)
													{
																cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="spearman",use="pairwise.complete.obs"),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																	cortest_pval=cortest_res$p.value
																	cortest_est=cortest_res$estimate
																	
																	if(cortest_est>=mergecorthresh){
																		
																		ttest_pval<-1
																		cortest_est<-1
																		cortest_pval<-0
																		
																	}
																	
																	if(cortest_est<mergecorthresh){
																		cortest_pval=1
																		
																					cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="pearson",use="pairwise.complete.obs"),silent=TRUE)
																				if (is(cortest_res, "try-error")){
																					
																				cortest_pval=1
																				cortest_est<-0
																				}else{
																					
																					cortest_pval=cortest_res$p.value
																					cortest_est=cortest_res$estimate
																					
																					if(cortest_est>=mergecorthresh){
																		
																		ttest_pval<-1
																		cortest_est<-1
																		cortest_pval<-0
																		
																	}
																					
																					if(cortest_est<mergecorthresh){
																						cortest_pval=1
																						cortest_est<-0
																					}
																					
																					}
																	}
																}
													
													
													
													if(diff_rt>30)
													{

													if(cortest_est>=mergecorthresh & ttest_pval>=merge.eval.pvalue)
                                                    {
													
													
															ttest_pval=1
													
													}else{
														ttest_pval=0
														}
													}else{
														
														if(cortest_est>=mergecorthresh | ttest_pval>=merge.eval.pvalue)
														{
													
																
																
																ttest_pval=1
														}
														
														
														
														if(diff_rt<2)
														{
															
															ttest_pval=1
															
														}
															
														
														}
													
													}
													else{
														
														ttest_pval<-ttest_pval
													}
													#ttest_pval<-max(c(ttest_pval,r6$p.value,r10$p.value),na.rm=TRUE)
													
													
													}
													else
													{
													
														ttest_pval=1
													}
													
													
													}
													else
													{
													
													ttest_pval=1
													}											
													
													
				                                  }
				                                 else
				                                 {
				                                 	
				                                 	
														ttest_pval=1
													
												}
												
				                                                                return(ttest_pval)
				                                                        })
return(ttest_sum)
}

#return(dataA)
  #Step 1 Group features by m/z
        mz_groups<-lapply(1:dim(dataA)[1],function(j){
                                commat={}
                                diffmz=new("list")
                                ppmb=(max.mz.diff)*(dataA$mz[j]/1000000)
                                getbind_same<-which(abs(dataA$mz-dataA$mz[j])<=ppmb)
                                diffmz[[diff_mz_num]]=getbind_same #dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
	})
	  
	#return(mz_groups)  
	
	del_list<-{}
	for(m in 1:length(mz_groups)){
	  
	   if((m%in%del_list)==FALSE){
		for(n in (m+1):length(mz_groups)){
		
			if(n>length(mz_groups)){
				break;
			}
			com1<-intersect(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
			if(length(com1)>0){
			
				mz_groups[[m]][[1]]<-c(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
				del_list<-c(del_list,n)
			}
		}
		
		mz_groups[[m]][[1]]<-unique(mz_groups[[m]][[1]])
		
		}
	  }
	  if(length(del_list)>0){
		 mz_groups<-mz_groups[-del_list]
		}
	 
	  
	mz_groups<-unique(mz_groups)
	
	diff_mz_num=1
	
	#cl1<-parallel::makeCluster(getOption("cl.cores", 2))
	#clusterExport(cl1,mz_groups)
	#clusterExport(cl1,dataA)
	
	mz_groups<-lapply(1:length(mz_groups),function(j){
                                commat={}
                                diffmz=new("list")
                                getbind_same=mz_groups[[j]][[1]]
                                diffmz[[diff_mz_num]]=dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })
          
	#  stopCluster(cl1)
	  
    fname<-paste("mz_groups",max.mz.diff,".Rda",sep="")
     
  # save(mz_groups,file=fname)

	if(mult.test.cor==TRUE){
		mz_group_size<-lapply(1:length(mz_groups),function(j){
		num_rows<-dim(mz_groups[[j]][[1]])
		n=num_rows[[1]][[1]]
		num_comp=(n*(n-1))/2
		})

	
		num_comparisons<-sum(unlist(mz_group_size))
	}else{
		num_comparisons=1
	}

	
	
        #Step 2 Sub-group features from Step 1 by Retention time
        #find the features with RT values within the defined range as compared to the query feature
         diffmat={} #length(mz_groups)
       # for(j in 1:length(mz_groups))
       diffmatres<-lapply(1:length(mz_groups),function(j)
        {
	
		
        	
    			
               temp_diff={}
                tempdata=mz_groups[[j]][[1]]
		
			
		#fname<-paste("mz_groups",round(tempdata$mz,4),tempdata$time,"_",j,".Rda",sep="")
     
             #  save(tempdata,file=fname)
		
                if(dim(tempdata)[1]>1)
                {
			
				
                                tempdata<-tempdata[order(tempdata$Qscore,decreasing=TRUE),]
				
				 tempdata<-unique(tempdata)
				rem_ind<-{}
				
				rem_ind<-which(rem_ind==TRUE)
								if(length(rem_ind)>0){
								
									
									 tempdata<-tempdata[-rem_ind,]
								}
				allrownames<-rownames(tempdata)
					
				dup_list<-{}
				temp_diff={}
				dup_ind<-{}
				curdupnames<-{}
				diffmat<-{}
                                for(d in 1:dim(tempdata)[1])
                                {
                                		
										if(d%in%dup_ind==FALSE & d%in%curdupnames==FALSE)
										{
										rowname<-rownames(tempdata[d,])                                      
										cur_feat=tempdata[d,]
										 other_feats=tempdata
									
										getbind_rtsame<-which(abs(other_feats$time-cur_feat$time)<=max.rt.diff)
					                                        commat={}
										 same_feat_ind={}
										 
											if(length(getbind_rtsame)>1)
											{
												other_feats<-other_feats[getbind_rtsame,]
		
												y=cur_feat[c(timeindex,sample.col.start:(dim(tempdata)[2]-9))]
						
												diff_rt<-abs(other_feats[,timeindex]-y[timeindex])

												ttest_sum<-compare_intensities_ttest(other_feats[,c(timeindex,sample.col.start:(dim(tempdata)[2]-9))],y,merge.eval.pvalue,mergecorthresh)
											
												same_feat_ind=which(ttest_sum>(merge.eval.pvalue/num_comparisons))
												same_feat_ind=c(same_feat_ind,which(is.na(ttest_sum)))
												dup_list<-c(dup_list,names(same_feat_ind))
												curnames<-names(same_feat_ind)
												curdupnames<-which(allrownames%in%curnames)
												
												curnames<-curnames[-which(curnames==rowname)]
												
												
													
												if(length(which(curdupnames%in%dup_ind==TRUE))<1)
													{
																dup_ind<-c(dup_ind,curdupnames)
														
					                                                        if(length(same_feat_ind)>1)
					                                                        {
					                                                                                        same_ind<-as.numeric(same_feat_ind)						
													other_feats<-as.data.frame(other_feats)
									
					                                                                commat=other_feats[c(same_ind),]
					                                                               maxint<-apply(commat,1,function(x){max(x,na.rm=TRUE)})
					                                                                
					                                                                maxint<-log((maxint+1),10)
					                                                                maxint<-maxint+0.001
					                                                                maxint<-maxint/max(maxint,na.rm=TRUE)
					                                                                
						                                           commat_qscore<-as.numeric(commat$Qscore)*maxint #cv_range/as.numeric(commat$numgoodsamples)*as.numeric(deltappm_res)
																	#inverse of the actual qscore
													best_level_index=which(as.numeric(commat_qscore)==max(as.numeric(commat_qscore),na.rm=TRUE))
					                                                                if(length(best_level_index)>0)
					                                                                {
					                                                                        best_level_index=best_level_index[1]
					                                                                }
					                                                                best_level_data=commat[best_level_index,]
					                                                        }else
					                                                        {
																					best_level_index=d
					                                                                best_level_data=cur_feat
					
					                                                        }
												
														
														
					                                                                       diffmat=rbind(diffmat, best_level_data)
					                                                                        #diffmat<-as.data.frame(best_level_data)

					                                                                        temp_diff=rbind(temp_diff,best_level_data)
														
																							dup_ind<-c(dup_ind,best_level_index)
														
													} else{
                                											if(d%in%dup_ind==FALSE){
				                                               					 diffmat<-rbind(diffmat,cur_feat)
				                                               					 #diffmat<-as.data.frame(cur_feat)

				                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
																				  dup_ind<-c(dup_ind,d)
				                                              					  }
                                				}
												
					
                                }
                                else{
                                			if(d%in%dup_ind==FALSE){
				                                               					 diffmat<-rbind(diffmat,cur_feat)
				                                               					 #diffmat<-as.data.frame(cur_feat)

				                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
																				  dup_ind<-c(dup_ind,d)
				                                              					  }
                                	}
                               }
                               
							}#for loop
                	
                }
                else{
								cur_feat<-as.matrix(mz_groups[[j]][[1]])
                	
                        
                                               					 diffmat<-rbind(diffmat,cur_feat)
                                               					 # diffmat<-as.data.frame(cur_feat)
                                              					 temp_diff=rbind(temp_diff,tempdata)
                        	}
                        
                        
                        
                
                diffmat<-unique(diffmat)
                 best_level_data={}

		return(diffmat)
        }
        )
       
        diffmat<-do.call(rbind,diffmatres)
        diffmat=unique(diffmat)
        return(diffmat)
        
}

merge.Results.child.ttestv32H<-function(dataA, max.mz.diff=15, max.rt.diff=300, merge.eval.pvalue=0.05,alignment.tool="apLCMS", mult.test.cor=FALSE,mergecorthresh=0.7,col.rm.index=NA,missingvalue=0)
{	
        diff_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}
        cnames=colnames(dataA)
        dataA<-as.data.frame(dataA)
        if(nrow(dataA)<2){
        	
        	return(dataA)
        }
        
        if(is.na(col.rm.index)==FALSE){
       	 dataA<-as.data.frame(dataA[,-c(col.rm.index)]) 
		}else{
			 dataA<-as.data.frame(dataA)
			}
			
		#dataA<-as.data.frame(dataA,2,as.numeric)
         if(alignment.tool=="apLCMS")
        {
              sample.col.start=6
	      timeindex=2
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    cnames[1]="mz"

                    cnames[4]="time"
                    colnames(dataA)=cnames
                    
		    timeindex=4
		    
		    dataA[,dim(dataA)[2]]<-1
              }
              else
              {
                    stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }

		if(is.na(missingvalue)==FALSE){
		
			dataA<-replace(as.matrix(dataA),which(dataA==missingvalue),NA)
			
		}
		
		dataA<-as.data.frame(dataA)
  
  compare_intensities_ttest<-function(other_feats,y,merge.eval.pvalue, mergecorthresh=0.6){
		
											ttest_sum=apply(other_feats,1,function(x)
				                                                        {
				                                                        	
												diff_rt<-abs(x[1]-y[1])
												
												x<-x[-c(1)]
												y<-y[-c(1)]
												x<-as.matrix(x)
												y<-as.matrix(y)
				                                
												if(diff_rt<2){
													
													return(1)
												}
												yind<-which(is.na(y)==TRUE)
												xind<-which(is.na(x)==TRUE)
												
												naind<-c(yind,xind)
												
												if(dim(x)[1]>dim(y)[1])
												{
													x<-t(as.numeric(x))
												}
												
												if(dim(y)[1]>dim(x)[1])
												{
													y<-t(as.numeric(y))
												}
												
												if(length(naind)>0){
												 x<-x[-naind]
												 y<-y[-naind]
				                                                                 
												 }
												
												 if(length(x)>1 & length(y)>1)
												 {
													num_com<-length(which(abs(x-y)<0.5))
													
													if(num_com<1)
													{
														
													
													
													mean_x<-mean((x+1),na.rm=TRUE)
													mean_y<-mean((y+1),na.rm=TRUE)
													
											
													log2fc<-abs(log2(mean_x/mean_y))
													
			
													
													if(log2fc>0.1)
													{
													
													
													ttest_res=try(t.test(x,y,paired=TRUE,na.omit=TRUE),silent=TRUE)
													
												
													if (is(ttest_res, "try-error"))
													{
														ttest_pval<-0
														ttest_pval2<-0

													}else
													{
														
														ttest_pval=ttest_res$p.value
														
														if(is.na(ttest_pval)==TRUE){
													
															ttest_pval=1
														}
														
														ttest_res2=try(t.test(x,y,paired=FALSE,na.omit=TRUE),silent=TRUE)
														ttest_pval2=ttest_res2$p.value
														if(is.na(ttest_pval2)==TRUE){
														
															ttest_pval2=1
														}
														
													}
													
													
													
													
													#if(ttest_pval<merge.eval.pvalue)
													
													
													#r6<-kruskal.test(list(x,y))
													
													#r10<-var.test(x,y)
													
													ttest_pval<-max(ttest_pval,ttest_pval2,na.rm=TRUE)
													
													if(length(x)>2 & length(y)>2)
													{
																cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="spearman",use="pairwise.complete.obs"),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																	cortest_pval=cortest_res$p.value
																	cortest_est=cortest_res$estimate
																	
																	if(cortest_est>=mergecorthresh){
																		
																		ttest_pval<-1
																		cortest_est<-1
																		cortest_pval<-0
																		
																	}
																	
																	if(cortest_est<mergecorthresh){
																		cortest_pval=1
																		
																					cortest_res=try(cor.test(as.numeric(x),as.numeric(y),method="pearson",use="pairwise.complete.obs"),silent=TRUE)
																				if (is(cortest_res, "try-error")){
																					
																				cortest_pval=1
																				cortest_est<-0
																				}else{
																					
																					cortest_pval=cortest_res$p.value
																					cortest_est=cortest_res$estimate
																					
																					if(cortest_est>=mergecorthresh){
																		
																		ttest_pval<-1
																		cortest_est<-1
																		cortest_pval<-0
																		
																	}
																					
																					if(cortest_est<mergecorthresh){
																						cortest_pval=1
																						cortest_est<-0
																					}
																					
																					}
																	}
																}
													
													
													
													if(diff_rt>30)
													{

														if(cortest_est>=mergecorthresh & ttest_pval>=merge.eval.pvalue)
														{
														
														
																ttest_pval=1
														
														}else{
															ttest_pval=0
														}
													}else{
														
														if(cortest_est>=mergecorthresh | ttest_pval>=merge.eval.pvalue)
														{
													
																
																
																ttest_pval=1
														}
														
														
														
														if(diff_rt<2)
														{
															
															ttest_pval=1
															
														}
															
														
														}
													
													}
													else{
														
														ttest_pval<-ttest_pval
													}
													#ttest_pval<-max(c(ttest_pval,r6$p.value,r10$p.value),na.rm=TRUE)
													
													
													}
													else
													{
													
														ttest_pval=1
													}
													
													
													}
													else
													{
													
													ttest_pval=1
													}											
													
													
				                                  }
				                                 else
				                                 {
				                                 	
				                                 	
														ttest_pval=1
													
												}
												
				                                                                return(ttest_pval)
				                                                        })
return(ttest_sum)
}

#return(dataA)
  #Step 1 Group features by m/z
        mz_groups<-lapply(1:dim(dataA)[1],function(j){
                                commat={}
                                diffmz=new("list")
                                ppmb=(max.mz.diff)*(dataA$mz[j]/1000000)
                                getbind_same<-which(abs(dataA$mz-dataA$mz[j])<=ppmb)
                                diffmz[[diff_mz_num]]=getbind_same #dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
	})
	  
	#return(mz_groups)  
	
	del_list<-{}
	for(m in 1:length(mz_groups)){
	  
	   if((m%in%del_list)==FALSE){
		for(n in (m+1):length(mz_groups)){
		
			if(n>length(mz_groups)){
				break;
			}
			com1<-intersect(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
			if(length(com1)>0){
			
				mz_groups[[m]][[1]]<-c(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
				del_list<-c(del_list,n)
			}
		}
		
		mz_groups[[m]][[1]]<-unique(mz_groups[[m]][[1]])
		
		}
	  }
	  if(length(del_list)>0){
		 mz_groups<-mz_groups[-del_list]
		}
	 
	  
	mz_groups<-unique(mz_groups)
	
	diff_mz_num=1
	
	
	mz_groups<-lapply(1:length(mz_groups),function(j){
                                commat={}
                                diffmz=new("list")
                                getbind_same=mz_groups[[j]][[1]]
                                diffmz[[diff_mz_num]]=dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })
          
    fname<-paste("mz_groups",max.mz.diff,".Rda",sep="")
     
  # save(mz_groups,file=fname)

	if(mult.test.cor==TRUE){
		mz_group_size<-lapply(1:length(mz_groups),function(j){
		num_rows<-dim(mz_groups[[j]][[1]])
		n=num_rows[[1]][[1]]
		num_comp=(n*(n-1))/2
		})

	
		num_comparisons<-sum(unlist(mz_group_size))
	}else{
		num_comparisons=1
	}

	
	
        #Step 2 Sub-group features from Step 1 by Retention time
        #find the features with RT values within the defined range as compared to the query feature
         diffmat={} #length(mz_groups)
       # for(j in 1:length(mz_groups))
       diffmatres<-lapply(1:length(mz_groups),function(j)
        {
	
		
        	
    			
               temp_diff={}
                tempdata=mz_groups[[j]][[1]]
		
			
		#fname<-paste("mz_groups",round(tempdata$mz,4),tempdata$time,"_",j,".Rda",sep="")
     
             #  save(tempdata,file=fname)
		
                if(dim(tempdata)[1]>1)
                {
			
				
                                tempdata<-tempdata[order(tempdata$Qscore,decreasing=TRUE),]
				
				 tempdata<-unique(tempdata)
				rem_ind<-{}
				
				rem_ind<-which(rem_ind==TRUE)
								if(length(rem_ind)>0){
								
									
									 tempdata<-tempdata[-rem_ind,]
								}
				allrownames<-rownames(tempdata)
					
				dup_list<-{}
				temp_diff={}
				dup_ind<-{}
				curdupnames<-{}
				diffmat<-{}
                                for(d in 1:dim(tempdata)[1])
                                {
                                		
										if(d%in%dup_ind==FALSE & d%in%curdupnames==FALSE)
										{
										rowname<-rownames(tempdata[d,])                                      
										cur_feat=tempdata[d,]
										 other_feats=tempdata
									
										getbind_rtsame<-which(abs(other_feats$time-cur_feat$time)<=max.rt.diff)
					                                        commat={}
										 same_feat_ind={}
										 
											if(length(getbind_rtsame)>1)
											{
												other_feats<-other_feats[getbind_rtsame,]
		
												y=cur_feat[c(timeindex,sample.col.start:(dim(tempdata)[2]-9))]
						
												diff_rt<-abs(other_feats[,timeindex]-y[timeindex])

												ttest_sum<-compare_intensities_ttest(other_feats[,c(timeindex,sample.col.start:(dim(tempdata)[2]-9))],y,merge.eval.pvalue,mergecorthresh)
											
												same_feat_ind=which(ttest_sum>(merge.eval.pvalue/num_comparisons))
												same_feat_ind=c(same_feat_ind,which(is.na(ttest_sum)))
												dup_list<-c(dup_list,names(same_feat_ind))
												curnames<-names(same_feat_ind)
												curdupnames<-which(allrownames%in%curnames)
												
												curnames<-curnames[-which(curnames==rowname)]
												
												
													
												if(length(which(curdupnames%in%dup_ind==TRUE))<1)
													{
																dup_ind<-c(dup_ind,curdupnames)
														
					                                                        if(length(same_feat_ind)>1)
					                                                        {
					                                                                                        same_ind<-as.numeric(same_feat_ind)						
													other_feats<-as.data.frame(other_feats)
									
					                                                                commat=other_feats[c(same_ind),]
					                                                               maxint<-apply(commat,1,function(x){max(x,na.rm=TRUE)})
					                                                                
					                                                                maxint<-log((maxint+1),10)
					                                                                maxint<-maxint+0.001
					                                                                maxint<-maxint/max(maxint,na.rm=TRUE)
					                                                                
						                                           commat_qscore<-as.numeric(commat$Qscore)*maxint #cv_range/as.numeric(commat$numgoodsamples)*as.numeric(deltappm_res)
																	#inverse of the actual qscore
													best_level_index=which(as.numeric(commat_qscore)==max(as.numeric(commat_qscore),na.rm=TRUE))
					                                                                if(length(best_level_index)>0)
					                                                                {
					                                                                        best_level_index=best_level_index[1]
					                                                                }
					                                                                best_level_data=commat[best_level_index,]
					                                                        }else
					                                                        {
																					best_level_index=d
					                                                                best_level_data=cur_feat
					
					                                                        }
												
														
														
					                                                                       diffmat=rbind(diffmat, best_level_data)
					                                                                        #diffmat<-as.data.frame(best_level_data)

					                                                                        temp_diff=rbind(temp_diff,best_level_data)
														
																							dup_ind<-c(dup_ind,best_level_index)
														
													} else{
                                											if(d%in%dup_ind==FALSE){
				                                               					 diffmat<-rbind(diffmat,cur_feat)
				                                               					 #diffmat<-as.data.frame(cur_feat)

				                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
																				  dup_ind<-c(dup_ind,d)
				                                              					  }
                                				}
												
					
                                }
                                else{
                                			if(d%in%dup_ind==FALSE){
				                                               					 diffmat<-rbind(diffmat,cur_feat)
				                                               					 #diffmat<-as.data.frame(cur_feat)

				                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
																				  dup_ind<-c(dup_ind,d)
				                                              					  }
                                	}
                               }
                               
							}#for loop
                	
                }
                else{
								cur_feat<-as.matrix(mz_groups[[j]][[1]])
                	
                        
                                               					 diffmat<-rbind(diffmat,cur_feat)
                                               					 # diffmat<-as.data.frame(cur_feat)
                                              					 temp_diff=rbind(temp_diff,tempdata)
                        	}
                        
                        
                        
                
                diffmat<-unique(diffmat)
                 best_level_data={}

		return(diffmat)
        }
        )
       
        diffmat<-do.call(rbind,diffmatres)
        diffmat=unique(diffmat)
        return(diffmat)
        
}


   	# 'expression_xls' is the expression index file (e.g. outputted by dChip); 'sample_info_file' is a tab-delimited text file containing the colums: Array  name, sample name, Batch, and any other covariates to be included in the modeling; 'type' currently supports two data file types 'txt' for a tab-delimited text file and 'csv' for an Excel .csv file (sometimes R handles the .csv file better, so use this if you have problems with a .txt file!); 'write' if 'T' ComBat writes adjusted data to a file, and if 'F' and ComBat outputs the adjusted data matrix if 'F' (so assign it to an object! i.e. NewData <- ComBat('my expression.xls','Sample info file.txt', write=F)); 'covariates=all' will use all of the columns in your sample info file in the modeling (except array/sample name), if you only want use a some of the columns in your sample info file, specify these columns here as a vector (you must include the Batch column in this list); 'par.prior' if 'T' uses the parametric adjustments, if 'F' uses the nonparametric adjustments--if you are unsure what to use, try the parametric adjustments (they run faster) and check the plots to see if these priors are reasonable; 'filter=value' filters the features with absent calls in > 1-value of the samples. The defaut here (as well as in dchip) is .8. Filter if you can as the EB adjustments work better after filtering. Filter must be numeric if your expression index file contains presence/absence calls (but you can set it >1 if you don't want to filter any features) and must be 'F' if your data doesn't have presence/absence calls; 'skip' is the number of columns that contain probe names and gene information, so 'skip=5' implies the first expression values are in column 6; 'prior.plots' if true will give prior plots with black as a kernal estimate of the empirical batch effect density and red as the parametric estimate. 
#type='txt'; write=T; covariates='all'; par.prior=T; filter=F; skip=0; prior.plots=T
MetabComBat <- function(dat, saminfo, type='txt', write=T, covariates='all', par.prior=T, filter=F, skip=0, prior.plots=T){
	#debug: expression_xls='exp.txt'; sample_info_file='sam.txt'; type='txt'; write=T; covariates='all'; par.prior=T; filter=F; skip=0; prior.plots=T
	#cat('Reading Sample Information File\n')
	#saminfo <- read.table(sample_info_file, header=T, sep='\t',comment.char='')
	if(sum(colnames(saminfo)=="Batch")!=1){stop('ERROR: Sample Information File does not have a Batch column!');
		return('ERROR: Sample Information File does not have a Batch column!')}
		
		if (skip>0){
              geneinfo <- as.matrix(dat[,1:skip])
              dat <- dat[,-c(1:skip)]
	}else{geneinfo=dat[,c(1:4)]; dat<-dat[,-c(1:4)]}
        #print(geneinfo[1:4])
        print(dat[1:2,1:3])
	
	print(dim(dat))
	
	print(dim(saminfo))
	
	
	if(filter){
		nfeatures <- nrow(dat)
		col <- ncol(dat)/2
		present <- apply(dat, 1, filter.absent, filter)
		dat <- dat[present, -(2*(1:col))]
		if (skip>0){geneinfo <- geneinfo[present,]}
		cat('Filtered features absent in more than',filter,'of samples. features remaining:',nrow(dat),'; features filtered:',nfeatures-nrow(dat),'\n')
		}

	if(any(apply(dat,2,mode)!='numeric')){stop('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option');return('ERROR: Array expression columns contain non-numeric values! (Check your .xls file for non-numeric values and if this is not the problem, make a .csv file and use the type=csv option)')}
	
	#cnames1<-colnames(dat)
	tmp <- match(colnames(dat),saminfo[,1])
	#tmp<-which(cnames1%in%saminfo[,1])
	#tmplen<-length(tmp)
	#if(tmplen!=length(cnames1)){return('ERROR: Sample Information File and Data Array Names are not the same!')}
	if(any(is.na(tmp))){stop('ERROR: Sample Information File and Data Array Names are not the same!');return('ERROR: Sample Information File and Data Array Names are not the same!')}
	
	if(ncol(dat)!=nrow(saminfo)){stop('ERROR: Sample Information File and Data Array Names are not the same!');return('ERROR: Sample Information File and Data Array Names are not the same!')}
	
	#tmp1 <- match(saminfo[,1],colnames(dat))
	#saminfo <- saminfo[tmp1[!is.na(tmp1)],]
	saminfo <- saminfo[tmp,]  ## Bug fixed 01/04/2011		

	if(any(covariates != 'all')){saminfo <- saminfo[,c(1:2,covariates)]}
	design1 <- design.mat(saminfo)	


	batches <- list.batch(saminfo)
	n.batch <- length(batches)
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)
	
	## Check for missing values
	NAs = any(is.na(dat))
	if(NAs){cat(c('Found',sum(is.na(dat)),'Missing Data Values\n'),sep=' ')}
        #print(dat[1:2,])
	##Standardize Data across features
	cat('Standardizing Data across features\n')
	
	#save(dat,file="dat.Rda")
	#save(design1,file="design1.Rda")
	
	if (!NAs){B.hat <- solve(t(design1)%*%design1)%*%t(design1)%*%t(as.matrix(dat))
	
	}else{B.hat=apply(dat,1,Beta.NA,design1)
	#print(length(B.hat))
	
	#print(dim(B.hat[[1]]))
	
	#save(B.hat,file="bhat.Rda")
	B.hat_orig<-B.hat
	
	
	numfeats<-length(B.hat)/length(n.batches)
	#print(numfeats)
	#print(n.batches)
	dim(B.hat)<-dim(matrix(0,nrow=length(n.batches),ncol=numfeats))
	
	} #Standarization Model
	
	#rows: number of batches; columns: number of metabs
	print("standardization done")
	print(dim(B.hat))
	
	
	grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
	if (!NAs){var.pooled <- ((dat-t(design1%*%B.hat))^2)%*%rep(1/n.array,n.array)}else{var.pooled <- apply(dat-t(design1%*%B.hat),1,var,na.rm=T)}

	stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
	if(!is.null(design1)){tmp <- design1;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
	s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))

	##Get regression batch effect parameters
	cat("Fitting L/S model and finding priors\n")
	batch.design1 <- design1[,1:n.batch]
	if (!NAs){gamma.hat <- solve(t(batch.design1)%*%batch.design1)%*%t(batch.design1)%*%t(as.matrix(s.data))}else{gamma.hat=apply(s.data,1,Beta.NA,batch.design1)}
	delta.hat <- NULL
	for (i in batches){
		delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
		}

	delta.hat<-replace(delta.hat,which(is.na(delta.hat)==TRUE),1)
	
			
	##Find Priors
	gamma.bar <- apply(gamma.hat, 1, mean)
	t2 <- apply(gamma.hat, 1, var)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)

	
	##Plot empirical and parametric priors

	if (prior.plots & par.prior){
		par(mfrow=c(2,2))
		tmp <- density(gamma.hat[1,])
		plot(tmp,  type='l', main="Density Plot")
		xx <- seq(min(tmp$x), max(tmp$x), length=100)
		lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
		qqnorm(gamma.hat[1,])	
		qqline(gamma.hat[1,], col=2)	
	
		tmp <- density(delta.hat[1,])
		invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
		tmp1 <- density(invgam)
		plot(tmp,  typ='l', main="Density Plot", ylim=c(0,max(tmp$y,tmp1$y)))
		lines(tmp1, col=2)
		qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
		lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
		title('Q-Q Plot')
	}
	
	##Find EB batch adjustments

	gamma.star <- delta.star <- NULL
	if(par.prior){
		cat("Finding parametric adjustments\n")
		for (i in 1:n.batch){
			#delta.hat[i,]<-replace(delta.hat[i,],which(is.na(delta.hat[i,])==TRUE),1)
			
			

			temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
	}else{
		cat("Finding nonparametric adjustments\n")
		for (i in 1:n.batch){
			bad_cols<-which(is.na(delta.hat[i,])==TRUE)
			#delta.hat[i,]<-replace(delta.hat[i,],which(is.na(delta.hat[i,])==TRUE),1)
			
			
			temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
		}


	### Normalize the Data ###
	cat("Adjusting the Data\n")

	bayesdata <- s.data
	j <- 1
	for (i in batches){
		bayesdata[,i] <- (bayesdata[,i]-t(batch.design1[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
		}

	bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean
	if(write){
		output_file <- paste('Adjusted_data_parprior',par.prior,'.xls',sep='_')
                 #print(geneinfo[1:2])
                 #print(bayesdata[1:2,1:4])
		 #cat(c(colnames(geneinfo),colnames(dat),'\n'),file=output_file,sep='\t')
		#suppressWarnings(write.table(cbind(geneinfo,formatC(as.matrix(bayesdata), format = "f")), file=output_file, sep="\t", quote=F,row.names=F,col.names=F,append=T))
                outdata <- cbind(ProbeID=geneinfo, bayesdata); write.table(outdata, file=output_file, sep="\t")
		cat("Adjusted data saved in file:",output_file,"\n")
		}else{return(cbind(geneinfo,bayesdata))}
	}

# filters data based on presence/absence call
filter.absent <- function(x,pct){
	present <- T
	col <- length(x)/2
	pct.absent <- (sum(x[2*(1:col)]=="A") + sum(x[2*(1:col)]=="M"))/col
	if(pct.absent > pct){present <- F}
	present
	}

# Next two functions make the design matrix (X) from the sample info file 
build.design <- function(vec, des=NULL, start=2){
	tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
	for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
	cbind(des,tmp)
	}

design.mat <- function(saminfo){
	tmp <- which(colnames(saminfo) == 'Batch')
	tmp1 <- as.factor(saminfo[,tmp])
	cat("Found",nlevels(tmp1),'batches\n')
	design <- build.design(tmp1,start=1)
	ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
	cat("Found",ncov,'covariate(s)\n')
	if(ncov>0){
		for (j in 1:ncov){
			tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
			design <- build.design(tmp1,des=design)
			}
		}
	design
	}

# Makes a list with elements pointing to which array belongs to which batch
list.batch <- function(saminfo){
	tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Batch')])
	batches <- NULL
	for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
	batches
	}

# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
trim.dat <- function(dat){
	tmp <- strsplit(colnames(dat),'\\.')
	tr <- NULL
	for (i in 1:length(tmp)){tr <- c(tr,tmp[[i]][1]!='X')}
	tr
	}

# Following four find empirical hyper-prior values
aprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (2*s2+m^2)/s2}
bprior <- function(gamma.hat){m=mean(gamma.hat); s2=var(gamma.hat); (m*s2+m^3)/s2}
postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}


# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
	n <- apply(!is.na(sdat),1,sum)
	g.old <- g.hat
	d.old <- d.hat
	change <- 1
	count <- 0
	#write.table(change,file="change.txt",sep="\t",row.names=FALSE)
	#write.table(conv,file="conv.txt",sep="\t",row.names=FALSE)
	
	while(change>conv){
		g.new <- postmean(g.hat,g.bar,n,d.old,t2)
		sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum,na.rm=T)
		d.new <- postvar(sum2,n,a,b)
		change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old,na.rm=TRUE)
		g.old <- g.new
		d.old <- d.new
		count <- count+1
		#if(is.na(change)==TRUE){
		#		change=conv
		#}
		
		}
	#cat("This batch took", count, "iterations until convergence\n")
	adjust <- rbind(g.new, d.new)
	rownames(adjust) <- c("g.star","d.star")
	adjust
	}

#likelihood function used below
L <- function(x,g.hat,d.hat){prod(dnorm(x,g.hat,sqrt(d.hat)))}

# Monte Carlo integration function to find the nonparametric adjustments
int.eprior <- function(sdat,g.hat,d.hat){
	g.star <- d.star <- NULL
	r <- nrow(sdat)
	for(i in 1:r){
		g <- g.hat[-i]
		d <- d.hat[-i]		
		x <- sdat[i,!is.na(sdat[i,])]
		n <- length(x)
		j <- numeric(n)+1
		dat <- matrix(as.numeric(x),length(g),n,byrow=T)
		resid2 <- (dat-g)^2
		sum2 <- resid2%*%j
		LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
		LH[LH=="NaN"]=0
		g.star <- c(g.star,sum(g*LH)/sum(LH))
		d.star <- c(d.star,sum(d*LH)/sum(LH))
		#if(i%%1000==0){cat(i,'\n')}
		}
	adjust <- rbind(g.star,d.star)
	rownames(adjust) <- c("g.star","d.star")
	adjust	
	} 

#fits the L/S model in the presence of missing data values

Beta.NA = function(y,X){
	des=X[!is.na(y),]
	#print(dim(des))
	sum_des<-apply(des,2,sum)
	
	#write.table(des,file="current_des1.txt",sep="\t",row.names=TRUE)
	
	bad_batches<-which(sum_des==0)
	if(length(bad_batches)>0){
	des<-des[,-bad_batches]
	}
	#write.table(des,file="current_des2.txt",sep="\t",row.names=TRUE)
	y1=y[!is.na(y)]
	B <- solve(t(des)%*%des)%*%t(des)%*%y1
	#print(dim(B))
	
	Ball<-matrix(0,nrow=length(sum_des),ncol=1)
	if(length(bad_batches)>0){
	Ball[-bad_batches,]<-B
	}else{
	Ball<-B
	}
	#print(Ball[1:10,])
	Ball
	}

Beta.NAv36 = function(y,X){
	des=X[!is.na(y),]
	#print(dim(des))
	
	#write.table(des,file="current_des1.txt",sep="\t",row.names=TRUE)
	#write.table(X,file="X.txt",sep="\t",row.names=TRUE)
	#write.table(y,file="y.txt",sep="\t",row.names=TRUE)
	
	#print(head(des))
	if(length(des)>0){
	des<-as.matrix(des)
	
	if(ncol(des)==1){
		des<-t(des)
		sum_des<-des
	}else{
			if(nrow(des)>1){
				sum_des<-apply(des,2,sum)
			}
	}
	
	}else{
		sum_des<-des
	}
	bad_batches<-which(sum_des==0)
	if(length(bad_batches)>0){
	des<-des[,-bad_batches]
	}
	#write.table(des,file="current_des2.txt",sep="\t",row.names=TRUE)
	y1=y[!is.na(y)]
	B <- solve(t(des)%*%des)%*%t(des)%*%y1
	#print(dim(B))
	
	Ball<-matrix(0,nrow=length(sum_des),ncol=1)
	if(length(bad_batches)>0){
	Ball[-bad_batches,]<-B
	}else{
	Ball<-B
	}
	#print(Ball[1:10,])
	Ball
	}
	
tic.eval<-function(dataA,outloc){
	
	dir.create(outloc,showWarnings=FALSE)
	setwd(outloc)
	
	tic<-apply(dataA,2,function(x){
		x<-replace(x,which(x==0),NA)
		return(sum(x,na.rm=TRUE))
	})
	
	
	mean_tic<-mean(tic)

	cv_tic<-100*sd(tic,na.rm=TRUE)/mean(tic,na.rm=TRUE)
	tic_res<-cbind(mean_tic,cv_tic)
	colnames(tic_res)<-c("Average_TIC","CV_TIC")



	main_lab<-paste("Total TIC using all features\n Average TIC=",round(mean_tic,2),"\n%CV TIC=",round(cv_tic,2),sep="")
	
	
	#tiff("barplot_TIC_using_all_features.tiff",width=2000,height=2000,res=300)
	#pdf("TIC_all_features.pdf")
    par(mfrow=c(1,1))
	barplot(tic,cex.names=0.4,cex.axis=1,main=main_lab,col="brown",cex.main=0.6)
		
	#boxplot(dataA,cex.names=0.35,cex.axis=1,main=main_lab)
	#dev.off()
	
	write.table(tic_res,"TIC_using_all_features.txt",sep="\t",quote=F,col.name=T,row.names=F)
	names(tic)<-c("sample_TIC")
	write.table(tic,"TIC_each_sample_using_all_features.txt",sep="\t",quote=F,col.name=T,row.names=T)
	
	#tiff("boxplot_sampleintensity_using_all_features.tiff",width=2000,height=2000,res=300)
	
}


msc.calib<-function (X, reference = NULL) 
{
  
    Z <- cbind(1, reference)
    B <- t(solve(crossprod(Z), t(X %*% Z)))
    #res <- (X - B[, 1])/B[, 2]
    #attr(res, "reference") <- reference
    #class(res) <- c("msc", "matrix")
    return(B)
}


eval.target.mz<-function(dataA,refMZ,feature.eval.result,mzthresh=10,timethresh=NA,outloc,folderheader=NA,alignment.tool=NA,mztype="raw",xMSanalyzer.outloc=NA){

feature09<-dataA
if(is.na(refMZ[1,1])==FALSE){
stddata<-refMZ
}else{
	Name<-paste("mz",seq(1,dim(feature09)[1]),sep="")
	
	if(is.na(timethresh)==TRUE){
	stddata<-cbind(feature09[,c(1)],Name)
	}else{
	stddata<-cbind(feature09[,c(1:2)],Name)
		
	}
}
qcresults<-feature.eval.result

rm(dataA)
outloc1<-paste(outloc,"/",folderheader,"targetedeval",mzthresh,"ppm/",sep="")
#dir.create(outloc1,showWarnings=FALSE)
#setwd(outloc1)

col.names.dataA<-colnames(feature09)

if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
		
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
              	
                    sample.col.start=9
                    col.names.dataA[1]="mz"
                    col.names.dataA[4]="time"
		    
                    colnames(data_a)=col.names.dataA
                   
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		   
		     col.names.dataA[2]="time"
		      sample.col.start=3
		     		      print("Using the 1st column as \"mz\" and 2nd columns as \"retention time\"")
		     
		    colnames(feature09)=col.names.dataA
                   
		}
		
par(mfrow=c(2,2))


overlapres5ppm<-getVenn(dataA=feature09, name_a=paste("Expdata",sep=""), dataB=stddata, name_b="Target", mz.thresh = mzthresh, time.thresh=timethresh,alignment.tool=alignment.tool,
xMSanalyzer.outloc=getwd(),plotvenn=FALSE,nearest.time.match=TRUE)

save(overlapres5ppm,qcresults,file="overlapRes.Rda")
match5ppmdata<-feature09[overlapres5ppm$common$index.A,]

name_mz<-{}
if(dim(match5ppmdata)[1]<2){
print("No matches found for targeted metabolites.")
return(name_mz);
}


qcresults5ppmdata<-qcresults[overlapres5ppm$common$index.A,]


fnames<-paste(xMSanalyzer.outloc,"/Targeted_feature_table_",mzthresh,"ppm_filtered.txt",sep="")

fnames<-paste("Boxplot_sampleintensity_usingtargetmatches",mzthresh,"ppm.tiff",sep="")


int_data<-log10(match5ppmdata[,-c(1:(sample.col.start-1))]+1)
#int_data<-match5ppmdata[,-c(1:(sample.col.start-1))]
main_lab<-paste("Intensity distribution (log10; all samples) \n of each m/z matching targets (+/- ",mzthresh,"ppm)",sep="")

mainlabtext=paste("Overlap with reference target list using ",mztype, " m/z",sep="")
vennDiagram(circle.col="red",overlapres5ppm$vennCounts,counts.col="blue",main=mainlabtext,cex.main=0.8)


cv5ppm<-apply(match5ppmdata[,-c(1:(sample.col.start-1))],1,function(x){
	x<-replace(x,which(x==0),NA)
	cvres<-100*sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)
	return(cvres)
	})

names(cv5ppm)<-round(match5ppmdata[,1],5)

fnames<-paste("TotalCVallsamples_refmz",mzthresh,"ppm.tiff",sep="")

#pdf("Targeted_mz_QC.pdf")
main_name<-paste("Total CV (across all samples) \n of each m/z matching targets (+/- ",mzthresh,"ppm)",sep="")

if(length(cv5ppm)>0){

try(barplot(cv5ppm,cex.names=0.4,cex.axis=1,main=main_name,col="brown",cex.main=0.7),silent=TRUE)
}

fnames<-paste("BarplotTIC_targetmatches",mzthresh,"ppm.tiff",sep="")

main_name<-paste("TIC per sample for matching targets (+/- ",mzthresh,"ppm)",sep="")

tic<-apply(match5ppmdata[,-c(1:(sample.col.start-1))],2,function(x){
	x<-replace(x,which(x==0),NA)
	res<-sum(x,na.rm=TRUE)
	return(res)
})

if(length(tic)>0){
try(barplot(tic,cex.names=0.4,cex.axis=1,main=main_name,col="brown",cex.main=0.7),silent=TRUE)
}

fnames<-paste("Pairwiseplot_overall",mzthresh,"ppm.tiff",sep="")
#tiff(fnames,width=2000,height=2000,res=300)
if(length(overlapres5ppm$common$index.A)>0){
try(plot(cbind(feature09[overlapres5ppm$common$index.A,c(1:2)],qcresults[overlapres5ppm$common$index.A,c(3,7,8)]),main="Pairwise plots of m/z, time, CV, Qscore",cex.main=0.7),silent=TRUE)
}
#dev.off()



#dev.off()


mean_tic<-mean(tic,na.rm=TRUE)


tic<-as.data.frame(tic)
names(tic)<-c("sample_TIC")

fnames<-paste("TICpersamp_targetmatches",mzthresh,"ppm.txt",sep="")
#write.table(tic,file=fnames,sep="\t",quote=F,col.name=T,row.names=TRUE)


#refMZ<-cbind(refMZ[overlapres5ppm$common$index.B,],match5ppmdata)

#colnames(refMZ)<-c("mz","Name")
if(is.na(alignment.tool)==FALSE){
if(alignment.tool=="apLCMS"){
name_mz<-cbind(feature09[overlapres5ppm$common$index.B,],stddata[overlapres5ppm$common$index.B,],feature09[overlapres5ppm$common$index.A,-c(1:4)])
}else{

	if(alignment.tool=="XCMS"){
		name_mz<-cbind(feature09[overlapres5ppm$common$index.A,],stddata[overlapres5ppm$common$index.B,],feature09[overlapres5ppm$common$index.A,-c(1:8)])
	}
}
}else{
		name_mz<-cbind(feature09[overlapres5ppm$common$index.A,c(1:2)],stddata[overlapres5ppm$common$index.B,],feature09[overlapres5ppm$common$index.A,-c(1:2)])
}

name_mz<-as.data.frame(name_mz)

#name_mz<-merge(refMZ,match5ppmdata,by="mz")

if(folderheader=="raw"){
fnames<-paste("../../Stage3b/",folderheader,"_target_featuretable_mz",mzthresh,"ppm_time",timethresh,"_RAW.txt",sep="")
}else{
	
	if(folderheader=="ComBat"){
        #fnames<-paste("../../Stage4b/",folderheader,"_target_featuretable",mzthresh,"ppm_Com.txt",sep="")
        fnames<-paste("../../Stage4b/",folderheader,"_target_mzcalibrated_featuretable_mz",mzthresh,"ppm_time",timethresh,"_ComBat.txt",sep="")

    }
}
fnames<-paste("Targetmz_matching_data",mzthresh,"ppmA.txt",sep="")
#write.table(name_mz,file=fnames,sep="\t",quote=F,col.name=T,row.names=F)




#outloc1<-paste(outloc1,"target_barplots/",sep="")
#dir.create(outloc1,showWarnings=FALSE)
#setwd(outloc1)

match5ppmdata<-name_mz
rm(name_mz)


delta_ppm_vec<-{}

match5ppmdata<-as.data.frame(match5ppmdata)
match5ppmdata_int<-match5ppmdata[,-c(1:(2+dim(refMZ)[2]))]
match5ppmdata_int<-apply(match5ppmdata_int,2,as.numeric)


if(length(overlapres5ppm$common$index.B)>0){
    
    delta_ppm_vec<-get_mz_error(match5ppmdata[,c(1,3)])
    
     delta_time_vec<-get_time_error(match5ppmdata[,c(2,4)])
    
    match5ppmdata<-cbind(match5ppmdata[,1:4],delta_ppm_vec, delta_time_vec,match5ppmdata[,-c(1:4)])
   
    fnames<-paste("Targetmz_matching_data",mzthresh,"ppm.txt",sep="")
   # write.table(match5ppmdata,file=fnames,sep="\t",quote=F,col.name=T,row.names=F)
    
for(i in 1:dim(match5ppmdata_int)[1]){
	fname<-paste("mz",round(match5ppmdata[i,1],5),".tiff",sep="")

	delta_ppm<-delta_ppm_vec[i] #10^6*(abs(match5ppmdata[i,1]-match5ppmdata[i,3]))/match5ppmdata[i,3]
	delta_ppm<-round(delta_ppm,2)


	if(length(overlapres5ppm$common$index.B)>1){
	cur_vec<-t(match5ppmdata_int[i,])
	}else{
	cur_vec<-t(match5ppmdata_int[i])
	}
	delta_ppm<-abs(delta_ppm)

	cv_val<-100*sd(cur_vec+1)/mean(cur_vec+1)
	cv_val<-round(cv_val,2)

	delta_s<- round(delta_time_vec[i],2)
    
    barplot(cur_vec,main=paste(mztype," mz: ",round(match5ppmdata[i,1],5),"; time:",round(match5ppmdata[i,2],0), "(s) ;\n",match5ppmdata$Name[i], " (delta m/z=",delta_ppm," ppm)\n delta time=",delta_s," s; %RSD or CV (across all samples): ",cv_val,sep=""),
    ylab="Intensity",xlab="Sample index",cex.axis=1,cex.main=0.7,cex.names=0.4,col="brown") #,silent=TRUE)
    
}

median_delta_ppm<-round(median(delta_ppm_vec,na.rm=TRUE),2)

median_delta_s<-round(median(delta_time_vec,na.rm=TRUE),2)


}else{
    median_delta_ppm<-mzthresh
    print("No target mz found.")
    
}
#dev.off()

return(list("targetdata"=match5ppmdata,"delta_ppm_error"=median_delta_ppm,"delta_time_error"=median_delta_s))
}





eval.target.mz.calib<-function(dataA,refMZ,feature.eval.result,mzthresh=10,timethresh=NA,outloc,folderheader=NA,alignment.tool=NA,mztype="raw",xMSanalyzer.outloc=NA){

feature09<-dataA
if(is.na(refMZ[1,1])==FALSE){
stddata<-refMZ
}else{
	Name<-paste("mz",seq(1,dim(feature09)[1]),sep="")
	
	if(is.na(timethresh)==TRUE){
	stddata<-cbind(feature09[,c(1)],Name)
	}else{
	stddata<-cbind(feature09[,c(1:2)],Name)
		
	}
}
qcresults<-feature.eval.result

rm(dataA)
outloc1<-paste(outloc,"/",folderheader,"targetedeval",mzthresh,"ppm/",sep="")
#dir.create(outloc1,showWarnings=FALSE)
#setwd(outloc1)

col.names.dataA<-colnames(feature09)

if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
		
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
              	
                    sample.col.start=9
                    col.names.dataA[1]="mz"
                    col.names.dataA[4]="time"
		    
                    colnames(data_a)=col.names.dataA
                   
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		   
		     col.names.dataA[2]="time"
		      sample.col.start=3
		     		      print("Using the 1st column as \"mz\" and 2nd columns as \"retention time\"")
		     
		    colnames(feature09)=col.names.dataA
                   
		}
		
par(mfrow=c(2,2))


overlapres5ppm<-getVenn(dataA=feature09, name_a=paste("Expdata",sep=""), dataB=stddata, name_b="Target", mz.thresh = mzthresh, time.thresh=timethresh,alignment.tool=alignment.tool,
xMSanalyzer.outloc=getwd(),plotvenn=FALSE,nearest.time.match=TRUE)


match5ppmdata<-feature09[overlapres5ppm$common$index.A,]

name_mz<-{}
if(dim(match5ppmdata)[1]<2){
print("No matches found for targeted metabolites.")
return(name_mz);
}


qcresults5ppmdata<-qcresults[overlapres5ppm$common$index.A,]


fnames<-paste(xMSanalyzer.outloc,"/Targeted_feature_table_",mzthresh,"ppm_filtered.txt",sep="")

fnames<-paste("Boxplot_sampleintensity_usingtargetmatches",mzthresh,"ppm.tiff",sep="")


int_data<-log10(match5ppmdata[,-c(1:(sample.col.start-1))]+1)
#int_data<-match5ppmdata[,-c(1:(sample.col.start-1))]
main_lab<-paste("Intensity distribution (log10; all samples) \n of each m/z matching targets (+/- ",mzthresh,"ppm)",sep="")

mainlabtext=paste("Overlap with reference target list for calibration using ",mztype, " m/z",sep="")
vennDiagram(circle.col="red",overlapres5ppm$vennCounts,counts.col="blue",main=mainlabtext,cex.main=0.8)


cv5ppm<-apply(match5ppmdata[,-c(1:(sample.col.start-1))],1,function(x){
	x<-replace(x,which(x==0),NA)
	cvres<-100*sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)
	return(cvres)
	})

names(cv5ppm)<-round(match5ppmdata[,1],5)

fnames<-paste("TotalCVallsamples_refmz",mzthresh,"ppm.tiff",sep="")

#pdf("Targeted_mz_QC.pdf")
main_name<-paste("Total CV (across all samples) \n of each m/z matching targets (+/- ",mzthresh,"ppm)",sep="")

if(length(cv5ppm)>0){

try(barplot(cv5ppm,cex.names=0.4,cex.axis=1,main=main_name,col="brown",cex.main=0.7),silent=TRUE)
}

fnames<-paste("BarplotTIC_targetmatches",mzthresh,"ppm.tiff",sep="")

main_name<-paste("TIC per sample for matching targets (+/- ",mzthresh,"ppm)",sep="")

tic<-apply(match5ppmdata[,-c(1:(sample.col.start-1))],2,function(x){
	x<-replace(x,which(x==0),NA)
	res<-sum(x,na.rm=TRUE)
	return(res)
})

if(length(tic)>0){
try(barplot(tic,cex.names=0.4,cex.axis=1,main=main_name,col="brown",cex.main=0.7),silent=TRUE)
}




#dev.off()


mean_tic<-mean(tic,na.rm=TRUE)


tic<-as.data.frame(tic)
names(tic)<-c("sample_TIC")

fnames<-paste("TICpersamp_targetmatches",mzthresh,"ppm.txt",sep="")
#write.table(tic,file=fnames,sep="\t",quote=F,col.name=T,row.names=TRUE)


#refMZ<-cbind(refMZ[overlapres5ppm$common$index.B,],match5ppmdata)

#colnames(refMZ)<-c("mz","Name")
if(is.na(alignment.tool)==FALSE){
if(alignment.tool=="apLCMS"){
name_mz<-cbind(feature09[overlapres5ppm$common$index.B,],stddata[overlapres5ppm$common$index.B,],feature09[overlapres5ppm$common$index.A,-c(1:4)])
}else{

	if(alignment.tool=="XCMS"){
		name_mz<-cbind(feature09[overlapres5ppm$common$index.A,],stddata[overlapres5ppm$common$index.B,],feature09[overlapres5ppm$common$index.A,-c(1:8)])
	}
}
}else{
		name_mz<-cbind(feature09[overlapres5ppm$common$index.A,c(1:2)],stddata[overlapres5ppm$common$index.B,],feature09[overlapres5ppm$common$index.A,-c(1:2)])
}

name_mz<-as.data.frame(name_mz)

#name_mz<-merge(refMZ,match5ppmdata,by="mz")

if(folderheader=="raw"){
fnames<-paste("../../Stage3b/",folderheader,"_target_featuretable_mz",mzthresh,"ppm_time",timethresh,"_RAW.txt",sep="")
}else{
	
	if(folderheader=="ComBat"){
        #fnames<-paste("../../Stage4b/",folderheader,"_target_featuretable",mzthresh,"ppm_Com.txt",sep="")
        fnames<-paste("../../Stage4b/",folderheader,"_target_mzcalibrated_featuretable_mz",mzthresh,"ppm_time",timethresh,"_ComBat.txt",sep="")

    }
}


match5ppmdata<-name_mz
rm(name_mz)


delta_ppm_vec<-{}

match5ppmdata<-as.data.frame(match5ppmdata)
match5ppmdata_int<-match5ppmdata[,-c(1:(2+dim(refMZ)[2]))]
match5ppmdata_int<-apply(match5ppmdata_int,2,as.numeric)


if(length(overlapres5ppm$common$index.B)>0){
    
    delta_ppm_vec<-get_mz_error(match5ppmdata[,c(1,3)])
    
     delta_time_vec<-get_time_error(match5ppmdata[,c(2,4)])
    
    match5ppmdata<-cbind(match5ppmdata[,1:4],delta_ppm_vec, delta_time_vec,match5ppmdata[,-c(1:4)])
   
    fnames<-paste("Targetmz_matching_data",mzthresh,"ppm.txt",sep="")
   # write.table(match5ppmdata,file=fnames,sep="\t",quote=F,col.name=T,row.names=F)
    


median_delta_ppm<-round(median(delta_ppm_vec,na.rm=TRUE),2)

median_delta_s<-round(median(delta_time_vec,na.rm=TRUE),2)


}else{
    median_delta_ppm<-mzthresh
    print("No target mz found.")
    
}
#dev.off()

return(list("targetdata"=match5ppmdata,"delta_ppm_error"=median_delta_ppm,"delta_time_error"=median_delta_s))
}



eval.target.calibratedmz<-function(dataA,refMZ,feature.eval.result,mzthresh=10,timethresh=NA,outloc,folderheader=NA,alignment.tool=NA,xMSanalyzer.outloc=NA){

feature09<-dataA
if(is.na(refMZ[1,1])==FALSE){
stddata<-refMZ
}else{
	Name<-paste("mz",seq(1,dim(feature09)[1]),sep="")
	stddata<-cbind(feature09[,c(1)],Name)
	
	}
qcresults<-feature.eval.result

rm(dataA)
outloc1<-paste(outloc,"/",folderheader,"calibratedtargetedeval",mzthresh,"ppm/",sep="")
#dir.create(outloc1,showWarnings=FALSE)
#setwd(outloc1)

col.names.dataA<-colnames(feature09)

if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
		
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
              	
                    sample.col.start=9
                    col.names.dataA[1]="mz"
                    col.names.dataA[4]="time"
		    
                    colnames(data_a)=col.names.dataA
                   
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		   
		     col.names.dataA[2]="time"
		      sample.col.start=3
		     		      print("Using the 1st column as \"mz\" and 2nd column as \"retention time\"")
		     
		    colnames(feature09)=col.names.dataA
                   
		}
		
par(mfrow=c(2,2))


overlapres5ppm<-getVenn(dataA=feature09, name_a=paste("Expdata",sep=""), dataB=stddata, name_b="Target", mz.thresh = mzthresh, time.thresh=timethresh,alignment.tool=alignment.tool,
xMSanalyzer.outloc=getwd(),plotvenn=FALSE,nearest.time.match=TRUE)

#save(overlapres5ppm,file="overlapRes.Rda")
match5ppmdata<-feature09[overlapres5ppm$common$index.A,]

name_mz<-{}
if(dim(match5ppmdata)[1]<2){
print("No matches found for targeted metabolites.")
return(name_mz);
}


qcresults5ppmdata<-qcresults[overlapres5ppm$common$index.A,]


#m1<-cbind(match5ppmdata[,c(1:(sample.col.start-1))],qcresults5ppmdata[,c(3,7,8)],match5ppmdata[,-c(1:(sample.col.start-1))])

fnames<-paste("../../Stage3b/Targeted_feature_table_",mzthresh,"ppm_filtered.txt",sep="")

vennDiagram(circle.col="red",overlapres5ppm$vennCounts,counts.col="blue",main="Overlap with reference target list using calibrated mz",cex.main=0.8)


if(is.na(alignment.tool)==FALSE){
if(alignment.tool=="apLCMS"){
name_mz<-cbind(feature09[overlapres5ppm$common$index.B,],stddata[overlapres5ppm$common$index.B,],feature09[overlapres5ppm$common$index.A,-c(1:4)])
}else{

	if(alignment.tool=="XCMS"){
		name_mz<-cbind(feature09[overlapres5ppm$common$index.A,],stddata[overlapres5ppm$common$index.B,],feature09[overlapres5ppm$common$index.A,-c(1:8)])
	}
}
}else{
		name_mz<-cbind(feature09[overlapres5ppm$common$index.A,c(1:2)],stddata[overlapres5ppm$common$index.B,],feature09[overlapres5ppm$common$index.A,-c(1:2)])
}

name_mz<-as.data.frame(name_mz)

#name_mz<-merge(refMZ,match5ppmdata,by="mz")

if(folderheader=="RAW"){
    #fnames<-paste("../../Stage3b/",folderheader,"_mzcalibratedtargetfeaturetable",mzthresh,"ppm.txt",sep="")
    fnames<-paste(xMSanalyzer.outloc,"/",folderheader,"_mzcalibrated_targeted_featuretable_mz",mzthresh,"ppm_time",timethresh,".txt",sep="")
    #fnames<-paste("../../Stage3b/",folderheader,"_mzcalibrated_targeted_featuretable_mz",mzthresh,"ppm_time",timethresh,".txt",sep="")


}else{
	
	if(folderheader=="ComBat"){
        #fnames<-paste("../../Stage4b/",folderheader,"_mzcalibratedtargetfeaturetable",mzthresh,"ppm.txt",sep="")
        #fnames<-paste("../../Stage4b/",folderheader,"_mzcalibrated_targeted_featuretable_mz",mzthresh,"ppm_time",timethresh,".txt",sep="")
        fnames<-paste(xMSanalyzer.outloc,"/",folderheader,"_mzcalibrated_targeted_featuretable_mz",mzthresh,"ppm_time",timethresh,".txt",sep="")
    }else{
        fnames<-paste(xMSanalyzer.outloc,"/",folderheader,"_mzcalibrated_targeted_featuretable_mz",mzthresh,"ppm_time",timethresh,".txt",sep="")
        
    }
}
#fnames<-paste("Targetmz_matching_data",mzthresh,"ppmA.txt",sep="")
#write.table(name_mz,file=fnames,sep="\t",quote=F,col.name=T,row.names=F)




#outloc1<-paste(outloc1,"target_barplots/",sep="")
#dir.create(outloc1,showWarnings=FALSE)
#setwd(outloc1)

match5ppmdata<-name_mz
rm(name_mz)


delta_ppm_vec<-{}

match5ppmdata<-as.data.frame(match5ppmdata)
match5ppmdata_int<-match5ppmdata[,-c(1:(2+dim(refMZ)[2]))]
match5ppmdata_int<-apply(match5ppmdata_int,2,as.numeric)


if(length(overlapres5ppm$common$index.B)>0){
    
    delta_ppm_vec<-get_mz_error(match5ppmdata[,c(1,3)])
    
     delta_time_vec<-get_time_error(match5ppmdata[,c(2,4)])
    
    match5ppmdata<-cbind(match5ppmdata[,1:4],delta_ppm_vec,match5ppmdata[,-c(1:4)])
   
   #  fnames<-paste("Target_calibratedmz_matching_data",mzthresh,"ppm.txt",sep="")
   
   if(folderheader=="raw"){
       
       folderheader="RAW"
   }
    
    
    #write.table(match5ppmdata,file=fnames,sep="\t",quote=F,col.name=T,row.names=F)
    
for(i in 1:dim(match5ppmdata_int)[1]){
	fname<-paste("mz",round(match5ppmdata[i,1],5),".tiff",sep="")

	delta_ppm<-delta_ppm_vec[i] #10^6*(abs(match5ppmdata[i,1]-match5ppmdata[i,3]))/match5ppmdata[i,3]
	delta_ppm<-round(delta_ppm,2)
	
	delta_s<-delta_time_vec[i]
	delta_s<-round(delta_s,2)
	

if(length(overlapres5ppm$common$index.B)>1){
cur_vec<-t(match5ppmdata_int[i,])
}else{
cur_vec<-t(match5ppmdata_int[i])
}
#delta_ppm<-abs(delta_ppm)
#delta_s<-abs(delta_s)

cv_val<-100*sd(cur_vec+1)/mean(cur_vec+1)
cv_val<-round(cv_val,2)

#try(
barplot(cur_vec,main=paste("Calibrated mz: ",round(match5ppmdata[i,1],5),"; time:",round(match5ppmdata[i,2],0), "(s) ;\n",match5ppmdata$Name[i], " (delta m/z=",delta_ppm," ppm)\n (delta time=",delta_s," s; %RSD or CV (across all samples): ",cv_val,sep=""),
ylab="Intensity",xlab="Sample index",cex.axis=1,cex.main=0.7,cex.names=0.4,col="brown") #,silent=TRUE)
	#dev.off()
    
    
}


median_delta_ppm<-round(median(delta_ppm_vec,na.rm=TRUE),2)
median_delta_s<-round(median(delta_time_vec,na.rm=TRUE),2)


}else{
    median_delta_ppm<-mzthresh
    print("No target mz found.")
    
}
#dev.off()

return(list("targetdata"=match5ppmdata,"delta_ppm_error"=median_delta_ppm,"delta_time_error"=median_delta_s))
}


runlmreg<-function(X,Y,fdrmethod="BH",fdrthresh=0.05,pvalue.thresh=0.05){
    
    data_m_fc<-X #[,-c(1:2)]
    classlabels_response_mat<-Y #[,-c(1)]
    data_m_fc_withfeats<-X
    logistic_reg<-FALSE
    rm(X)
    fileheader="lmreg"
    
    #save(data_m_fc,file="data_m_fc.Rda")
    
   
    res1<-apply(data_m_fc,1,function(x){
        
        xvec<-x
        
        if(dim(classlabels_response_mat)[2]>1){
            
            for(cnum in 2:dim(classlabels_response_mat)[2]){
                
                classlabels_response_mat[,cnum]<-as.numeric(classlabels_response_mat[,cnum])
                
            }
        }else{
            
            classlabels_response_mat[,1]<-as.numeric(classlabels_response_mat[,1])
        }
        
        data_mat_anova<-cbind(xvec,classlabels_response_mat)
        
        cnames<-colnames(data_mat_anova)
        cnames[1]<-"Response"
        
        colnames(data_mat_anova)<-cnames
        
        
        
        anova_res<-diffexplmreg(dataA=data_mat_anova,logistic_reg)
        
        return(anova_res)
    })
    
    #save(res1,file="lmregres.Rda")
    main_pval_mat<-{}
    
    posthoc_pval_mat<-{}
    pvalues<-{}
    
    all_inf_mat<-{}
    
    for(i in 1:length(res1)){
        
        main_pval_mat<-rbind(main_pval_mat,res1[[i]]$mainpvalues)
        pvalues<-c(pvalues,res1[[i]]$mainpvalues[1])
        
        cur_pvals<-t(res1[[i]]$mainpvalues)
        cur_est<-t(res1[[i]]$estimates)
        cur_stderr<-t(res1[[i]]$stderr)
        cur_tstat<-t(res1[[i]]$statistic)
        cur_res<-cbind(cur_pvals,cur_est,cur_stderr,cur_tstat)
        
        all_inf_mat<-rbind(all_inf_mat,cur_res)
        
        
        
    }
    if(fdrmethod=="BH"){
        fdr_adjust_pvalue<-p.adjust(pvalues,method="BH")
    }else{
        if(fdrmethod=="ST"){
            fdr_adjust_pvalue<-qvalue(pvalues)
            fdr_adjust_pvalue<-fdr_adjust_pvalue$qvalues
        }else{
            if(fdrmethod=="Strimmer"){
                pdf("fdrtool.pdf")
                
                fdr_adjust_pvalue<-fdrtool(as.vector(pvalues),statistic="pvalue",verbose=FALSE)
                fdr_adjust_pvalue<-fdr_adjust_pvalue$qval
                dev.off()
            }else{
                if(fdrmethod=="none"){
                    #fdr_adjust_pvalue<-pvalues
                    fdr_adjust_pvalue<-p.adjust(pvalues,method="none")
                }else{
                    if(fdrmethod=="BY"){
                        fdr_adjust_pvalue<-p.adjust(pvalues,method="BY")
                    }else{
                        if(fdrmethod=="bonferroni"){
                            fdr_adjust_pvalue<-p.adjust(pvalues,method="bonferroni")
                            
                        }
                    }
																}
            }
        }
        
        
        
    }
    
    
    if(fdrmethod=="none"){
        filename<-paste(fileheader,"_pvalall_withfeats.txt",sep="")
        
    }else{
        filename<-paste(fileheader,"_fdrall_withfeats.txt",sep="")
    }
    cnames_tab<-colnames(data_m_fc_withfeats)
    cnames_tab<-c("P.value","adjusted.P.value",cnames_tab)
    
    pvalues<-as.data.frame(pvalues)
    
    final.pvalues<-pvalues
    sel.diffdrthresh<-fdr_adjust_pvalue<fdrthresh & final.pvalues<pvalue.thresh
    
    data_limma_fdrall_withfeats<-cbind(pvalues,fdr_adjust_pvalue,data_m_fc_withfeats)
    
    #colnames(data_limma_fdrall_withfeats)<-as.character(cnames_tab)
    
    
    filename<-paste(fileheader,"_pval_coef_stderr.txt",sep="")
    
    data_allinf_withfeats<-cbind(all_inf_mat,data_m_fc_withfeats)
    

    cnames_tab<-colnames(data_m_fc_withfeats)
    
    class_column_names<-colnames(classlabels_response_mat)
    
    
    
    cnames_tab<-c(paste("P.value_",class_column_names,sep=""),
    paste("Estimate_",class_column_names,sep=""), paste("StdError_",class_column_names,sep=""),
    paste("statistic_",class_column_names,sep=""),cnames_tab)
    
    colnames(data_allinf_withfeats)<-as.character(cnames_tab)
    data_allinf_withfeats<-data_allinf_withfeats[,c(1:4)]
 
    return(data_allinf_withfeats)
    
}

get_mz_error<-function(dataA){
    
    dataA<-as.data.frame(dataA)
    
    delta_vec<-apply(dataA,1,function(x){
        x[1]<-as.numeric(as.character(x[1]));
        x[2]<-as.numeric(as.character(x[2]));
        
        ppmerror=10^6*(x[2]-x[1])/(x[2]);
        return(ppmerror);
    })
    
    return(delta_vec)
}

get_time_error<-function(dataA){
    
    dataA<-as.data.frame(dataA)
    
    delta_vec<-apply(dataA,1,function(x){
        x[1]<-as.numeric(as.character(x[1]));
        x[2]<-as.numeric(as.character(x[2]));
        
       timeerror=(x[2]-x[1])
        return(timeerror);
    })
    
    return(delta_vec)
}


get_calibrated_mz_data<-function(dataA,delta_ppm_error=NA,delta_time_error=NA,refMZ=NA,mzthresh=10,timethresh=NA,outloc=NA,feature.eval.result=NA,calibration.method=c("median.adjustment","multiplicative.signal.correction")){
    
    dataA<-as.data.frame(dataA)
    
    original_mz<-dataA$mz
    
   check_name<-grep(x=colnames(refMZ),pattern="Calibration")
   
   if(length(check_name)>0){
    
	refMZ=refMZ[which(refMZ$Calibration==1),]
    }
    
if(calibration.method[1]=="median.adjustment"){

 print("Calibrating m/z")
	if(is.na(delta_ppm_error)==TRUE){
	targeted_feat_raw<-eval.target.mz.calib(dataA=dataA,refMZ=refMZ,feature.eval.result=feature.eval.result,mzthresh=mzthresh,
	timethresh=timethresh,outloc=outloc,folderheader="",xMSanalyzer.outloc=outloc)
            
            
	delta_ppm_error<-median(targeted_feat_raw$delta_ppm_error,na.rm=TRUE)
	
	delta_time_error<-median(targeted_feat_raw$delta_time_error,na.rm=TRUE)
	}
    
    
    adjusted_mz_vec<-lapply(1:length(original_mz),function(j){
        original_mz[j]<-as.numeric(as.character(original_mz[j]));
       
       
        correction_factor<-(delta_ppm_error*original_mz[j])/(10^6)
        
        #adjust mz
        adjusted_mz<-original_mz[j]+correction_factor
        
        return(adjusted_mz);
    })
	adjusted_mz_vec<-unlist(adjusted_mz_vec)
    
    dataA$mz<-adjusted_mz_vec
    
    if(is.na(timethresh)==FALSE){
    
    print("Calibrating time")
     original_time<-dataA$time
	
if(FALSE){	
		 adjusted_time_vec<-lapply(1:length(original_time),function(j){
original_time[j]<-as.numeric(as.character(original_time[j]));


correction_factor<-(delta_time_error)

#adjust time
adjusted_time<-original_time[j]+correction_factor

return(adjusted_time);
})
    }
    
     correction_factor<-(delta_time_error)
        
	 original_time<-as.numeric(as.character(original_time));
	 
        #adjust time
        adjusted_time_vec<-original_time+correction_factor
        
	if(min(adjusted_time_vec,na.rm=TRUE)<0){
	
		adjusted_time_vec[which(adjusted_time_vec<0)]<-original_time[which(adjusted_time_vec<0)]
	}
	
    
    dataA$time<-adjusted_time_vec
    

	}
}else{

		msc.mz.B<-msc.calib(dataA$mz,reference=refMZ$mz)
		
		
		
		adjusted_mz_vec<-(dataA$mz-msc.mz.B[,1])/msc.mz.B[,2]
		
		 if(is.na(timethresh)==FALSE){
			msc.time.B<-msc.calib(dataA$time,reference=refMZ$time)
			adjusted_time_vec<-(dataA$time-msc.time.B[,1])/msc.time.B[,2]
			
			
				if(min(adjusted_time_vec,na.rm=TRUE)<0){
	
						adjusted_time_vec[which(adjusted_time_vec<0)]<-original_time[which(adjusted_time_vec<0)]
				}
		}
		
		  dataA$mz<-adjusted_mz_vec
		 dataA$time<-adjusted_time_vec
		 
}
	
	
    
    return(dataA)
}



#Nointeraction
diffexplmreg<-function(dataA,logistic_reg=FALSE){
    
    dataA<-as.data.frame(dataA)
    
    # save(dataA,file="lmreg_func.Rda")
    
    
    if(logistic_reg==TRUE){
        cnames1<-colnames(dataA)
        cnames1[2]<-"Class"
        colnames(dataA)<-cnames1
        
        labels_1<-levels(as.factor(dataA$Class))
        
        dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[1]),0)
        dataA$Class<-replace(dataA$Class,which(dataA$Class==labels_1[2]),1)
        
        
        
        a1 <- glm(dataA$Class ~ .,family=binomial(logit),data=dataA)
    }else{
        a1 <- lm(dataA$Response ~ .,data=dataA) # aov(dataA$Response ~ .,data=dataA) # + chocolate$Factor1*chocolate$Factor2)
    }
    s1<-summary(a1)
    
    #save(s1,file="s1.Rda")
    
    
    if(FALSE){
        anova_res<-anova(a1)
        num_rows<-dim(anova_res)
        pvalues_factors<-data.frame(t(anova_res["Pr(>F)"][-c(num_rows),]))
    }
    
    if(logistic_reg==FALSE){
        
        r2<-s1$adj.r.squared
    }else{
        r2<-NA
    }
    s1<-s1$coefficients
    
    
    s1<-s1[-c(1),]
    
 
    if(dim(dataA)[2]<3){ # && dim(dataA)[1]<3){
        #s1<-as.data.frame(s1)
        s1<-t(s1)
        
    }
    
    
    confint_lower<-s1[,1]-(1.96*s1[,2])
    confint_upper<-s1[,1]+(1.96*s1[,2])
    
    
    return(list("mainpvalues"=s1[,4],"estimates"=s1[,1],"statistic"=s1[,3],"stderr"=s1[,2],"r2"=r2,"confint"=c(confint_lower,confint_upper)))
    
    
}

pca.eval<-function(X,samplelabels,filename=NA,ncomp=5,center=TRUE,scale=TRUE,legendlocation="topright",legendcex=0.5,outloc=getwd()){
	
	
	X<-as.matrix(X)
	print("Performing PCA")
	batchlabels<-samplelabels
	

    if(length(samplelabels)<10){

		ncomp=2
	}
    metabpcaresultlog2allmetabs5pcs<-mixOmics::pca(X,ncomp=ncomp,center=center,scale=scale)
    
    #metabpcaresultnotransform10pcsallmetabs<-metabpcaresult
    result<-metabpcaresultlog2allmetabs5pcs
    
    r1<-100*result$explained_variance #sdev/sum(result$sdev)
    r1<-round(r1,2)
    
    t1<-table(samplelabels)
    
    if(length(t1)<30){
    col_vec<-c("mediumpurple4","mediumpurple1","blueviolet","darkblue","blue","cornflowerblue","cyan4","skyblue",
    "darkgreen", "seagreen1", "green","yellow","orange","pink", "coral1", "palevioletred2",
    "red","saddlebrown","brown","brown3","white","darkgray","aliceblue",
    "aquamarine","aquamarine3","bisque","burlywood1","lavender","khaki3","black")
    }else{
        
        col_vec<-topo.colors(length(t1))
    }
    
    #col_vec<-col_vec[sample(1:length(col_vec),length(col_vec))]
    
    
    l1<-levels(samplelabels)
    col_all=topo.colors(256)
    
    dir.create(outloc,showWarnings=FALSE)
    setwd(outloc)
    print(paste("Generating PCA plots to evaluate batch-effect in ",filename," data",sep=""))
    
    fname<-paste("PCA_batcheffecteval",filename,".tiff",sep="")
    ## 1) raw data
    #tiff(fname,res=300, width=2000,height=2000)
    col <- rep(col_vec[1:length(t1)], t1)
    
    #col<-rep(col_all[1:length(l1)],t1)
    ## Choose different size of points
    cex <- rep(1.5, dim(X)[1])
    
    ## Choose the form of the points (square, circle, triangle and diamond-shaped
    #pch <- replace(pch,which(pch==2),21)
    #pch <- replace(pch,which(pch==3),17)
    
    pch <- rep(15,dim(X)[1])
    
    
    ## comp is the two dimensions you want to display
    ## ind.names is whether you want labels on each points
    ## rep.space determines the subspace to project the individuals ("X-variate",
    ## "Y-variate" or "XY-variate")
    #plotIndiv(result) #, comp = c(1,2), ind.names = FALSE, rep.space = "X-variate", col = col, cex = cex, pch = pch, X.label="PC1",Y.label="PC2")
   
    save(result,file="pcares.Rda")
   save(col,file="col.Rda")
   save(samplelabels,file="samplelabels.Rda")
    
    #classlabels_mat<-cbind(result$names$sample,samplelabels)
    classlabels_mat<-as.data.frame(samplelabels)
    
    colnames(classlabels_mat)<-c("Batch")
    
    if(length(l1)>1){
    reg_res<-runlmreg(X=t(result$x),Y=classlabels_mat,fdrmethod="BH",fdrthresh=0.05)
    }else{
    
	reg_res<-matrix(NA,nrow=ncol(result$x),ncol=1)
    }
    
    
    #if(FALSE)
    {
        
        main_lab<-paste("Batch-effect evaluation; \n Batch ~ PC1; p=",round(reg_res[1,1],2),"\n Batch ~ PC2; p=",round(reg_res[2,1],2),sep="")
       
       #print(main_lab)
       #save(main_lab,file="main_lab.Rda")
       #plotIndiv(result,col.per.group=unique(col),group=samplelabels,legend=TRUE,ind.names=FALSE,cex.legend=1,pch=c(21),size.legend=rel(0.8),title=main_lab)
       
         pca_res<-suppressWarnings(try(plotIndiv(result,col.per.group=unique(col),group=samplelabels,legend=TRUE,ind.names=FALSE,cex.legend=1,pch=c(21),size.legend=rel(0.8),title=main_lab,size.title = rel(1)),silent=TRUE))
        
        #try(print(pca_res),silent=TRUE)
   
   main_lab<-paste("Batch-effect evaluation; \n Batch ~ PC2; p=",round(reg_res[2,1],2),"\n Batch ~ PC3; p=",round(reg_res[3,1],2),sep="")
   
   
   pca_res<-suppressWarnings(try(plotIndiv(result,comp=c(2,3),col.per.group=unique(col),group=samplelabels,legend=TRUE,ind.names=FALSE,cex.legend=1,pch=c(21),size.legend=rel(0.8),title=main_lab,size.title = rel(1)),silent=TRUE))
   
   # try(print(pca_res),silent=TRUE)
   
      main_lab<-paste("Batch-effect evaluation; \n Batch ~ PC1; p=",round(reg_res[1,1],2),"\n Batch ~ PC3; p=",round(reg_res[3,1],2),sep="")
   
   
   pca_res<-suppressWarnings(try(plotIndiv(result,comp=c(1,3),col.per.group=unique(col),group=samplelabels,legend=TRUE,ind.names=FALSE,cex.legend=1,pch=c(21),size.legend=rel(0.8),title=main_lab,size.title = rel(1)),silent=TRUE))
   
   #print(pca_res)
   
	if(is(pca_res, "try-error")){
		
        #print("PCA plot could not be generated.")
        
        #try(plotIndiv(result,col=col,ind.names=FALSE),silent=TRUE)
        
        pca_res<-suppressWarnings(try(plotIndiv(result,col.per.group=unique(col),group=samplelabels,add.legend=TRUE,ind.names=FALSE,pch=c(21),cex=0.8,main=main_lab),silent=TRUE))
        
        #print(pca_res)
	}else{
        

    #print("PCA plot could not be generated.")

    }
}
return(result)

}

merge.Results.child.cortest<-function(dataA, max.mz.diff=15, max.rt.diff=300, merge.eval.pvalue=0.05,alignment.tool="apLCMS", mult.test.cor=FALSE)
{
    
        diff_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}
        cnames=colnames(dataA)
        dataA<-as.data.frame(dataA)
	
		#dataA<-as.data.frame(dataA,2,as.numeric)
         if(alignment.tool=="apLCMS")
        {
              sample.col.start=6
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    cnames[1]="mz"

                    cnames[4]="time"
                    colnames(dataA)=cnames
                    
              }
              else
              {
                    stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
              }
        }

compare_intensities_cortest<-function(other_feats,y){
											cortest_sum=apply(other_feats,1,function(x)
				                                                        {
				                                                        	x<-as.matrix(x)
				                                                        	y<-as.matrix(y)
												yind<-which(y==0)
												xind<-which(x==0)
												naind<-c(yind,xind)
												
												
				                                                                if(dim(x)[1]>dim(y)[1])
												{
													x<-t(x)
												}
												if(dim(y)[1]>dim(x)[1])
												{
													y<-t(y)
												}
				                                                                      if(length(naind)>0){
														 x<-x[-naind]
														 y<-y[-naind]
				                                                                 
													}
												
													#if(max(x)!=0 & max(y)!=0)
													if(length(x)>2 & length(y)>2)
													{
														
														
														cortest_res=try(cor.test(x,y),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																cortest_pval=cortest_res$p.value
																cortest_est=cortest_res$estimate
																if(cortest_est<0){
																	cortest_pval=1
																}
																}
													}else{
													cortest_res=try(cor.test(x,y),silent=TRUE)
																if (is(cortest_res, "try-error")){
																	cortest_pval<-1
																	cortest_est<-0
																}else{
																cortest_pval=cortest_res$p.value
																cortest_est=cortest_res$estimate
																if(cortest_est<0){
																	cortest_pval=1
																}
																}
													}
												
												
				                                                                return(cortest_pval)
				                                                        })
return(cortest_sum)
}


        #Step 1 Group features by m/z
        mz_groups<-lapply(1:dim(dataA)[1],function(j){
                                commat={}
                                diffmz=new("list")
                                ppmb=(max.mz.diff)*(dataA$mz[j]/1000000)
                                getbind_same<-which(abs(dataA$mz-dataA$mz[j])<=ppmb)
                                diffmz[[diff_mz_num]]=getbind_same #dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })
	  
	  del_list<-{}
	  for(m in 1:length(mz_groups)){
	  
	   if((m%in%del_list)==FALSE){
		for(n in (m+1):length(mz_groups)){
		
			if(n>length(mz_groups)){
				break;
			}
			com1<-intersect(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
			if(length(com1)>0){
			
				mz_groups[[m]][[1]]<-c(mz_groups[[m]][[1]],mz_groups[[n]][[1]])
				del_list<-c(del_list,n)
			}
		}
		
		mz_groups[[m]][[1]]<-unique(mz_groups[[m]][[1]])
		
		}
	  }
	  if(length(del_list)>0){
		 mz_groups<-mz_groups[-del_list]
		}
	 
	  
	mz_groups<-unique(mz_groups)
	
	diff_mz_num=1
	mz_groups<-lapply(1:length(mz_groups),function(j){
                                commat={}
                                diffmz=new("list")
                                getbind_same=mz_groups[[j]][[1]]
                                diffmz[[diff_mz_num]]=dataA[getbind_same,]
                                diff_mz_num=diff_mz_num+1
                                return(diffmz)
          })
	  
	if(mult.test.cor==TRUE){
		mz_group_size<-lapply(1:length(mz_groups),function(j){
		num_rows<-dim(mz_groups[[j]][[1]])
		n=num_rows[[1]][[1]]
		num_comp=(n*(n-1))/2
		})

	
		num_comparisons<-sum(unlist(mz_group_size))
	}else{
		num_comparisons=1
	}
	print("num comparisons")
	print(num_comparisons)
	
        #Step 2 Sub-group features from Step 1 by Retention time
        #find the features with RT values within the defined range as compared to the query feature
        diffmat={}
        for(j in 1:length(mz_groups))
        {
                temp_diff={}
                tempdata=mz_groups[[j]][[1]]
		
                if(dim(tempdata)[1]>1)
                {
		
                                tempdata<-tempdata[order(tempdata$median),]
				
				rem_ind<-apply(tempdata,1,function(x){as.numeric(x[(dim(tempdata)[2]-6)])==as.numeric(x[(dim(tempdata)[2]-1)])})
				
				rem_ind<-which(rem_ind==TRUE)
								if(length(rem_ind)>0){
								
									
									 tempdata<-tempdata[-rem_ind,]
								}
				if(dim(tempdata)[1]>0){
				dup_list<-{}
				temp_diff={}
                                for(d in 1:dim(tempdata)[1])
                                {
					rowname<-rownames(tempdata[d,])
					print(rowname)
					print(tempdata[,1:5])
					print(length(tempdata))
					if((rowname%in%dup_list)==FALSE)
					{
                                        #tempdata=mz_groups[[j]][[1]]
                                        cur_feat=tempdata[d,]
                                        other_feats=tempdata
                                        getbind_rtsame<-which(abs(other_feats$time-cur_feat$time)<=max.rt.diff)
                                        commat={}
					 same_feat_ind={}
                                        if(length(getbind_rtsame)>1)
                                        {
                                                other_feats<-other_feats[getbind_rtsame,]
		
					y=cur_feat[sample.col.start:(dim(tempdata)[2]-7)]
		
										 
										ttest_sum<-compare_intensities_cortest(other_feats[,sample.col.start:(dim(tempdata)[2]-7)],y)

                                                        same_feat_ind=which(ttest_sum<(merge.eval.pvalue/num_comparisons))
                                                        same_feat_ind=c(same_feat_ind,which(is.na(ttest_sum)))
							print(ttest_sum)
							dup_list<-c(dup_list,names(same_feat_ind))
									 				
                                                        if(length(same_feat_ind)>0)
                                                        {
                                                                commat=other_feats[same_feat_ind,]
								#commat_qscore<-as.numeric(commat$median)/apply(commat[,sample.col.start:(dim(tempdata)[2]-7)],1,countpeaks)
                                                                #best_level_index=which(as.numeric(commat$median)==min(as.numeric(commat$median)))
								
								 commat_qscore<-as.numeric(commat$median)/as.numeric(commat$numgoodsamples)
								 
								best_level_index=which(as.numeric(commat_qscore)==min(as.numeric(commat_qscore)))
                                                                if(length(best_level_index)>0)
                                                                {
                                                                        best_level_index=best_level_index[1]
                                                                }
                                                                best_level_data=commat[best_level_index,]
                                                        }else
                                                        {
                                                                best_level_data=cur_feat

                                                        }
                          
                                                                        diffmat=rbind(diffmat, best_level_data)
                                                                        temp_diff=rbind(temp_diff,best_level_data)
					
                                        }
                                        else
                                        {
          

                                               					 diffmat<-rbind(diffmat,cur_feat)
                                              					  temp_diff=rbind(temp_diff,tempdata[d,])
                                              					  
                                        		
                                        		
                                        }
					}
                                }
				}
                }
                else
                {
                	cur_feat<-as.matrix(mz_groups[[j]][[1]])
                	
                        									#tempdata=mz_groups[[j]][[1]]
                         if(tempdata$min!=tempdata$max)
			 {
                                               					 diffmat<-rbind(diffmat,cur_feat)
                                              					 temp_diff=rbind(temp_diff,tempdata)
                        	}
                        
                        
                        
                }
                diffmat<-unique(diffmat)
                 best_level_data={}
                

        }
        diffmat=unique(diffmat)
        return(diffmat)
}
 countpeaks<-function(intensity_vec, missing_val=0){
                
						if(is.na(missing_val)==TRUE){
							zeroind<-which(is.na(intensity_vec)==FALSE)
						}else{
							if(missing_val==0){
								zeroind<-which(intensity_vec>0)
							}
							else{
							 stop(paste("Invalid value for \"missing_val\". Please use either \"NA\" or \"0\"", sep=""))
							}
						}
						

						return(length(zeroind))
					}

#################################################################
#Function: feat.batch.annotation.KEGG
#Description: search m/z in Metlin and KEGG
#mz: single m/z value in ppm (eg: 121.945)
#max.mz.diff: +/- m/z match tolerance in ppm for database matching
#           (eg: 5)
##################################################################
feat.batch.annotation.KEGG<-function(dataA,max.mz.diff=10, queryadductlist=c("M+H"), xMSanalyzer.outloc, numnodes=1,syssleep=1)
{
	data_a<-as.data.frame(dataA)
	 
	print("Using the 1st column as \"mz\" for annotation.")
		    
	mzlist<-data_a[,1]
	
	dir.create(xMSanalyzer.outloc,showWarnings=FALSE)
	setwd(xMSanalyzer.outloc)
	
    adductlist=c(1.007276,22.989218,38.963158,-35.012676,-17.0027,0.0227,7.01597,18.033823,
    33.033486,42.033826,44.97116,64.015768,2.014552,23.996494,45.978436,3.021828,
    25.00377,46.985712,-19.01839,-1.007276,18.998371,20.974666,34.969402,36.948606,
    44.998194,59.013851,78.918885,-2.014552,-3.021828)
    alladducts<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M-H2O+NH4", "M+Li","M+NH4",
    "M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
    "M+2H+Na","M+2Na+H","M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
    "M+FA-H","M+CH3COO-H","M+Br","M-2H","M-3H")
    names(adductlist)<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M-H2O+NH4", "M+Li","M+NH4",
    "M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
    "M+2H+Na","M+2Na+H","M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
    "M+FA-H","M+CH3COO-H","M+Br","M-2H","M-3H")
    
    mult_charge<-c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3,3,3,1,1,1,1,1,1,1,1,1,2,3)
    names(mult_charge)<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M-H2O+NH4", "M+Li","M+NH4",
    "M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
    "M+2H+Na","M+2Na+H","M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
    "M+FA-H","M+CH3COO-H","M+Br","M-2H","M-3H")
    
	if(queryadductlist[1]=="positive")
	{
		queryadductlist<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M+Li","M+NH4",
		"M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
		"M+2H+Na","M+2Na+H")
	}else{
		if(queryadductlist[1]=="negative")
		{
			queryadductlist<-c("M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
			"M+FA-H","M+CH3COO","M+Br","M-2H","M-3H")
		}else{
		if(queryadductlist[1]=="all"){
		
			
		queryadductlist<-alladducts
		
		
		}else{
			if(length(which(queryadductlist%in%alladducts==FALSE))>0){
			
				errormsg<-paste("Adduct should be one of:",sep="")
				for(i in alladducts){errormsg<-paste(errormsg, i,sep=" ; ")}
				stop(errormsg, "\n\nUsage: feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"M+H\", \"M+Na\"), xMSanalyzer.outloc, numnodes=1)", 
				"\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"positive\"), xMSanalyzer.outloc, numnodes=1)",
				"\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"negative\"), xMSanalyzer.outloc, numnodes=1)",
				"\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"all\"), xMSanalyzer.outloc, numnodes=1)"
				)
			}
		
		}
	}	
		}
	parentres={}
	#cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "Metlin", "Compound.Name", "CASID", "KEGG.Compound.ID", "KEGG.Pathway.name", "KEGG.Pathway.ID", "HMDB.ID", "PubChem.Substance.ID", "PubChem.Compound.ID", "ChEBI.ID", "LIPID.MAPS.ID")
  cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)", "Metlin", "Compound.Name",  "Chemical.Formula", "Exact Mass", "CASID", "KEGG.Compound.ID",
  "KEGG.Pathway.ID", "KEGG.Pathway.name","HMDB.ID", "PubChem.Substance.ID", "ChEBI.ID", "LIPID.MAPS.ID","KEGG.Brite.Category")
 
	
	for(adnum in 1:length(queryadductlist))
	{
		adductname=queryadductlist[adnum]
		adductmass=adductlist[as.character(adductname)]
		adductcharge=mult_charge[as.character(adductname)]
		
		cl<-parallel::makeCluster(numnodes)
		
		clusterEvalQ(cl, library(XML))
		clusterEvalQ(cl, library(RCurl))
		clusterEvalQ(cl, "feat.batch.annotation.child")
		
		print(paste("Query adduct: ",adductname,sep=""))
		
		mz.annot.res<-new("list")
		min_mz<-min(mzlist)
		max_mz<-max(mzlist)

		mz_group<-ceiling(max_mz/min_mz)
		
		#mz_group<-round(max_mz/10)
		#length(mzlist)
		num_mz<-1
		
		for(mind in seq(1,length(mzlist),mz_group)){
		
			stopind<-mind+mz_group
			if(stopind>length(mzlist)){
			
			stopind<-length(mzlist)
			}
			s1<-mzlist[mind:stopind]
			s1<-unique(s1)
			num_mz<-num_mz+length(s1)
			
			if(num_mz%%50>0){
			Sys.sleep((syssleep/2))	
			}else{
			Sys.sleep(syssleep)	
			}
			if(length(s1)>1){
				repeat{
				cur.annot.res<-parLapply(cl,s1,feat.batch.annotation.child,max.mz.diff=max.mz.diff,adductname=adductname,adductmass=adductmass,adductcharge=adductcharge, syssleep=syssleep)
				if(is(cur.annot.res,"try-error")){
					Sys.sleep(10)
					cur.annot.res<-parLapply(cl,s1,feat.batch.annotation.child,max.mz.diff=max.mz.diff,adductname=adductname,adductmass=adductmass,adductcharge=adductcharge, syssleep=syssleep)
				
				}else{
				break
				}
				}
				mz.annot.res<-c(mz.annot.res,cur.annot.res)
				
			}else{
				for(i in 1:length(s1)){
					rescur<-feat.batch.annotation.child(mz.val=s1[i],max.mz.diff=max.mz.diff,adductname=adductname,adductmass=adductmass,adductcharge=adductcharge, syssleep=syssleep)
					#print(length(rescur))
					if(length(rescur)>0){
						rescur<-as.matrix(rescur)
						#print(dim(rescur))
						if(dim(rescur)[2]==1){
							rescur<-t(rescur)
							rescur<-as.data.frame(rescur)
						} 
						rescur<-as.data.frame(rescur)
						#print(dim(rescur))
					mz.annot.res<-c(mz.annot.res,rescur)
					}
				}
			}
			if(mind%%10>0){
				Sys.sleep((syssleep/10))			
			}else{
				Sys.sleep(syssleep)
				stopCluster(cl)
				cl<-parallel::makeCluster(numnodes)
				
			
				clusterEvalQ(cl, library(XML))
				clusterEvalQ(cl, library(RCurl))
				clusterEvalQ(cl, "feat.batch.annotation.child")
			}
			
		}
		
		stopCluster(cl)
		res={}
		#print(adductname)
			if(length(mz.annot.res)>0){
		for(mzl in 1:length(mz.annot.res))
		{
			res=rbind(res,mz.annot.res[[mzl]])
			
		}
		}
		res<-unique(res)
		
		text_res<-{}
		
		if(length(res)>0){
			
		adductname=c(rep(adductname,dim(res)[1]))
		
		temp_res<-cbind(adductname,res)
		
		temp_res<-as.matrix(temp_res)
		
		
		
		if(dim(temp_res)[2]==1){
		
			temp_res<-t(temp_res)
			temp_res<-as.data.frame(temp_res)
		
		}
		
		bad_rows<-which(temp_res[,2]=="1")
		
		if(length(bad_rows)>0){
		temp_res<-temp_res[-bad_rows,]
		temp_res<-as.matrix(temp_res)
		
		#temp_res<-t(temp_res)
			if(dim(temp_res)[2]==1){
		
			temp_res<-t(temp_res)
			
		
		}

		
		}
		#temp_res<-as.data.frame(temp_res)
		colnames(temp_res)=NULL
		
			
		#text_resindex<-c(1,2,5,6,7,8,11,10,13,15,17,19,21)
		#text_resindex<-c(1,2,5,6,7,8,9,12,11,14,16,18,20,22)
		
		#text_resindex<-c(1,2,5,6:8,9,11,13,14,16,18,20,21)
		#text_resindex<-c(1,2,5,6:7,4,8,10,12,13,17,19,21)
		
		text_resindex<-c(1,2,5,6:7,4,8,9,11,13,14,16,18,20,22)
		text_resindex<-text_resindex+1
		#print(dim(temp_res))
		text_res<-temp_res[,c(1,text_resindex)]
		
		text_res<-as.matrix(text_res)
		
		if(dim(text_res)[2]==1){
		
			text_res<-t(text_res)
			
		
		}
		text_res<-as.data.frame(text_res)
		bad_rows<-which(text_res[,2]=="1")
		
		
		if(length(bad_rows)>0){
		text_res<-text_res[-bad_rows,]
		text_res<-as.matrix(text_res)
		text_res<-t(text_res)
		}
		text_res<-as.data.frame(text_res)
		
		sernum=seq(1,dim(text_res)[1])
		text_res<-cbind(sernum,text_res)
		colnames(text_res)=cnames
		text_res<-text_res[,-c(5)]
		
		parentres=rbind(parentres,temp_res)
		rm(temp_res)
		colnames(parentres)=NULL
		
		}
		#num_cols<-dim(text_res)[2]
		#text_res<-cbind(text_res[,c(1:10)],text_res[,c(num_cols)],text_res[,c(11:(num_cols-1))])
		fname=paste(xMSanalyzer.outloc,"/KEGG_annotation_results_",queryadductlist[adnum],".txt",sep="")
		write.table(text_res,file=fname,sep="\t",row.names=FALSE)
		
		
		
		Sys.sleep(syssleep)
	
	}
	
	html_res<-{}
	text_res<-{}
	
		
	if(length(parentres)>0){
		
		res<-parentres[order(parentres[,2]),]

	#html_resindex<-c(1,2,5,4,6:7,9,11:12,14,16,18,20,22)
	res<-as.matrix(res)


	
	if(dim(res)[2]==1){res<-t(res)}
	
	
	
	#html_resindex<-c(1,2,5,6:7,9,11:12,14,16,18,20,22)
	
	#html_resindex<-c(1,2,5,6:7,4,8,10,12,13,15,17,19,21)
	html_resindex<-c(1,2,5,6:7,4,8,10,12,13,15,17,19,21,22)
	
	html_resindex<-html_resindex+1
	html_res<-res[,c(1,html_resindex)]
		html_res<-as.matrix(html_res)
		
	
		if(dim(html_res)[2]==1){html_res<-t(html_res)}
		
		sernum=seq(1,dim(html_res)[1])
		
		#html_res<-cbind(sernum,html_res)
		html_res<-as.data.frame(html_res)
		
		cnames<-c("Adduct","Query.m/z", "Search mass \n tolerance range (+/-)","Metlin", "Compound.Name",  "Chemical.Formula",  "Exact Mass", "CASID", "KEGG.Compound.ID",
		"KEGG.Pathway.ID", "KEGG.Pathway.Name", "HMDB.ID", "PubChem.Substance.ID", "ChEBI.ID", "LIPID.MAPS.ID","KEGG.Brite.Category")
 
		#"KEGG.Gene.ID",
		
		colnames(html_res)<-cnames
		
		html_res<-html_res[,-c(4)]
		#fname=paste("Annotation_results",sep="")
		#num_cols<-dim(html_res)[2]
		#text_res<-cbind(html_res[,c(1:10)],html_res[,c(num_cols)],html_res[,c(11:(num_cols-1))])
		
		fname=paste("KEGG_annotation_results",sep="")
		unlink(fname)
		#fname=paste(xMSanalyzer.outloc,"/KEGG_annotation_results.html",sep="")
		HTMLInitFile(filename=fname,Title="KEGG annotation results", outdir=xMSanalyzer.outloc)
		fname=paste(xMSanalyzer.outloc,"/KEGG_annotation_results.html",sep="")
		HTML(html_res,file=fname,Border=1,innerBorder=1,useCSS=TRUE)
		HTMLEndFile(file=fname)
		
		
		
		cnames<-c("","Adduct","Query.m/z", "Search mass \n tolerance range (+/-)","Metlin", "Compound.Name",  "Chemical.Formula",  "Exact Mass", "CASID", "KEGG.Compound.ID",
		"KEGG.Pathway.ID", "KEGG.Pathway.Name", "HMDB.ID", "PubChem.Substance.ID", "ChEBI.ID", "LIPID.MAPS.ID","KEGG.Brite.Category")
		
		
		text_resindex<-c(1,2,5,6:7,4,8,9,11,13,14,16,18,20,22)
		text_resindex<-text_resindex+1
		text_res<-res[,c(1,text_resindex)]
		
		text_res<-as.matrix(text_res)
		if(dim(text_res)[2]==1){text_res<-t(text_res)}
		
		if(length(text_res)>0){
		sernum=seq(1,dim(text_res)[1])
		}else{
			sernum={}
			}
		text_res<-cbind(sernum,text_res)
		
		text_res<-as.data.frame(text_res)
		colnames(text_res)=cnames
		text_res<-text_res[,-c(5)]
		#num_cols<-dim(text_res)[2]
		#text_res<-cbind(text_res[,c(1:10)],text_res[,c(num_cols)],text_res[,c(11:(num_cols-1))])
		fname=paste(xMSanalyzer.outloc,"/KEGG_annotation_results_alladducts.txt",sep="")
		write.table(text_res,file=fname,sep="\t",row.names=FALSE)
		}
	return(list("text.res"=text_res,"html.res"=html_res))
} 

#################################################################
#Function: feat.batch.annotation.child
#Description: search m/z in Metlin with links to KEGG, HMDB, PubChem, LipidMaps, ChEBI, and CAS 
#mzorig: single m/z value in ppm (eg: 121.945)
#max.mz.diff: +/- m/z match tolerance in ppm for database matching
#           (eg: 5)
#adductmass: mass of the adduct to be selected (eg: 1.00727 for M+H, or -35.012729 for M+H-2H2O)
##################################################################
##################################################################
feat.batch.annotation.child<-function(mz.val,max.mz.diff, adductname, adductmass=NA, adductcharge=NA, syssleep)
{
	
	adductlist=c(1.007276,22.989218,38.963158,-35.012676,-17.0027,0.0227,7.01597,18.033823,
	33.033486,42.033826,44.97116,64.015768,2.014552,23.996494,45.978436,3.021828,
	25.00377,46.985712,-19.01839,-1.007276,18.998371,20.974666,34.969402,36.948606,
	44.998194,59.013851,78.918885,-2.014552,-3.021828)
	alladducts<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M-H2O+NH4", "M+Li","M+NH4",
	"M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
	"M+2H+Na","M+2Na+H","M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
	"M+FA-H","M+CH3COO-H","M+Br","M-2H","M-3H")
	names(adductlist)<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M-H2O+NH4", "M+Li","M+NH4",
	"M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
	"M+2H+Na","M+2Na+H","M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
	"M+FA-H","M+CH3COO-H","M+Br","M-2H","M-3H")
	
	mult_charge<-c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3,3,3,1,1,1,1,1,1,1,1,1,2,3)
	names(mult_charge)<-c("M+H","M+Na","M+K","M+H-2H2O","M+H-H2O", "M-H2O+NH4", "M+Li","M+NH4",
	"M+CH3OH+H","M+ACN+H","M+2Na-H","M+ACN+Na","M+2H", "M+H+Na","M+2Na","M+3H",
	"M+2H+Na","M+2Na+H","M-H2O-H", "M-H", "M+F","M+Na-2H","M+Cl","M+K-2H",
	"M+FA-H","M+CH3COO-H","M+Br","M-2H","M-3H")
	
    
    adductmass=adductlist[as.character(adductname)]
    adductcharge=mult_charge[as.character(adductname)]
    
    print(adductmass)
    print(adductcharge)
	#print(mz.val)
	
	
        delta_ppm=(max.mz.diff)*(mz.val/1000000)
      
        min_mz=round((mz.val-delta_ppm),5)
        max_mz=round((mz.val+delta_ppm),5)
        
        print(mz.val)
        print(min_mz)
        print(max_mz)
        
        #convert to neutral mass
        min_mz=(min_mz*adductcharge)-adductmass
        max_mz=(max_mz*adductcharge)-adductmass
        
        
	
	res={} #c("-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-")
	mzorig=round(mz.val,5)
	delta_ppm=round(delta_ppm,5)
	
	syssleep1<-(syssleep/5)
	Sys.sleep(syssleep1)
	
    #write.table(mz.val,file="mzval.txt",sep="\t",row.names=FALSE)

	html_link<-"-"
	search_link=paste("http://rest.kegg.jp/find/compound/",min_mz,"-",max_mz,"/exact_mass",sep="")
	
	#print(search_link)
	d1<-try(readLines(search_link),silent=TRUE)

	if(is(d1,"try-error")){

	res<-c(mz.val,rep("NC",27))
	
	 write.table(mz.val,file="kegg_bad_mzs.txt",sep="\t",row.names=FALSE,append=TRUE)

	}else{
	
	
	
	#print(dim(d1))
	cnames<-c("ENTRY","NAME","FORMULA","EXACT_MASS","REACTION","PATHWAY","ENZYME","PubChem","ChEBI","PDB")
	
	#pattern_list<-c("C[0-9]{3,5}","[:blank:]{2,}[0-9|A-Z|:punct:|(|:print:][[:punct:]|[:alnum:]]*{3,}", "FORMULA","EXACT_MASS","CAS:","PubChem:","KNApSAcK:","PDB-CCD")
	
	#pattern_list<-c("C[0-9]{3,5}","NAME", "FORMULA","EXACT_MASS","ko[0-9]{5}","CAS:","ChEBI:","LIPIDMAPS:","PubChem:", "KNApSAcK:","PDB-CCD:")
	
	#pattern_list<-c("C[0-9]{3,5}","NAME", "FORMULA","EXACT_MASS","ko[0-9]{5}","CAS:","ChEBI:","LIPIDMAPS:","PubChem:", "KNApSAcK:","PDB-CCD:", "map")
	
	pattern_list<-c("EXACT_MASS","NAME", "FORMULA","CAS:","PubChem:","ChEBI:","LIPIDMAPS:", "BRITE", "map")
	
	pattern_keggid<-"C[0-9]{3,5}"
	
	#if(dim(d1)[1]>0)
	id_list<-"-"
		CName<-"-"
		mass<-"-"
		casID<-"-"
		keggID<-"-"
		kegglink<-"-"
		keggpathid<-"-"
		keggpathname<-"-"
		keggpathlink<-"-"
		hmdbID<-"-"
		hmdblink<-"-"
		pubchemsid<-"-"
		pubchemslink<-"-"
		pubchemcid<-"-"
		pubchemclink<-"-"
		chebiid<-"-"
		chebilink<-"-"
		lipidmapsid<-"-"
		lipidmapslink<-"-"
		chemformula<-"-"
		
	if(length(d1)>0){ 
	for(i in 1:length(d1))
	{
		if(i%%5>0){
		syssleep1<-(syssleep/5)
		Sys.sleep(syssleep1)
		}else{
		syssleep1<-(syssleep/3)
		Sys.sleep(syssleep1)
		}
		id_list<-"-"
		CName<-"-"
		mass<-"-"
		casID<-"-"
		keggID<-"-"
		kegglink<-"-"
		keggpathid<-"-"
		keggpathname<-"-"
		keggpathlink<-"-"
		hmdbID<-"-"
		hmdblink<-"-"
		pubchemsid<-"-"
		pubchemslink<-"-"
		pubchemcid<-"-"
		pubchemclink<-"-"
		chebiid<-"-"
		chebilink<-"-"
		lipidmapsid<-"-"
		lipidmapslink<-"-"
		chemformula<-"-"
		keggpathinf<-{}
		
		#l1<-grep(d1[i],pattern=pattern_list[5])
		str_text=d1[i]
		t2<-gregexpr(pattern=pattern_keggid,perl=FALSE,text=str_text)
			if(t2[[1]][1]>0)
			{
				t3=t2[[1]]
				strlength=attr(t3,"match.length")-1
				t4=strsplit(as.character(str_text),"")
				
				keggID<-t4[[1]][t3[1]:(t3[1]+strlength)]
			
				
				keggID<-paste(keggID,collapse="")
				kegglink<-paste("<a href=http://www.genome.jp/dbget-bin/www_bget?cpd:",keggID,">",keggID,"</a>",sep="")
				
				
				#html_res=readHTMLTable(kegglink)			
				search_link1=paste("http://rest.genome.jp/link/cpd:",keggID,"+-e",sep="")
				
				#dlink<-readLines(search_link1)
				dlink<-getURL(search_link1)
	
				if(dlink!=""){
					dlink<-read.delim(search_link1,header=FALSE)
					dlink2<-as.data.frame(dlink)
					if(dim(dlink2)[2]>0){
					for(l in 1:dim(dlink2)[1]){
						link_text=dlink2[l,2]
						t2<-gregexpr(pattern="HMDB[0-9]{2,}",perl=FALSE,text=link_text)
						t3=t2[[1]]
						strlength=attr(t3,"match.length")-1
						t4=strsplit(as.character(link_text),"")
						if(strlength>0)
						{
							hmdbID<-t4[[1]][t3[1]:(t3[1]+strlength)]
							hmdbID<-paste(hmdbID,collapse="")
							hmdblink<-paste("<a href=http://www.hmdb.ca/metabolites/",hmdbID,">",hmdbID,"</a>",sep="")
						}
					 }
					}
				}
				
				#keggID<-"C00392"
				
				#keggID<-"C00157"
		#keggID<-"C00082"
		search_link=paste("http://rest.kegg.jp/get/cpd:",as.character(keggID),sep="")
		d2<-read.delim(search_link,header=FALSE)
		d3<-as.data.frame(d2)
		
		
		
		
		if(length(d3)>0){
		pat.res<-{}
	

		url_vec<-{}
			url_strs<-c("-","-","-","-","<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=","<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:",
			"<a href=http://www.lipidmaps.org/data/get_lm_lipids_dbgif.php?LM_ID=","-")

		url_vec<-{}
		
		if(length(d3)>0){
		for(j in 1:(length(pattern_list)-1))
			{
				l1<-grep(d3[,1],pattern=pattern_list[j])
				if(length(l1)>0){
				for(ind1 in 1:length(l1)){
				
				str1<-gsub(as.character(d3[l1[ind1],1]),pattern=" ",replacement="_")
				s1<-strsplit(str1," ")
				#print(s1)
				if(j==8){
				p1<-paste("(DBLINKS)|[_]{2,}|;*",pattern_list[j],sep="")
				s2<-gsub(as.character(s1[[1]]),pattern=p1,replacement="")
				s2<-gsub(s2,pattern="_",replacement=" ")
				}else{
				p1<-paste("(DBLINKS)|[_]*|:*|;*",pattern_list[j],sep="")
				s2<-gsub(as.character(s1[[1]]),pattern=p1,replacement="")
				
				}
				#p1<-paste("([DBLINKS])|[_]|:|;",pattern_list[j],sep="")
				
				
				
				
				#paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:",chebiid,">",chebiid,"</a>",sep="")
				
				url_str_cur<-paste(url_strs[j],s2,">",s2,"</a>",sep="")
				
				#if(ind1>1){s2<-paste(s2,";",sep="")}
				
				pat.res<-c(pat.res,s2)
				pat.res<-c(pat.res,url_str_cur)		
				}
				}else{
				pat.res<-c(pat.res,rep("-",2))
				}
			}
				l1<-grep(d3[,1],pattern=pattern_list[length(pattern_list)])
				if(length(l1)>0){
				
				keggpathid<-""
				keggpathname<-""
				keggpathlink<-""
				
				for(ind1 in 1:length(l1)){
				temp.pat.res<-{}
				str1<-gsub(as.character(d3[l1[ind1],1]),pattern=" ",replacement="_")
				s1<-strsplit(str1," ")
				
				
				p1<-paste("(DBLINKS)|[_]{3,}|:*|;*|PATHWAY",sep="")
				#p1<-paste("([DBLINKS])|[_]|:|;",pattern_list[j],sep="")
				s2<-gsub(as.character(s1[[1]]),pattern=p1,replacement="")
				
				s3<-strsplit(s2,"__")
				#s2<-gsub(s2,"__",";",sep="")
				
				#temp.pat.res<-c(temp.pat.res,s3[[1]][1])
				keggpathurl<-paste("<a href=http://www.genome.jp/kegg-bin/show_pathway?",s3[[1]][1],"+",keggID,">",s3[[1]][1],"</a>",sep="")
				#temp.pat.res<-c(temp.pat.res,keggpathlink)				
				s4<-gsub(as.character(s3[[1]][2]),pattern="_",replacement=" ")
				
				#temp.pat.res<-c(temp.pat.res,s4)
				keggpathid<-paste(keggpathid,paste(s3[[1]][1],";",sep=""),sep="")
				keggpathlink<-paste(keggpathlink,paste(keggpathurl,";",sep=""),sep="<br>")
				keggpathname<-paste(keggpathname,paste(s4,";",sep=""),sep="<br>")
				}
				
				pat.res<-c(pat.res,keggpathid,keggpathlink,keggpathname)
				
				}else{
				pat.res<-c(pat.res,"-","-","-")
				}
				
				
		
				
			
		
		}
		

	#pattern_list<-c("EXACT_MASS","NAME", "FORMULA","CAS:","C[0-9]{3,5}","PubChem:","ChEBI:","LIPIDMAPS:", "map")
	
		#res<-rbind(res,c(mzorig,delta_ppm,as.character(id_list), mass, html_link, CName,chemformula,casID,keggID,kegglink,keggpathid,keggpathname,keggpathlink,hmdbID,hmdblink,pubchemsid,pubchemslink, pubchemcid,pubchemclink,chebiid,chebilink, lipidmapsid, lipidmapslink))
				
	res<-rbind(res,c(mzorig,delta_ppm,as.character(id_list), pat.res[1], html_link, pat.res[3],pat.res[5],pat.res[7],keggID,kegglink,pat.res[17],pat.res[18],pat.res[19],hmdbID,hmdblink,pat.res[9],pat.res[10],pat.res[11],pat.res[12],pat.res[13],pat.res[14],pat.res[15]))
			
	
	
     
	
			
		
		}
		}
		
		}
	}
	metres<-html_link
	#write.table(res,file="kegg_cur_res.txt",sep="\t",append=TRUE,row.names=FALSE)
	}
	syssleep1<-(syssleep/5)
	Sys.sleep(syssleep1)

	return(res)
}


#########################################################################################
###Function to find metabolic characteristics of individuals
#curdata: output matrix from apLCMS or XCMS
#min.samps: minimum number of samples in which a feature signal 
#	   should be detected in at least min.reps replicates
#min.reps: minimum proportion of replicates in which a signal is present (eg: 0.5 or 1)
#num_replicats: number of replicats for each sample
#alignment.tool: name of feature alignment tool eg: "apLCMS" or "XCMS"
########################################################################################

check.mz.in.replicates<-function(dataA,min.samps=2,min.reps=2,num_replicates=3)
{
        mean_replicate_difference<-{}
        sd_range_duplicate_pairs<-{}
	
	curdata<-dataA
	rm(dataA)

	#curdata<-curdata[1:10,]
        numfeats=dim(curdata)[1]
        numsamp=dim(curdata)[2]
        rnames<-colnames(curdata)
        rnames<-gsub(".cdf", "", rnames, ignore.case=TRUE)
	quantcolnames=c("min", "first_quartile", "median", "mean", "third_quartile", "max")
	

        
                newrow={}
                finalmat={}
		
		
		cl<-parallel::makeCluster(2)
			
		
		
		clusterExport(cl, "check.mz.in.replicates.child") 
		
		cv.res<-parApply(cl,curdata,1,check.mz.in.replicates.child,min.samps=min.samps,min.reps=min.reps,num_replicates=num_replicates)
		print("done")
	dim(cv.res)=dim(matrix(nrow=1,ncol=numfeats))	
	print("done")
	stopCluster(cl)
        #final_set<-as.data.frame(cv.res)
	
        #rownames(final_set)=NULL
	#final_set<-apply(final_set,2,as.numeric)
	#final_set<-as.data.frame(t(final_set))
	#colnames(final_set)<-quantcolnames
        final_set<-curdata[which(cv.res==1),]
	#final_set<-cbind(curdata_mz_rt_info,final_set)
	return(final_set)
}


check.mz.in.replicates.child<-function(curdata,min.samps,min.reps,num_replicates)
{
        
        numsamp=length(curdata)
        textp1=""
        t=1

        finalmat={}
      
	
                replicate_check=0
		intvec={}
		num.samps.check=0
                for(samp in seq(1,(numsamp),num_replicates))
                {
                        newrow={}
                	intvec={}
		        for(replicate in 1:num_replicates)
			{ 
				i=samp+replicate-1
				intvec=c(intvec,curdata[i])
				
			}
			if(length(which(intvec>0))>=(min.reps))
			{
                                	replicate_check=1
				 	num.samps.check=num.samps.check+1
			}
                
		}
		if(num.samps.check>=min.samps)
		{	
			replicate_check=1
		}
		else
		{
			replicate_check=0
		}
               
	return(replicate_check)

}



find.Overlapping<-function (dataA, dataB, mz.thresh = 10, time.thresh = 30, compare.intensities=FALSE, numnodes=2,nearest.time.match=FALSE) 
{

  data_a <- as.data.frame(dataA)
  data_b <- as.data.frame(dataB)
  rm(dataA)
  rm(dataB)
  
  colnames(data_a)[1] = "mz"
  colnames(data_b)[1] = "mz"
  
  library(parallel)

  cl <- parallel::makeCluster(getOption("cl.cores", numnodes))
      clusterExport(cl, "find.overlapping.single")
  
    
  mz_groups <- parLapply(cl,1:dim(data_a)[1], find.overlapping.single,time.thresh=time.thresh,data_a=data_a,data_b=data_b,mz.thresh=mz.thresh,nearest.time.match=nearest.time.match)
  
  stopCluster(cl)
  
  commat = data.frame()
  if (length(mz_groups) > 0) {
    
    for(j in 1:length(mz_groups)){
      if(dim(mz_groups[[j]])[1]>0){
        commat <- rbind(commat,mz_groups[[j]])
      }
    }
    if(dim(commat)[1]>0){
      rownames(commat)=1:dim(commat)[1]
    }else{
      
    }
    return(commat)
  }else{
    return(commat)
  }
  
}

find.overlapping.single<-function(j,time.thresh,data_a,data_b,mz.thresh,nearest.time.match) {
    commat = {}
      if (is.na(time.thresh) == FALSE) {
    colnames(data_a)[2] = "time"
    colnames(data_b)[2] = "time"
    mznames = c("index.A", "mz.data.A", "time.data.A", "index.B", 
                "mz.data.B", "time.data.B", 'mz.difference', "time.difference")
   # print("Using the 1st columns as 'mz' and 2nd column as 'retention time'")
  }else {
    mznames = c("index.A", "mz.data.A", "index.B", "mz.data.B", 'mz.difference')
    #print("Using the 1st columns as 'mz'")
  }
  
    ppmb = (mz.thresh) * (data_a$mz[j]/1e+06)
    getbind_same <- which(abs(data_b$mz - data_a$mz[j]) <= ppmb)
    
  
    if (is.na(time.thresh) == FALSE) {
      if (length(getbind_same) > 0) {
        
        tmp = suppressWarnings(cbind(cbind(j, data_a[j, c(1, 2)]),cbind(getbind_same, data_b[getbind_same, c(1, 2)])))
        tmp[,"mzdiff"]=(abs(tmp[,2]-tmp[,5])/tmp[,2])*10^6
        tmp[,"timediff"]=abs(tmp[,3]-tmp[,6])
        colnames(tmp)=mznames
        tmp = tmp[tmp$time.difference<=time.thresh,]
	
	print(tmp)
	if(nearest.time.match==TRUE){
		if(dim(tmp)[1]>0){
		  commat=tmp[which(tmp$time.difference==min(tmp$time.difference)),]
		}
	      
	}else{
		commat=tmp
	}
      
      }
      
    } else {
      if (length(getbind_same) > 0) {
        
        tmp = suppressWarnings(cbind(as.data.frame(cbind(j, data_a[j, c(1)])),as.data.frame(cbind(getbind_same, data_b[getbind_same, 1]))))
        tmp[,"mzdiff"]=(abs(tmp[,2]-tmp[,4])/tmp[,2])*10^6
        colnames(tmp)=mznames
        commat = tmp
        
      }
    }
    return(as.data.frame(commat))
  }

#Function:find.Overlapping.mzs
#Description: This function matches features between two or more datasets using the
#following user defined criteria:
#1) Maximum m/z difference (+/-) ppm
#2) Maximum retention time difference in seconds
#Input:
#data_a->apLCMS output for dataset A,
#data_b->apLCMS output for dataset B,
#max.mz.diff->Maximum m/z difference (+/-) ppm
#max.rt.diff->Maximum retention time difference in seconds
#
#Output:
#Data frame that includes mz and retention time of common features
#
#Usage:
#common_features<-matchFeaturesmulti(data_a, data_b, max.mz.diff=10, max.rt.diff=300)
############################################
find.Overlapping.mzs<-function(dataA, dataB, mz.thresh=10, time.thresh=NA, alignment.tool=NA,nearest.time.match=FALSE,numnodes=2)
{

        data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	#data_a<-unique(data_a)
	rm(dataA)
	rm(dataB)
        
     #   data_b<-unique(data_b)
        
        com_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}

        commat={}

	col.names.dataA=colnames(data_a)
	col.names.dataB=colnames(data_b)


	if(is.na(alignment.tool)==FALSE){
		 if(alignment.tool=="apLCMS")
		{
		      sample.col.start=5
			data_a<-data_a[,c(1,2)]
			    data_b<-data_b[,c(1,2)]
			    col.names.dataA<-col.names.dataA[1:2]
			    col.names.dataB<-col.names.dataB[1:2]
			    col.names.dataA[1]="mz"
			    col.names.dataA[2]="time"
			    col.names.dataB[1]="mz"
			    col.names.dataB[2]="time"
			    colnames(data_a)=col.names.dataA
			    colnames(data_b)=col.names.dataB
		}
		else
		{
		      if(alignment.tool=="XCMS")
		      {
			    sample.col.start=9
			
			    data_a<-data_a[,c(1,4)]
			    data_b<-data_b[,c(1,4)]
			    col.names.dataA<-col.names.dataA[1:2]
			    col.names.dataB<-col.names.dataB[1:2]
			    col.names.dataA[1]="mz"
			    col.names.dataA[2]="time"
			    col.names.dataB[1]="mz"
			    col.names.dataB[2]="time"
			    colnames(data_a)=col.names.dataA
			    colnames(data_b)=col.names.dataB
		      }
		      
		}
	
	}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		sample.col.start=3	   
		    col.names.dataA[1]="mz"
		   
		    col.names.dataB[1]="mz"
		    if(is.na(time.thresh)==FALSE){
		     col.names.dataA[2]="time"
		     col.names.dataB[2]="time"
		      print("Using the 1st column as \"mz\" and 2nd column as \"retention time\"")
		     }else{
		      print("Using the 1st column as \"mz\"")
		     }
		    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
	}

       #data_a<-data_a[order(data_a$mz),]
       #data_b<-data_b[order(data_b$mz),]
       data_a<-as.data.frame(data_a)
	data_b<-as.data.frame(data_b)
	colnames(data_a)=col.names.dataA
	colnames(data_b)=col.names.dataB
        #create header for the matrix with common features
	if(is.na(time.thresh)==FALSE){
	mznames=c("index.A","mz.data.A", "time.data.A", "index.B","mz.data.B","time.data.B", "time.difference") 
        }else{
	mznames=c("index.A","mz.data.A", "index.B","mz.data.B") 
	}
	
	commat<-find.Overlapping(dataA=data_a, dataB=data_b, mz.thresh = mz.thresh, time.thresh = time.thresh,numnodes=numnodes,nearest.time.match=nearest.time.match) 

        return(commat)
}


find.Unique.mzs.sameset<-function(dataA, dataB, mz.thresh=10, time.thresh=NA, alignment.tool=NA)
{

#print(dim(dataA))
        data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	data_a<-unique(data_a)
	rm(dataA)
	rm(dataB)
        
        data_b<-unique(data_b)
        
        com_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}

        commat={}

	col.names.dataA=colnames(data_a)
	col.names.dataB=colnames(data_b)

	if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
		data_a<-data_a[,c(1,2)]
		    data_b<-data_b[,c(1,2)]
		    col.names.dataA<-col.names.dataA[1:2]
		    col.names.dataB<-col.names.dataB[1:2]
                    col.names.dataA[1]="mz"
                    col.names.dataA[2]="time"
		    col.names.dataB[1]="mz"
                    col.names.dataB[2]="time"
                    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                    data_a<-data_a[,c(1,4)]
		    data_b<-data_b[,c(1,4)]
		    col.names.dataA<-col.names.dataA[1:2]
		    col.names.dataB<-col.names.dataB[1:2]
                    col.names.dataA[1]="mz"
                    col.names.dataA[2]="time"
		    col.names.dataB[1]="mz"
                    col.names.dataB[2]="time"
                    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		    col.names.dataB[1]="mz"
		    if(is.na(time.thresh)==FALSE){
		     col.names.dataA[2]="time"
		     col.names.dataB[2]="time"
		      print("Using the 1st column as \"mz\" and 2nd column as \"retention time\"")
		     }else{
		      print("Using the 1st column as \"mz\"")
		     }
		    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
	}

       #data_a<-data_a[order(data_a$mz),]
       #data_b<-data_b[order(data_b$mz),]
       data_a<-as.data.frame(data_a)
	data_b<-as.data.frame(data_b)
	colnames(data_a)=col.names.dataA
	colnames(data_b)=col.names.dataB
	
        #create header for the matrix with common features
	if(is.na(time.thresh)==FALSE){
	mznames=c("index.A","mz.data.A", "time.data.A", "index.B","mz.data.B","time.data.B", "time.difference") 
        }else{
	mznames=c("index.A","mz.data.A", "index.B","mz.data.B") 
	}
	
	overlap_res<-find.Overlapping.mzs(dataA=data_a, dataB=data_b, mz.thresh, time.thresh, alignment.tool)
	
    #print(overlap_res)

	filt_ind<-which(overlap_res$index.A!=overlap_res$index.B)
	
	if(length(filt_ind)>0){
	uniqueA<-data_a[-c(filt_ind),]
	}else{
		uniqueA<-data_a
	}
	
        return(list("uniqueA"=uniqueA))
}


find.Unique.mzs<-function(dataA, dataB, mz.thresh=10, time.thresh=NA, alignment.tool=NA)
{

	
        data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	data_a<-unique(data_a)
	rm(dataA)
	rm(dataB)
        
        data_b<-unique(data_b)
        
        com_mz_num=1
        unique_mz={}
        ppm_v={}
        rt_v={}

        commat={}

	col.names.dataA=colnames(data_a)
	col.names.dataB=colnames(data_b)

	if(is.na(alignment.tool)==FALSE){
	 if(alignment.tool=="apLCMS")
        {
              sample.col.start=5
		data_a<-data_a[,c(1,2)]
		    data_b<-data_b[,c(1,2)]
		    col.names.dataA<-col.names.dataA[1:2]
		    col.names.dataB<-col.names.dataB[1:2]
                    col.names.dataA[1]="mz"
                    col.names.dataA[2]="time"
		    col.names.dataB[1]="mz"
                    col.names.dataB[2]="time"
                    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
        }
        else
        {
              if(alignment.tool=="XCMS")
              {
                    sample.col.start=9
                   data_a<-data_a[,c(1,4)]
		    data_b<-data_b[,c(1,4)]
		    col.names.dataA<-col.names.dataA[1:2]
		    col.names.dataB<-col.names.dataB[1:2]
                    col.names.dataA[1]="mz"
                    col.names.dataA[2]="time"
		    col.names.dataB[1]="mz"
                    col.names.dataB[2]="time"
                    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
              }
	      
	}}else{
                    #stop(paste("Invalid value for alignment.tool. Please use either \"apLCMS\" or \"XCMS\"", sep=""))
		   
		    col.names.dataA[1]="mz"
		   
		    col.names.dataB[1]="mz"
		    if(is.na(time.thresh)==FALSE){
		     col.names.dataA[2]="time"
		     col.names.dataB[2]="time"
		      print("Using the 1st column as \"mz\" and 2nd column as \"retention time\"")
		     }else{
		      print("Using the 1st column as \"mz\"")
		     }
		    colnames(data_a)=col.names.dataA
                    colnames(data_b)=col.names.dataB
	}

       #data_a<-data_a[order(data_a$mz),]
       #data_b<-data_b[order(data_b$mz),]
       data_a<-as.data.frame(data_a)
	data_b<-as.data.frame(data_b)
	colnames(data_a)=col.names.dataA
	colnames(data_b)=col.names.dataB
	
        #create header for the matrix with common features
	if(is.na(time.thresh)==FALSE){
	mznames=c("index.A","mz.data.A", "time.data.A", "index.B","mz.data.B","time.data.B", "time.difference") 
        }else{
	mznames=c("index.A","mz.data.A", "index.B","mz.data.B") 
	}
	
	overlap_res<-find.Overlapping.mzs(dataA=data_a, dataB=data_b, mz.thresh, time.thresh, alignment.tool)
	
	

	if(length(overlap_res$index.A)>0){
	uniqueA<-data_a[-c(overlap_res$index.A),]
	}else{
		uniqueA<-data_a
	}
	
	if(length(overlap_res$index.B)>0){
	uniqueB<-data_b[-c(overlap_res$index.B),]
	}else{
		uniqueB<-data_b
	}
        return(list("uniqueA"=uniqueA,"uniqueB"=uniqueB))
}


#########################################################
#
#
#
#
#########################################################
getVenn<-function(dataA,name_a, dataB,name_b,mz.thresh=10,time.thresh=30,alignment.tool=NA, xMSanalyzer.outloc=NA,use.unique.mz=FALSE,plotvenn=TRUE,nearest.time.match=FALSE,numnodes=2)
{


	
	if(is.na(xMSanalyzer.outloc)==FALSE){
		dir.create(xMSanalyzer.outloc,showWarnings=FALSE)
	}
	data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	rm(dataA)
	rm(dataB)

	############################################

	
	if(use.unique.mz==TRUE){
		data_a<-find.Unique.mzs.sameset(dataA=data_a,dataB=data_a,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
		data_a<-data_a$uniqueA
		
		data_b<-find.Unique.mzs.sameset(dataA=data_b,dataB=data_b,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
		data_b<-data_b$uniqueA
	}
	common<-find.Overlapping.mzs(data_a,data_b,mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool,nearest.time.match=nearest.time.match,numnodes=numnodes)
	
	if(length(common)>0){
		commonA<-data_a[c(common$index.A),]
		commonA<-unique(commonA)
		
		commonB<-data_b[c(common$index.B),]
		commonB<-unique(commonB)
		
		data_a<-data_a[-c(common$index.A),]
		data_b<-data_b[-c(common$index.B),]
		
		

		rm_index<-which(data_a$mz%in%common$mz.data.A)
		
		if(length(rm_index)>0){
			uniqueA<-data_a[-rm_index,]
		}else{
			uniqueA<-data_a
		}
		
		rm_index<-which(data_b$mz%in%common$mz.data.B)
		
		if(length(rm_index)>0){
			uniqueB<-data_b[-rm_index,]
		}else{
			uniqueB<-data_b
		}
		
	
		num_commonA<-length(unique(common$index.A))
		num_commonB<-length(unique(common$index.B))
		
		num_common<-min(num_commonA,num_commonB)[1]
	}else{
		uniqueA<-data_a
		uniqueB<-data_b
		num_common<-0
		commonA<-{}
		commonB<-{}
	}
	num_commonA<-num_common
	num_commonB<-num_common
	num_uniqueA<-dim(uniqueA)[1]
	num_uniqueB<-dim(uniqueB)[1]
	

	g1 <-c(seq(1,(num_commonA+num_uniqueA)))
	g2<-c(seq(1,(num_commonB+num_uniqueB)))

	g1[1:num_commonA]=paste("x_",g1[1:num_commonA],sep="")
	g2[1:num_commonB]=paste("x_",g2[1:num_commonB],sep="")
	
	if(num_uniqueA>0){
	g1[(num_commonA+1):(num_commonA+num_uniqueA)]=paste("y_",g1[(num_commonA+1):(num_commonA+num_uniqueA)],sep="")
	}
	if(num_uniqueA>0){
		g2[(num_commonB+1):(num_commonB+num_uniqueB)]=paste("z_",g2[(num_commonB+1):(num_commonB+num_uniqueB)],sep="")
	}
	set1=as.character(g1)
	set2=as.character(g2)
	universe <- sort(unique(c(set1,set2)))
	
	Counts <- matrix(0, nrow=length(universe), ncol=2)
	colnames(Counts) <- c(name_a, name_b)
	for (i in 1:length(universe))
	{
		Counts[i,1] <- universe[i] %in% set1
		Counts[i,2] <- universe[i] %in% set2
	}
	
	if(is.na(xMSanalyzer.outloc)==TRUE){
		venn_counts<-try(vennCounts(Counts),silent=TRUE)
		if(plotvenn==TRUE){
			try(vennDiagram(venn_counts),silent=TRUE)
		}
		
	}else{
		fname<-paste(xMSanalyzer.outloc,"/Venn", name_a,"_",name_b,"_",mz.thresh,"ppm",time.thresh,"s.pdf",sep="")
		venn_counts<-try(vennCounts(Counts),silent=TRUE)
		if(plotvenn==TRUE){
			
			pdf(fname)
			try(vennDiagram(venn_counts),silent=TRUE)
			dev.off()
		}
		
	}
	return(list("common"=common,"commonA"=commonA,"uniqueA"=uniqueA,"commonB"=commonB,"uniqueB"=uniqueB,"vennCounts"=venn_counts))
	
}

getVennmultiple<-function(dataA,name_a, dataB,name_b,dataC,name_c,mz.thresh=10,time.thresh=30,alignment.tool=NA, xMSanalyzer.outloc, use.unique.mz=FALSE,plotvenn=TRUE)
{
	dir.create(xMSanalyzer.outloc,showWarnings=FALSE)
	
	 data_a<-as.data.frame(dataA)
        data_b<-as.data.frame(dataB)
	 data_c<-as.data.frame(dataC)
	
	rm(dataA)
	rm(dataB)
	rm(dataC)
	
	data_a<-as.data.frame(data_a)
	data_b<-as.data.frame(data_b)
	data_c<-as.data.frame(data_c)

	if(use.unique.mz==TRUE){
	    data_a<-find.Unique.mzs.sameset(data_a,data_a,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
	    data_a<-data_a$uniqueA
	    
	data_b<-find.Unique.mzs.sameset(data_b,data_b,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
	 data_b<-data_b$uniqueA
	    
	    data_c<-find.Unique.mzs.sameset(data_c,data_c,mz.thresh=mz.thresh,time.thresh=time.thresh,alignment.tool=alignment.tool)
	     data_c<-data_c$uniqueA
	}
	commonAB<-find.Overlapping.mzs(data_a,data_b,mz.thresh,time.thresh,alignment.tool=alignment.tool)
	commonAC<-find.Overlapping.mzs(data_a,data_c,mz.thresh,time.thresh,alignment.tool=alignment.tool)
	commonBC<-find.Overlapping.mzs(data_b,data_c,mz.thresh,time.thresh,alignment.tool=alignment.tool)
	
	#if(length(commonAB)>0)
	{
	data_ab<-data_a[commonAB$index.A,]
	
	data_ab<-as.data.frame(data_ab)
	
	commonABC<-find.Overlapping.mzs(data_ab,data_c,mz.thresh,time.thresh,alignment.tool=alignment.tool)

	}
	data_ba<-data_b[commonAB$index.B,]
	data_ba<-as.data.frame(data_ba)
	
	commonBAC<-find.Overlapping.mzs(data_ba,data_c,mz.thresh,time.thresh=NA,alignment.tool=alignment.tool)
	
	#get unique A
	rm_index<-which(data_a$mz%in%commonAB$mz.data.A)
	
	rm_index<-c(rm_index,which(data_a$mz%in%commonAC$mz.data.A))
	
	if(length(rm_index)>0){
	uniqueA<-data_a[-rm_index,]
	uniqueA<-as.data.frame(uniqueA)
	}
	

	rm_index<-which(data_b$mz%in%commonAB$mz.data.B)
	
	rm_index<-c(rm_index,which(data_b$mz%in%commonBC$mz.data.A))
	
	#rm_index<-which(uniqueB$mz%in%commonBC$mz.data.A)

	if(length(rm_index)>0){
	uniqueB<-data_b[-rm_index,]
	uniqueB<-as.data.frame(uniqueB)
	}	
	
	
	rm_index<-which(data_c$mz%in%commonAC$mz.data.B)

	rm_index<-c(rm_index,which(data_c$mz%in%commonBC$mz.data.B))

	if(length(rm_index)>0){
	uniqueC<-data_c[-rm_index,]
	uniqueC<-as.data.frame(uniqueC)
	}
	
	
	
	num_commonAB<-length(unique(commonAB$mz.data.A))-length(which(unique(commonAB$mz.data.A)%in%unique(commonABC$mz.data.A)))
	
	num_commonBC<-length(unique(commonBC$mz.data.A))-length(which(unique(commonBC$mz.data.A)%in%unique(commonBAC$mz.data.A)))
	
	num_commonAC<-length(unique(commonAC$mz.data.A))-length(which(unique(commonAC$mz.data.A)%in%unique(commonABC$mz.data.A)))
	
	num_commonCB<-length(unique(commonBC$mz.data.B))-length(which(unique(commonBC$mz.data.B)%in%unique(commonBAC$mz.data.B)))
	
	num_commonCA<-length(unique(commonAC$mz.data.B))-length(which(unique(commonAC$mz.data.B)%in%unique(commonABC$mz.data.B)))
	
	
	num_commonABC<-min(length(unique(commonABC$mz.data.A)),length(unique(commonABC$mz.data.B)))
	
	num_uniqueA<-dim(uniqueA)[1]
	num_uniqueB<-dim(uniqueB)[1]
	num_uniqueC<-dim(uniqueC)[1]
	
	g1 <-paste("a",seq(num_commonAB+num_commonAC+num_commonABC+num_uniqueA),sep="")

	g2<-paste("b",seq(num_commonAB+num_commonBC+num_commonABC+num_uniqueB),sep="")
	
	g3<-paste("c",seq(num_commonCA+num_commonCB+num_commonABC+num_uniqueC),sep="")

	#x: AB; w:AC; v:BC;u:ABC
	if(num_commonAB>0)
	{
		g1[1:num_commonAB]=paste("x_",seq(1,num_commonAB),sep="")
		
		g2[1:num_commonAB]=paste("x_",seq(1,num_commonAB),sep="")
	}
	
	if(num_commonAC>0)
	{
		g1[(num_commonAB+1):(num_commonAB+num_commonAC)]=paste("w_",seq(1,num_commonAC),sep="")
		
		g3[1:num_commonAC]=paste("w_",seq(1,num_commonAC),sep="")
		
	}
	
	if(num_commonBC>0)
	{
		g2[(num_commonAB+1):(num_commonAB+num_commonBC)]=paste("v_",seq(1,num_commonBC),sep="")
	
		g3[(num_commonAC+1):(num_commonAC+num_commonBC)]=paste("v_",seq(1,num_commonBC),sep="")
	}
	
	g1[(num_commonAB+num_commonAC+1):(num_commonAB+num_commonAC+num_commonABC)]=paste("u_",seq(1,num_commonABC),sep="")
	g2[(num_commonAB+num_commonBC+1):(num_commonAB+num_commonBC+num_commonABC)]=paste("u_",seq(1,num_commonABC),sep="")
	g3[(num_commonAC+num_commonBC+1):(num_commonAC+num_commonBC+num_commonABC)]=paste("u_",seq(1,num_commonABC),sep="")
	

	set1=as.character(g1)
	set2=as.character(g2)
	set3=as.character(g3)
	universe <- sort(unique( c(set1,set2,set3)))
	Counts <- matrix(0, nrow=length(universe), ncol=3)
	colnames(Counts) <- c(name_a, name_b,name_c)
	for (i in 1:length(universe))
	{
		Counts[i,1] <- universe[i] %in% set1
		Counts[i,2] <- universe[i] %in% set2
		Counts[i,3] <- universe[i] %in% set3
	}
	venn_counts<-vennCounts(Counts)
	fname<-paste(xMSanalyzer.outloc,"/Venn", name_a,"_",name_b,"_",name_c,"_",mz.thresh,"ppm",time.thresh,"s.pdf",sep="")
	
	if(plotvenn==TRUE){
		pdf(fname)
		vennDiagram(venn_counts)
		dev.off()
	}
	return(list("commonABC"=commonABC,"uniqueA"=uniqueA,"uniqueB"=uniqueB,"uniqueC"=uniqueC,"commonAB"=commonAB, 
	"commonBC"=commonBC,"commonAC"=commonAC,"vennCounts"=venn_counts))

}




simpleAnnotation<-function(dataA,max.mz.diff=10,num_nodes=2,queryadductlist=c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H","M+NH4","M+Na","M+ACN+H","M+ACN+Na","M+2ACN+H","2M+H","2M+Na","2M+ACN+H"),
gradienttype="Acetonitrile",mode="pos",outloc,db_name="KEGG")
{
    
    
    if(db_name=="KEGG"){
        
        
        data(keggCompMZ)
        chemCompMZ<-keggCompMZ
        suppressWarnings(rm(keggCompMZ))
    }else{
        if(db_name=="HMDB"){
            data(hmdbCompMZ)
            chemCompMZ<-hmdbCompMZ
            suppressWarnings(rm(hmdbCompMZ))
        }else{
            if(db_name=="T3DB"){
                data(t3dbCompMZ)
                chemCompMZ<-t3dbCompMZ
                suppressWarnings(rm(t3dbCompMZ))
            }else{
                
                if(db_name=="LipidMaps"){
                    data(lipidmapsCompMZ)
                    chemCompMZ<-lipidmapsCompMZ
                    suppressWarnings(rm(lipidmapsCompMZ))
                }
            }
            
        }
    }
    NOPS_check=TRUE
    chemCompMZ<-as.data.frame(chemCompMZ)

    chemCompMZ$mz<-as.numeric(as.character(chemCompMZ$mz))
   
    chemCompMZ$mz<-round(chemCompMZ$mz,4)
    
    dataA<-as.data.frame(dataA)

    dataA$mz<-round(dataA$mz,4) 
    data(adduct_table)
    allowWGCNAThreads(nThreads=num_nodes)
    #rm(adduct_table)
    #data(adduct_table)
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="2M\\+2H",replacement="M+H")
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="3M\\+3H",replacement="M+H")
    
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="2M\\+2Na",replacement="M+Na")
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="3M\\+3Na",replacement="M+Na")
    
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="2M\\+2K",replacement="M+K")
    adduct_table$Adduct<-gsub(adduct_table$Adduct,pattern="3M\\+3K",replacement="M+K")
    
    adduct_table<-unique(adduct_table)
    
    suppressWarnings(dir.create(outloc))
    
    if(queryadductlist=="all" & mode=="pos"){
        
        adduct_names<-adduct_table$Adduct[(adduct_table$Type=="S" & adduct_table$Mode=="positive") | (adduct_table$Type==gradienttype & adduct_table$Mode=="positive")]
        
        adduct_table<-adduct_table[which(adduct_table$Adduct%in%adduct_names),]
        
    }else{
        if(queryadductlist=="all" & mode=="neg"){
            
            adduct_names<-adduct_table$Adduct[(adduct_table$Type=="S" & adduct_table$Mode=="negative") | (adduct_table$Type==gradienttype & adduct_table$Mode=="negative")]
            adduct_table<-adduct_table[which(adduct_table$Adduct%in%adduct_names),]
        }else{
            
            adduct_names<-adduct_table$Adduct[which(adduct_table$Adduct%in%queryadductlist)]
            
            adduct_table<-adduct_table[which(adduct_table$Adduct%in%adduct_names),]
            
        }
    }
    
    adduct_names<-unique(adduct_names)
    
    
    
    
    chemCompMZ<-chemCompMZ[which(chemCompMZ$Adduct%in%adduct_names),]
    
    print("name of DB")
    print(db_name)
    print("dim of Db")
    print(dim(chemCompMZ))
   library(parallel)
   
    cl<-parallel::makeCluster(num_nodes)
    
    clusterEvalQ(cl, library(XML))
    clusterEvalQ(cl, library(R2HTML))
    clusterEvalQ(cl, library(RCurl))
    clusterEvalQ(cl, library(SSOAP))
    clusterEvalQ(cl, library(limma))
    
    clusterEvalQ(cl, library(plyr))
    
        clusterEvalQ(cl, library(parallel))
    
    clusterEvalQ(cl, "processWSDL")
    clusterEvalQ(cl, library(png))
    clusterExport(cl, "Annotationbychemical_IDschild")
    
    clusterExport(cl, "find.Overlapping.mzs")
    clusterExport(cl, "find.Overlapping")
     clusterExport(cl, "find.overlapping.single")
	   clusterExport(cl, "getVenn")
       
       
       s1<-seq(1,length(adduct_names))
       print("Mapping m/z to metabolites:")
       l2<-parLapply(cl,s1,Annotationbychemical_IDschild,dataA=dataA, queryadductlist=c(adduct_names),adduct_type=c("S",gradienttype),
       adduct_table=adduct_table,max.mz.diff=max.mz.diff,outloc=outloc,chemCompMZ=chemCompMZ,otherdbs=FALSE,otherinfo=FALSE)
       
       stopCluster(cl)
    
       
       levelB_res<-{};
       for(j in 1:length(l2)){
           if(length(l2[[j]])>1){
               levelB_res<-rbind(levelB_res,l2[[j]])
           }
       }
       
       rm(l2)
       
       
       levelB_res$mz<-as.numeric(as.character(levelB_res$mz))
       
       levelB_res$time<-as.numeric(as.character(levelB_res$time))
       
       levelB_res<-as.data.frame(levelB_res)
      
      #save(levelB_res,file="levelB_res.Rda")
 
       uniq_formula<-as.character(unique(levelB_res$Formula))
       
       bad_formula<-which(is.na(uniq_formula)==TRUE)
       if(length(bad_formula)>0){
           uniq_formula<-uniq_formula[-c(bad_formula)]
       }
       
       cl<-parallel::makeCluster(num_nodes)
       
       clusterExport(cl, "check_golden_rules")
       clusterExport(cl, "check_element")
       #clusterExport(cl, "uniq_formula")
       #clusterExport(cl, "NOPS_check")
       
       levelB_res_check<-parLapply(cl,1:length(uniq_formula),function(j,uniq_formula,NOPS_check){
           
           curformula<-as.character(uniq_formula[j])
           return(check_golden_rules(curformula,NOPS_check=NOPS_check))
           
       },uniq_formula=uniq_formula,NOPS_check=NOPS_check)
       stopCluster(cl)
       
       #save(levelB_res_check,file="xMSannotator_levelB_check.Rda")
       levelB_res_check2<-ldply(levelB_res_check,rbind)
       
       levelB_res_check3<-levelB_res_check2[which(levelB_res_check2[,2]==1),]
       
       
       levelB_res<-levelB_res[which(levelB_res$Formula%in%as.character(levelB_res_check3[,1])),]
       
       water_adducts<-c("M+H-H2O","M+H-2H2O","M-H2O-H")
       
       water_adduct_ind<-which(levelB_res$Adduct%in%water_adducts)
       
       cl<-parallel::makeCluster(num_nodes)
       
       
       clusterExport(cl, "check_element")
       
       
       
       if(length(water_adduct_ind)>0){
           levelB_res2<-levelB_res[c(water_adduct_ind),]
           
           levelB_res<-levelB_res[-c(water_adduct_ind),]
           
           sind1<-seq(1:dim(levelB_res2)[1])
           
           levelB_res_check3<-parLapply(cl,sind1,function(j){
               
               adduct<-as.character(levelB_res2$Adduct[j])
               curformula<-as.character(levelB_res2$Formula[j])
               
               numoxygens<-check_element(curformula,"O")
               
               if(numoxygens>0){
                   bool_check<-1
               }else{
                   bool_check<-0
               }
               
               res<-cbind(curformula,bool_check)
               res<-as.data.frame(res)
               return(res)
               
               
           })
           
           levelB_res_check4<-ldply(levelB_res_check3,rbind)
           
           valid_form<-{}
           
           if(length(which(levelB_res_check4[,2]==1))>0){
               levelB_res_check4<-levelB_res_check4[which(levelB_res_check4[,2]==1),]
               
               
               valid_form<-which(levelB_res2$Formula%in%as.character(levelB_res_check4[,1]))
           }
           if(length(valid_form)>0){
               levelB_res2<-levelB_res2[valid_form,]
               levelB_res<-rbind(levelB_res,levelB_res2)
           }
           
       }
       multiresmat<-levelB_res
       
	rm(levelB_res)
       dupmz<-multiresmat$mz[which(duplicated(multiresmat$mz)==TRUE)]
       
       
       MatchCategory<-rep("Multiple",dim(multiresmat)[1])
      
	if(length(dupmz)>0){ 
      	 MatchCategory[-which(multiresmat$mz%in%dupmz)]<-"Unique"
       }
       levelB_res<-cbind(MatchCategory,multiresmat)
       
       #save(levelB_res,file="xMSannotator_levelB.Rda")
       return(levelB_res)
}




Annotationbychemical_IDschild<-function(adduct_index=NA,dataA,queryadductlist=c("M+H"),adduct_type=c("S","Acetonitrile"),adduct_table,max.mz.diff=10,outloc, otherdbs=FALSE,otherinfo=FALSE,chemCompMZ){
    
    #load("~/Documents/Emory/JonesLab/Projects/xMSannotator/keggCompMZ.Rda")
    dataA<-as.data.frame(dataA)
    adduct_names<-as.character(adduct_table[,1])
    adductlist<-adduct_table[,4]
    mult_charge<-adduct_table[,3]
    num_mol<-adduct_table[,2]
    names(adductlist)<-as.character(adduct_names)
    names(mult_charge)<-as.character(adduct_names)
    names(num_mol)<-as.character(adduct_names)
    alladducts<-adduct_names
    
    
    keggCompMZ<-chemCompMZ
    
    rm(chemCompMZ)
    
    if(is.na(adduct_index)==FALSE){
        
        queryadductlist=queryadductlist[adduct_index]
    }
    
    #load("~/Documents/Emory/JonesLab/Projects/xMSannotator/keggCompMZ.Rda")
    
    alladducts<-adduct_names
    #print(queryadductlist)
    #print(alladducts)
    if(queryadductlist[1]=="positive")
    {
        queryadductlist<-adduct_names[which(adduct_table[,5]=="positive")]
        
    }else{
        if(queryadductlist[1]=="negative")
        {
            
            queryadductlist<-adduct_names[which(adduct_table[,5]=="negative")]
            
        }else{
            if(queryadductlist[1]=="all"){
                
                
                queryadductlist<-alladducts
                
                
            }else{
                if(length(which(queryadductlist%in%alladducts==FALSE))>0){
                    
                    errormsg<-paste("Adduct should be one of:",sep="")
                    for(i in alladducts){errormsg<-paste(errormsg, i,sep=" ; ")}
                    stop(errormsg, "\n\nUsage: feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"M+H\", \"M+Na\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"positive\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"negative\"), xMSannotator.outloc, numnodes=1)",
                    "\n\n OR  feat.batch.annotation.KEGG(dataA,max.mz.diff=10, queryadductlist=c(\"all\"), xMSannotator.outloc, numnodes=1)"
                    )
                }
                
            }
        }
    }
    adduct_table<-as.data.frame(adduct_table)
    
    adduct_table<-adduct_table[which(adduct_table$Type%in%adduct_type),] #=="S" | adduct_table$Type=="Acetonitrile",]
    
    #adduct_table<-adduct_table[which(adduct_table$Mode%in%adduct_mode),] #=="positive",]
    
    adduct_table<-adduct_table[which(adduct_table$Adduct%in%queryadductlist),] #adduct_table$Adduct=="M+H" | adduct_table$Adduct=="M+Na",]
    
    
    suppressWarnings(dir.create(outloc))
    
    setwd(outloc)
    #keggres<-KEGG.annotation(dataA=mz_search_list,queryadductlist = c("positive"),xMSannotator.outloc)
    
    
    #cur_fname<-"/Users/karanuppal/Documents/Emory/JonesLab/Projects/MaHPIC/Exp2/c18/apLCMS_with_xMSanalyzer_merged_data/apLCMS_feature_list_at_p1_U_p2cor0.7_CV100.txt"
    #dataA<-read.table(cur_fname,sep="\t",header=TRUE)
    
    mz_search_list_1<-as.data.frame(keggCompMZ[which(keggCompMZ$Adduct%in%adduct_table$Adduct),c(1,7)])
    
    
    
    mz_search_list_1<-apply(mz_search_list_1,2,as.numeric)
   
    gcur<-getVenn(dataA=dataA,name_a="Experimental",name_b="DB",dataB=mz_search_list_1,mz.thresh=max.mz.diff,time.thresh=NA,
    xMSanalyzer.outloc=outloc,alignment.tool=NA)

    #print("here 1")
    #print(length(gcur$common))
    #print("here 2")
     
    if(length(gcur$common)>0){
    
	keggCompMZ$mz<-round(keggCompMZ$mz,4)
	gcur$common$mz.data.B<-round(gcur$common$mz.data.B,4)
	gcur$common$mz.data.A<-round(gcur$common$mz.data.A,4)
	    
	dataA$mz<-round(dataA$mz,4)

        mcur<-merge(keggCompMZ[which(keggCompMZ$Adduct%in%adduct_table$Adduct),],gcur$common,by.x="mz",by.y="mz.data.B")
      
        #print(mcur[1,1:4])
        #print(dataA[1,])
       
        mcur_ind<-{} #paste(mcur$HMDBID,"_",mcur$index.A,sep="")

	if(length(which(duplicated(mcur_ind)==TRUE))>0){
	mcur<-mcur[-which(duplicated(mcur_ind)==TRUE),] 	 
        
	}

        mcur_2<-merge(mcur,dataA,by.x="mz.data.A",by.y="mz")
        
        mcur_3<-mcur_2[order(mcur_2$Name),]
        
        mcur_3<-mcur_3[,-c(9,10)]
        
	
        cnames<-colnames(mcur_3)
        cnames[1]<-"mz"
        cnames[2]<-"theoretical.mz"
        colnames(mcur_3)<-as.character(cnames)
        
        mcur_4<-as.data.frame(mcur_3)
        rm(keggCompMZ)
        rm(dataA)
        rm(mcur_2)
        rm(mcur_3)
        rm(mz_search_list_1)
        
        
        if(dim(mcur)[1]>1){
            
            h1<-table(mcur_4$mz) #Adduct
            
            if(length(h1)>0){
                #u1<-c(u1,which(h1<=1))
                u1<-which(h1<=1)
            }
            match_status<-rep("Multiple",dim(h1)[1])
            
            
            uniq_kegg_matches<-names(u1)
            
            
            
            match_status[u1]<-"Unique"
            
            
            
            match_status_mat<-cbind(rownames(h1),match_status)
            
            
            
            colnames(match_status_mat)<-c("mz","MatchCategory")
            match_status_mat<-as.data.frame(match_status_mat)
            
            mcur_5<-merge(match_status_mat,mcur_4,by="mz")
            
            
            rm(mcur_4)
            
            #h1<-table(mcur_5$chemical_ID,mcur_5$mz)
            
            #s2<-apply(h1,2,sum)
            
            
            mcur_5<-as.data.frame(mcur_5)
            
            #mcur_5$mz<-as.numeric(mcur_5$mz)
            
            mcur_5<-mcur_5[order(mcur_5$mz),]
            
        }else{
            MatchCategory<-"Unique"
            cnames1<-colnames(mcur_4)
            cnames1<-c(cnames1[1],"MatchCategory",cnames[-c(1)])
            mcur_5<-cbind(mcur_4[1,1],MatchCategory,mcur_4[,-c(1)])
            mcur_5<-as.data.frame(mcur_5)
            colnames(mcur_5)<-as.character(cnames1)
            
            
            #mcur_5$mz<-as.numeric(mcur_5$mz)
            
            mcur_5<-mcur_5[order(mcur_5$mz),]
            
        }
        if(otherinfo==TRUE){
            info_mat<-sapply(1:dim(mcur_5)[1],function(j){
                
                
                b1<-keggGet(paste("cpd:",mcur_5[j,1],sep=""))
                brite_inf<-paste(b1[[1]]$BRITE,collapse=";")
                path_inf<-paste(b1[[1]]$PATHWAYS,collapse=";")
                otherdb_inf<-paste(b1[[1]]$DBLINKS,collapse=";")
                r1<-c(as.character(mcur_5[j,1]),as.character(brite_inf),as.character(path_inf),as.character(otherdb_inf))
                
                return(r1)
            })
            
            
            info_mat_1<-as.data.frame(t(info_mat))
            colnames(info_mat_1)<-c("chemical_ID","BriteCategory","Pathways","ExternalLinks")
            
            
            mcur_6<-merge(info_mat_1,mcur_5,by="chemical_ID")
            
            mcur_7<-unique(mcur_6)
            
            rm(mcur_6)
            
            if(otherdbs==TRUE){
                info_mat_2<-sapply(1:dim(mcur_7)[1],function(j){
                    
                    b1<-keggLink(paste("cpd:",mcur_7[j,1],"+-e",sep=""))
                    hmdbID<-"-"
                    lipidmapsID<-"-"
                    
                    link_text<-b1[,2]
                    
                    t2<-gregexpr(pattern="hmdb:",perl=FALSE,text=link_text)
                    
                    if(length(t2)>1){
                        g_ind<-which(t2==1)
                        
                        if(length(g_ind)>0){
                            if(length(g_ind)>1){
                                for(g in g_ind){
                                    t3=t2[[g]]
                                    
                                    hmdbID<-c(hmdbID,gsub(b1[g,2],pattern="hmdb:",replacement=""))
                                }
                                if(length(g_ind)>1){hmdbID<-paste(hmdbID,collapse=";")}
                            }else{
                                
                                hmdbID<-gsub(b1[g_ind,2],pattern="hmdb:",replacement="")
                            }
                        }
                    }
                    
                    
                    t2<-gregexpr(pattern="lipidmaps:",perl=FALSE,text=link_text)
                    
                    if(length(t2)>1){
                        g_ind<-which(t2==1)
                        
                        if(length(g_ind)>0){
                            
                            if(length(g_ind)>1){
                                for(g in g_ind){
                                    t3=t2[[g]]
                                    
                                    lipidmapsID<-c(lipidmapsID,gsub(b1[g,2],pattern="lipidmaps:",replacement=""))
                                    
                                    
                                }
                                lipidmapsID<-paste(lipidmapsID,collapse=";")
                                
                            }else{lipidmapsID<-gsub(b1[g_ind,2],pattern="lipidmaps:",replacement="")}
                            
                        }
                        
                    }
                    
                    
                    return(list(keggid=as.character(mcur_7[j,1]),hmdb=hmdbID,lipidmaps=lipidmapsID))
                    c1<-c(as.character(mcur_7[j,1]),hmdbID,lipidmapsID)
                    c1<-as.data.frame(c1)
                    return(c1)
                })
                
                info_mat_3<-{}
                #for(i in 1:dim(info_mat_2)[1]){
                
                cdata<-rbind(info_mat_2[1,],info_mat_2[2,])
                cdata<-rbind(cdata,info_mat_2[3,])
                cdata<-as.data.frame(cdata)
                info_mat_3<-rbind(info_mat_3,cdata)
                
                
                #}
                
                #info_mat_3<-as.data.frame(t(info_mat_2))
                info_mat_3<-t(info_mat_3)
                colnames(info_mat_3)<-c("chemical_ID","HMDBID","LIPIDMAPS")
                
                mcur_7<-as.data.frame(mcur_7)
                
                mcur_8<-cbind(info_mat_3,mcur_7) #,by="chemical_ID")
                mcur_8<-unique(mcur_8)
                rownames(mcur_8)<-NULL
                return(mcur_8)
            }else{
                mcur_7<-as.data.frame(mcur_7)
                
                rownames(mcur_7)<-NULL
                return(mcur_7)
            }
        }else{
            mcur_5<-unique(mcur_5)
            return(mcur_5)
            
        }
        
    }else{return("no match found.")}
    #}else{return("no match found.")}
}



##################################################################
#Function:
#Description: wrapper function based on apLCMS.align,evaluate.Features,
#            evaluate.Samples,merge.Results,search.Metlin, and
#            search.KEGG
################################################################
xMSwrapper.apLCMS<-function(cdfloc, apLCMS.outloc, xMSanalyzer.outloc, min.run.list=c(12,3), min.pres.list=c(0.5,0.8), minexp.pct=0.1, mztol=1e-5, alignmztol=1e-5,
alignchrtol=NA,numnodes=NA,run.order.file=NA, apLCMSmode="untargeted", known_table=NA, match_tol_ppm=5, max.mz.diff=10,max.rt.diff=300,merge.eval.pvalue=0.2,
mergecorthresh=0.7,deltamzminmax.tol=100, subs=NA ,num_replicates=3,
mz.tolerance.dbmatch=10, adduct.list=c("M+H"), samp.filt.thresh=0.70,feat.filt.thresh=50,cormethod="pearson", mult.test.cor=TRUE,
missingvalue=0,ignore.missing=TRUE,filepattern=".cdf",sample_info_file=NA,refMZ=NA,refMZ.mz.diff=10,refMZ.time.diff=NA,void.vol.timethresh=30,
replacezeroswithNA=TRUE,
scoreweight=30,charge_type="pos",syssleep=0.5,plotEICs="target",rawprofileloc=NA,peak.score.thresh=0,baseline.correct.noise.percentile = 0.25, 

    shape.model = "bi-Gaussian", baseline.correct = NA, peak.estim.method = "moment", 

    min.bw = NA, max.bw = NA, sd.cut = c(1, 60), sigma.ratio.lim = c(0.33, 

        3), max.align.mz.diff = 0.01, pre.process = FALSE, recover.mz.range = NA, 

recover.chr.range = NA, use.observed.range = TRUE, recover.min.count = 3, new_feature_min_count=4, component.eliminate = 0.01, moment.power = 2,
reference_sample_index=NA,merge.pairwise=FALSE,min.samp.percent=0.6,impute.bool=TRUE,summarize.replicates=TRUE,
summary.method="median",max.number.of.replicates.with.missingvalue=1,summary.na.replacement="zeros",db_name=c("KEGG","HMDB","LipidMaps"),qc_label=NA,data.norm.pipeline="AC",calibration.method=c("median.adjustment","multiplicative.signal.correction"))
{

	suppressWarnings(suppressWarnings(sink(file=NULL)))
    
    dir.create(xMSanalyzer.outloc,showWarnings=FALSE)

	library(parallel)
	
    rep.num.max.missing.thresh=max.number.of.replicates.with.missingvalue
    #variables not used
    impute.method.MAR="KNN"
    impute.method.MNAR="QRILC"
    # summarize.replicates=FALSE
    #summary.method="median"
    #rep.num.max.missing.thresh=1
    #summary.na.replacement="zeros"
	#cdfloc=NA
	
	targeted_feat_raw<-{}
	
    
    x<-date()
    x<-strsplit(x,split=" ")
	x1<-unlist(x)
	
    inp_params<-{}
    inp_params<-rbind(inp_params,cbind("cdfloc: ",cdfloc))
    inp_params<-rbind(inp_params,cbind("apLCMS.outloc: ",apLCMS.outloc))
    inp_params<-rbind(inp_params,cbind("xMSanalyzer.outloc: ",xMSanalyzer.outloc))
    inp_params<-rbind(inp_params,cbind("num_replicates: ",num_replicates))
    inp_params<-rbind(inp_params,cbind("max.mz.diff: ",max.mz.diff))
    inp_params<-rbind(inp_params,cbind("max.rt.diff: ",max.rt.diff))
    inp_params<-rbind(inp_params,cbind("merge.eval.pvalue: ",merge.eval.pvalue))
    inp_params<-rbind(inp_params,cbind("mergecorthresh: ",mergecorthresh))
    inp_params<-rbind(inp_params,cbind("deltamzminmax.tol: ",deltamzminmax.tol))
    inp_params<-rbind(inp_params,cbind("mz.tolerance.dbmatch: ",mz.tolerance.dbmatch))
    inp_params<-rbind(inp_params,cbind("adduct.list: ",paste(adduct.list,collapse=";")))
    inp_params<-rbind(inp_params,cbind("samp.filt.thresh: ",samp.filt.thresh))
    inp_params<-rbind(inp_params,cbind("feat.filt.thresh: ",feat.filt.thresh))
    inp_params<-rbind(inp_params,cbind("cormethod: ",cormethod))
    inp_params<-rbind(inp_params,cbind("charge_type: ",charge_type))
    inp_params<-rbind(inp_params,cbind("db_name: ",paste(db_name,collapse=";")))
    inp_params<-rbind(inp_params,cbind("impute.bool: ",impute.bool))
    inp_params<-rbind(inp_params,cbind("merge.pairwise: ",merge.pairwise))

    inp_params<-rbind(inp_params,cbind("summarize.replicates: ",summarize.replicates))
    inp_params<-rbind(inp_params,cbind("summary.method: ",summary.method))
    inp_params<-rbind(inp_params,cbind("max.number.of.replicates.with.missingvalue: ",rep.num.max.missing.thresh))
    inp_params<-rbind(inp_params,cbind("summary.na.replacement: ",summary.na.replacement))
     inp_params<-rbind(inp_params,cbind("data.norm.pipeline: ", data.norm.pipeline))
   
	    inp_params<-rbind(inp_params,cbind("calibration.method: ",calibration.method[1]))
   
    fname<-paste(xMSanalyzer.outloc,"/Input_parameters.csv",sep="")
    
    colnames(inp_params)<-c("Parameter","Value")
    write.csv(inp_params,file=fname,row.names=FALSE)
    
    
    x1<-gsub(x1,pattern=":",replacement="_")
	fname<-paste(x1[2:4],collapse="")
	fname<-gsub(fname,pattern=":",replacement="_")
	#fname<-paste(fname,x1[5],sep="")
	x1[5]<-gsub(x1[5],pattern=":",replacement="_")
	fname<-paste(fname,x1[5],sep="_")

	if(is.na(rawprofileloc)==TRUE){

		setwd(apLCMS.outloc)
		setwd("../")
		rawprofileloc=getwd() #apLCMS.outloc
	}	
	
	fname<-paste(xMSanalyzer.outloc,"/Log_",fname,".txt",sep="")
	print(paste("Program is running. Please check the logfile for runtime status: ",fname,sep=""))

	 data_rpd_all=new("list")
     peak_score_list=new("list")
	data_rpd=new("list")
	union_list=new("list")
	feateval_list<-new("list")
	rsqres_list<-new("list")
	annot.res<-{}
	annotres_list<-new("list")
	minexp<-2
	
	if(is.na(cdfloc)==FALSE){
	
		if(cdfloc=="NA"){
			cdfloc=NA
		}
	}
	
	if(is.na(sample_info_file)==FALSE){
	
		if(sample_info_file=="NA"){
			sample_info_file=NA
		}
	}
    #print(is.na(cdfloc))
    #print(cdfloc)
    sampleid_mapping<-{}
    
	filepattern=".cdf|.mzxml|mXML"
	 if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
    {
    	
    	
		sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)

		if(is.na(cdfloc)==FALSE){
			l1<-list.files(cdfloc,filepattern,ignore.case=TRUE)

			if(is.na(minexp.pct)==FALSE){	
			minexp<-round(minexp.pct*length(l1))
			}else{
			minexp<-2

			}
		}else{
				l1<-rep(1,num_replicates)	
			}
		if(length(l1)!=dim(sampleid_mapping)[1] & (is.na(cdfloc)==FALSE))
		{
			num_mis_files<-dim(sampleid_mapping)[1]-length(l1)
			
			if(is.na(subs)==TRUE){
				stop(paste("ERROR: Only ",length(l1)," spectral files were found. ",num_mis_files," files are missing.",sep=""))
			}
		}
	}else{
	
	
		if(is.na(cdfloc)==FALSE){
			l1<-list.files(cdfloc,filepattern,ignore.case=TRUE)
			
			 if(is.na(minexp.pct)==FALSE){
                        minexp<-round(minexp.pct*length(l1))
                        }else{
                        minexp<-2

                        }
			
		
		}else{
			l1<-rep(1,num_replicates)	
		}
	}
if(length(l1)%%num_replicates>0)
{stop(paste("ERROR: Not all samples have ",num_replicates," replicates.",sep=""))
}


dir.create(apLCMS.outloc,showWarnings=FALSE)
	dir.create(xMSanalyzer.outloc,showWarnings=FALSE)
	setwd(xMSanalyzer.outloc)	
	
	sink(fname)
	print(sessionInfo())
	    if(is.na(refMZ)==FALSE){
			stddata<-read.table(refMZ,sep="\t",header=TRUE)
			print(refMZ)
			print(head(stddata))
			
		}else{
			if(charge_type=="pos"){
			data(example_target_list_pos)
			stddata<-example_target_list_pos
			}else{
				
				if(charge_type=="neg"){
						data(example_target_list_neg)
						stddata<-example_target_list_neg
					}else{
						stop("Invalid option. \'charge_type\' should be \'pos\' or \'neg\'.")
						}
				}
		}
        ############################################
        #1) Align profiles using the cdf.to.ftr wrapper function in apLCMS
        if(is.na(apLCMS.outloc)==TRUE)
	{
		stop("Undefined value for parameter, apLCMS.outloc. Please define the apLCMS output location.")

	}
	 if(is.na(xMSanalyzer.outloc)==TRUE)
        {
                stop("Undefined value for parameter, xMSanalyzer.outloc. Please define the xMSanalyzer output location.")

        }
	
	

	if(is.na(cdfloc)==FALSE)
        {
                setwd(cdfloc)
                if(is.na(apLCMS.outloc)==FALSE)
                {
			setwd(apLCMS.outloc)
			
			fname<-paste("apLCMSalign.Rda")
			test1<-suppressWarnings(try(load(file=fname),silent=TRUE))
			
			if(is(test1, "try-error")){
						
			data_rpd_all<-apLCMS.align(cdfloc=cdfloc, apLCMS.outloc=apLCMS.outloc,min.run.list=min.run.list, min.pres.list=min.pres.list, minexp=minexp, mztol=mztol, 
			alignmztol=alignmztol, alignchrtol=alignchrtol, 
			numnodes=numnodes,run.order.file=run.order.file,subs=subs,filepattern=filepattern,apLCMSmode=apLCMSmode, 
			known_table=known_table, match_tol_ppm=match_tol_ppm,
			refMZ.mz.diff=refMZ.mz.diff,refMZ.time.diff=refMZ.time.diff,target.mz.list = stddata,plotEICs=plotEICs,peak.score.thresh=peak.score.thresh,
			num_replicates=num_replicates,
			baseline.correct.noise.percentile = baseline.correct.noise.percentile, 

			    shape.model = shape.model, baseline.correct = baseline.correct, peak.estim.method = peak.estim.method, 

			    min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, sigma.ratio.lim = sigma.ratio.lim,  

			    max.align.mz.diff = max.align.mz.diff, pre.process = pre.process, recover.mz.range = recover.mz.range, 

			    recover.chr.range = recover.chr.range, use.observed.range = use.observed.range, 
			    recover.min.count = recover.min.count,reference_sample_index=reference_sample_index,
			    cvthresh=feat.filt.thresh,xMSanalyzer.outloc=xMSanalyzer.outloc,component.eliminate =component.eliminate, moment.power = moment.power)
    
			setwd(apLCMS.outloc)
			fname<-paste("apLCMSalign.Rda")
			save(data_rpd_all,file=fname)
   
   
			data_rpd<-new("list")
			
			#return(data_rpd_all)
			alignmentresults<-data_rpd_all$filenames_list
			peak_score_list<-data_rpd_all$peak_score_list
			data_rpd_all<-data_rpd_all$aligned_data_list
		
			cnames<-colnames(data_rpd_all[[i]])

                cnames[1]<-"mz"
                cnames[2]<-"time"
                colnames(data_rpd_all[[i]])<-cnames
	
			}else{
			alignmentresults<-data_rpd_all$filenames_list
			peak_score_list<-data_rpd_all$peak_score_list
			data_rpd_all<-data_rpd_all$aligned_data_list
				cnames<-colnames(data_rpd_all[[i]])

                cnames[1]<-"mz"
                cnames[2]<-"time"
                colnames(data_rpd_all[[i]])<-cnames


			}
			
                }
                else
                {
                        stop("Undefined value for parameter, apLCMS.outloc. Please define the output location.")
                }
        }else{
	
			setwd(apLCMS.outloc)
			fname<-paste("apLCMSalign.Rda")
			test1<-suppressWarnings(try(load(file=fname),silent=TRUE))
            
            #print("here")
			if(is(test1, "try-error")){
			
            print("apLCMSalign.Rda not found. Using the following files as apLCMS output:")
            
            alignmentresults<-list.files(".","*.txt")
            
	    print(alignmentresults)
	    
            for(i in 1:length(alignmentresults)){
	    
	   
                
                data_rpd_all[[i]]<-read.table(paste(apLCMS.outloc,"/",alignmentresults[i],sep=""),header=TRUE)
                #data_rpd_all[[i]]<-data_rpd_all[[i]]
                
                #data_rpd_all[[i]]<-data_rpd_all[[i]][1:1000,]
                data_rpd_all[[i]]<-unique(data_rpd_all[[i]])
                
                peak_score_list[[i]]<-rep(1,nrow(data_rpd_all[[i]]))
                
                plotEICs="none"

		cnames<-colnames(data_rpd_all[[i]])

		cnames[1]<-"mz"
		cnames[2]<-"time"
		colnames(data_rpd_all[[i]])<-cnames
            
	    }
            
			}else{
			alignmentresults<-data_rpd_all$filenames_list
			peak_score_list<-data_rpd_all$peak_score_list
			data_rpd_all<-data_rpd_all$aligned_data_list
			cnames<-colnames(data_rpd_all[[i]])

	                cnames[1]<-"mz"
        	        cnames[2]<-"time"
                	colnames(data_rpd_all[[i]])<-cnames


			}
	}
	
		setwd(apLCMS.outloc)
                
		if(is.na(sample_info_file)==FALSE)
		{
			cnames_cur=colnames(data_rpd_all[[1]][,-c(1:4)])
			cnames_cur<-gsub(cnames_cur,pattern=".mzXML|.cdf|.mzML",replacement="",ignore.case=T)
			
			
			 match_names_check<-match(sampleid_mapping[,1],cnames_cur)
			 
			
			 
			 if(length(which(is.na(match_names_check)==TRUE))>0){
				stop("Sample names do not match between sequence file and feature table.")
			 }
	
		}
		
		
		
		print("Filenames for apLCMS feature tables:")
		print(alignmentresults)
		for(i in 1:length(alignmentresults)){
		
			#data_rpd_all[[i]]<-read.table(paste(apLCMS.outloc,"/",alignmentresults[i],sep=""),header=TRUE)
			#data_rpd_all[[i]]<-data_rpd_all[[i]]
	
			data_rpd_all[[i]]<-unique(data_rpd_all[[i]])		
			#if(is.na(cdfloc)==TRUE)
			{
			l1<-dim(data_rpd_all[[i]])[2]-4
	
			}
			
			
			numfeats<-apply(data_rpd_all[[i]][,-c(1:4)],1,function(x){
					return(countpeaks(x, missingvalue))
					})
					
					numfeats<-unlist(numfeats)
					 if(is.na(minexp.pct)==FALSE)
					{
					minexp<-round(minexp.pct*l1)
					
					if(length(which(numfeats>=minexp))>0){
					#print(paste("Removing ",length(which(numfeats<minexp))," features based on missing value criteria",sep=""))
					#data_rpd_all[[i]]<-data_rpd_all[[i]][which(numfeats>=minexp),]
					}else{
						stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
						}
					
					
					}
			peak_score_list[[i]]<-unlist(peak_score_list[[i]])
		}
	
	
	
	#subdir1<-paste(xMSanalyzer.outloc,"/Quality_assessment_files",sep="")
	#subdir2<-paste(xMSanalyzer.outloc,"/apLCMS_filtered_data",sep="")
	#subdir3<-paste(xMSanalyzer.outloc,"/apLCMS_with_xMSanalyzer_merged_data",sep="")
	
	#dir.create(subdir1,showWarnings=FALSE)
	#dir.create(subdir2,showWarnings=FALSE)
	#dir.create(subdir3,showWarnings=FALSE)
        
        subdir1<-paste(xMSanalyzer.outloc,"/Stage1",sep="")  #QC individual parameter settings
        subdir2<-paste(xMSanalyzer.outloc,"/Stage2",sep="")  #Data filtering
        subdir3<-paste(xMSanalyzer.outloc,"/Stage3a",sep="")  #Data merger/parameter optimization
        subdir3b<-paste(xMSanalyzer.outloc,"/Stage3b",sep="")
        subdir4a<-paste(xMSanalyzer.outloc,"/Stage4a",sep="")	 #Raw QC: batch effect eval, TIC, etc
        
	
	dir.create(subdir1,showWarnings=FALSE)
	dir.create(subdir2,showWarnings=FALSE)
	dir.create(subdir3,showWarnings=FALSE)
	dir.create(subdir3b,showWarnings=FALSE)
	dir.create(subdir4a,showWarnings=FALSE)
	
	  #if(is.na(sample_info_file)==FALSE)
	   if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
    {
		subdir4b<-paste(xMSanalyzer.outloc,"/Stage4b",sep="")	 #Batch-effect corrected QC: batch effect eval, TIC, etc
		dir.create(subdir4b,showWarnings=FALSE)
		
	}
	
	if(is.na(adduct.list)==FALSE){
	subdir5<-paste(xMSanalyzer.outloc,"/Stage5",sep="")	 #Putative unprocessed annotations;
	

	dir.create(subdir5,showWarnings=FALSE)
	}
	
	bestscore<-(-1000000)
	
                #stop("Undefined value for parameter, cdfloc. Please enter path of the folder where the CDF files to be processed are located.")
                #change location to the output folder
                setwd(apLCMS.outloc)
              
                if(length(alignmentresults)>0)
                {
                          curdata_dim={}
                         
			  
			  parent_bad_list<-{}
			  
                          if(num_replicates==2)
                          {
                                  fileroot="_PID"
                          }
                          else
                          {
                                  if(num_replicates>2)
                                  {
                                          fileroot="_CV"
                                  }
                                  else
                                  {
                                          fileroot=""
                                            #stop("Need at least 2 technical replicates per sample.")
			        }
                          }
			  
			  pfilenames<-{}
			  
			   
			  
			   print("*******xMSanalyzer Stage 1: QC evaluation of invidual parameters*******")
			   
			   cat("\n")

			rep_cor_mat<-{}
            parent_bad_list<-{}
                          for(i in 1:length(alignmentresults))
                          {
                          	
				  print(paste("*Evaluating apLCMS results from parameter setting ",i,"*",sep=""))				  
                                  ############################################
                                  #2)Calculate pairwise correlation coefficients
                                  file_name=sapply(strsplit(alignmentresults[[i]],".txt"),head)
				  
				  print(paste("*File name: ",alignmentresults[[i]],"*",sep=""))	
				  
				  crow<-c(paste("p",i,sep=""),file_name)
				  pfilenames<-rbind(pfilenames,crow)

                                  curdata=data_rpd_all[[i]]
                                  #############################################
                                  ############################################
                                  #3) Calculate Percent Intensity Difference
deltamzminmax.tol=NA
                             #  if(num_replicates>1)
                               {   
				  
				peak_score_vec=unlist(peak_score_list[[i]])
					
				print("Start")
                                feat.eval.result=evaluate.Features(curdata, numreplicates=num_replicates,min.samp.percent=min.samp.percent,alignment.tool="apLCMS",impute.bool=impute.bool,
				peak_scores=peak_score_vec,numnodes=numnodes)
                                
				print("End")
                                cnames=colnames(feat.eval.result)
                                feat.eval.result<-apply(feat.eval.result,2,as.numeric)
                                feat.eval.result<-as.data.frame(feat.eval.result)
                                feat.eval.result.mat=cbind(curdata[,c(1:4)],feat.eval.result)
				
                                feat.eval.outfile=paste(subdir1,"/",file_name,fileroot,"featureassessment.txt",sep="")
				
				print(feat.eval.outfile)
				print(head(feat.eval.result.mat))
                                #write results
                                write.table(feat.eval.result.mat, feat.eval.outfile,sep="\t", row.names=FALSE)
                                
                                curdata<-curdata[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),]
                                
                                feat.eval.result.mat<-feat.eval.result.mat[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),]
                                
                                 
                                curdata<-as.data.frame(curdata)
                                curdata<-replace(as.matrix(curdata),which(is.na(curdata)==TRUE),0)
                                
				
                          
                                
                                feateval_list[[i]]<-feat.eval.result.mat
                                data_rpd_all[[i]]<-curdata
								
								  if(num_replicates>1)
								  {
									  print(paste("**calculating pairwise ",cormethod," correlation**",sep=""))

									  
									  
											rsqres_list<-evaluate.Samples(curdata, num_replicates, alignment.tool="apLCMS", cormethod,missingvalue,ignore.missing)

											rsqres<-as.data.frame(rsqres_list$cor.matrix)
											
											curdata<-as.data.frame(rsqres_list$feature.table)
											snames<-colnames(curdata[,-c(1:4)])
											snames_1<-snames[seq(1,length(snames),num_replicates)]
											rownames(rsqres)<-snames_1
											pcor_outfile=paste(subdir1,"/",file_name,"_sampleassessment_usinggoodfeatures.txt",sep="")
											write.table(rsqres, pcor_outfile,sep="\t",row.names=TRUE)
											rsqres_list[[i]]<-rsqres
											
                                            t1_fname<-paste(subdir1,"/","rsqres.Rda",sep="")
                                            # save(rsqres,file=t1_fname)
									
								if(num_replicates>2)
								{
									rep_cor_mat<-cbind(rep_cor_mat,rsqres$meanCorrelation)
									bad_samples<-which(rsqres$meanCorrelation<samp.filt.thresh)
								}else
								{
									bad_samples<-which(rsqres<samp.filt.thresh)
									rep_cor_mat<-rsqres
								}
								
								if(length(bad_samples)>0){
									bad_sample_names<-snames_1[bad_samples]
									
									feat.eval.outfile=paste(subdir1,"/",file_name,"_badsamples_at_cor",samp.filt.thresh,".txt",sep="")
									bad_sample_names<-as.data.frame(bad_sample_names)
									colnames(bad_sample_names)<-paste("Samples with correlation between technical replicates <", samp.filt.thresh,sep="")
									write.table(bad_sample_names, file=feat.eval.outfile,sep="\t", row.names=FALSE)
								}
								
								bad_list={}
								if(length(bad_samples)>0)
								{
									for(n1 in 1:length(bad_samples))
									{	
										if(bad_samples[n1]>1)
										{
											bad_samples[n1]=bad_samples[n1]+(bad_samples[n1]-1)*(num_replicates-1)
										}
											
									}
									for(n1 in 1:num_replicates)
									{
										bad_list<-c(bad_list,(bad_samples+n1-1))
									}
									bad_list<-bad_list[order(bad_list)]
									
								}
								if(i>1){
										parent_bad_list<-intersect(parent_bad_list,bad_list)
									}
									else{
									    
									    parent_bad_list<-bad_list

									}
									  }
								  else
								  {
									  print("**skipping sample evaluataion as only one replicate is present**")
								  }
								
								#curdata<-cbind(curdata,feat.eval.result[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),])
								
								
								
							
                                  }
				#  else
				#  {
					#  print("*********skipping feature evaluataion as only one replicate is present******")
				#  }
                                
	       	}
	       	
		file_name="parameter_filenames.txt"
		pfile_name=paste(subdir3,"/",file_name,sep="")
		write.table(pfilenames,file=pfile_name,sep="\t",row.names=FALSE)
		
		cat("\n")
		 print("*********Stage 2: Quality based filtering of results from each parameter setting based on sample and feature checks********")
		 cat("\n")

		if(num_replicates>2)
		{
			max_rep_mean_cor<-apply(rep_cor_mat,1,max)
		}else{
			max_rep_mean_cor<-rep_cor_mat
		}
        #print(head(max_rep_mean_cor))
        
        #save(max_rep_mean_cor,file="max_rep_mean_cor.Rda")

		 if(length(parent_bad_list)>0){
             
             
             if(is.na(sample_info_file)==FALSE && sample_info_file!="NA"){
					sampleid_mapping<-sampleid_mapping[-c(parent_bad_list),]
					sampleidfilteredfile=paste(subdir2,"/","filtered_sampleid_mapping.txt",sep="")
					write.table(sampleid_mapping,file=sampleidfilteredfile,sep="\t",row.names=FALSE)
             }

		}

		  for(i in 1:length(alignmentresults))
                          {
                                  file_name=sapply(strsplit(alignmentresults[[i]],".txt"),head)
                                  feat.eval.file=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,deltamzminmax.tol,"ppmmzrangefiltereddata.txt",sep="")	
                                 
				  
				  curdata<-data_rpd_all[[i]]
				  feat.eval.result.mat<-feateval_list[[i]]
				  if(length(parent_bad_list)>0){
					  curdata<-curdata[,-c(parent_bad_list+4)]
					  #maxint<-apply(curdata[,-c(1:4,((dim(curdata)[2]-6):dim(curdata)[2]))],1,max)
					  maxint<-apply(curdata[,-c(1:4)],1,max)
					  badfeats<-which(maxint==0)
					  if(length(badfeats)>0){
						curdata<-curdata[-c(badfeats),]
						
						
						
						  feat.eval.result.mat<- feat.eval.result.mat[-c(badfeats),]
						  feateval_list[[i]]<-feat.eval.result.mat
						  peak_score_list[[i]]<-peak_score_list[[i]][-c(badfeats),]
					  }
					
					}
					
					
                    #print(dim(curdata[which(as.numeric(feat.eval.result.mat$median)<=feat.filt.thresh),]))
					data_rpd_all[[i]]<-curdata #[which(as.numeric(feat.eval.result.mat$median)<=feat.filt.thresh),]
				 
				
				 feat.eval.outfile=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
								
				 #write results
                 write.table(data_rpd_all[[i]], feat.eval.outfile,sep="\t", row.names=FALSE)
				 
				
			 }
                          ###########################################
                          #4) Merge two or more parameter settings
                          cat("\n")
                          print("*******Stage 3a: Merging features detected at different parameter settings**********")
							cat("\n")                         
                          num_pairs=1
                          finalres={}
                          rnames={}
                          
			 if(merge.pairwise==TRUE){
				for(i in 1:length(alignmentresults))
				{
                                  file_name=sapply(strsplit(alignmentresults[[i]],".txt"),head)
                                  feat.eval.file=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
                                 
				 

				  a1=sapply(strsplit(as.character(alignmentresults[i]), "\\_pres"), head, n=2)[2]
                                  minpres=sapply(strsplit(as.character(a1), "\\_run"), head, n=2)[1]
                                  minpres=as.numeric(minpres)
                                  a2=sapply(strsplit(as.character(alignmentresults[i]), "\\_run"), head, n=2)[2]
                                  minrun=sapply(strsplit(as.character(a2), "\\_"), head, n=2)[1]
                                  minrun=as.numeric(minrun)
                                  p1=paste(minrun,"_",minpres,sep="")
					    			
	
					for(j in i:length(alignmentresults))
					{
								bool_num<-1
								if(i==j){
									if(length(alignmentresults)>1){
										bool_num<-0
									}
									else{
										bool_num<-1
									}
								}
				    
								if(bool_num==1)
								{
								  file_name=sapply(strsplit(alignmentresults[[j]],".txt"),head)
								  feat.eval.file=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
							      

								  a1=sapply(strsplit(as.character(alignmentresults[j]), "\\_pres"), head, n=2)[2]
								  minpres=sapply(strsplit(as.character(a1), "\\_run"), head, n=2)[1]
								  minpres=as.numeric(minpres)
								  a2=sapply(strsplit(as.character(alignmentresults[j]), "\\_run"), head, n=2)[2]
								  minrun=sapply(strsplit(as.character(a2), "\\_"), head, n=2)[1]
								  minrun=as.numeric(minrun)
								  p2=paste(minrun,"_",minpres,sep="")
							
								  p1_p2=paste("p",i,"_U_","p",j,sep="")
								  
								  
								cnames1<-colnames(feateval_list[[i]])
								
								
								
								 feat.eval.A<-feateval_list[[i]] #cbind(feateval_list[[i]],peak_score_list[[i]])
								 feat.eval.B<-feateval_list[[j]] #cbind(feateval_list[[j]],peak_score_list[j]])
								 
								 feat.eval.A<-as.data.frame(feat.eval.A)
								 feat.eval.B<-as.data.frame(feat.eval.B)
								 

								 
								 feat.eval.A<-feat.eval.A[which(as.numeric(feat.eval.A$median)<=feat.filt.thresh),]
								 feat.eval.B<-feat.eval.B[which(as.numeric(feat.eval.B$median)<=feat.filt.thresh),]
								 
									print(paste("Number of good quality features from setting ",i,":", dim(data_rpd_all[[i]])[1],sep=": "))
									
									if(i!=j){
									print(paste("Number of good quality features from setting ",j,":",dim(data_rpd_all[[j]])[1],sep=": "))
									}
								 
								data_m=merge.Results(data_rpd_all[[i]],data_rpd_all[[j]],feat.eval.A,feat.eval.B,max.mz.diff,max.rt.diff,merge.eval.pvalue,alignment.tool="apLCMS",
								 numnodes=numnodes,mult.test.cor,mergecorthresh,missingvalue)
								
								numcols<-dim(data_m)[2]

								data_m[,5]<-data_m[,5]-1
								 data_m_int<-data_m[,-c(1:5,(numcols-8):numcols)]
								 numsamps<-dim(data_m_int)[2]/num_replicates
								 
								
                                 
                                 if(is.na(minexp)==FALSE){
                                     
                                    
                                     if(length(which(data_m[,5]>=minexp))>0){
                                         data_m<-data_m[which(data_m[,5]>=minexp),]
                                         
                                         print(paste("Number of features after minexp ",minexp," filtering",sep=""))
                                         print(dim(data_m)[1])
                                         
                                     }else{
                                         stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
                                     }
                                     
                                     
                                 }else{
                                     
                                     if(is.na(minexp.pct)==FALSE){
                                         
                                         minexp<-round(minexp.pct*dim(data_m_int)[2])
                                         
                                         if(length(which(data_m[,5]>=minexp))>0){
                                             data_m<-data_m[which(data_m[,5]>=minexp),]
                                         }else{
                                             stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
                                         }
                                         
                                         
                                     }
                                 }
							
									
								
								 numcols<-dim(data_m)[2]

								 data_m_int<-data_m[,-c(1:5,(numcols-8):numcols)]
								 numsamps<-dim(data_m_int)[2]/num_replicates
								 
								 
								
								 
								maxint<-apply(data_m_int,1,function(x){max(x,na.rm=TRUE)})
								
				
									
								
								
								union_list[[num_pairs]]<-data_m[,c(1:5)]
								
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$numgoodsamples)
								
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$median)
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$peak_score)
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$Qscore)
								
								
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],maxint)
								
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m_int)
								
								featinfo<-colnames(data_m[,c(1:4)])
								cnames<-colnames(data_m_int)
								
								merge.res.colnames<-c(featinfo,"NumPres.All.Samples","NumPres.Biological.Samples",paste("median",fileroot,sep=""),"PeakScore","Qscore","Max.Intensity",cnames)
								colnames(union_list[[num_pairs]])<-as.character(merge.res.colnames)
								
								finalname=paste("apLCMSfeaturetable_at_", p1_p2,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
								  
								  #Output merge results
                                   write.table(union_list[[num_pairs]],file=paste(subdir3,"/",finalname,sep=""), sep="\t",row.names=FALSE)
							
								
								  curres={}
								  curres=cbind(curres, p1_p2)
								  curres=cbind(curres, dim(union_list[[num_pairs]])[1])
								  curres=cbind(curres, mean(as.numeric(union_list[[num_pairs]][,7])))
								  curres=cbind(curres, mean(as.numeric(data_m$Qscore)))
								  
								  curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(union_list[[num_pairs]][,9]))))
								  
								  if(curscore>bestscore){
									
									bestscore<-curscore
									best_i<-num_pairs
									best_pair<-p1_p2
								  }
								  
								  
								  curres=cbind(curres,curscore)
								  
								  curres<-as.data.frame(curres)
								  
								  finalres=rbind(finalres,curres)
								  

										 num_pairs=num_pairs+1
								}
					}				  
                          
				}
			
			}else{
			
							data_rpd_all_parameters<-{}
							 feat.eval.all<-{}
							 p1_p2<-"pall"
								for(i in 1:length(alignmentresults))
								{
									 
										
													
									data_rpd_all_parameters<-rbind(data_rpd_all_parameters,data_rpd_all[[i]])
									
									feat.eval.A<-feateval_list[[i]]
									 feat.eval.A<-feat.eval.A[which(as.numeric(feat.eval.A$median)<=feat.filt.thresh),]
								
									 feat.eval.all<-rbind(feat.eval.all,feat.eval.A)
									
									print(dim(data_rpd_all[[i]]))
									print(dim(feateval_list[[i]]))
									
								}
								
								  
								  
								cnames1<-colnames(feateval_list[[i]])
								
								
								
								#feat.eval.all<-unique(feat.eval.all)
								 #data_rpd_all_parameters<-unique(data_rpd_all_parameters)
								 
								 feat.eval.all<-as.data.frame(feat.eval.all)
								 data_rpd_all_parameters<-as.data.frame(data_rpd_all_parameters)
								 

								 
								
								 
						
								 
								data_m=merge.Results(data_rpd_all_parameters,data_rpd_all_parameters, feat.eval.all,feat.eval.all,max.mz.diff,max.rt.diff,merge.eval.pvalue,alignment.tool="apLCMS",
								 numnodes=numnodes,mult.test.cor,mergecorthresh,missingvalue)
								
								numcols<-dim(data_m)[2]

								data_m[,5]<-data_m[,5]-1
								 data_m_int<-data_m[,-c(1:5,(numcols-8):numcols)]
								 numsamps<-dim(data_m_int)[2]/num_replicates
								 
								 if(is.na(minexp.pct)==FALSE){
								 
									minexp<-round(minexp.pct*dim(data_m_int)[2])
								
								if(length(which(data_m[,5]>=minexp))>0){
									data_m<-data_m[which(data_m[,5]>=minexp),]
								}else{
									stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
									}
								
								
								}
							
									
								
								 numcols<-dim(data_m)[2]

								 data_m_int<-data_m[,-c(1:5,(numcols-8):numcols)]
								 numsamps<-dim(data_m_int)[2]/num_replicates
								 
								 
								
								 
								maxint<-apply(data_m_int,1,function(x){max(x,na.rm=TRUE)})
								
				
									
								
								
								union_list[[num_pairs]]<-data_m[,c(1:5)]
								
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$numgoodsamples)
								
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$median)
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$peak_score)
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$Qscore)
								
								
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],maxint)
								
								union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m_int)
								
								featinfo<-colnames(data_m[,c(1:4)])
								cnames<-colnames(data_m_int)
								
								merge.res.colnames<-c(featinfo,"NumPres.All.Samples","NumPres.Biological.Samples",paste("median",fileroot,sep=""),"PeakScore","Qscore","Max.Intensity",cnames)
								colnames(union_list[[num_pairs]])<-as.character(merge.res.colnames)
								
								finalname=paste("apLCMSfeaturetable_at_", p1_p2,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
								  
								  #Output merge results
                                  write.table(union_list[[num_pairs]],file=paste(subdir3,"/",finalname,sep=""), sep="\t",row.names=FALSE)
							
								
								  curres={}
								   curres=cbind(curres, p1_p2)
								  curres=cbind(curres, dim(union_list[[num_pairs]])[1])
								  curres=cbind(curres, mean(as.numeric(union_list[[num_pairs]][,7])))
								  curres=cbind(curres, mean(as.numeric(data_m$Qscore)))
								  
								  curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(union_list[[num_pairs]][,9]))))
								  
								  if(curscore>bestscore){
									
									bestscore<-curscore
									best_i<-num_pairs
									best_pair<-p1_p2
								  }
								  

								  
								  curres=cbind(curres,curscore)
								  
								  curres<-as.data.frame(curres)
								  
								  finalres=rbind(finalres,curres)
								  
								  num_pairs<-num_pairs+1
			
			}
                          finalres<-as.data.frame(finalres)
                         


            if(merge.pairwise==TRUE){
			 
			  colnames(finalres)<-c("Parameter Combination", "Number of Features", "median PID/CV between sample replicates","mean Qscore (Quality score)", "Parameter score")
              write.table(finalres,file=paste(xMSanalyzer.outloc,"/apLCMS_with_xMSanalyzer_merge_summary.txt",sep=""), sep="\t", row.names=FALSE)
		          
		          }
		          


    
    print("Most optimal feature setting:")
    print(best_pair)
      cat("\n")
    #########################################################################
    
    
    cat("\n")
      print(paste("********Stage 3b: Generating final (pre-batcheffect correction) untargeted and targeted feature tables using ",best_pair," results******",sep=""))
      cat("\n")
      #pdf(file=paste("subdir4a/Stage4a_QE_plots.pdf",sep=""))
      
      
 	d1<-union_list[[best_i]]
    d1<-unique(d1)

	#d1<-d1[1:50,]	
	num_samples<-dim(d1)[2]-10
	uniq_samp_ind<-seq(1,num_samples,num_replicates)
	num_EIC<-min(30,num_samples)
	num_EIC<-(num_EIC)*(1/num_replicates)
 	samp_ind<-sample(uniq_samp_ind,size=num_EIC)
	samp_ind<-samp_ind[order(samp_ind)]
	samp_ind1<-samp_ind+1
	samp_ind2<-{}
	if(num_replicates>=3){
		samp_ind2<-samp_ind+2
	}
	samp_ind<-c(samp_ind,samp_ind1,samp_ind2)
	samp_ind<-samp_ind[order(samp_ind)]

temp_qres<-{}
	#if(FALSE)
	
	if(plotEICs=="target" | plotEICs=="all")
	{	
	pdfname<-paste(subdir4a,"/EIC_postfiltering",".pdf",sep="")
	pdf(pdfname)
	#par(mfrow=c(3,3))
	
	plotEICs_parent<-plotEICs
	
	for(i in 1:length(min.run.list)){
	
		
		plotEICs<-plotEICs_parent
		
		presval=min.pres.list[i]
		runval=min.run.list[i]

		setwd(apLCMS.outloc)
		#minexp=208
                fname<-paste(apLCMS.outloc, "/apLCMS",apLCMSmode,"_aligned", "_pres",  presval, "_run", runval,"_", minexp, "exppostrecovery.Rda", sep="")

		load_fname<-try(load(fname),silent=TRUE)
		
		if(is(load_fname, "try-error")){
		
			plotEICs=="none"
		}
		
		
		
	if(plotEICs=="all")
	{
		Name<-paste("mz",seq(1,dim(d1)[1]),sep="")
	        stddata_temp<-cbind(d1[,c(1:2)],Name)
		stddata_temp<-as.data.frame(stddata)
			
		stddata_temp$mz<-as.numeric(as.character(stddata_temp$mz))
		stddata_temp$time<-as.numeric(as.character(stddata_temp$time))
			
			
		

	
		

		finalfeatmat=aligned$final.ftrs
		finalfeatmat<-as.data.frame(finalfeatmat)

		 rows_ind=seq(1,dim(finalfeatmat)[1])
		files <- colnames(aligned$final.ftrs)
	
		if(length(files)<1){
			aligned$final.ftrs<-aligned$aligned.ftrs
			files <- colnames(aligned$final.ftrs)
			aligned$aligned.ftrs<-NULL
			finalfeatmat=aligned$final.ftrs
			finalfeatmat<-as.data.frame(finalfeatmat)
		}
	

#	apLCMS.EIC.plot(aligned, rows = rows_ind, colors = NA, transform = "none",
  #  subset = NA, mz.list=NA,time.list=NA,minrt=NA, maxrt=NA, min.run=runval,min.pres=presval, max.spline.time.points = 1000,
 #   chem.names=NA,rawprofileloc=cdfloc,cvvec=NA,plotEIC=FALSE,numnodes=NA)
    
    

		stddata_temp<-as.data.frame(stddata_temp)
		
		stddata_temp$mz<-as.numeric(as.character(stddata_temp$mz))
		stddata_temp$time<-as.numeric(as.character(stddata_temp$time))

             try(apLCMS.target.EIC(aligned=aligned,finalfeatmat=aligned$aligned.ftrs[,c(1:2)],refMZ.mz.diff=1,refMZ.time.diff=1,
	     target.mz.list = stddata_temp,subs=samp_ind,colors=NA,apLCMS.outloc,runval=min.run.list[i],presval=min.pres.list[i],
	     minexp,cdfloc=rawprofileloc,cvvec=d1[,7]),silent=TRUE) 

		}else{
		
		if(plotEICs=="target"){

		stddata_temp<-stddata

		stddata_temp<-as.data.frame(stddata_temp)
		
		stddata_temp$mz<-as.numeric(as.character(stddata_temp$mz))
		#stddata_temp$time<-as.numeric(as.character(stddata_temp$time))

             try(apLCMS.target.EIC(aligned=aligned,finalfeatmat=aligned$aligned.ftrs[,c(1:2)],refMZ.mz.diff=refMZ.mz.diff,refMZ.time.diff=refMZ.time.diff,
	     target.mz.list = stddata,subs=samp_ind,colors=NA,apLCMS.outloc,runval=min.run.list[i],presval=min.pres.list[i],
	     minexp,cdfloc=rawprofileloc),silent=TRUE) 
	       }
	       
	       
		}

	}

	dev.off()		
		
	
	plotEICs<-plotEICs_parent
}
	
		cnames<-colnames(d1)
		cnames[10]<-"Max.Intensity"
		colnames(d1)<-cnames
		d1<-as.data.frame(d1)
		d1<-unique(d1)
	
		if(is.na(peak.score.thresh)==FALSE){	
		

			d1<-d1[which(d1$PeakScore>peak.score.thresh),]
		}
        #print(d1[1:3,1:15])
		
        #pdf(file=paste(subdir4a,"/Stage4a_QE_plots.pdf",sep=""))
    
    pdf(file=paste(subdir4a,"/QE_plots_All_samples_RAW.pdf",sep=""))
	time_thresh<-NA
	
	if(is.na(void.vol.timethresh)==FALSE){
        
        
        #dfirst15<-d1[which(d1[,2]<void.vol.timethresh),]
        
        dfirst15<-d1[which(d1$time<void.vol.timethresh),]
        
        
        
        if(nrow(dfirst15)>1){
            
            
            ind1<-which(dfirst15$Max.Intensity==max(dfirst15$Max.Intensity))[1]
            
            time_thresh<-dfirst15$time[ind1]
            
            time_thresh<-time_thresh-(0.30*time_thresh)
            
            time_thresh<-round(time_thresh,1)
            
            plot(dfirst15$time,dfirst15$Max.Intensity,xlab="Time (s)", col="brown",ylab="Max intensity \nacross all samples", main=paste("Estimated void volume time: ",time_thresh," s",sep=""))
            abline(v=time_thresh,col=4,lty=3)
            
            
            
            d1<-d1[which(d1$time>=time_thresh),]
            
            print("Estimated void volume time cutoff")
            print(time_thresh)
            
            
        }else{

			print("No features eluting before void volume time cutoff")
		}

			
		
		#sfname<-paste(xMSanalyzer.outloc,"/stage5/feature_table_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtered.txt",sep="")
		finalname<-paste("featuretable_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,"_voidtimefilt.txt",sep="")
		
        #write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
	   

	}else{
		
		
		finalname<-paste("featuretable_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
		
        #write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
	   

		
		}

d1_int<-round(d1[,-c(1:10)],0)

if(num_replicates>2){
rsqres_list<-evaluate.Samples(d1_int, num_replicates, alignment.tool=NA, cormethod,missingvalue,ignore.missing)

rsqres<-as.data.frame(rsqres_list$cor.matrix)

rsqres<-as.data.frame(rsqres)
snames<-colnames(d1_int)
snames_1<-snames[seq(1,length(snames),num_replicates)]
rownames(rsqres)<-snames_1
pcor_outfile=paste(subdir4a,"/Pairwise_Pearson_correlation_technical_replicates.txt",sep="")

write.table(rsqres, pcor_outfile,sep="\t",row.names=TRUE)
}

#if(is.na(sample_info_file)==FALSE)
if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
{
    #sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)
    
    sampleid_mapping[,1]<-gsub(sampleid_mapping[,1],pattern=filepattern,replacement="",ignore.case=TRUE)
    cnames<-colnames(data_m)
    
    sampleid_mapping[,2]<-tolower(sampleid_mapping[,2])
    
    qc_label<-tolower(qc_label)
    
    cnames<-gsub(cnames,pattern=filepattern,replacement="",ignore.case=TRUE)
    colnames(data_m)<-cnames
    
    class_labels<-as.factor(sampleid_mapping[,2])
    
    if(is.na(qc_label)==FALSE){
    
        qc_label<-tolower(qc_label)
    qc_label_check<-gregexpr(pattern=qc_label,text=sampleid_mapping[,2])
    
    qc_index<-which(qc_label_check>0)
    
    if(length(qc_index)>0){
  
    print("QC file index")
    print(qc_index)
    
	d1_int_qc<-d1_int[,qc_index]

    
    rsqres_listQC<-evaluate.Samples(d1_int_qc, numreplicates=length(qc_index), alignment.tool=NA, cormethod,missingvalue,ignore.missing)
   
    rsqres_listQC$cor.matrix<-round(rsqres_listQC$cor.matrix,2)
    
    r1<-matrix(0,nrow=length(qc_index),ncol=length(qc_index))
    
    #r1[upper.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
    diag(r1)<-1
    #r1[lower.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
    
    start_count=0
    for(r1_row in 1:(length(qc_index)-1))
    {
        
          
            r1[r1_row,c((r1_row+1):length(qc_index))]<-rsqres_listQC$cor.matrix[,(start_count+1):(start_count+length(qc_index)-r1_row)]
           
           
            start_count=(start_count+length(qc_index)-r1_row) #(r1_row*length(qc_index))-r1_row-1 #(length(qc_index)-1)


    }
    
    r1[lower.tri(r1)]<-NA #r1[upper.tri(r1,diag=TRUE)] #[-length(rsqres_listQC$cor.matrix)]
    
    
   
   TIC <- vector("list",length(qc_index))
   
   for(qc_i in 1:length(qc_index)){
       
       
       #plot( d1$time, d1_int_qc[,qc_i],xlab="Retention Time",ylab="TIC",main=mainlab1)
       TIC[[qc_i]] <-cbind(d1$time, d1_int_qc[,qc_i])
   }
    }else{
        
        print("No matches found for the QC samples.")
        qc_label=NA
    }
   
    }
   
  
}



if(num_replicates>2)
{
    max_rep_mean_cor<-apply(rsqres,1,max)
    
    max_rep_mean_cor<-as.data.frame(max_rep_mean_cor)
}else{

	if(num_replicates==2){
	max_rep_mean_cor<-rsqres
	}else{
	max_rep_mean_cor<-1
	}
    
}

d1<-cbind(d1[,c(1:10)],d1_int)
rm(d1_int)


 cat("\n")
   print(paste("********Stage 4a: Performing QC checks using ",best_pair," results******",sep=""))
         
      cat("\n")
      Sys.sleep(1)
	
    
		
		print("Dim data after void time filtering")
		print(dim(d1))
        
        max_ylim<-nrow(d1)+100
        
        num_features<-nrow(d1)
        
        
        par(mfrow=c(1,2))
        if(num_replicates>2){
            
            rep_cv<-d1$median_CV
            
            rep_cv_ylab<-"CV"
            
            h1<-hist(d1$median_CV,breaks=seq(0,max(d1$median_CV,na.rm=TRUE)+10,10),main=paste("Histogram median CV \n (using all ",num_features," features)",sep=""),col="brown",xlab="median CV%", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
            lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
        }else{
            
	    if(num_replicates==2){
            rep_cv<-d1$median_PID
            
            rep_cv_ylab<-"PID"
            
            h1<-hist(d1$median_PID,breaks=seq(0,max(d1$median_PID,na.rm=TRUE)+10,10),main=paste("Histogram median PID \n (using all ",num_features," features)",sep=""),col="brown",xlab="median PID%", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
            lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
	    
	    }else{
		h1<-{}
		rep_cv<-{}
	    }
	    
        }
        #lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
        #pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main="Pie chart of median CVs (%) \n using all features",cex.main=0.7)
        
        if(num_replicates>2){
            pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main=paste("Pie chart of median CVs (%) \n using all ",num_features," features\n; average=",round(mean(d1$median_CV),2),sep=""),cex.main=0.7)
        }else{
	
		if(num_replicates==2){
		pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main=paste("Pie chart of median PIDs (%) \n using all ",num_features," features\n; average=",round(mean(d1$median_PID),2),sep=""),cex.main=0.7)
            }
        }
        
        par(mfrow=c(1,1))

hist(d1$NumPres.All.Samples,main="Histogram NumPres.All.Samples",col="brown",xlab="Number of samples (including replicates) \n with non-zero intensity values", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
			
			hist(d1$NumPres.Biological.Samples,main="Histogram NumPres.Biological.Samples",col="brown",xlab="Number of biological samples \n with non-zero intensity values", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
			#h1<-hist(d1$median_CV,seq(0,max(d1$median_CV,na.rm=TRUE),10),main="Histogram median CV \n (using all data)",col="brown",xlab="median CV%", ylab="Number of features",cex.main=0.7)
			
			#h1<-hist(d1$median_CV,breaks=seq(0,max(d1$median_CV,na.rm=TRUE)+10,10),main="Histogram median CV \n (using all data)",col="brown",xlab="median CV%", ylab="Number of features",cex.main=0.7)
			
			par(mfrow=c(1,1))
			
            if(num_replicates>1){
                
                min_rep_pcor<-round(min(max_rep_mean_cor[,1],na.rm=TRUE),2)
                max_rep_pcor<-round(max(max_rep_mean_cor[,1],na.rm=TRUE),2)
                
                if(is.na(samp.filt.thresh)==FALSE){
                    
                
                parent_bad_list<-which(max_rep_mean_cor[,1]<samp.filt.thresh)
                
                if(length(parent_bad_list)>0){
                    
                    if(is.na(sample_info_file)==FALSE && sample_info_file!="NA"){
                    
                    filt_samp_names<-sampleid_mapping[parent_bad_list,1]
                    filt_samp_names<-paste(filt_samp_names,collapse=";")
                    
                    filt_samp_names<-length(parent_bad_list)
                    }
                }else{
                    filt_samp_names<-"None"
                    
                }
                }else{
                    filt_samp_names<-"None"
                    
                }
                
                
                hist(max_rep_mean_cor[,1],breaks=seq(0,1,0.1),main=paste("Histogram for mean Pearson correlation (min: ",min_rep_pcor,"; max: ",max_rep_pcor,") \n within technical replicates after filtering \n # of files filtered at threshold ",samp.filt.thresh,": ",filt_samp_names,sep=""),col="brown",xlab="mean replicate Correlation", ylab="Number of samples",cex.main=0.7)
                
                
            }

				
			
            par(mfrow=c(2,2))
            #5ppm
            data_a<-find.Unique.mzs.sameset(dataA=d1[,c("mz","time")],dataB=d1[,c("mz","time")],mz.thresh=5,time.thresh=NA,alignment.tool=NA)
           num_unique_features<-nrow(data_a$uniqueA)
			total_features<-nrow(d1)
            
            xlab2<-paste("overlapping m/z features \n based on +/- ",5,"ppm overlap criteria",sep="")
            
            temp_df<-cbind(total_features,num_unique_features)
            
            pie_v1<-round(temp_df[1,2]/temp_df[1,1],2)
            pie_v2<-1-pie_v1
            
            temp_dfA<-cbind(temp_df[1,1],temp_df[1,2], (temp_df[1,1]-temp_df[1,2]))#cbind(100*pie_v1,100*pie_v2)
            
            colnames(temp_dfA)<-c("Total","Unique","Overlapping")
            
            #barplot(temp_dfA,col="brown",ylab="Percentage (%) of total features", main=paste("Overlapping vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,100))
            
             barplot(temp_dfA,col="brown",ylab="Number of features", main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,temp_df[1,1]+100))
            #non-overlapping m/z features")
            
            # pie(c(pie_v2,pie_v1),labels=c("Overlapping","Unique"),col=c("orange","brown"),main=paste("Total vs ",xlab2,sep=""),cex.main=0.8)
            
            pie_label1<-paste("Overlapping \n(",100*pie_v2,"%)",sep="")
            pie_label2<-paste("Unique \n(",100*pie_v1,"%)",sep="")
            
            pie(c(pie_v2,pie_v1),labels=c(pie_label1,pie_label2),col=c("orange","brown"),main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8)
            
            #10ppm
            data_a<-find.Unique.mzs.sameset(dataA=d1[,c("mz","time")],dataB=d1[,c("mz","time")],mz.thresh=10,time.thresh=NA,alignment.tool=NA)
            num_unique_features<-nrow(data_a$uniqueA)
            total_features<-nrow(d1)
            
            #xlab2<-paste("Unique (non-overlapping) m/z features \n based on +/- ",10,"ppm overlap criteria",sep="")
            
            xlab2<-paste("overlapping m/z features \n based on +/- ",10,"ppm overlap criteria",sep="")
           
           
            temp_df<-cbind(total_features,num_unique_features)
            #colnames(temp_df)<-c("Total","Unique (non-overlapping)")
            
            pie_v1<-round(temp_df[1,2]/temp_df[1,1],2)
            pie_v2<-1-pie_v1
            
            temp_dfB<-cbind(temp_df[1,1],temp_df[1,2], (temp_df[1,1]-temp_df[1,2])) #cbind(100*pie_v1,100*pie_v2)
            
            colnames(temp_dfB)<-c("Total","Unique","Overlapping")
            
            #barplot(temp_dfB,col="brown",ylab="Percentage (%) of total features", main=paste("Total vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,100)) #non-overlapping m/z features")
            
            barplot(temp_dfB,col="brown",ylab="Number of features", main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,temp_df[1,1]+100))
            
        
        pie_label1<-paste("Overlapping \n(",100*pie_v2,"%)",sep="")
        pie_label2<-paste("Unique \n(",100*pie_v1,"%)",sep="")
            pie(c(pie_v2,pie_v1),labels=c(pie_label1,pie_label2),col=c("orange","brown"),main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8)
            
            par(mfrow=c(1,1))
            
            hist(d1$Qscore,main="Histogram Qscore",col="brown",xlab="Quality score", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
            
            num_features<-nrow(d1)

    #pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main="Pie chart of median CVs (%) \n using all features",cex.main=0.7)
			
			
            d2<-d1 #[order(d1$time),]
			d2<-as.data.frame(d2)
			#d2<-apply(d2,1,as.numeric)
			
            plot(d2$time,d2$mz,main="m/z vs Time",col="brown",xlab="Time (s)",ylab="m/z",cex.main=0.7)
            
            
            par(mfrow=c(1,1))
            
            plot(d2$mz,d2$Max.Intensity,main="Intensity vs m/z",col="brown",xlab="m/z",ylab="Max intensity \nacross all samples/profiles",cex.main=0.7)
            plot(d2$time,d2$Max.Intensity,main="Intensity vs time",col="brown",xlab="Time (s)",ylab="Max intensity \nacross all samples/profiles",cex.main=0.7)
            
	    if(num_replicates>1){
            plot(rep_cv,d2$Max.Intensity,main=paste("Intensity vs ",rep_cv_ylab,sep=""),col="brown",xlab=rep_cv_ylab,ylab="Max intensity \nacross all samples/profiles",cex.main=0.7)
            }
           
            par(mfrow=c(1,1))
            
            plot(d2$mz,d2$Qscore,main="Qscore vs m/z",col="brown",xlab="m/z",ylab="Qscore",cex.main=0.7)
            plot(d2$time,d2$Qscore,main="Qscore vs Time",col="brown",xlab="Time (s)",ylab="Qscore",cex.main=0.7)
            
	     if(num_replicates>1){
	    plot(rep_cv,d2$Qscore,main=paste("Qscore vs ",rep_cv_ylab,sep=""),col="brown",xlab=rep_cv_ylab,ylab="Qscore",cex.main=0.7)
            }
            
            par(mfrow=c(1,1))
		
			
            targeted_feat_raw<-eval.target.mz(dataA=d1[,-c(3:10)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="raw",xMSanalyzer.outloc=subdir4a)
            
            
            median_error<-median(targeted_feat_raw$delta_ppm_error,na.rm=TRUE)
	    
	     median_time_error<-median(targeted_feat_raw$delta_time_error,na.rm=TRUE)
	    
	    raw_mz_time<-d1[,c("mz","time")]
	    
        
	    write.table(raw_mz_time,file=paste(subdir3b,"/Precalibration_mz_time.txt",sep=""),sep="\t",row.names=FALSE)
            
	    save(d1,targeted_feat_raw,median_error,stddata,feat.eval.result,refMZ.mz.diff,refMZ.time.diff,file="test.Rda")
	    
	    #here1
            d1<-get_calibrated_mz_data(dataA=d1,refMZ=stddata,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,feature.eval.result=feat.eval.result,calibration.method=calibration.method[1])

		#delta_ppm_error=median_error,delta_time_error=median_time_error,r
            
            finalname<-paste("RAW_mzcalibrated_untargeted_featuretable.txt",sep="")
            
    #        write.table(d1,file=finalname,sep="\t",row.names=FALSE)
    
            write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
    
    #write.table(d1,file=paste(xMSanalyzer.outloc,finalname,sep=""),sep="\t",row.names=FALSE)

    	    targeted_feat_calibratedmz<-eval.target.calibratedmz(dataA=d1[,-c(3:10)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="RAW",xMSanalyzer.outloc=subdir4a)
            
            fnamesT<-paste(subdir3b,"/RAW_mzcalibrated_targeted_featuretable_mz",refMZ.mz.diff,"ppm_time",refMZ.time.diff,".txt",sep="")
            write.table(targeted_feat_calibratedmz$targetdata,file=fnamesT,sep="\t",row.names=FALSE)
            
	      
	data_m<-d1[,-c(1:10)]
   
    if(replacezeroswithNA==TRUE){
	data_m<-replace(as.matrix(data_m),which(data_m==0),NA)
	d1<-cbind(d1[,c(1:10)],data_m)
	
	}
	
	
	counna<-apply(data_m,1,function(x){length(which(is.na(x)==TRUE))})

	maxzeros<-1*dim(data_m)[2]
	
	if(length(which(counna<maxzeros))){
		data_m<-data_m[which(counna<maxzeros),]
	}

    maxint<-apply(data_m,1,function(x){max(x,na.rm=TRUE)})
    
   	maxint_ord<-order(maxint,decreasing=TRUE)
   	#[maxint_ord[1:5000]
    X<-t(data_m) #[maxint_ord[1:2000],])
    X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
    
    tic.eval(d1[,-c(1:10)],outloc=subdir4a)
    
    #save(d1,file="d1.Rda")
    feat.eval.result=evaluate.Features(d1[,-c(5:10)], numreplicates=num_replicates,min.samp.percent=min.samp.percent,
    alignment.tool="apLCMS",impute.bool=impute.bool,peak_scores=d1$PeakScore,numnodes=numnodes)
    
    
										cnames=colnames(feat.eval.result)

				
    Sys.sleep(1)
			




    Sys.sleep(1)
    
s1<-"Stage 1 results: QC evaluation of invidual parameters from apLCMS"
s2<-"Stage 2 results: filtered results from each paramter setting based on sample and feature quality (CV within replicates) checks"
s3<-"Stage 3a results: merged results using stage 2 filtered files"
s3b<-"Stage 3b results: RAW m/z calibrated untargeted and targeted feature tables (averaged and non-averaged)"
s4a<-"Stage 4a results: QC evaluation of targeted and untargeted data before batch-effect correction"

sm<-rbind(s1,s2,s3,s3b,s4a)

qc_res_mat<-{}
rnames_qc_results<-{}
    #if(is.na(sample_info_file)==FALSE)
     if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
    {
    	#sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)
    	
    	sampleid_mapping[,1]<-gsub(sampleid_mapping[,1],pattern=filepattern,replacement="",ignore.case=TRUE)
	
	cnames<-colnames(d1)
		
		cnames<-gsub(cnames,pattern=filepattern,replacement="",ignore.case=TRUE)
		colnames(d1)<-cnames
		
		
		cnames<-colnames(data_m)
		
		cnames<-gsub(cnames,pattern=filepattern,replacement="",ignore.case=TRUE)
		colnames(data_m)<-cnames

	
		
        class_labels<-as.factor(sampleid_mapping[,2])
    	batch_inf<-sampleid_mapping[,3]
    	
    	batch_labels<-as.factor(sampleid_mapping[,3])
    
    	l1<-levels(batch_labels)

	batch_levels<-levels(batch_labels)
		    try(pca.eval(X=X,samplelabels=batch_labels,filename="raw",ncomp=5,center=TRUE,scale=TRUE,legendlocation="bottomleft",legendcex=0.5,outloc=subdir4a),silent=TRUE)
			 
			Sys.sleep(1)
			
			

			
			
			   		    
    		if(dim(sampleid_mapping)[2]<4){
   		    mod<-rep(1,dim(data_m)[2])
   		    }else{
   		    		mod<-sampleid_mapping[,-c(1:3)]
   		    	}
    

  	
    		dev.off()
    		try(dev.off(),silent=TRUE)





    		Sys.sleep(2)
            
            
            if(is.na(qc_label)==FALSE){
            pdfname="QE_plots_QC_samples_RAW.pdf"
            #save(TIC,file="TIC.Rda")
            
            
            pdf(pdfname,w=12,h=12) #,w=16,h=10)
            
            tic<-apply(d1_int_qc,2,function(x){
                x<-replace(x,which(x==0),NA)
                return(sum(x,na.rm=TRUE))
            })
            
            
            mean_tic<-mean(tic)
            
            cv_tic<-100*sd(tic,na.rm=TRUE)/mean(tic,na.rm=TRUE)
            
            
            
            main_lab<-paste("Total ion intensity using all features\n for QC samples\nmean=",round(mean_tic,2),"\n%CV total ion intensity=",round(cv_tic,2),sep="")
            
            
            #tiff("barplot_TIC_using_all_features.tiff",width=2000,height=2000,res=300)
            #pdf("TIC_all_features.pdf")
            par(mfrow=c(1,1))
            barplot(tic,cex.names=0.4,cex.axis=1,main=main_lab,col="brown",cex.main=0.6)
            

            N=length(qc_index)
            cols <- rainbow(length(qc_index))
            lty = 1:N
            pch = rep(15,N)
            xlim = range(sapply(TIC, function(x) range(x[,1])))
            ylim = range(sapply(TIC, function(x) range(x[,2])))
            plot(0,0, type="n", xlim = xlim, ylim = ylim, main = "Feature intensity vs Time (QC samples)", xlab = "Retention Time", ylab = "Intensity")
            for (i in 1:N) {
                
                tic <- TIC[[i]]
                tic<-tic[order(tic[,1]),]
                
                points(tic[,1], (tic[,2]), col = cols[i], pch = pch[i], type="l")
            }
            legend("topleft",paste(class_labels[qc_index]), col = cols, lty = lty, pch = pch,cex=0.3)
            
            
            
            heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
            heatmap_cols<-rev(heatmap_cols)
            
            mean_all_rep_pcor<-rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)]
            mean_all_rep_cv<-round(cv_tic,2)
            
            qc_res_mat<-rbind(qc_res_mat,cbind(mean_all_rep_pcor,mean_all_rep_cv))
            
            rnames_qc_results<-c(rnames_qc_results,"RAW (all replicates)")
            
            cor_mainlab<-paste("Pairwise Pearson correlation \n within QC samples (all replicates)\n mean: ",rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)],"; range: ",min(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)])," to ",max(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)]),sep="")
              heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
            heatmap_cols<-rev(heatmap_cols)
            
            h73<-suppressWarnings(heatmap.2(r1, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8))
            
            #h73<-heatmap.2(r1, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8)
            
            
            d1_qc<-cbind(d1[,c(1:10)],d1_int_qc)
            
            targeted_feat_calibratedmzQC<-eval.target.calibratedmz(dataA=d1_qc[,-c(3:10)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="rawQC",xMSanalyzer.outloc=subdir4a)
            try(dev.off(),silent=TRUE)
            
            }
            
            if(summarize.replicates==TRUE){
                
                Ymat=rep(1,dim(d1[,-c(1:10)])[1])
                
                
                if(num_replicates>1){
                    
                    if(data.norm.pipeline=="AC"){
			minexp=minexp/num_replicates
                    sampleid_mapping<-sampleid_mapping[seq(1,dim(d1[,-c(1:10)])[2],num_replicates),]
                    
                    
                    batch_inf<-batch_inf[seq(1,dim(d1[,-c(1:10)])[2],num_replicates)]
                    
                    mod<-mod[seq(1,dim(d1[,-c(1:10)])[2],num_replicates)]
                    
		    batch_labels<-batch_labels[seq(1,dim(d1[,-c(1:10)])[2],num_replicates)]
		    
		    l1<-levels(batch_labels)
    
                    data_matrix<-data_summarize(Xmat=d1[,- c(3:10)],Ymat=NA,feature_table_file=NA,parentoutput_dir=subdir3b,class_labels_file=NA,num_replicates=num_replicates,summarize.replicates=summarize.replicates,
                    summary.method=summary.method,missing.val=missingvalue, rep.num.max.missing.thresh=rep.num.max.missing.thresh,summary.na.replacement=summary.na.replacement)
                    
                    check_bad_feats_allzeros<-apply(data_matrix[,-c(1:2)],1,sum)
		    bad_feats_allzeros_index<-which(check_bad_feats_allzeros==0)
		    
		    if(length(bad_feats_allzeros_index)>0){
		    
			data_matrix<-data_matrix[-bad_feats_allzeros_index,]
			
			d1<-d1[-bad_feats_allzeros_index,]
		    }
		    
                    d1<-cbind(d1[,c(1:10)],data_matrix[,-c(1:2)])
                    
                    d1_int<-d1[,-c(1:10)]
		    
		    
		    
                    X<-t(data_matrix[,-c(1:2)])
                    
                    data_m<-data_matrix[,-c(1:2)]
                   
                   
                    num_replicates=1
                    
                    sampleid_mapping[,1]<-gsub(sampleid_mapping[,1],pattern=filepattern,replacement="",ignore.case=TRUE)
                    cnames<-colnames(data_m)
                    
                    cnames<-gsub(cnames,pattern=filepattern,replacement="",ignore.case=TRUE)
                    colnames(data_m)<-cnames
                    
                    class_labels<-as.factor(sampleid_mapping[,2])
		    
		     pdfname=paste(subdir4a,"/QE_plots_All_samples_RAW_replicate ",summary.method," summarization.pdf",sep="")
                        #pdf(pdfname)
                        pdf(pdfname,w=12,h=12) #,w=16,h=10)
                        
                        d1_int_qc<-d1_int
                        all_index<-seq(1,dim(d1_int)[2])
			
                        TIC <- vector("list",length(all_index))
                        
                        for(qc_i in 1:length(all_index)){
                            
                            
                            #plot( d1$time, d1_int_qc[,qc_i],xlab="Retention Time",ylab="TIC",main=mainlab1)
                            TIC[[qc_i]] <-cbind(d1$time, d1_int_qc[,qc_i])
                        }
                        
                        N=length(all_index)
                        cols <- rainbow(length(all_index))
                        lty = 1:N
                        pch = rep(15,N)
                        xlim = range(sapply(TIC, function(x) range(x[,1])))
                        ylim = range(sapply(TIC, function(x) range(x[,2])))
                        plot(0,0, type="n", xlim = xlim, ylim = ylim, main = "Feature intensity vs Time \n(all samples replicates summarized)", xlab = "Retention Time", ylab = "Intensity")
                        for (i in 1:N) {
                            
                            tic <- TIC[[i]]
                            tic<-tic[order(tic[,1]),]
                            
                            points(tic[,1], (tic[,2]), col = cols[i], pch = pch[i], type="l")
                        }
                        legend("topleft",paste(class_labels[all_index]), col = cols, lty = lty, pch = pch,cex=0.3)
                        
            
            
                        
                        
                        tic<-apply(d1_int_qc,2,function(x){
                            x<-replace(x,which(x==0),NA)
                            return(sum(x,na.rm=TRUE))
                        })
                        
                        
                        mean_tic<-mean(tic)
                        
                        cv_tic<-100*sd(tic,na.rm=TRUE)/mean(tic,na.rm=TRUE)
                        
                        
                        
                       
                        main_lab<-paste("Total ion intensity using all features\n for all samples (replicate ",summary.method," summarization)\nmean=",round(mean_tic,2),"\n%CV total ion intensity=",round(cv_tic,2),sep="")
                        
                        #tiff("barplot_TIC_using_all_features.tiff",width=2000,height=2000,res=300)
                        #pdf("TIC_all_features.pdf")
                        par(mfrow=c(1,1))
                        barplot(tic,cex.names=0.4,cex.axis=1,main=main_lab,col="brown",cex.main=0.6)
                        
                            
                        rsqres_listQC<-evaluate.Samples(d1_int_qc, numreplicates=length(all_index), alignment.tool=NA, cormethod,missingvalue,ignore.missing)
                        
                        rsqres_listQC$cor.matrix<-round(rsqres_listQC$cor.matrix,2)
                        
                        mean_avg_raw_pcor<-rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)]
                        mean_avg_raw_cv<-round(cv_tic,2)
                        
                       # qc_res_mat<-rbind(qc_res_mat,cbind(mean_avg_raw_pcor,mean_avg_raw_cv))
                        
                      #  rnames_qc_results<-c(rnames_qc_results,"RAW post-replicate summarization")
                        
                        r2<-matrix(0,nrow=length(all_index),ncol=length(all_index))
                        
                        #r1[upper.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
                        diag(r2)<-1
                        #r1[lower.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
                        
                        start_count=0
                        for(r1_row in 1:(length(all_index)-1))
                        {
                            
                            
                            r2[r1_row,c((r1_row+1):length(all_index))]<-rsqres_listQC$cor.matrix[,(start_count+1):(start_count+length(all_index)-r1_row)]
                            
                            
                            start_count=(start_count+length(all_index)-r1_row) #(r1_row*length(all_index))-r1_row-1 #(length(all_index)-1)
                            
                            
                        }
			
			  heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
			heatmap_cols<-rev(heatmap_cols)
            
                        
                        r2[lower.tri(r2)]<-NA #r1[upper.tri(r1,diag=TRUE)] #[-length(rsqres_listQC$cor.matrix)]
                        
                        cor_mainlab<-paste("Pairwise Pearson correlation \n within all samples (post replicate ",summary.method," summarization) \n mean: ",rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)],"; range: ",min(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)])," to ",max(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)]),sep="")
                        
                        h73<-suppressWarnings(heatmap.2(r2, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8))
                        
                        #h73<-heatmap.2(r2, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8)
                        
                        
                        d1_qc<-cbind(d1[,c(1:10)],d1_int_qc)
                        
                        targeted_feat_calibratedmzAvg<-eval.target.calibratedmz(dataA=d1_qc[,-c(3:10)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="rawQCavg",xMSanalyzer.outloc=subdir4a)
                        try(dev.off(),silent=TRUE)
			
			rm(d1_qc)
			rm(d1_int_qc)
                    
                    if(is.na(qc_label)==FALSE){
			
				    qc_label<-tolower(qc_label)
                    qc_label_check<-gregexpr(pattern=qc_label,text=sampleid_mapping[,2])
                    
                    qc_index<-which(qc_label_check>0)
                    
                    if(length(qc_index)>0){
                        
                        pdfname=paste(subdir4a,"/QE_plots_QCsamples_RAW_replicate_",summary.method," summarization.pdf",sep="")
                        #pdf(pdfname)
                        pdf(pdfname,w=12,h=12) #,w=16,h=10)
                        
                        d1_int_qc<-d1_int[,qc_index]
                        
                        TIC <- vector("list",length(qc_index))
                        
                        for(qc_i in 1:length(qc_index)){
                            
                            
                            #plot( d1$time, d1_int_qc[,qc_i],xlab="Retention Time",ylab="TIC",main=mainlab1)
                            TIC[[qc_i]] <-cbind(d1$time, d1_int_qc[,qc_i])
                        }
                        
                        N=length(qc_index)
                        cols <- rainbow(length(qc_index))
                        lty = 1:N
                        pch = rep(15,N)
                        xlim = range(sapply(TIC, function(x) range(x[,1])))
                        ylim = range(sapply(TIC, function(x) range(x[,2])))
                        plot(0,0, type="n", xlim = xlim, ylim = ylim, main = paste("Feature intensity vs Time \n(QC samples; replicate ",summary.method," summarization)",sep=""), xlab = "Retention Time", ylab = "Intensity")
                        for (i in 1:N) {
                            
                            tic <- TIC[[i]]
                            tic<-tic[order(tic[,1]),]
                            
                            points(tic[,1], (tic[,2]), col = cols[i], pch = pch[i], type="l")
                        }
                        legend("topleft",paste(class_labels[qc_index]), col = cols, lty = lty, pch = pch,cex=0.3)
                        
            
            
                        print("QC file index")
                        print(qc_index)
                        
                        tic<-apply(d1_int_qc,2,function(x){
                            x<-replace(x,which(x==0),NA)
                            return(sum(x,na.rm=TRUE))
                        })
                        
                        
                        mean_tic<-mean(tic)
                        
                        cv_tic<-100*sd(tic,na.rm=TRUE)/mean(tic,na.rm=TRUE)
                        
                        
                        
                       
                        main_lab<-paste("Total ion intensity using all features\n for QC samples (replicate ",summary.method," summarization)\nmean=",round(mean_tic,2),"\n%CV total ion intensity=",round(cv_tic,2),sep="")
                        
                        #tiff("barplot_TIC_using_all_features.tiff",width=2000,height=2000,res=300)
                        #pdf("TIC_all_features.pdf")
                        par(mfrow=c(1,1))
                        try(barplot(tic,cex.names=0.4,cex.axis=1,main=main_lab,col="brown",cex.main=0.6),silent=TRUE)
                        
                            
                        rsqres_listQC<-evaluate.Samples(d1_int_qc, numreplicates=length(qc_index), alignment.tool=NA, cormethod,missingvalue,ignore.missing)
                        
                        rsqres_listQC$cor.matrix<-round(rsqres_listQC$cor.matrix,2)
                        
                        mean_avg_raw_pcor<-rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)]
                        mean_avg_raw_cv<-round(cv_tic,2)
                        
                        qc_res_mat<-rbind(qc_res_mat,cbind(mean_avg_raw_pcor,mean_avg_raw_cv))
                        
			if(summary.method=="mean"){
                        rnames_qc_results<-c(rnames_qc_results,"RAW after averaging")
                        }else{
			rnames_qc_results<-c(rnames_qc_results,"RAW after median summarization")
			}
                        r2<-matrix(0,nrow=length(qc_index),ncol=length(qc_index))
                        
                        #r1[upper.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
                        diag(r2)<-1
                        #r1[lower.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
                        
                        start_count=0
                        for(r1_row in 1:(length(qc_index)-1))
                        {
                            
                            
                            r2[r1_row,c((r1_row+1):length(qc_index))]<-rsqres_listQC$cor.matrix[,(start_count+1):(start_count+length(qc_index)-r1_row)]
                            
                            
                            start_count=(start_count+length(qc_index)-r1_row) #(r1_row*length(qc_index))-r1_row-1 #(length(qc_index)-1)
                            
                            
                        }
                          heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
				heatmap_cols<-rev(heatmap_cols)
            
                        r2[lower.tri(r2)]<-NA #r1[upper.tri(r1,diag=TRUE)] #[-length(rsqres_listQC$cor.matrix)]
                        
                        cor_mainlab<-paste("Pairwise Pearson correlation \n within QC samples (replicate ",summary.method," summarization) \n mean: ",rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)],"; range: ",min(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)])," to ",max(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)]),sep="")
                        
                        h73<-suppressWarnings(try(heatmap.2(r2, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8),silent=TRUE))
                        
                        #h73<-heatmap.2(r2, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8)
                        
                        
                        d1_qc<-cbind(d1[,c(1:10)],d1_int_qc)
                        
                        targeted_feat_calibratedmzAvg<-eval.target.calibratedmz(dataA=d1_qc[,-c(3:10)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="rawQCavg",xMSanalyzer.outloc=subdir4a)
                        try(dev.off(),silent=TRUE)
                        }
                    
                    }
                    }
                    
                    
                }
                
            }
            
            
            
            
            try(dev.off(),silent=TRUE)
            
		#print(batch_levels)
#aplcms
	if(length(batch_levels)>1){ 
		     #######################################################################
		    cat("\n")
		     print(paste("Stage 4b: ComBat processing and post-ComBat QC evaluation using ",best_pair," results",sep=""))
		    cat("\n")
		    #pdf("Stage4b_QE_plots.pdf")
		   # pdf(file=paste("subdir4b/Stage4b_QE_plots.pdf",sep=""))
           #pdf(file=paste(subdir4b,"/Stage4b_QE_plots.pdf",sep=""))
           
           if(summarize.replicates==FALSE){
            pdf(file=paste(subdir4b,"/QE_plots_All_samples_ComBat.pdf",sep=""))
           }else{
               
	        if(summarize.replicates==TRUE && data.norm.pipeline=="AC"){
		
			if(summary.method=="mean"){
			pdf(file=paste(subdir4b,"/QE_plots_All_samples_1averaged_2ComBat.pdf",sep=""))
			}else{
				pdf(file=paste(subdir4b,"/QE_plots_All_samples_1mediansummarized_2ComBat.pdf",sep=""))
			}
		}else{
		
			pdf(file=paste(subdir4b,"/QE_plots_All_samples_ComBat.pdf",sep=""))
		
		}
		
		
               
           }
		    ##################################################################
		
		
	
	
	
	
	 if(is.na(missingvalue)==FALSE){
			
			  if(replacezeroswithNA==TRUE){
					#data_m<-replace(as.matrix(data_m),which(data_m==0),NA)
					data_m_na=replace(as.matrix(data_m),which(data_m==missingvalue),NA)
					d1<-cbind(d1[,c(1:10)],data_m_na)
					
					sum_check<-apply(data_m_na,1,function(x){length(which(is.na(x)==FALSE))})
					
					
					if(length(which(sum_check<minexp))>0){
					data_m_na<-data_m_na[-which(sum_check<minexp),]
					d1<-d1[-which(sum_check<minexp),]
					
					}
					
					if(summarize.replicates==TRUE){
					 print(paste("Number of features after minexp ",minexp," filtering post-averaging",sep=""))
                                         print(dim(data_m)[1])
					 }
					
					}else{
						data_m_na<-data_m
					}

			}else{
					data_m_na<-data_m
					sum_check<-apply(data_m_na,1,function(x){length(which(is.na(x)==FALSE))})
					
					if(length(which(sum_check<minexp))>0){
						data_m_na<-data_m_na[-which(sum_check<minexp),]
						d1<-d1[-which(sum_check<minexp),]
					}
					
					 if(summarize.replicates==TRUE){
					 print(paste("Number of features after minexp ",minexp," filtering post-averaging",sep=""))
                                         print(dim(data_m)[1])
					 }
		   }
		   
		  
		   
		    adjdata<-try(sva::ComBat(dat=data_m_na,batch=batch_inf,mod=mod,par.prior=TRUE),silent=TRUE)
		    
		    if(is(adjdata,"try-error")){
		 
				print("Dim input")
				print(dim(data_m_na))
				print(dim(d1)[1])
				print(dim(sampleid_mapping))
				data_m1<-cbind(d1[,c(1:4)],data_m_na)
				
				
				adjdata<-MetabComBat(dat=data_m1,saminfo=sampleid_mapping,par.prior=T,filter=F,write=F,prior.plots=F)
				print(adjdata[1:3,1:5])
				adjdata<-adjdata[,-c(1:4)]
				
         
		    }
            
            #save(adjdata,file="adjdata.Rda")
    				  
                      #adjdata<-as.data.frame(adjdata)
		    maxint<-apply(adjdata,1,function(x){max(x,na.rm=TRUE)})
		    
		    adjdata2<-replace(as.matrix(adjdata),which(is.na(adjdata)==TRUE),missingvalue)
		    adjdata2<-replace(as.matrix(adjdata2),which(adjdata2<0),0)
		    
		    adjdata2<-round(adjdata2,1)
		    adjdata2<-cbind(d1[,c(1:4)],adjdata2)
            
           
            
		    adjdata2$time<-round(adjdata2$time,1) 
		    adjdata2$mz<-round(adjdata2$mz,5) 
            #adjdata2$min.mz<-round(adjdata2$min.mz,5)
            #adjdata2$max.mz<-round(adjdata2$max.mz,5)
	
    #expression_xls<-paste(subdir4b,"/ComBat_corrected_",best_pair,"_feature_table.txt",sep="")
	
    if(impute.bool=="TRUEBfsfs1" & summarize.replicates==TRUE){
        
        expression_xls<-paste(xMSanalyzer.outloc,"/ComBat_mzcalibrated_untargeted_imputed_averaged_featuretable.txt",sep="")
    }else{
        
        if(impute.bool=="TRUEBfsfs1" & summarize.replicates==FALSE){
            expression_xls<-paste(xMSanalyzer.outloc,"/ComBat_mzcalibrated_untargeted_imputed_featuretable.txt",sep="")
        }else{
            
            if(impute.bool==FALSE & summarize.replicates==TRUE){
                expression_xls<-paste(xMSanalyzer.outloc,"/ComBat_mzcalibrated_untargeted_averaged_featuretable.txt",sep="")
            }else{
                
                if(impute.bool==FALSE & summarize.replicates==FALSE){
                    expression_xls<-paste(xMSanalyzer.outloc,"/ComBat_mzcalibrated_untargeted_featuretable.txt",sep="")
                }

            }

        }
    }
    
    if(summarize.replicates==TRUE && data.norm.pipeline=="AC"){
    
	if(summary.method=="mean"){
        expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_averaged_featuretable.txt",sep="")
	}else{
		if(summary.method=="median"){
			expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt",sep="")
		}
	
	}
    }else{
    expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_featuretable.txt",sep="")
    }
    
    #write.table(adjdata2,file=expression_xls,sep="\t",row.names=FALSE)
		   							
			cnames=colnames(feat.eval.result)
										
			feat.eval.result<-apply(feat.eval.result,2,as.numeric)
			feat.eval.result<-as.data.frame(feat.eval.result)
			#feat.eval.result.mat=cbind(adjdata2[,c(1:4)],feat.eval.result)  
			
		    numpeaks<-apply(adjdata2[,-c(1:4)],1,countpeaks)						
		    
		    maxint<-apply(adjdata2[,-c(1:4)],1,function(x){max(x,na.rm=TRUE)})

    
		    feat.eval.result=evaluate.Features(adjdata2, numreplicates=num_replicates,min.samp.percent=min.samp.percent,
		    alignment.tool="apLCMS",impute.bool=impute.bool,peak_scores=d1$PeakScore,numnodes=numnodes)
			
		    adjdata2<-cbind(adjdata2[,c(1:4)],numpeaks,feat.eval.result$numgoodsamples,feat.eval.result$median,
		    d1$PeakScore,feat.eval.result$Qscore,maxint,adjdata2[,-c(1:4)])
		    
		    colnames(adjdata2)<-colnames(d1)
		    setwd(subdir4b)
		    write.table(adjdata2,file=expression_xls,sep="\t",row.names=FALSE)
		    
		    
		    	    
			
		   # X<-t(adjdata2[maxint_ord[1:2000],-c(1:9)])
		   
		    
		    
		     tic.eval(adjdata2[,-c(1:10)],outloc=subdir4b)
    
		   
		     X<-t(adjdata2[,-c(1:10)])
		   
		     
		    X<-as.matrix(X)
		    X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
		    try(pca.eval(X,samplelabels=batch_labels,filename="ComBat",ncomp=5,center=TRUE,scale=TRUE,legendlocation="topright",legendcex=0.5,outloc=subdir4b),silent=TRUE)
			    
		    
		    #eval.reference.metabs(dataA=adjdata2,stdData=stddata,mzthresh=10,outloc=getwd(),folderheader="ComBat")
		    
		targeted_feat_combat<-eval.target.calibratedmz(dataA=adjdata2[,-c(3:10)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4b,folderheader="ComBat",xMSanalyzer.outloc=subdir4b)
			
            
            fnamesT<-paste(subdir4b,"/ComBat_mzcalibrated_targeted_featuretable_mz",refMZ.mz.diff,"ppm_time",refMZ.time.diff,".txt",sep="")
            write.table(targeted_feat_combat$targetdata,file=fnamesT,sep="\t",row.names=FALSE)
            
            
            
		s4b<-"Stage 4b results: Batch-effect evaluation post ComBat and QC evaluation of targeted data post batch-effect correction"
		sm<-rbind(sm,s4b)
	
		}else{
			print("Only one batch found. Skipping Stage 4b.")
			adjdata2=d1
			
			adjdata2<-replace(as.matrix(adjdata2),which(is.na(adjdata2)==TRUE),missingvalue)
			adjdata2<-replace(as.matrix(adjdata2),which(adjdata2<0),0)
		    
			
                        targeted_feat_combat<-targeted_feat_raw
		}
		

		}else{
			
			adjdata2=d1
			targeted_feat_combat<-targeted_feat_raw
			
			adjdata2<-replace(as.matrix(adjdata2),which(is.na(adjdata2)==TRUE),missingvalue)
			adjdata2<-replace(as.matrix(adjdata2),which(adjdata2<0),0)
		    
			
			
			}
    
			try(dev.off(),silent=TRUE)
            
	    if(length(batch_levels)>1){
            if(is.na(qc_label)==FALSE){
                   qc_label<-tolower(qc_label)
                d1_int<-adjdata2[,-c(1:10)]
                qc_label_check<-gregexpr(pattern=qc_label,text=sampleid_mapping[,2])
                
                qc_index<-which(qc_label_check>0)
                
		print(qc_index)
                if(length(qc_index)>0){
                    
                    
                    if(summarize.replicates==TRUE && data.norm.pipeline=="AC"){
                        if(summary.method=="mean"){
				mainlab_C<-"Replicate averaging+ComBat"
				pdfname=paste(subdir4b,"/QE_plots_QCsamples_1averaged_2ComBat.pdf",sep="")
		    
			}else{
				 if(summary.method=="median"){
				mainlab_C<-"Replicate median summarization+ComBat"
				pdfname=paste(subdir4b,"/QE_plots_QCsamples_1mediansummarized_2ComBat.pdf",sep="")
		    
				}
				
			}
                    }else{
                        
                        mainlab_C<-"ComBat"
                        pdfname=paste(subdir4b,"/QE_plots_QCsamples_ComBat.pdf",sep="")
                    }
                    #pdf(pdfname)
                    pdf(pdfname,w=12,h=12) #,w=16,h=10)
                    
                    d1_int_qc<-d1_int[,qc_index]
                    
                    TIC <- vector("list",length(qc_index))
                    
                    for(qc_i in 1:length(qc_index)){
                        
                        
                        #plot( d1$time, d1_int_qc[,qc_i],xlab="Retention Time",ylab="TIC",main=mainlab1)
                        TIC[[qc_i]] <-cbind(d1$time, d1_int_qc[,qc_i])
                    }
                    
		    #save(TIC,file="TIC.Rda")
		    
                    N=length(qc_index)
                    cols <- rainbow(length(qc_index))
                    lty = 1:N
                    pch = rep(15,N)
                    xlim = range(sapply(TIC, function(x) range(x[,1])))
                    ylim = range(sapply(TIC, function(x) range(x[,2])))
                    plot(0,0, type="n", xlim = xlim, ylim = ylim, main = paste("Feature intensity vs Time \n(QC samples ",mainlab_C,")",sep=""), xlab = "Retention Time", ylab = "Intensity")
                    for (i in 1:N) {
                        
                        tic <- TIC[[i]]
                        tic<-tic[order(tic[,1]),]
                        
                        points(tic[,1], (tic[,2]), col = cols[i], pch = pch[i], type="l")
                    }
                    legend("topleft",paste(class_labels[qc_index]), col = cols, lty = lty, pch = pch,cex=0.3)
                    
                    
                    
                    print("QC file index")
                    print(qc_index)
                    
                    tic<-apply(d1_int_qc,2,function(x){
                        x<-replace(x,which(x==0),NA)
                        return(sum(x,na.rm=TRUE))
                    })
                    
                    
                    mean_tic<-mean(tic)
                    
                    cv_tic<-100*sd(tic,na.rm=TRUE)/mean(tic,na.rm=TRUE)
                    
                    
                    
                    
                    main_lab<-paste("Total ion intensity using all features\n for QC samples after ",mainlab_C,")","\nmean=",round(mean_tic,2),"\n%CV total ion intensity=",round(cv_tic,2),sep="")
                    
                    #tiff("barplot_TIC_using_all_features.tiff",width=2000,height=2000,res=300)
                    #pdf("TIC_all_features.pdf")
                    par(mfrow=c(1,1))
                    barplot(tic,cex.names=0.4,cex.axis=1,main=main_lab,col="brown",cex.main=0.6)
                    



                    rsqres_listQC<-evaluate.Samples(d1_int_qc, numreplicates=length(qc_index), alignment.tool=NA, cormethod,missingvalue,ignore.missing)
                    
                    rsqres_listQC$cor.matrix<-round(rsqres_listQC$cor.matrix,2)
                    
                    
                    mean_combat_pcor<-rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)]
                    mean_combat_cv<-round(cv_tic,2)
                    
                    qc_res_mat<-rbind(qc_res_mat,cbind(mean_combat_pcor,mean_combat_cv))
                    
                    if(summarize.replicates==TRUE && data.norm.pipeline=="AC"){
		    
		    
                    rnames_qc_results<-c(rnames_qc_results,paste(summary.method, " summarization+ComBat (AC)",sep=""))
                    }else{
                        rnames_qc_results<-c(rnames_qc_results,"ComBat")
                        
                    }
                    
                    r2<-matrix(0,nrow=length(qc_index),ncol=length(qc_index))
                    
                    #r1[upper.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
                    diag(r2)<-1
                    #r1[lower.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
                    
                    start_count=0
                    for(r1_row in 1:(length(qc_index)-1))
                    {
                        
                        
                        r2[r1_row,c((r1_row+1):length(qc_index))]<-rsqres_listQC$cor.matrix[,(start_count+1):(start_count+length(qc_index)-r1_row)]
                        
                        
                        start_count=(start_count+length(qc_index)-r1_row) #(r1_row*length(qc_index))-r1_row-1 #(length(qc_index)-1)
                        
                        
                    }
                    
                    r2[lower.tri(r2)]<-NA #r1[upper.tri(r1,diag=TRUE)] #[-length(rsqres_listQC$cor.matrix)]
                    
                    cor_mainlab<-paste("Pairwise Pearson correlation \n within QC samples (after ComBat) \n mean: ",rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)],"; range: ",min(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)])," to ",max(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)]),sep="")
                      heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
			heatmap_cols<-rev(heatmap_cols)
            
                    h73<-suppressWarnings(heatmap.2(r2, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8))
                    
                    #h73<-heatmap.2(r2, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8)
                    
                    
                    d1_qc<-cbind(adjdata2[,c(1:10)],d1_int_qc)
                    
                    targeted_feat_calibratedmzComBat<-eval.target.calibratedmz(dataA=d1_qc[,-c(3:10)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="ComBat",xMSanalyzer.outloc=subdir4b)
                    try(dev.off(),silent=TRUE)
                    
                    }
                    
                    
                    if(data.norm.pipeline=="CA" && summarize.replicates==TRUE){
                        sampleid_mapping<-sampleid_mapping[seq(1,dim(adjdata2[,-c(1:10)])[2],num_replicates),]
                        
                        
                        batch_inf<-batch_inf[seq(1,dim(adjdata2[,-c(1:10)])[2],num_replicates)]
                        
                        mod<-mod[seq(1,dim(adjdata2[,-c(1:10)])[2],num_replicates)]
                        
			batch_labels<-batch_labels[seq(1,dim(adjdata2[,-c(1:10)])[2],num_replicates)]
			
			l1<-levels(batch_labels)
                       
			 minexp=minexp/num_replicates 
                        data_matrix<-data_summarize(Xmat=adjdata2[,- c(3:10)],Ymat=NA,feature_table_file=NA,parentoutput_dir=subdir4b,class_labels_file=NA,num_replicates=num_replicates,summarize.replicates=summarize.replicates,
                        summary.method=summary.method,missing.val=missingvalue, rep.num.max.missing.thresh=rep.num.max.missing.thresh,summary.na.replacement=summary.na.replacement,fileheader="ComBat")
                        
                          check_bad_feats_allzeros<-apply(data_matrix[,-c(1:2)],1,sum)
		    bad_feats_allzeros_index<-which(check_bad_feats_allzeros==0)
		    
		    if(length(bad_feats_allzeros_index)>0){
		    
			data_matrix<-data_matrix[-bad_feats_allzeros_index,]
			
			adjdata2<-adjdata2[-bad_feats_allzeros_index,]
		    }
		    
                        adjdata2<-cbind(adjdata2[,c(1:10)],data_matrix[,-c(1:2)])
                        
                        d1_int<-adjdata2[,-c(1:10)]
                        X<-t(data_matrix[,-c(1:2)])
                        
                        data_m<-data_matrix[,-c(1:2)]
                        
                        
                        num_replicates=1
                        
                        sampleid_mapping[,1]<-gsub(sampleid_mapping[,1],pattern=filepattern,replacement="",ignore.case=TRUE)
                        cnames<-colnames(data_m)
                        
                        cnames<-gsub(cnames,pattern=filepattern,replacement="",ignore.case=TRUE)
                        colnames(data_m)<-cnames
                        
                        class_labels<-as.factor(sampleid_mapping[,2])
                        
                        if(is.na(qc_label)==FALSE){
                                qc_label<-tolower(qc_label)
                            qc_label_check<-gregexpr(pattern=qc_label,text=sampleid_mapping[,2])
                            
                            qc_index<-which(qc_label_check>0)
                            
                            if(length(qc_index)>0){
                                
				
				 if(summary.method=="mean"){
				
				pdfname=paste(subdir4b,"/QE_plots_QCsamples_1ComBat_2averaging.pdf",sep="")
		    
				}else{
					if(summary.method=="median"){
				
					pdfname=paste(subdir4b,"/QE_plots_QCsamples_1ComBat_2mediansummarization.pdf",sep="")
		    
					}
				
				}
			
                                #pdfname=paste(subdir4b,"/QE_plots_QCsamples_ComBat_averaged.pdf",sep="")
                                #pdf(pdfname)
                                pdf(pdfname,w=12,h=12) #,w=16,h=10)
                                
                                d1_int_qc<-d1_int[,qc_index]
                                
                                TIC <- vector("list",length(qc_index))
                                
                                for(qc_i in 1:length(qc_index)){
                                    
                                    
                                    #plot( d1$time, d1_int_qc[,qc_i],xlab="Retention Time",ylab="TIC",main=mainlab1)
                                    TIC[[qc_i]] <-cbind(d1$time, d1_int_qc[,qc_i])
                                }
                                
                                N=length(qc_index)
                                cols <- rainbow(length(qc_index))
                                lty = 1:N
                                pch = rep(15,N)
                                xlim = range(sapply(TIC, function(x) range(x[,1])))
                                ylim = range(sapply(TIC, function(x) range(x[,2])))
				main_text<-paste("Feature intensity vs Time \n(QC samples ComBat+replicate ",summary.method," summarization)",sep="")
                                plot(0,0, type="n", xlim = xlim, ylim = ylim, main = main_text, xlab = "Retention Time", ylab = "Intensity")
                                for (i in 1:N) {
                                    
                                    tic <- TIC[[i]]
                                    tic<-tic[order(tic[,1]),]
                                    
                                    points(tic[,1], (tic[,2]), col = cols[i], pch = pch[i], type="l")
                                }
                                legend("topleft",paste(class_labels[qc_index]), col = cols, lty = lty, pch = pch,cex=0.3)
                                
                                
                                
                                print("QC file index")
                                print(qc_index)
                                
                                tic<-apply(d1_int_qc,2,function(x){
                                    x<-replace(x,which(x==0),NA)
                                    return(sum(x,na.rm=TRUE))
                                })
                                
                                
                                mean_tic<-mean(tic)
                                
                                cv_tic<-100*sd(tic,na.rm=TRUE)/mean(tic,na.rm=TRUE)
                                
                                
                                
                                
                                main_lab<-paste("Total ion intensity using all features\n for QC samples (ComBat+replicate ",summary.method," summarization)\nmean=",round(mean_tic,2),"\n%CV total ion intensity=",round(cv_tic,2),sep="")
                                
                                #tiff("barplot_TIC_using_all_features.tiff",width=2000,height=2000,res=300)
                                #pdf("TIC_all_features.pdf")
                                par(mfrow=c(1,1))
                                barplot(tic,cex.names=0.4,cex.axis=1,main=main_lab,col="brown",cex.main=0.6)
                                
                                
                                rsqres_listQC<-evaluate.Samples(d1_int_qc, numreplicates=length(qc_index), alignment.tool=NA, cormethod,missingvalue,ignore.missing)
                                
                                rsqres_listQC$cor.matrix<-round(rsqres_listQC$cor.matrix,2)
                                
                                mean_avg_raw_pcor<-rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)]
                                mean_avg_raw_cv<-round(cv_tic,2)
                                
                                qc_res_mat<-rbind(qc_res_mat,cbind(mean_avg_raw_pcor,mean_avg_raw_cv))
                                
                                rnames_qc_results<-c(rnames_qc_results,paste("ComBat+",summary.method," summarization (CA)",sep=""))
                                
                                r2<-matrix(0,nrow=length(qc_index),ncol=length(qc_index))
                                
                                #r1[upper.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
                                diag(r2)<-1
                                #r1[lower.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
                                
                                start_count=0
                                for(r1_row in 1:(length(qc_index)-1))
                                {
                                    
                                    
                                    r2[r1_row,c((r1_row+1):length(qc_index))]<-rsqres_listQC$cor.matrix[,(start_count+1):(start_count+length(qc_index)-r1_row)]
                                    
                                    
                                    start_count=(start_count+length(qc_index)-r1_row) #(r1_row*length(qc_index))-r1_row-1 #(length(qc_index)-1)
                                    
                                    
                                }
                                
                                r2[lower.tri(r2)]<-NA #r1[upper.tri(r1,diag=TRUE)] #[-length(rsqres_listQC$cor.matrix)]
                                
                                cor_mainlab<-paste("Pairwise Pearson correlation \n within QC samples (ComBat+ replicate ",summary.method," summarization) \n mean: ",rsqres_listQC$cor.matrix[length(rsqres_listQC$cor.matrix)],"; range: ",min(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)])," to ",max(rsqres_listQC$cor.matrix[-length(rsqres_listQC$cor.matrix)]),sep="")
                                  heatmap_cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
					heatmap_cols<-rev(heatmap_cols)
            
                                h73<-suppressWarnings(heatmap.2(r2, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8))
                                
                                #h73<-heatmap.2(r2, Rowv=NULL, Colv=NULL,  col=heatmap_cols, scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.7, cexCol=0.7,xlab="",ylab="", main=cor_mainlab,cex.main=0.8)
                                
                                
                                d1_qc<-cbind(d1[,c(1:10)],d1_int_qc)
                                
                                targeted_feat_calibratedmzAvg<-eval.target.calibratedmz(dataA=d1_qc[,-c(3:10)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="ComBatQCavg",xMSanalyzer.outloc=subdir4a)
                                try(dev.off(),silent=TRUE)
                            }
                            
                        }
                    }
                }
                
            }
			
            
            if(is.na(qc_label)==FALSE){
             
             if(length(qc_res_mat)>0){
                 colnames(qc_res_mat)<-c("mean pairwise Pearson correlation within QC samples","%CV for total ion intensity for QC samples")
                 rownames(qc_res_mat)<-rnames_qc_results
                 write.csv(qc_res_mat,file=paste(xMSanalyzer.outloc,"/Evaluation_report_",qc_label,"_QCsamples.csv",sep=""))
                
             }
            }
					 metlin.res={}
                      kegg.res={}
                      
					 #length(union_list[[num_pairs]]$mz
					 if(is.na(adduct.list)==FALSE){
					
					cat("\n")
					  print("*********Stage 5: Mapping m/z values to known metabolites*********")
					 cat("\n")
					 
                     #annot.res<-feat.batch.annotation.KEGG(adjdata2,mz.tolerance.dbmatch,adduct.list,subdir5, numnodes=numnodes,syssleep=syssleep)
                     
                     for(db in db_name){
                     annot.res<-simpleAnnotation(dataA=adjdata2,max.mz.diff=mz.tolerance.dbmatch,num_nodes=numnodes,queryadductlist=adduct.list,
                     gradienttype="Acetonitrile",mode=charge_type,outloc=subdir5,db_name=db)
                     
                      fname<-paste(subdir5,"/DBmatches_",db,".txt",sep="")
                     
                     #fname<-paste(xMSanalyzer.outloc,"/DBmatches_",db,".txt",sep="")
                     write.table(annot.res,file=fname,sep="\t",row.names=FALSE)
                     
                     }
                     
                                    
					    annotres_list<-annot.res
					    s5<-"Stage 5 results: Annotation of features"
					    sm<-rbind(sm,s5)
                     			}
				

		          
               		  print("*************Processing complete**********") 

                }
                else
                {
                        stop(paste("No files exist in",apLCMS.outloc, "Please check the input value for cdfloc", sep=""))
                }
                
                suppressWarnings(sink(file=NULL))
                #sm<-"The results include pairwise Pearson correlation within technical replicates; 2) file with raw m/z and time for each feature; 3) raw (non batch-effect corrected) feature table after m/z calibration; 4) ComBat processed feature table after m/z calibration;  5) Simple (m/z based) database searches"
sm<-as.data.frame(sm)
colnames(sm)<-"Output_description"
setwd(xMSanalyzer.outloc)
write.table(sm,file="Readme.txt",sep="\t",row.names=FALSE)
                                print(paste("Processing is complete. Program results can be found at: ",xMSanalyzer.outloc,sep=""))
		
        #return(list("mergeresult"=union_list, "apLCMSres"=data_rpd_all,"apLCMSresfilt"=data_rpd))
	return(list("apLCMS.merged.res"=union_list, "apLCMS.ind.res"=data_rpd_all,"apLCMS.ind.res.filtered"=data_rpd,"final.feat.table.annot"=annotres_list, "feat.eval.ind"=feateval_list, "sample.eval.ind"=rsqres_list,"final.feat.table.raw"=d1,"final.feat.table.combat"=adjdata2,
	"final.targeted.feat.table.raw"=targeted_feat_raw,"final.targeted.feat.table.combat"=targeted_feat_combat))
	
	
}

##################################################################
#Function:xMSwrapper.XCMS
#Description: wrapper function based on apLCMS.align,evaluate.Features,
#            evaluate.Samples,merge.Results,search.Metlin, and 
#            search.KEGG
################################################################
xMSwrapper.XCMS.matchedFilter<-function(cdfloc, XCMS.outloc,xMSanalyzer.outloc,step.list=c(0.001,0.01,0.1), mz.diff.list=c(0.001,0.01,0.1), sn.thresh.list=c(10),
max=50, bw=c(10), minfrac.val=0.5, minsamp=2, mzwid=0.25, sleep=0, retcor.family = "symmetric", retcor.plottype = "mdevden", groupval.method = "medret",
numnodes=2,run.order.file=NA,max.mz.diff=15,max.rt.diff=300,
merge.eval.pvalue=0.2,mergecorthresh=0.7,deltamzminmax.tol=10,num_replicates=3,subs=NA, mz.tolerance.dbmatch=10,
adduct.list=c("M+H"), samp.filt.thresh=0.70,feat.filt.thresh=50,cormethod="pearson", mult.test.cor=TRUE,
missingvalue=0,ignore.missing=TRUE,sample_info_file=NA,refMZ=NA,refMZ.mz.diff=10,refMZ.time.diff=NA,
void.vol.timethresh=30,replacezeroswithNA=TRUE,scoreweight=30,
filepattern=".cdf",charge_type="pos", minexp.pct=0.1,syssleep=0.5,merge.pairwise=FALSE,
min.samp.percent=0.6,impute.bool=TRUE,summarize.replicates=TRUE,
summary.method="median",max.number.of.replicates.with.missingvalue=1,summary.na.replacement="zeros",db_name=c("KEGG","HMDB","LipidMaps"),data.norm.pipeline="AC",calibration.method=c("median.adjustment","multiplicative.signal.correction"))
{
suppressWarnings(sink(file=NULL))
dir.create(xMSanalyzer.outloc,showWarnings=FALSE)

targeted_feat_raw<-{}
library(parallel)
x<-date()
x<-strsplit(x,split=" ")
x1<-unlist(x)

rep.num.max.missing.thresh=max.number.of.replicates.with.missingvalue
#cdfloc=NA

inp_params<-{}
inp_params<-rbind(inp_params,cbind("cdfloc: ",cdfloc))
inp_params<-rbind(inp_params,cbind("XCMS.matchedFilter.outloc: ",XCMS.outloc))
inp_params<-rbind(inp_params,cbind("xMSanalyzer.outloc: ",xMSanalyzer.outloc))
inp_params<-rbind(inp_params,cbind("num_replicates: ",num_replicates))
inp_params<-rbind(inp_params,cbind("max.mz.diff: ",max.mz.diff))
inp_params<-rbind(inp_params,cbind("max.rt.diff: ",max.rt.diff))
inp_params<-rbind(inp_params,cbind("merge.eval.pvalue: ",merge.eval.pvalue))
inp_params<-rbind(inp_params,cbind("mergecorthresh: ",mergecorthresh))
inp_params<-rbind(inp_params,cbind("deltamzminmax.tol: ",deltamzminmax.tol))
inp_params<-rbind(inp_params,cbind("mz.tolerance.dbmatch: ",mz.tolerance.dbmatch))
inp_params<-rbind(inp_params,cbind("adduct.list: ",paste(adduct.list,collapse=";")))
inp_params<-rbind(inp_params,cbind("samp.filt.thresh: ",samp.filt.thresh))
inp_params<-rbind(inp_params,cbind("feat.filt.thresh: ",feat.filt.thresh))
inp_params<-rbind(inp_params,cbind("cormethod: ",cormethod))
inp_params<-rbind(inp_params,cbind("charge_type: ",charge_type))
inp_params<-rbind(inp_params,cbind("db_name: ",paste(db_name,collapse=";")))
inp_params<-rbind(inp_params,cbind("impute.bool: ",impute.bool))
inp_params<-rbind(inp_params,cbind("merge.pairwise: ",merge.pairwise))
inp_params<-rbind(inp_params,cbind("summarize.replicates: ",summarize.replicates))
inp_params<-rbind(inp_params,cbind("summary.method: ",summary.method))
inp_params<-rbind(inp_params,cbind("max.number.of.replicates.with.missingvalue: ",rep.num.max.missing.thresh))
inp_params<-rbind(inp_params,cbind("summary.na.replacement: ",summary.na.replacement))



fname<-paste(xMSanalyzer.outloc,"/Input_parameters.csv",sep="")

colnames(inp_params)<-c("Parameter","Value")
write.csv(inp_params,file=fname,row.names=FALSE)

#variables not used
#summarize.replicates=FALSE
#summary.method="median"
#rep.num.max.missing.thresh=1
#summary.na.replacement="zeros"


x1<-gsub(x1,pattern=":",replacement="_")
fname<-paste(x1[2:4],collapse="")
fname<-gsub(fname,pattern=":",replacement="_")
#fname<-paste(fname,x1[5],sep="")
x1[5]<-gsub(x1[5],pattern=":",replacement="_")
fname<-paste(fname,x1[5],sep="_")

	

	fname<-paste(xMSanalyzer.outloc,"/Log_",fname,".txt",sep="")
	print(paste("Program is running. Please check the logfile for runtime status: ",fname,sep=""))

	data_rpd_all=new("list")
	data_rpd=new("list")
	union_list=new("list")
	feateval_list<-new("list")
	rsqres_list<-new("list")
	annot.res<-{}
	annotres_list<-new("list")
	best_score<-(-1000000)
		if(is.na(cdfloc)==FALSE){
	
		if(cdfloc=="NA"){
			cdfloc=NA
		}
	}
	
	if(is.na(sample_info_file)==FALSE){
	
		if(sample_info_file=="NA"){
			sample_info_file=NA
		}
	}

filepattern=".cdf|.mzxml|mXML"	
	
	if(is.na(sample_info_file)==FALSE)
    {
    	
    	
		sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)

		if(is.na(cdfloc)==FALSE){
			l1<-list.files(cdfloc,filepattern,ignore.case=TRUE)
	
			minexp<-round(minfrac.val*length(l1))
			
			}else{
				l1<-rep(1,num_replicates)	
			}
		if(length(l1)!=dim(sampleid_mapping)[1] & (is.na(cdfloc)==FALSE))
		{
			num_mis_files<-dim(sampleid_mapping)[1]-length(l1)
			if(is.na(subs)==TRUE){
				stop(paste("ERROR: Only ",length(l1)," spectral files were found. ",num_mis_files," files are missing.",sep=""))
			}
		}
	}else{
	
	
		if(is.na(cdfloc)==FALSE){
			l1<-list.files(cdfloc,filepattern,ignore.case=TRUE)
			
			minexp<-round(minfrac.val*length(l1))
		}else{
			l1<-rep(1,num_replicates)	
		}
	}

if(length(l1)%%num_replicates>0)
{stop(paste("ERROR: Not all samples have ",num_replicates," replicates.",sep=""))
}

dir.create(XCMS.outloc,showWarnings=FALSE)
	dir.create(xMSanalyzer.outloc,showWarnings=FALSE)
	#subdir1<-paste(xMSanalyzer.outloc,"/Quality_assessment_files",sep="")
	#subdir2<-paste(xMSanalyzer.outloc,"/XCMS_filtered_data",sep="")
	#subdir3<-paste(xMSanalyzer.outloc,"/XCMS_with_xMSanalyzer_merged_data",sep="")
	
	
	sink(fname)
	print(sessionInfo())

	if(is.na(refMZ)==FALSE){
                        stddata<-read.table(refMZ,sep="\t",header=TRUE)
                        print(refMZ)
                        print(head(stddata))

                }else{
                        if(charge_type=="pos"){
                        data(example_target_list_pos)
                        stddata<-example_target_list_pos
                        }else{

                                if(charge_type=="neg"){
                                                data(example_target_list_neg)
                                                stddata<-example_target_list_neg
                                        }else{
                                                stop("Invalid option. \'charge_type\' should be \'pos\' or \'neg\'.")
                                                }
                                }
                }


        ############################################
        #1) Align profiles using the cdf.to.ftr wrapper function in apLCMS
	 if(is.na(XCMS.outloc)==TRUE)
        {
                stop("Undefined value for parameter, XCMS.outloc. Please define the XCMS output location.")

        }
         if(is.na(xMSanalyzer.outloc)==TRUE)
        {
                stop("Undefined value for parameter, xMSanalyzer.outloc. Please define the xMSanalyzer output location.")

        }


        if(is.na(cdfloc)==FALSE)
        {
                setwd(cdfloc)
                
                if(is.na(XCMS.outloc)==FALSE)
                {
                    
                    
                        data_rpd_all=XCMS.align.matchedFilter(cdfloc, XCMS.outloc, step.list = step.list, mz.diff.list = mz.diff.list,
                         sn.thresh.list = sn.thresh.list, max = max, bw.val = bw, minfrac.val = minfrac.val,
                         minsamp.val = minsamp, mzwid.val = mzwid, sleep.val = sleep, run.order.file = run.order.file,
                         subs = subs, retcor.family = retcor.family, retcor.plottype = retcor.plottype,
                         groupval.method = groupval.method,target.mz.list = stddata,xMSanalyzer.outloc=xMSanalyzer.outloc)
                    
                    
                    
                        
                        #data_rpd_all=XCMS.align.matchedFilter(cdfloc, XCMS.outloc,step.list, mz.diff.list, sn.thresh.list, max, bw, minfrac.val, minsamp, mzwid, sleep,
                        #run.order.file,subs, retcor.family, retcor.plottype,groupval.method)
			
                       
                }
                else
                {
                        stop("Undefined value for parameter, XCMS.outloc. Please define the output location.")
                }
        }
	
		setwd(XCMS.outloc)
                alignmentresults<-list.files(XCMS.outloc, "*.txt")
		print("Files found in XCMS output location:")
                print(alignmentresults)
		for(i in 1:length(alignmentresults)){
	
			data_rpd_all[[i]]<-read.table(paste(XCMS.outloc,"/",alignmentresults[i],sep=""),header=TRUE)
			data_rpd_all[[i]]<-unique(data_rpd_all[[i]]) 	
			
		}
		if(is.na(sample_info_file)==FALSE)
		{
			 match_names_check<-match(sampleid_mapping[,1],colnames(data_rpd_all[[1]][,-c(1:4)]))
			 
			 if(length(which(is.na(match_names_check)==TRUE))>0){
				stop("Sample names do not match between sequence file and feature table.")
			 }
	
		}
	
    subdir1<-paste(xMSanalyzer.outloc,"/Stage1",sep="")  #QC individual parameter settings
    subdir2<-paste(xMSanalyzer.outloc,"/Stage2",sep="")  #Data filtering
    subdir3<-paste(xMSanalyzer.outloc,"/Stage3a",sep="")  #Data merger/parameter optimization
    subdir3b<-paste(xMSanalyzer.outloc,"/Stage3b",sep="")
    subdir4a<-paste(xMSanalyzer.outloc,"/Stage4a",sep="")	 #Raw QC: batch effect eval, TIC, etc
    
    
    dir.create(subdir1,showWarnings=FALSE)
    dir.create(subdir2,showWarnings=FALSE)
    dir.create(subdir3,showWarnings=FALSE)
    dir.create(subdir3b,showWarnings=FALSE)
    dir.create(subdir4a,showWarnings=FALSE)
    
    #if(is.na(sample_info_file)==FALSE)
	   if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
       {
           subdir4b<-paste(xMSanalyzer.outloc,"/Stage4b",sep="")	 #Batch-effect corrected QC: batch effect eval, TIC, etc
           dir.create(subdir4b,showWarnings=FALSE)
           
       }
       
       if(is.na(adduct.list)==FALSE){
           subdir5<-paste(xMSanalyzer.outloc,"/Stage5",sep="")	 #Putative unprocessed annotations;
           
           
           dir.create(subdir5,showWarnings=FALSE)
       }
		bestscore<-(-1000000)
	
        {
                #stop("Undefined value for parameter, cdfloc. Please enter path of the folder where the CDF files to be processed are located.")
                #change location to the output folder
                setwd(XCMS.outloc)
                alignmentresults<-list.files(XCMS.outloc, "*.txt")
                
                if(length(data_rpd_all)>0)
                {
                          curdata_dim={}
                          if(num_replicates==2)
                          {
                                  fileroot="_PID"
                          }
                          else
                          {
                                  if(num_replicates>2)
                                  {
                                          fileroot="_CV"
                                  }
                                  else
                                  {
                                          fileroot=""
                             
					#stop("Need at least 2 technical replicates per sample.")
				     }
                          }
                          #for(i in 1:length(alignmentresults))
                          
                          cat("\n")
                          print("*******xMSanalyzer Stage 1: QC evaluation of invidual parameters*******")
                          cat("\n")
			  
			rep_cor_mat<-{}
            rep_cor_mat<-{}
            parent_bad_list<-{}
            
			       for(i in 1:length(data_rpd_all))
                          {
                              
                              print("dim of data is ")
                              print(dim(data_rpd_all[[i]]))
                          	
				  print(paste("**Evaluating XCMS results from parameter setting ",i,"**",sep=""))				  
                                  ############################################
                                  #2)Calculate pairwise correlation coefficients
                                  file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
                                  #curdata=read.table(paste(XCMS.outloc,"/",alignmentresults[i],sep=""),header=TRUE)
                                  #curdata=check.mz.in.replicates(curdata)
                                  curdata=data_rpd_all[[i]]
                                  #############################################
                                  ############################################
                                  #3) Calculate Percent Intensity Difference
			                   if(num_replicates>1)
                                  {
							
								
								
								
								feat.eval.result=evaluate.Features(curdata, numreplicates=num_replicates,min.samp.percent=min.samp.percent,alignment.tool="XCMS",impute.bool=impute.bool,numnodes=numnodes)
								cnames=colnames(feat.eval.result)
								feat.eval.result<-apply(feat.eval.result,2,as.numeric)
								feat.eval.result<-as.data.frame(feat.eval.result)
								feat.eval.result.mat=cbind(curdata[,c(1:8)],feat.eval.result)  
								feat.eval.outfile=paste(subdir1,"/",file_name,fileroot,"featureassessment.txt",sep="")					  
								#write results
								write.table(feat.eval.result.mat, feat.eval.outfile,sep="\t", row.names=FALSE)
								
                                curdata<-curdata[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),]
                                
                                feat.eval.result.mat<-feat.eval.result.mat[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),]
								
								#curdata<-cbind(curdata,feat.eval.result[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),])
								
								curdata<-as.data.frame(curdata)
								curdata<-replace(as.matrix(curdata),which(is.na(curdata)==TRUE),0)
								
                                if(is.na(deltamzminmax.tol)==FALSE){
                                    print("filtering by delta m/z")
                                    mz_min_max<-cbind(curdata[,2],curdata[,3])
                                    mz_min_max<-as.data.frame(mz_min_max)
                                    
                                    deltappm_res<-apply(mz_min_max,1,get_deltappm)
                                    
                                    curdata<-curdata[which(as.numeric(deltappm_res)<=deltamzminmax.tol),]
                                    feat.eval.result.mat<-feat.eval.result.mat[which(as.numeric(deltappm_res)<=deltamzminmax.tol),]
                                }
                                
                                feateval_list[[i]]<-feat.eval.result.mat
                                data_rpd_all[[i]]<-curdata
								
								  if(num_replicates>1)
								  {
									  print(paste("**calculating pairwise ",cormethod," correlation**",sep=""))

									  
									  
											rsqres_list<-evaluate.Samples(curdata, num_replicates, alignment.tool="XCMS", cormethod,missingvalue,ignore.missing)
											
											rsqres<-as.data.frame(rsqres_list$cor.matrix)
											
											curdata<-as.data.frame(rsqres_list$feature.table)
											rsqres<-as.data.frame(rsqres)
											snames<-colnames(curdata[,-c(1:8)])
											snames_1<-snames[seq(1,length(snames),num_replicates)]
											rownames(rsqres)<-snames_1
											pcor_outfile=paste(subdir1,"/",file_name,"_sampleassessment_usinggoodfeatures.txt",sep="")
											write.table(rsqres, pcor_outfile,sep="\t",row.names=TRUE)
											rsqres_list[[i]]<-rsqres
								  }
								  else
								  {
									  print("**skipping sample evaluataion as only one replicate is present**")
								  }
									
                                    #save(rsqres,file="rsqres.Rda")
								if(num_replicates>2)
								{
									rep_cor_mat<-cbind(rep_cor_mat,rsqres$meanCorrelation)
									bad_samples<-which(rsqres$meanCorrelation<samp.filt.thresh)
								}else
								{
									bad_samples<-which(rsqres<samp.filt.thresh)
                                    #rep_cor_mat<-cbind(rep_cor_mat,rsqres$Correlation)
                                    rep_cor_mat<-rsqres
                                }
								
								if(length(bad_samples)>0){
									bad_sample_names<-snames_1[bad_samples]
									
									feat.eval.outfile=paste(subdir1,"/",file_name,"_badsamples_at_cor",samp.filt.thresh,".txt",sep="")
									bad_sample_names<-as.data.frame(bad_sample_names)
									colnames(bad_sample_names)<-paste("Samples with correlation between technical replicates <", samp.filt.thresh,sep="")
									write.table(bad_sample_names, file=feat.eval.outfile,sep="\t", row.names=FALSE)
								}
								
								bad_list={}
								if(length(bad_samples)>0)
								{
									for(n1 in 1:length(bad_samples))
									{	
										if(bad_samples[n1]>1)
										{
											bad_samples[n1]=bad_samples[n1]+(bad_samples[n1]-1)*(num_replicates-1)
										}
											
									}
									for(n1 in 1:num_replicates)
									{
										bad_list<-c(bad_list,(bad_samples+n1-1))
									}
									bad_list<-bad_list[order(bad_list)]
									
								}
								if(i>1){
										parent_bad_list<-intersect(parent_bad_list,bad_list)
									}
									else{
									    
									    parent_bad_list<-bad_list

									}
								
								
							
                                  }
				  else
				  {
					  print("**skipping feature evaluataion as only one replicate is present**")
				  }
                                
	       	}
		
		cat("\n")
		 print("********Stage 2: Filtering results from each parameter setting based on sample and feature quality checks*******")
		 cat("\n")

		if(num_replicates>2)
		{
			max_rep_mean_cor<-apply(rep_cor_mat,1,max)
		}else{
			max_rep_mean_cor<-rep_cor_mat
		}
        
      
		 if(length(parent_bad_list)>0){
		
        
        if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
        {
					  sampleid_mapping<-sampleid_mapping[-c(parent_bad_list),]

				

					sampleidfilteredfile=paste(subdir2,"/","filtered_sampleid_mapping.txt",sep="")
					write.table(sampleid_mapping,file=sampleidfilteredfile,sep="\t",row.names=FALSE)
                    
        }

		}
		  for(i in 1:length(alignmentresults))
                          {
                                
				 
				 
				    file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
                                  feat.eval.file=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
                                  #data_rpd_all[[i]]=read.table(feat.eval.file,header=TRUE)
				  
				  curdata<-data_rpd_all[[i]]
				 
				  feat.eval.result.mat<-feateval_list[[i]]
				  if(length(parent_bad_list)>0){
					  
					 
					  curdata<-curdata[,-c(parent_bad_list+8)]
					 
					 
					  #maxint<-apply(curdata[,-c(1:4,((dim(curdata)[2]-6):dim(curdata)[2]))],1,max)
					  maxint<-apply(curdata[,-c(1:8)],1,max)
					  badfeats<-which(maxint==0)
					  if(length(badfeats)>0){
						curdata<-curdata[-c(badfeats),]
						
						
						
						  feat.eval.result.mat<- feat.eval.result.mat[-c(badfeats),]
						  feateval_list[[i]]<-feat.eval.result.mat
					  }
					
					}
				 data_rpd_all[[i]]<-curdata # [which(as.numeric(feat.eval.result.mat$median)<=feat.filt.thresh),]
				 
				 
				 
				 feat.eval.outfile=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
								
				 #write results
                 write.table(data_rpd_all[[i]], feat.eval.outfile,sep="\t", row.names=FALSE)
				 
				 feat.eval.outfile=paste(subdir1,"/",file_name,fileroot,"featureassessment.txt",sep="")					  
								#write results
								write.table(feateval_list[[i]], feat.eval.outfile,sep="\t", row.names=FALSE)
							
				 
			 }
                          ###########################################
                          #4) Merge two or more parameter settings
                          cat("\n")
                          print("*************Stage 3a: Merging features detected at different parameter settings********************")
                          cat("\n")
                          num_pairs=1
                          finalres={}
                          rnames={}
		         
			  if(merge.pairwise==TRUE){
                          for(i in 1:length(alignmentresults))
                          {
                                  file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
                                
                                  #feat.eval.file=paste(xMSanalyzer.outloc,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
				  #data_rpd_all[[i]]=read.table(feat.eval.file,header=TRUE)
                                  a1=sapply(strsplit(as.character(alignmentresults[i]), "\\_thresh"), head, n=2)[2]
                                  a2=sapply(strsplit(as.character(a1), "\\_"), head, n=2)
                                  threshval=as.numeric(a2[1])
                                  stepval=sapply(strsplit(as.character(a2[2]), "\\."), head, n=2)[2]
                                  stepval=paste(".",stepval,sep="")
                                  stepval=as.numeric(stepval)
                                  a1=sapply(strsplit(as.character(alignmentresults[i]), "\\_mzdiff"), head, n=2)[2]
                                  a2=sapply(strsplit(as.character(a1[1]), "\\_max"), head, n=2)
                                  mzdiff=as.numeric(a2[1])
                                  a2=sapply(strsplit(as.character(a2[2]), "\\.txt"), head, n=2)
                                  max=as.numeric(a2[1])
                                  p1=paste(threshval,"_",stepval,"_",mzdiff,"_",max,sep="")

								bool_num<-1
								
                                  for(j in i:length(alignmentresults))
                                  {
                                  	bool_num<-1
                                  	if(i==j){
                                  	if(length(alignmentresults)>1){
                                  		bool_num<-0
                                  	}
                                  	else{
                                  		bool_num<-1
                                  		}
                                  	}
				                	#if(i!=j)
				                	 if(bool_num==1)
				    				  {
				          file_name=sapply(strsplit(alignmentresults[j],".txt"),head)
					
                                          
                                          
                                         
                                         # if(i!=j)
                                          {
                                                  p1_p2=paste("p",i,"_U_","p",j,sep="")
                                          }
                                          #else
                                          #{
                                       #           p1_p2="p1"
                                          #}
					  
					   feat.eval.A<-feateval_list[[i]]
					 feat.eval.B<-feateval_list[[j]]
					 feat.eval.A<-feat.eval.A[which(as.numeric(feat.eval.A$median)<=feat.filt.thresh),]
					 feat.eval.B<-feat.eval.B[which(as.numeric(feat.eval.B$median)<=feat.filt.thresh),]
					 
					 #data_m=merge.Results(data_rpd_all[[i]],data_rpd_all[[j]],feateval[[i]],feateval[[j]],max.mz.diff,max.rt.diff,merge.eval.method,merge.eval.pvalue,alignment.tool="XCMS",
					
					print(paste("Number of good quality features from setting ",i,":", dim(data_rpd_all[[i]])[1],sep=": "))
					if(i!=j)
                                          {
					print(paste("Number of good quality features from setting ",j,":",dim(data_rpd_all[[j]])[1],sep=": "))
						}
                                         data_m=merge.Results(data_rpd_all[[i]],data_rpd_all[[j]],feat.eval.A,feat.eval.B,max.mz.diff,max.rt.diff, merge.eval.pvalue,alignment.tool="XCMS",
					numnodes=numnodes, mult.test.cor,mergecorthresh,missingvalue)
					data_m<-unique(data_m)
					
					numcols<-dim(data_m)[2]



					 data_m_int<-data_m[,-c(1:8,(numcols-8):numcols)]
					numsamps<-dim(data_m_int)[2]/num_replicates
					
                                           numcols<-dim(data_m)[2]

if(is.na(minexp.pct)==FALSE)
{
    minexp<-round(minexp.pct*dim(data_m_int)[2])
    
    if(length(which(data_m[,8]>=minexp))>0){
        data_m<-data_m[which(data_m[,8]>=minexp),]
    }else{
        stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
    }
    
    
}

					 data_m_int<-data_m[,-c(1:8,(numcols-8):numcols)]
					numsamps<-dim(data_m_int)[2]/num_replicates
					 
					maxint<-apply(data_m_int,1,function(x){max(x,na.rm=TRUE)})
					#numpeaks<-lapply(1:dim(data_m_int)[1],function(j){length(which(is.na(data_m_int[j,])==FALSE))})
				
					 
					 
					numpeaks<-apply(data_m_int,1,countpeaks)		
					 
					  union_list[[num_pairs]]<-data_m[,c(1:7)]
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],numpeaks)
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$numgoodsamples)
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$median)
					#Qscore<-(as.numeric(data_m$numgoodsamples)/as.numeric(data_m$median+0.1))
					#Qscore<-100*(Qscore/numsamps)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$Qscore)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],maxint)
					
			
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m_int)
					
					featinfo<-colnames(data_m[,c(1:7)])
					cnames<-colnames(data_m_int)
					merge.res.colnames<-c(featinfo,"NumPres.All.Samples","NumPres.Biological.Samples",paste("median",fileroot,sep=""),"Qscore","Max.Intensity",cnames)
					colnames(union_list[[num_pairs]])<-as.character(merge.res.colnames)
					
					
					
					  curres={}
					   curres=cbind(curres, p1_p2)
                                         
                      curres=cbind(curres, dim(union_list[[num_pairs]])[1])
                      curres=cbind(curres, mean(as.numeric(data_m$median)))
                      curres=cbind(curres, mean(as.numeric(data_m$Qscore)))
                      
                       #  curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$Qscore))))
                         curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$median))))
                                                           
                                          if(curscore>bestscore){
                                          	
                                          	bestscore<-curscore
                                          	best_i<-num_pairs
                                          	best_pair<-p1_p2
                                          }
                                          
                                          
                                          curres=cbind(curres,curscore)
                      
					  curres<-as.data.frame(curres)
                                          
                                          finalres=rbind(finalres,curres)
					  
					
                                             finalname=paste("XCMS_feature_list_at_", p1_p2,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
					  
                                          #Output merge results
                                          write.table(union_list[[num_pairs]],file=paste(subdir3,"/",finalname,sep=""), sep="\t",row.names=FALSE)
					  
							
								  num_pairs=num_pairs+1
                                     }
				  }
                          }
			  
			  }else{
			  
									data_rpd_all_parameters<-{}
							 feat.eval.all<-{}
							 p1_p2<-"pall"
								for(i in 1:length(alignmentresults))
								{
									 
										
													
									data_rpd_all_parameters<-rbind(data_rpd_all_parameters,data_rpd_all[[i]])
									
									feat.eval.A<-feateval_list[[i]]
									 feat.eval.A<-feat.eval.A[which(as.numeric(feat.eval.A$median)<=feat.filt.thresh),]
								
									 feat.eval.all<-rbind(feat.eval.all,feat.eval.A)
									
								}
								
								  
                                  
								cnames1<-colnames(feateval_list[[i]])

								
								
								
								feat.eval.all<-unique(feat.eval.all)
								 data_rpd_all_parameters<-unique(data_rpd_all_parameters)
								 
								 feat.eval.all<-as.data.frame(feat.eval.all)
								 data_rpd_all_parameters<-as.data.frame(data_rpd_all_parameters)
								 
                                 
                                 
								 
								data_m=merge.Results(data_rpd_all_parameters,data_rpd_all_parameters, feat.eval.all,feat.eval.all,max.mz.diff,max.rt.diff,merge.eval.pvalue,alignment.tool="XCMS",
								 numnodes=numnodes,mult.test.cor,mergecorthresh,missingvalue)
								 
								 
                                
                                 
								 numcols<-dim(data_m)[2]



					 data_m_int<-data_m[,-c(1:8,(numcols-8):numcols)]
					numsamps<-dim(data_m_int)[2]/num_replicates
					
                                           numcols<-dim(data_m)[2]

if(is.na(minexp.pct)==FALSE)
{
    minexp<-round(minexp.pct*dim(data_m_int)[2])
    
    if(length(which(data_m[,8]>=minexp))>0){
        data_m<-data_m[which(data_m[,8]>=minexp),]
    }else{
        stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
    }
    
    
}

					 data_m_int<-data_m[,-c(1:8,(numcols-8):numcols)]
					numsamps<-dim(data_m_int)[2]/num_replicates
					 
					maxint<-apply(data_m_int,1,function(x){max(x,na.rm=TRUE)})
					#numpeaks<-lapply(1:dim(data_m_int)[1],function(j){length(which(is.na(data_m_int[j,])==FALSE))})
				
					 
					 #numpeaks<-lapply(1:dim(data_m_int)[1],function(j){length(which(data_m_int[j,]>0))})
					
					 numpeaks<-apply(data_m_int,1,countpeaks)		
					  union_list[[num_pairs]]<-data_m[,c(1:7)]
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],numpeaks)
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$numgoodsamples)
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$median)
					#Qscore<-(as.numeric(data_m$numgoodsamples)/as.numeric(data_m$median+0.1))
					#Qscore<-100*(Qscore/numsamps)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$Qscore)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],maxint)
					
			
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m_int)
					
					featinfo<-colnames(data_m[,c(1:7)])
					cnames<-colnames(data_m_int)
					merge.res.colnames<-c(featinfo,"NumPres.All.Samples","NumPres.Biological.Samples",paste("median",fileroot,sep=""),"Qscore","Max.Intensity",cnames)
					colnames(union_list[[num_pairs]])<-as.character(merge.res.colnames)
					
					
					
					  curres={}
					   curres=cbind(curres, p1_p2)
                                         
                      curres=cbind(curres, dim(union_list[[num_pairs]])[1])
                      curres=cbind(curres, mean(as.numeric(data_m$median)))
                      curres=cbind(curres, mean(as.numeric(data_m$Qscore)))
                      
                       #  curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$Qscore))))
                         curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$median))))
                                                           
                                          if(curscore>bestscore){
                                          	
                                          	bestscore<-curscore
                                          	best_i<-num_pairs
                                          	best_pair<-p1_p2
                                          }
                                          
                                          
                                          curres=cbind(curres,curscore)
                      
					  curres<-as.data.frame(curres)
                                          
                                          finalres=rbind(finalres,curres)
					  
					
                                             finalname=paste("XCMS_feature_list_at_", p1_p2,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
					  
                                          #Output merge results
                                          write.table(union_list[[num_pairs]],file=paste(subdir3,"/",finalname,sep=""), sep="\t",row.names=FALSE)
					  
								
			  
			  }
                          finalres<-as.data.frame(finalres)
                        
			
            if(merge.pairwise==TRUE){
			  colnames(finalres)<-c("Parameter Combination", "Number of Features", "median PID/CV between sample replicates","mean Qscore (Quality score)","Parameter score")
                          write.table(finalres,file=paste(xMSanalyzer.outloc,"/xcms_with_xMSanalyzer_merge_summary.txt",sep=""), sep="\t", row.names=FALSE)
            }
			 
	print("Most optimal feature setting:")
    print(best_pair)
     cat("\n")
    #########################################################################
    
     cat("\n")
      print(paste("********Stage 3b: Generating final (pre-batcheffect correction) untargeted and targeted feature tables using ",best_pair," results******",sep=""))
      cat("\n")
    
    #rawQCeval/
   
  # pdf("Stage4a_QE_plots.pdf")
    #pdf(file=paste("subdir4a/Stage4a_QE_plots.pdf",sep=""))
    
    pdf(file=paste(subdir4a,"/QE_plots_All_samples_RAW.pdf",sep=""))
    
    #pdf(file=paste(xMSanalyzer.outloc,"/QE_plots_RAW.pdf",sep=""))
 	d1<-union_list[[best_i]]
    
    d1<-unique(d1)
    
  if(is.na(refMZ)==FALSE){
			stddata<-read.table(refMZ,sep="\t",header=TRUE)
			print(refMZ)
			print(head(stddata))
			
		}else{
			if(charge_type=="pos"){
			data(example_target_list_pos)
			stddata<-example_target_list_pos
			}else{
				
				if(charge_type=="neg"){
						data(example_target_list_neg)
						stddata<-example_target_list_neg
					}else{
						stop("Invalid option. \'charge_type\' should be \'pos\' or \'neg\'.")
						}
				}
		}
  
  d1_int<-round(d1[,-c(1:12)],0)
  
  
  rsqres_list<-evaluate.Samples(d1_int, num_replicates, alignment.tool=NA, cormethod,missingvalue,ignore.missing)
  
  rsqres<-as.data.frame(rsqres_list$cor.matrix)
  
  rsqres<-as.data.frame(rsqres)
  snames<-colnames(d1_int)
  snames_1<-snames[seq(1,length(snames),num_replicates)]
  rownames(rsqres)<-snames_1
  pcor_outfile=paste(subdir4a,"/Pairwise_Pearson_correlation_technical_replicates.txt",sep="")
  
  write.table(rsqres, pcor_outfile,sep="\t",row.names=TRUE)
  
  
  if(num_replicates>2)
  {
      max_rep_mean_cor<-apply(rsqres,1,max)
      
      max_rep_mean_cor<-as.data.frame(max_rep_mean_cor)
  }else{
      max_rep_mean_cor<-rsqres
      
      
  }
  
  
  d1<-cbind(d1[,c(1:12)],d1_int)
  rm(d1_int)

  
			Sys.sleep(1)
            
            max_ylim<-nrow(d1)+100
            
           
            
            num_features<-nrow(d1)
            
            
            par(mfrow=c(1,2))
            
            if(num_replicates>2){
                
                rep_cv<-d1$median_CV
                
                rep_cv_ylab<-"CV"
                
                h1<-hist(d1$median_CV,breaks=seq(0,max(d1$median_CV,na.rm=TRUE)+10,10),main=paste("Histogram median CV \n (using all ",num_features," features)",sep=""),col="brown",xlab="median CV%", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
                lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
            }else{
                
                rep_cv<-d1$median_PID
                
                rep_cv_ylab<-"PID"
                
                h1<-hist(d1$median_PID,breaks=seq(0,max(d1$median_PID,na.rm=TRUE)+10,10),main=paste("Histogram median CV \n (using all ",num_features," features)",sep=""),col="brown",xlab="median PID%", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
                lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
            }
            
            lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
            #pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main="Pie chart of median CVs (%) \n using all features",cex.main=0.7)
            
            if(num_replicates>2){
                pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main=paste("Pie chart of median CVs (%) \n using all ",num_features," features\n; average=",round(mean(d1$median_CV),2),sep=""),cex.main=0.7)
            }else{
                pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main=paste("Pie chart of median PIDs (%) \n using all ",num_features," features\n; average=",round(mean(d1$median_PID),2),sep=""),cex.main=0.7)
                
            }
            
            par(mfrow=c(1,1))
			
		hist(d1$NumPres.All.Samples,main="Histogram NumPres.All.Samples",col="brown",xlab="Number of samples (including replicates) \n with non-zero intensity values", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
			
			hist(d1$NumPres.Biological.Samples,main="Histogram NumPres.Biological.Samples",col="brown",xlab="Number of biological samples \n with non-zero intensity values", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
			#h1<-hist(d1$median_CV,main="Histogram median CV \n (using all data)",col="brown",xlab="median CV%", ylab="Number of features",cex.main=0.7)
			
			#h1<-hist(d1$median_CV,breaks=seq(0,max(d1$median_CV,na.rm=TRUE)+10,10),main="Histogram median CV \n (using all data)",col="brown",xlab="median CV%", ylab="Number of features",cex.main=0.7)
			
			
			par(mfrow=c(1,1))
            if(num_replicates>1){
                
                min_rep_pcor<-round(min(max_rep_mean_cor[,1],na.rm=TRUE),2)
                max_rep_pcor<-round(max(max_rep_mean_cor[,1],na.rm=TRUE),2)
                
                if(is.na(samp.filt.thresh)==FALSE){
                    
                    
                    parent_bad_list<-which(max_rep_mean_cor[,1]<samp.filt.thresh)
                    
                    if(length(parent_bad_list)>0){
                        
                        if(is.na(sample_info_file)==FALSE && sample_info_file!="NA"){
                            
                            filt_samp_names<-sampleid_mapping[parent_bad_list,1]
                            filt_samp_names<-paste(filt_samp_names,collapse=";")
                            
                            filt_samp_names<-length(parent_bad_list)
                        }
                    }else{
                        filt_samp_names<-"None"
                        
                    }
                }else{
                    filt_samp_names<-"None"
                    
                }

                
                
                hist(max_rep_mean_cor[,1],breaks=seq(0,1,0.1),main=paste("Histogram for mean Pearson correlation (min: ",min_rep_pcor,"; max: ",max_rep_pcor,") \n within technical replicates after filtering \n # of files filtered at threshold ",samp.filt.thresh,": ",filt_samp_names,sep=""),col="brown",xlab="mean replicate Correlation", ylab="Number of samples",cex.main=0.7)
                
                
            }
            

            
            

				
			
            
            par(mfrow=c(2,2))
            #5ppm
            data_a<-find.Unique.mzs.sameset(dataA=d1[,c("mz","time")],dataB=d1[,c("mz","time")],mz.thresh=5,time.thresh=NA,alignment.tool=NA)
            num_unique_features<-nrow(data_a$uniqueA)
            total_features<-nrow(d1)
            
            xlab2<-paste("overlapping m/z features \n based on +/- ",5,"ppm overlap criteria",sep="")
            
            temp_df<-cbind(total_features,num_unique_features)
            
            pie_v1<-round(temp_df[1,2]/temp_df[1,1],2)
            pie_v2<-1-pie_v1
            
            temp_dfA<-cbind(temp_df[1,1],temp_df[1,2], (temp_df[1,1]-temp_df[1,2]))#cbind(100*pie_v1,100*pie_v2)
            
            colnames(temp_dfA)<-c("Total","Unique","Overlapping")
            
            #barplot(temp_dfA,col="brown",ylab="Percentage (%) of total features", main=paste("Overlapping vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,100))
            
            barplot(temp_dfA,col="brown",ylab="Number of features", main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,temp_df[1,1]+100))
            #non-overlapping m/z features")
            
            # pie(c(pie_v2,pie_v1),labels=c("Overlapping","Unique"),col=c("orange","brown"),main=paste("Total vs ",xlab2,sep=""),cex.main=0.8)
            
            pie_label1<-paste("Overlapping \n(",100*pie_v2,"%)",sep="")
            pie_label2<-paste("Unique \n(",100*pie_v1,"%)",sep="")
            
            pie(c(pie_v2,pie_v1),labels=c(pie_label1,pie_label2),col=c("orange","brown"),main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8)
            
            #10ppm
            data_a<-find.Unique.mzs.sameset(dataA=d1[,c("mz","time")],dataB=d1[,c("mz","time")],mz.thresh=10,time.thresh=NA,alignment.tool=NA)
            num_unique_features<-nrow(data_a$uniqueA)
            total_features<-nrow(d1)
            
            #xlab2<-paste("Unique (non-overlapping) m/z features \n based on +/- ",10,"ppm overlap criteria",sep="")
            
            xlab2<-paste("overlapping m/z features \n based on +/- ",10,"ppm overlap criteria",sep="")
            
            
            temp_df<-cbind(total_features,num_unique_features)
            #colnames(temp_df)<-c("Total","Unique (non-overlapping)")
            
            pie_v1<-round(temp_df[1,2]/temp_df[1,1],2)
            pie_v2<-1-pie_v1
            
            temp_dfB<-cbind(temp_df[1,1],temp_df[1,2], (temp_df[1,1]-temp_df[1,2])) #cbind(100*pie_v1,100*pie_v2)
            
            colnames(temp_dfB)<-c("Total","Unique","Overlapping")
            
            #barplot(temp_dfB,col="brown",ylab="Percentage (%) of total features", main=paste("Total vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,100)) #non-overlapping m/z features")
            
            barplot(temp_dfB,col="brown",ylab="Number of features", main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,temp_df[1,1]+100))
            
            
            pie_label1<-paste("Overlapping \n(",100*pie_v2,"%)",sep="")
            pie_label2<-paste("Unique \n(",100*pie_v1,"%)",sep="")
            pie(c(pie_v2,pie_v1),labels=c(pie_label1,pie_label2),col=c("orange","brown"),main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8)
            
            par(mfrow=c(1,1))
            
            
            
            hist(d1$Qscore,main="Histogram Qscore",col="brown",xlab="Quality score", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
            
            
			#pie(h1$counts,col=rainbow(length(h1$counts)),labels=h1$breaks,main="Pie chart of median CVs (%) \n using all features",cex.main=0.7)
			
			#pie(h1$counts,col=rainbow(length(h1$counts)),labels=h1$breaks,main="Pie chart of median CVs (%) \n using all features",cex.main=0.7)
			
            
			
			d2<-d1[order(d1$time),]
			
            plot(d2$time,d2$mz,main="m/z vs Time",col="brown",xlab="Time (s)",ylab="m/z",cex.main=0.7)
            
            
            par(mfrow=c(1,1))
            
            plot(d2$mz,d2$Max.Intensity,main="Intensity vs m/z",col="brown",xlab="m/z",ylab="Max intensity \nacross all samples/profiles",cex.main=0.7)
            plot(d2$time,d2$Max.Intensity,main="Intensity vs time",col="brown",xlab="Time (s)",ylab="Max intensity \nacross all samples/profiles",cex.main=0.7)
            
            plot(rep_cv,d2$Max.Intensity,main=paste("Intensity vs ",rep_cv_ylab,sep=""),col="brown",xlab=rep_cv_ylab,ylab="Max intensity \nacross all samples/profiles",cex.main=0.7)
            
            
            par(mfrow=c(1,1))
            
            plot(d2$mz,d2$Qscore,main="Qscore vs m/z",col="brown",xlab="m/z",ylab="Qscore",cex.main=0.7)
            plot(d2$time,d2$Qscore,main="Qscore vs Time",col="brown",xlab="Time (s)",ylab="Qscore",cex.main=0.7)
            plot(rep_cv,d2$Qscore,main=paste("Qscore vs ",rep_cv_ylab,sep=""),col="brown",xlab=rep_cv_ylab,ylab="Qscore",cex.main=0.7)
            
            
			
            par(mfrow=c(1,1))
            
	max_numzeros<-dim(d1)[2]*1

	if(is.na(void.vol.timethresh)==FALSE){
		
        dfirst15<-d1[which(d1$time<void.vol.timethresh),]
        
        
        
        if(nrow(dfirst15)>1){
            
            
            ind1<-which(dfirst15$Max.Intensity==max(dfirst15$Max.Intensity))[1]
            
            time_thresh<-dfirst15$time[ind1]
            
            time_thresh<-time_thresh-(0.30*time_thresh)
            
            time_thresh<-round(time_thresh,1)
            
            plot(dfirst15$time,dfirst15$Max.Intensity,xlab="Time (s)", col="brown",ylab="Max intensity \nacross all samples", main=paste("Estimated void volume time: ",time_thresh," s",sep=""))
            abline(v=time_thresh,col=4,lty=3)
            
            
            
            d1<-d1[which(d1$time>=time_thresh),]
            
            print("Estimated void volume time cutoff")
            print(time_thresh)
            
            
        }else{
            
            print("No features eluting before void volume time cutoff")
        }
        

        
        
		
	
		#sfname<-paste(xMSanalyzer.outloc,"/stage5/feature_table_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtered.txt",sep="")
		finalname<-paste("featuretable_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,"_voidtimefilt.txt",sep="")
		
        #write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
	   

	}else{
		
		
		finalname<-paste("feature_table_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
		
        #write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
	   

		
		}



		cat("\n")
     print(paste("*********Stage 4a: Search for target features, performing QC evaluation using ",best_pair," results********",sep=""))
     cat("\n")
     
     Sys.sleep(1)
	
		
		print("Dim data after void time filtering")
		print(dim(d1))
        
        
        targeted_feat_raw<-eval.target.mz(dataA=d1[,-c(2:3,5:12)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="raw",xMSanalyzer.outloc=subdir4a)
        
	   raw_mz_time<-d1[,c("mz","time")]
	    
        write.table(raw_mz_time,file=paste(subdir3b,"/Precalibration_mz_time.txt",sep=""),sep="\t",row.names=FALSE)
        
        
        median_error<-median(targeted_feat_raw$delta_ppm_error,na.rm=TRUE)
        
	  median_time_error<-median(targeted_feat_raw$delta_time_error,na.rm=TRUE)
	
        #d1<-get_calibrated_mz_data(dataA=d1,delta_ppm_error=median_error) #,refMZ=stddata,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff)
        
	#  d1<-get_calibrated_mz_data(dataA=d1,delta_ppm_error=median_error,delta_time_error=median_time_error)
	    #here1
           # d1<-get_calibrated_mz_data(dataA=d1,delta_ppm_error=median_error,delta_time_error=median_time_error,refMZ=stddata,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff)
	   
	   
	     d1<-get_calibrated_mz_data(dataA=d1,refMZ=stddata,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,feature.eval.result=feat.eval.result,calibration.method=calibration.method[1])
            
        finalname<-paste("RAW_mzcalibrated_untargeted_featuretable.txt",sep="")
        
#        write.table(d1,file=finalname,sep="\t",row.names=FALSE)
write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)

#write.table(d1,file=paste(xMSanalyzer.outloc,finalname,sep=""),sep="\t",row.names=FALSE)



 targeted_calibfeat_raw<-eval.target.calibratedmz(dataA=d1[,-c(2:3,5:12)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="raw",xMSanalyzer.outloc=subdir4a)
        

	data_m<-d1[,-c(1:12)]
   
    if(replacezeroswithNA==TRUE){
		data_m<-replace(as.matrix(data_m),which(data_m==0),NA)
		d1<-cbind(d1[,c(1:12)],data_m)
	}
	
	
	counna<-apply(data_m,1,function(x){length(which(is.na(x)==TRUE))})

	maxzeros<-1*dim(data_m)[2]
	
	if(length(which(counna<maxzeros))){
		data_m<-data_m[which(counna<maxzeros),]
	}

    maxint<-apply(data_m,1,function(x){max(x,na.rm=TRUE)})
   	maxint_ord<-order(maxint,decreasing=TRUE)
   	#[maxint_ord[1:5000]
   	
    X<-t(data_m) #[maxint_ord[1:2000],])
    X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
    
    tic.eval(d1[,-c(1:12)],outloc=subdir4a)
    
    #save(d1,file="d1.Rda")
    feat.eval.result=evaluate.Features(d1[,-c(9:12)], numreplicates=num_replicates,min.samp.percent=min.samp.percent,alignment.tool="XCMS",impute.bool=impute.bool,numnodes=numnodes)
			
    		
    Sys.sleep(1)
    s1<-"Stage 1 results: QC evaluation of invidual parameters from XCMS"
s2<-"Stage 2 results: filtered results from each paramter setting based on sample and feature quality (CV within replicates) checks"
s3<-"Stage 3a results: merged results using stage 2 filtered files"
s3b<-"Stage 3b results: RAW m/z calibrated untargeted and targeted feature tables (averaged and non-averaged)"
s4a<-"Stage 4a results: QC evaluation of targeted and untargeted data before batch-effect correction"

sm<-rbind(s1,s2,s3,s3b,s4a)

 batch_levels<-{}
#targeted_feat_raw<-eval.target.mz(dataA=d1[,-c(3:12)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="raw",xMSanalyzer.outloc=xMSanalyzer.outloc)


if(summarize.replicates==TRUE && data.norm.pipeline=="AC"){
    
    if(num_replicates>1){
    sampleid_mapping<-sampleid_mapping[seq(1,dim(d1[,-c(1:12)])[2],num_replicates),]
    
    data_matrix<-data_summarize(Xmat=d1[,- c(3:12)],Ymat=NA,feature_table_file=NA,parentoutput_dir=subdir3b,class_labels_file=NA,num_replicates=num_replicates,summarize.replicates=summarize.replicates,
    summary.method=summary.method,missing.val=missingvalue, rep.num.max.missing.thresh=rep.num.max.missing.thresh,summary.na.replacement=summary.na.replacement)
    
    print(dim(data_matrix))
    
    d1<-cbind(d1[,c(1:12)],data_matrix[,-c(1:2)])
   	 minexp=minexp/num_replicates 
   
    X<-t(data_matrix[,-c(1:2)])
    data_m<-data_matrix[,-c(1:2)]
    }
    num_replicates=1
}

   # if(is.na(sample_info_file)==FALSE)
   if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
    {
        #sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)
    	
    	sampleid_mapping[,1]<-gsub(sampleid_mapping[,1],pattern=filepattern,replacement="",ignore.case=TRUE)
		cnames<-colnames(data_m)
		
		cnames<-gsub(cnames,pattern=filepattern,replacement="",ignore.case=TRUE)
		colnames(data_m)<-cnames


    	batch_inf<-sampleid_mapping[,3]
    	
    	batch_labels<-as.factor(sampleid_mapping[,3])
    
    	l1<-levels(batch_labels)
     batch_levels<-levels(batch_labels)
    	    	   

		    try(pca.eval(X=X,samplelabels=batch_labels,filename="raw",ncomp=5,center=TRUE,scale=TRUE,legendlocation="bottomleft",legendcex=0.5,outloc=subdir4a),silent=TRUE)
			 
			Sys.sleep(1)
			
			
										

			cnames=colnames(feat.eval.result)
		
			

			if(dim(sampleid_mapping)[2]<4){
   		    mod<-rep(1,dim(data_m)[2])
   		    }else{
   		    		mod<-sampleid_mapping[,-c(1:3)]
   		    	}
    
  	
    		try(dev.off(),silent=TRUE)

    		Sys.sleep(1)
    #matchedfilter
	if(length(batch_levels)>1){
		     #######################################################################
		     cat("\n")
		    print(paste("**********Stage 4b: Processing data using ComBat and post-correction QC evaluation using ",best_pair," results*******",sep=""))
		      cat("\n")
     #pdf("Stage4b_QE_plots.pdf")
     #pdf(file=paste("subdir4b/Stage4b_QE_plots.pdf",sep=""))
     
    
     
                if(summarize.replicates==FALSE){
            pdf(file=paste(subdir4b,"/QE_plots_All_samples_ComBat.pdf",sep=""))
           }else{
               
	        if(summarize.replicates==TRUE && data.norm.pipeline=="AC"){
			if(summary.method=="mean"){
			pdf(file=paste(subdir4b,"/QE_plots_All_samples_1averaged_2ComBat.pdf",sep=""))
			}else{
				pdf(file=paste(subdir4b,"/QE_plots_All_samples_1mediansummarized_2ComBat.pdf",sep=""))
			}
		}else{
		
			pdf(file=paste(subdir4b,"/QE_plots_All_samples_ComBat.pdf",sep=""))
		
		}
		
		
               
           }
     
     # pdf(file=paste(xMSanalyzer.outloc,"/QE_plots_ComBat.pdf",sep=""))
		    ##################################################################
		
		    #adjdata<-try(sva::ComBat(dat=data_m,batch=batch_inf,mod=mod,par.prior=TRUE),silent=TRUE)
		    
		#     if(is.na(missingvalue)==FALSE){
		#	data_m_na=replace(as.matrix(data_m),which(data_m==missingvalue),NA)
		#   }else{
		   #data_m_na<-data_m
		#   }
		    if(is.na(missingvalue)==FALSE){
			
			  if(replacezeroswithNA==TRUE){
					#data_m<-replace(as.matrix(data_m),which(data_m==0),NA)
					data_m_na=replace(as.matrix(data_m),which(data_m==missingvalue),NA)
					d1<-cbind(d1[,c(1:12)],data_m_na)
					
					sum_check<-apply(data_m_na,1,function(x){length(which(is.na(x)==FALSE))})
					
					
					if(length(which(sum_check<minexp))>0){
					data_m_na<-data_m_na[-which(sum_check<minexp),]
					d1<-d1[-which(sum_check<minexp),]
					
					}
					
					 if(summarize.replicates==TRUE){
					 print(paste("Number of features after minexp ",minexp," filtering post-averaging",sep=""))
                                         print(dim(data_m)[1])
					 }
					
					}else{
						data_m_na<-data_m
					}

			}else{
					data_m_na<-data_m
					sum_check<-apply(data_m_na,1,function(x){length(which(is.na(x)==FALSE))})
					
					if(length(which(sum_check<minexp))>0){
						data_m_na<-data_m_na[-which(sum_check<minexp),]
						d1<-d1[-which(sum_check<minexp),]
					}
					
					if(summarize.replicates==TRUE){
					 print(paste("Number of features after minexp ",minexp," filtering post-averaging",sep=""))
                                         print(dim(data_m)[1])
					 }
					
		   }
		   
		    adjdata<-try(sva::ComBat(dat=data_m_na,batch=batch_inf,mod=mod,par.prior=TRUE),silent=TRUE)
		    
		    
		    if(is(adjdata,"try-error")){
		 
		
				data_m1<-cbind(d1[,c(1:4)],data_m_na)
				
			
		    	adjdata<-MetabComBat(dat=data_m1,saminfo=sampleid_mapping,par.prior=T,filter=F,write=F,prior.plots=F)
		    	adjdata<-adjdata[,-c(1:4)]

#	print("Done with MetabComBat")
#print(dim(adjdata))
#				print(dim(d1))
		    }
    
    				    
		    maxint<-apply(adjdata,1,function(x){max(x,na.rm=TRUE)})
		    
		    adjdata2<-replace(as.matrix(adjdata),which(is.na(adjdata)==TRUE),missingvalue)
		    adjdata2<-replace(as.matrix(adjdata2),which(adjdata2<0),0)
		    
		    adjdata2<-cbind(d1[,c(1:8)],adjdata2)
		    
            #save(adjdata2,file="adjdata2.Rda")
		    
            #expression_xls<-paste(subdir4b,"/ComBat_corrected_",best_pair,"_feature_table.txt",sep="")
		    
            expression_xls<-paste(xMSanalyzer.outloc,"/ComBat_mzcalibrated_untargeted_featuretable.txt",sep="")
 
		 if(summarize.replicates==TRUE){
    
     if(summary.method=="mean"){
        expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_averaged_featuretable.txt",sep="")
	}else{
		if(summary.method=="median"){
			expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt",sep="")
		}
	
	}
 }else{
     expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_featuretable.txt",sep="")
 }
	
		    feat.eval.result=evaluate.Features(adjdata2, numreplicates=num_replicates,min.samp.percent=min.samp.percent,alignment.tool="XCMS",impute.bool=impute.bool,numnodes=numnodes)
		    
			cnames=colnames(feat.eval.result)
										
			feat.eval.result<-apply(feat.eval.result,2,as.numeric)
			feat.eval.result<-as.data.frame(feat.eval.result)
			#feat.eval.result.mat=cbind(adjdata2[,c(1:4)],feat.eval.result)  
			
			numpeaks<-apply(adjdata2[,-c(1:8)],1,countpeaks)						
		    
		    maxint<-apply(adjdata2[,-c(1:8)],1,function(x){max(x,na.rm=TRUE)})
		    
		    adjdata2<-cbind(adjdata2[,c(1:7)],numpeaks,feat.eval.result$numgoodsamples,feat.eval.result$median,feat.eval.result$Qscore,maxint,adjdata2[,-c(1:8)])
		    
		    colnames(adjdata2)<-colnames(d1)
		 
		    adjdata2$time<-round(adjdata2$time,1)
                    adjdata2$mz<-round(adjdata2$mz,4)  
                    adjdata2$mzmin<-round(adjdata2$mzmin,5)
                    adjdata2$mzmax<-round(adjdata2$mzmax,5)
		    adjdata2$rtmin<-round(adjdata2$rtmin,1)
                    adjdata2$rtmax<-round(adjdata2$rtmax,1)
 
 
 if(summarize.replicates==TRUE){
    
     if(summary.method=="mean"){
        expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_averaged_featuretable.txt",sep="")
	}else{
		if(summary.method=="median"){
			expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt",sep="")
		}
	
	}
 }else{
     expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_featuretable.txt",sep="")
 }
 
 #expression_xls<-paste(xMSanalyzer.outloc,"/ComBat_mzcalibrated_untargeted_featuretable.txt",sep="")
 
		    write.table(adjdata2,file=expression_xls,sep="\t",row.names=FALSE)
		    
		   # X<-t(adjdata2[maxint_ord[1:2000],-c(1:9)])
		   
		    tic.eval(adjdata2[,-c(1:12)],outloc=subdir4b)
    

		    X<-t(adjdata2[,-c(1:12)])
		   
		     
		    X<-as.matrix(X)
		    
		    X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
		    
		    try(pca.eval(X,samplelabels=batch_labels,filename="ComBat",ncomp=5,center=TRUE,scale=TRUE,legendlocation="topright",legendcex=0.5,outloc=subdir4b),silent=TRUE)
			    
		    
		    #eval.reference.metabs(dataA=adjdata2,stdData=stddata,mzthresh=10,outloc=getwd(),folderheader="ComBat")
			
			targeted_feat_combat<-eval.target.calibratedmz(dataA=adjdata2[,-c(2:3,5:12)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4b,folderheader="ComBat",xMSanalyzer.outloc=subdir4b)
			
			s4b<-"Stage 4b results: Batch-effect evaluation post ComBat and QC evaluation of targeted data post batch-effect correction"
		sm<-rbind(sm,s4b)
			
				
	dev.off()
	
		}else{
			
			adjdata2=d1
			targeted_feat_combat<-targeted_feat_raw
			
			adjdata2<-replace(as.matrix(adjdata2),which(is.na(adjdata2)==TRUE),missingvalue)
			adjdata2<-replace(as.matrix(adjdata2),which(adjdata2<0),0)
		    
			
			
			}
    
}
			 
			 
			 
			  	 metlin.res={}
                 kegg.res={}
                      
					 #length(union_list[[num_pairs]]$mz
					 if(is.na(adduct.list)==FALSE){
					
					cat("\n")
					  print("*********Stage 5: Mapping m/z values to known metabolites*********")
					 cat("\n")
					 
                     #annot.res<-feat.batch.annotation.KEGG(d1,mz.tolerance.dbmatch,adduct.list,subdir5, numnodes=numnodes,syssleep=syssleep)
                        
                        
                        for(db in db_name){
                            annot.res<-simpleAnnotation(dataA=adjdata2,max.mz.diff=mz.tolerance.dbmatch,num_nodes=numnodes,queryadductlist=adduct.list,
                            gradienttype="Acetonitrile",mode=charge_type,outloc=subdir5,db_name=db)
                            
                            #fname<-paste("DBmatches_",db,".txt",sep="")

                                fname<-paste(subdir5,"/DBmatches_",db,".txt",sep="")
                            write.table(annot.res,file=fname,sep="\t",row.names=FALSE)
                            
                        }
                        
					    annotres_list<-annot.res
					    
					     s5<-"Stage 5 results: Annotation of features"
					    sm<-rbind(sm,s5)

                     }
				
			 
			 
			  print("*************Processing complete**********")

			  #print("*********Characterizing metabolites*********")
		}             
                else
                {
                        stop(paste("No files exist in",XCMS.outloc, "Please check the input value for cdfloc", sep=""))
                }
                    
                
                
        }
        
        suppressWarnings(sink(file=NULL))
        #sm<-"The results include pairwise Pearson correlation within technical replicates; 2) file with raw m/z and time for each feature; 3) raw (non batch-effect corrected) feature table after m/z calibration; 4) ComBat processed feature table after m/z calibration;  5) Simple (m/z based) database searches"

sm<-as.data.frame(sm)
colnames(sm)<-"Output_description"
setwd(xMSanalyzer.outloc)
write.table(sm,file="Readme.txt",sep="\t",row.names=FALSE)

		print(paste("Processing is complete. Program results can be found at: ",xMSanalyzer.outloc,sep=""))
		
	 return(list("XCMS.merged.res"=union_list, "XCMS.ind.res"=data_rpd_all,"XCMS.ind.res.filtered"=data_rpd, "final.feat.table.annot"=annotres_list, "feat.eval.ind"=feateval_list, "sample.eval.ind"=rsqres_list,"final.feat.table.raw"=d1,"final.feat.table.combat"=adjdata2,
	 "final.targeted.feat.table.raw"=targeted_feat_raw,"final.targeted.feat.table.combat"=targeted_feat_combat))
}
                                           
xMSwrapper.XCMS.centWave<-function(cdfloc, XCMS.outloc,xMSanalyzer.outloc,ppm.list=c(2.5), mz.diff.list=c(-0.00005), sn.thresh.list=c(3), prefilter.list=c(3,1000), bw.val=c(5,10),groupval.method="medret",
step.list=c(0.1,1),max=50,minfrac.val=0.1, minsamp.val=2, mzwid.val=0.015, sleep.val=0,retcor.method="obiwarp",retcor.family="symmetric", retcor.plottype="deviation", peakwidth=c(10,60),
numnodes=2,run.order.file=NA,max.mz.diff=5,max.rt.diff=300, merge.eval.pvalue=0.05,mergecorthresh=0.7,deltamzminmax.tol=100,
num_replicates=3,subs=NA, mz.tolerance.dbmatch=10, adduct.list=c("M+H"), samp.filt.thresh=0.70,
feat.filt.thresh=50,cormethod="pearson",mult.test.cor=TRUE,missingvalue=0,ignore.missing=TRUE,
sample_info_file=NA,refMZ=NA,refMZ.mz.diff=10,refMZ.time.diff=NA,void.vol.timethresh=30,
replacezeroswithNA=TRUE,scoreweight=30,filepattern=".cdf",charge_type="pos", minexp.pct=0.1,syssleep=0.5,merge.pairwise=FALSE,min.samp.percent=0.6,impute.bool=TRUE,summarize.replicates=TRUE,
summary.method="median",max.number.of.replicates.with.missingvalue=1,summary.na.replacement="zeros",db_name=c("KEGG","HMDB","LipidMaps"),qc_label=NA,data.norm.pipeline="AC",calibration.method=c("median.adjustment","multiplicative.signal.correction"))
{
	suppressWarnings(sink(file=NULL))
    dir.create(xMSanalyzer.outloc,showWarnings=FALSE)

    targeted_feat_raw<-{}
	library(parallel)
    x<-date()
    x<-strsplit(x,split=" ")
    x1<-unlist(x)
    rep.num.max.missing.thresh=max.number.of.replicates.with.missingvalue
    
    #cdfloc<-NA
    
    inp_params<-{}
    inp_params<-rbind(inp_params,cbind("cdfloc: ",cdfloc))
    inp_params<-rbind(inp_params,cbind("XCMS.centWave.outloc: ",XCMS.outloc))
    inp_params<-rbind(inp_params,cbind("xMSanalyzer.outloc: ",xMSanalyzer.outloc))
    inp_params<-rbind(inp_params,cbind("num_replicates: ",num_replicates))
    inp_params<-rbind(inp_params,cbind("max.mz.diff: ",max.mz.diff))
    inp_params<-rbind(inp_params,cbind("max.rt.diff: ",max.rt.diff))
    inp_params<-rbind(inp_params,cbind("merge.eval.pvalue: ",merge.eval.pvalue))
    inp_params<-rbind(inp_params,cbind("mergecorthresh: ",mergecorthresh))
    inp_params<-rbind(inp_params,cbind("deltamzminmax.tol: ",deltamzminmax.tol))
    inp_params<-rbind(inp_params,cbind("mz.tolerance.dbmatch: ",mz.tolerance.dbmatch))
    inp_params<-rbind(inp_params,cbind("adduct.list: ",paste(adduct.list,collapse=";")))
    inp_params<-rbind(inp_params,cbind("samp.filt.thresh: ",samp.filt.thresh))
    inp_params<-rbind(inp_params,cbind("feat.filt.thresh: ",feat.filt.thresh))
    inp_params<-rbind(inp_params,cbind("cormethod: ",cormethod))
    inp_params<-rbind(inp_params,cbind("charge_type: ",charge_type))
    inp_params<-rbind(inp_params,cbind("db_name: ",paste(db_name,collapse=";")))
    inp_params<-rbind(inp_params,cbind("impute.bool: ",impute.bool))
    inp_params<-rbind(inp_params,cbind("merge.pairwise: ",merge.pairwise))
    
    inp_params<-rbind(inp_params,cbind("summarize.replicates: ",summarize.replicates))
    inp_params<-rbind(inp_params,cbind("summary.method: ",summary.method))
    inp_params<-rbind(inp_params,cbind("max.number.of.replicates.with.missingvalue: ",rep.num.max.missing.thresh))
    inp_params<-rbind(inp_params,cbind("summary.na.replacement: ",summary.na.replacement))

    
    fname<-paste(xMSanalyzer.outloc,"/Input_parameters.csv",sep="")
    
    colnames(inp_params)<-c("Parameter","Value")
    write.csv(inp_params,file=fname,row.names=FALSE)
    
    #variables not used
	#summarize.replicates=FALSE
    #summary.method="median"
    #rep.num.max.missing.thresh=1
    #summary.na.replacement="zeros"

    
    x1<-gsub(x1,pattern=":",replacement="_")
    fname<-paste(x1[2:4],collapse="")
    fname<-gsub(fname,pattern=":",replacement="_")
    #fname<-paste(fname,x1[5],sep="")
    x1[5]<-gsub(x1[5],pattern=":",replacement="_")
    fname<-paste(fname,x1[5],sep="_")
        
        
       
	

	fname<-paste(xMSanalyzer.outloc,"/Log_",fname,".txt",sep="")
	
	fname<-paste(xMSanalyzer.outloc,"/Log.txt",sep="")
	
	#print(paste("Program running. Please check the logfile for runtime status: ",fname,sep=""))
print(paste("Program is running. Please check the logfile for runtime status: ",fname,sep=""))


	data_rpd_all=new("list")
	data_rpd=new("list")
	union_list=new("list")
	feateval_list<-new("list")
	rsqres_list<-new("list")
	annot.res<-{}
	annotres_list<-new("list")
	if(is.na(cdfloc)==FALSE){
	
		if(cdfloc=="NA"){
			cdfloc=NA
		}
	}
	
	if(is.na(sample_info_file)==FALSE){
	
		if(sample_info_file=="NA"){
			sample_info_file=NA
		}
	}
	print(is.na(cdfloc))
	print(cdfloc)
	
	print(length(subs))
	
	
filepattern=".cdf$|.mzxml$|mXML$"
if(is.na(sample_info_file)==FALSE)
    {
    	
    	
		sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)

		if(is.na(cdfloc)==FALSE){
			l1<-list.files(cdfloc,filepattern,ignore.case=TRUE)
	
			#minexp<-round(minfrac.val*length(l1))
			
			if(is.na(subs)==TRUE){
			
			
			minexp<-round(minfrac.val*length(l1))
			}else{
			
				l1<-l1[subs]
				
				print(l1)
				print(length(l1))
				
				minexp<-round(minfrac.val*length(l1))
			}
			
			
			}else{
				l1<-rep(1,num_replicates)	
			}
		if(length(l1)!=dim(sampleid_mapping)[1] & (is.na(cdfloc)==FALSE))
		{
			num_mis_files<-dim(sampleid_mapping)[1]-length(l1)
			if(is.na(subs)==TRUE){
				stop(paste("ERROR: Only ",length(l1)," spectral files were found. ",num_mis_files," files are missing.",sep=""))
			}
		}
	}else{
	
	
		if(is.na(cdfloc)==FALSE){
			l1<-list.files(cdfloc,filepattern,ignore.case=TRUE)
			if(is.na(subs)==TRUE){
			
			
			minexp<-round(minfrac.val*length(l1))
			}else{
			
				l1<-l1[subs]
				minexp<-round(minfrac.val*length(l1))
			}
		}else{
			l1<-rep(1,num_replicates)	
		}
	}
if(length(l1)%%num_replicates>0)
{stop(paste("ERROR: Not all samples have ",num_replicates," replicates.",sep=""))
}







        ############################################
        #1) Align profiles using the cdf.to.ftr wrapper function in apLCMS
	 if(is.na(XCMS.outloc)==TRUE)
        {
                stop("Undefined value for parameter, XCMS.outloc. Please define the XCMS output location.")

        }
         if(is.na(xMSanalyzer.outloc)==TRUE)
        {
                stop("Undefined value for parameter, xMSanalyzer.outloc. Please define the xMSanalyzer output location.")

        }

dir.create(XCMS.outloc,showWarnings=FALSE)
	dir.create(xMSanalyzer.outloc,showWarnings=FALSE)
	

	sink(fname)
	print(sessionInfo())

	if(is.na(refMZ)==FALSE){
                        stddata<-read.table(refMZ,sep="\t",header=TRUE)
                        print(refMZ)
                        print(head(stddata))

                }else{
                        if(charge_type=="pos"){
                        data(example_target_list_pos)
                        stddata<-example_target_list_pos
                        }else{

                                if(charge_type=="neg"){
                                                data(example_target_list_neg)
                                                stddata<-example_target_list_neg
                                        }else{
                                                stop("Invalid option. \'charge_type\' should be \'pos\' or \'neg\'.")
                                                }
                                }
                }	
        if(is.na(cdfloc)==FALSE)
        {
                setwd(cdfloc)
                if(is.na(XCMS.outloc)==FALSE)
                {
                   #     data_rpd_all=XCMS.align.centWave(cdfloc, XCMS.outloc,ppm.list, mz.diff.list, sn.thresh.list, prefilter.list, bw.val,groupval.method, 
			#step.list,max,minfrac.val, minsamp.val, mzwid.val, sleep.val, run.order.file,subs, retcor.method,retcor.family, retcor.plottype, peakwidth)


data_rpd_all<-XCMS.align.centWave(cdfloc, XCMS.outloc,ppm.list=ppm.list, mz.diff.list=mz.diff.list, sn.thresh.list=sn.thresh.list, prefilter.list=prefilter.list,
bw.val=bw.val,groupval.method=groupval.method, 
step.list=step.list,max=max,minfrac.val=minfrac.val, minsamp.val=minsamp.val, mzwid.val=mzwid.val, sleep.val=sleep.val, run.order.file,subs, retcor.method=retcor.method,
retcor.family=retcor.family, retcor.plottype=retcor.plottype, peakwidth=peakwidth,target.mz.list=stddata,nSlaves=numnodes,xMSanalyzer.outloc=xMSanalyzer.outloc)

                       
                }
                else
                {
                        stop("Undefined value for parameter, XCMS.outloc. Please define the output location.")
                }
        }
	


		setwd(XCMS.outloc)
                alignmentresults<-list.files(XCMS.outloc, "*.txt")
		print("Files found in XCMS output location:")
                print(alignmentresults)
		for(i in 1:length(alignmentresults))
		{
		
			data_rpd_all[[i]]<-read.table(paste(XCMS.outloc,"/",alignmentresults[i],sep=""),header=TRUE)
			data_rpd_all[[i]]<-unique(data_rpd_all[[i]]) 	
			
			
		}
	
		if(is.na(sample_info_file)==FALSE)
		{
			 match_names_check<-match(sampleid_mapping[,1],colnames(data_rpd_all[[1]][,-c(1:8)]))
			 
			 if(length(which(is.na(match_names_check)==TRUE))>0){
				stop("Sample names do not match between sequence file and feature table.")
			 }
	
		}
	#subdir1<-paste(xMSanalyzer.outloc,"/Quality_assessment_files",sep="")
	#subdir2<-paste(xMSanalyzer.outloc,"/XCMS_filtered_data",sep="")
	#subdir3<-paste(xMSanalyzer.outloc,"/XCMS_with_xMSanalyzer_merged_data",sep="")
	#dir.create(subdir1,showWarnings=FALSE)
	#dir.create(subdir2,showWarnings=FALSE)
	#dir.create(subdir3,showWarnings=FALSE)
	
	
    subdir1<-paste(xMSanalyzer.outloc,"/Stage1",sep="")  #QC individual parameter settings
    subdir2<-paste(xMSanalyzer.outloc,"/Stage2",sep="")  #Data filtering
    subdir3<-paste(xMSanalyzer.outloc,"/Stage3a",sep="")  #Data merger/parameter optimization
    subdir3b<-paste(xMSanalyzer.outloc,"/Stage3b",sep="")
    subdir4a<-paste(xMSanalyzer.outloc,"/Stage4a",sep="")	 #Raw QC: batch effect eval, TIC, etc
    
    
    dir.create(subdir1,showWarnings=FALSE)
    dir.create(subdir2,showWarnings=FALSE)
    dir.create(subdir3,showWarnings=FALSE)
    dir.create(subdir3b,showWarnings=FALSE)
    dir.create(subdir4a,showWarnings=FALSE)
    
    #if(is.na(sample_info_file)==FALSE)
	   if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
       {
           subdir4b<-paste(xMSanalyzer.outloc,"/Stage4b",sep="")	 #Batch-effect corrected QC: batch effect eval, TIC, etc
           dir.create(subdir4b,showWarnings=FALSE)
           
       }
       
       if(is.na(adduct.list)==FALSE){
           subdir5<-paste(xMSanalyzer.outloc,"/Stage5",sep="")	 #Putative unprocessed annotations;
           
           
           dir.create(subdir5,showWarnings=FALSE)
       }
	
	bestscore<-(-1000000)
	
        {
                #stop("Undefined value for parameter, cdfloc. Please enter path of the folder where the CDF files to be processed are located.")
                #change location to the output folder
                setwd(XCMS.outloc)
                alignmentresults<-list.files(XCMS.outloc, "*.txt")
                #alignmentresults<-list.files(XCMS.outloc, pattern="(XCMS).*(bw).*\\.txt")

                if(length(data_rpd_all)>0)
                {
                          curdata_dim={}
                          if(num_replicates==2)
                          {
                                  fileroot="_PID"
                          }
                          else
                          {
                                  if(num_replicates>2)
                                  {
                                          fileroot="_CV"
                                  }
                                  else
                                  {
                                          fileroot=""
                                  
					#  stop("Need at least 2 technical replicates per sample.")
				   }
                          }
                          #for(i in 1:length(alignmentresults))
                          cat("\n")
                           print("*******xMSanalyzer Stage 1: QC evaluation of invidual parameters*******")
                           cat("\n")
				rep_cor_mat<-{}
                parent_bad_list<-{}
			         for(i in 1:length(data_rpd_all))
                          {
				 				  print(paste("******Evaluating XCMS results from parameter setting ",i,"*******",sep=""))				  
                                  ############################################
                                  #2)Calculate pairwise correlation coefficients
                                  file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
                                  #curdata=read.table(paste(XCMS.outloc,"/",alignmentresults[i],sep=""),header=TRUE)
                                  #curdata=check.mz.in.replicates(curdata)
                                  curdata=data_rpd_all[[i]]
                            	  
			          #############################################
                                  ############################################
                                  #3) Calculate Percent Intensity Difference
			                   if(num_replicates>1)
                                  {
							
                            feat.eval.result=evaluate.Features(curdata, numreplicates=num_replicates,min.samp.percent=min.samp.percent,alignment.tool="XCMS",impute.bool=impute.bool,numnodes=numnodes)
                            cnames=colnames(feat.eval.result)
                            feat.eval.result<-apply(feat.eval.result,2,as.numeric)
                            feat.eval.result<-as.data.frame(feat.eval.result)
                            feat.eval.result.mat=cbind(curdata[,c(1:8)],feat.eval.result)
                            feat.eval.outfile=paste(subdir1,"/",file_name,fileroot,"featureassessment.txt",sep="")
                            #write results
                            write.table(feat.eval.result.mat, feat.eval.outfile,sep="\t", row.names=FALSE)
                            
                            curdata<-curdata[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),]
                            
                            feat.eval.result.mat<-feat.eval.result.mat[which(as.numeric(feat.eval.result$median)<=feat.filt.thresh),]
                            
			
                            curdata<-as.data.frame(curdata)
                            curdata<-replace(as.matrix(curdata),which(is.na(curdata)==TRUE),0)
                            
                            if(is.na(deltamzminmax.tol)==FALSE){
                                print("filtering by delta m/z")
                                mz_min_max<-cbind(curdata[,2],curdata[,3])
                                mz_min_max<-as.data.frame(mz_min_max)
                                
                                deltappm_res<-apply(mz_min_max,1,get_deltappm)
                                
                                curdata<-curdata[which(as.numeric(deltappm_res)<=deltamzminmax.tol),]
                                feat.eval.result.mat<-feat.eval.result.mat[which(as.numeric(deltappm_res)<=deltamzminmax.tol),]
                            }
                            
                            feateval_list[[i]]<-feat.eval.result.mat
                            data_rpd_all[[i]]<-curdata
								
								  if(num_replicates>1)
								  {
									  print(paste("**calculating pairwise ",cormethod," correlation**",sep=""))

									  
									  
											rsqres_1<-evaluate.Samples(curdata, num_replicates, alignment.tool="XCMS", cormethod,missingvalue,ignore.missing)
											
									
											rsqres<-as.data.frame(rsqres_1$cor.matrix)
											
											curdata<-as.data.frame(rsqres_1$feature.table)
											rsqres<-as.data.frame(rsqres)
											snames<-colnames(curdata[,-c(1:8)])
											snames_1<-snames[seq(1,length(snames),num_replicates)]
											rownames(rsqres)<-snames_1
											pcor_outfile=paste(subdir1,"/",file_name,"_sampleassessment_usinggoodfeatures.txt",sep="")
											write.table(rsqres, pcor_outfile,sep="\t",row.names=TRUE)
											rsqres_list[[i]]<-rsqres
								  }
								  else
								  {
									  print("**skipping sample evaluation as only one replicate is present**")
								  }
									
                                    
                                    
								if(num_replicates>2)
								{
									rep_cor_mat<-cbind(rep_cor_mat,rsqres$meanCorrelation)
									bad_samples<-which(rsqres$meanCorrelation<samp.filt.thresh)
								}else
								{
									bad_samples<-which(rsqres<samp.filt.thresh)
                                    #rep_cor_mat<-cbind(rep_cor_mat,rsqres[,1])
                                    rep_cor_mat<-rsqres
								}
								
								if(length(bad_samples)>0){
									bad_sample_names<-snames_1[bad_samples]
									
									feat.eval.outfile=paste(subdir1,"/",file_name,"_badsamples_at_cor",samp.filt.thresh,".txt",sep="")
									bad_sample_names<-as.data.frame(bad_sample_names)
									colnames(bad_sample_names)<-paste("Samples with correlation between technical replicates <", samp.filt.thresh,sep="")
									write.table(bad_sample_names, file=feat.eval.outfile,sep="\t", row.names=FALSE)
								}
								
								bad_list={}
								if(length(bad_samples)>0)
								{
									for(n1 in 1:length(bad_samples))
									{	
										if(bad_samples[n1]>1)
										{
											bad_samples[n1]=bad_samples[n1]+(bad_samples[n1]-1)*(num_replicates-1)
										}
											
									}
									for(n1 in 1:num_replicates)
									{
										bad_list<-c(bad_list,(bad_samples+n1-1))
									}
									bad_list<-bad_list[order(bad_list)]
									
								}
								if(i>1){
										parent_bad_list<-intersect(parent_bad_list,bad_list)
									}
									else{
									    
									    parent_bad_list<-bad_list

									}
								
								
							
                                  }
				  else
				  {
					  print("*********skipping feature evaluataion as only one replicate is present******")
				  }
                                
	       	}
		cat("\n")
		 print("********Stage 2: Filtering results from each paramter setting based on sample and feature quality checks*********")
		 cat("\n")
		 if(length(parent_bad_list)>0){
		

				if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
				{
				
						print("filtering bad sample id names")
									  sampleid_mapping<-sampleid_mapping[-c(parent_bad_list),]

									sampleidfilteredfile=paste(subdir2,"/","filtered_sampleid_mapping.txt",sep="")
									write.table(sampleid_mapping,file=sampleidfilteredfile,sep="\t",row.names=FALSE)
				}
		}

			if(num_replicates>2)
		{
			max_rep_mean_cor<-apply(rep_cor_mat,1,max)
		}else{
			max_rep_mean_cor<-rep_cor_mat
		}
        
        


rm(curdata)
		  for(i in 1:length(alignmentresults))
                          {
                                 
				 
				
				 
				    file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
                                  feat.eval.file=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
                                  #data_rpd_all[[i]]=read.table(feat.eval.file,header=TRUE)
				  
				  curdata<-data_rpd_all[[i]]
				 
				  feat.eval.result.mat<-feateval_list[[i]]
				  
				  
				
				  if(length(parent_bad_list)>0){
					  
					 print("filtering bad sample id names from feature table")
					  curdata<-curdata[,-c(parent_bad_list+8)]
					  
					  
					  
					  #sampleid_mapping<-sampleid_mapping[-c(parent_bad_list),]
					  #maxint<-apply(curdata[,-c(1:4,((dim(curdata)[2]-6):dim(curdata)[2]))],1,max)
					  maxint<-apply(curdata[,-c(1:8)],1,max)
					
					  badfeats<-which(maxint==0)
					  if(length(badfeats)>0){
						curdata<-curdata[-c(badfeats),]
						
						
						
						  feat.eval.result.mat<- feat.eval.result.mat[-c(badfeats),]
						  feateval_list[[i]]<-feat.eval.result.mat
						   #[which(as.numeric(feat.eval.result.mat$median)<=feat.filt.thresh),]
					  }
					  
				
					}
					
					data_rpd_all[[i]]<-curdata 
				 
				
				 feat.eval.outfile=paste(subdir2,"/",file_name,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtereddata.txt",sep="")	
								
				 #write results
                  write.table(data_rpd_all[[i]], feat.eval.outfile,sep="\t", row.names=FALSE)
				 
				 feat.eval.outfile=paste(subdir1,"/",file_name,fileroot,"featureassessment.txt",sep="")					  
								#write results
								write.table(feateval_list[[i]], feat.eval.outfile,sep="\t", row.names=FALSE)
							
			 }
                          ###########################################
                          #4) Merge two or more parameter settings
                          cat("\n")
                          print("*************Stage 3: merging features detected at different parameter settings********************")
                          cat("\n")
                          num_pairs=1
                          finalres={}
                          rnames={}
		         
		          if(merge.pairwise==TRUE){
		         			if(length(alignmentresults)>1){
                          for(i in 1:length(alignmentresults))
                          {
                                  file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
                                
                                 bool_num<-1
                                 
                                  for(j in i:length(alignmentresults))
                                  {
                                  	bool_num<-1
                                      #if(i!=j)
					if(i==j){
                                  	if(length(alignmentresults)>1){
                                  		bool_num<-0
                                  	}
                                  	else{
                                  		bool_num<-1
                                  		}
                                  	}
                	#if(i!=j)
                	 if(bool_num==1)                
				      {
				          file_name=sapply(strsplit(alignmentresults[j],".txt"),head)
					
                                         
                                         # if(i!=j)
                                          {
                                                  p1_p2=paste("p",i,"_U_","p",j,sep="")
                                          }
                                         # else
                                          #{
                                           #       p1_p2="p1"
                                          #}
					  
					    feat.eval.A<-feateval_list[[i]]
					 feat.eval.B<-feateval_list[[j]]
					 feat.eval.A<-feat.eval.A[which(as.numeric(feat.eval.A$median)<=feat.filt.thresh),]
					 feat.eval.B<-feat.eval.B[which(as.numeric(feat.eval.B$median)<=feat.filt.thresh),]
					 
					
					#print(p1_p2)
					 
					 print(paste("Number of good quality features from setting ",i,":", dim(data_rpd_all[[i]])[1],sep=": "))
						print(paste("Number of good quality features from setting ",j,":",dim(data_rpd_all[[j]])[1],sep=": "))
	
					  data_m=merge.Results(data_rpd_all[[i]],data_rpd_all[[j]],feat.eval.A,feat.eval.B,max.mz.diff,max.rt.diff, merge.eval.pvalue,alignment.tool="XCMS",
					numnodes=numnodes, mult.test.cor,mergecorthresh,missingvalue)
                                      
                      numcols<-dim(data_m)[2]
					 



					 data_m_int<-data_m[,-c(1:8,(numcols-8):numcols)]


					 numsamps<-dim(data_m_int)[2]/num_replicates
                     maxint<-apply(data_m_int,1,function(x){max(x,na.rm=TRUE)})
                     
                     
                     
                    # numpeaks<-lapply(1:dim(data_m_int)[1],function(j){length(which(is.na(data_m_int[j,])==FALSE))})
			numpeaks<-apply(data_m_int,1,countpeaks)		
			
                     
                     data_m[,8]<-numpeaks
                                      
                      if(is.na(minexp.pct)==FALSE)
                      {
					minexp<-round(minexp.pct*dim(data_m_int)[2])
					
					if(length(which(data_m[,8]>=minexp))>0){
					data_m<-data_m[which(data_m[,8]>=minexp),]
					}else{
						stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
						}
					
					
					}                 
                                     
                    
					
					
					
					
					union_list[[num_pairs]]<-data_m[,c(1:8)]
					
					#union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],numpeaks)
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$numgoodsamples)
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$median)
					#Qscore<-(as.numeric(data_m$numgoodsamples)/as.numeric(data_m$median+0.1))
					#Qscore<-100*(Qscore/numsamps)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$Qscore)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],maxint)
					
			
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m_int)
					
					featinfo<-colnames(data_m[,c(1:7)])
					cnames<-colnames(data_m_int)
					merge.res.colnames<-c(featinfo,"NumPres.All.Samples","NumPres.Biological.Samples",paste("median",fileroot,sep=""),"Qscore","Max.Intensity",cnames)
					colnames(union_list[[num_pairs]])<-as.character(merge.res.colnames)
					
					
					
					 
					
					  	  curres={}
					   curres=cbind(curres, p1_p2)
                                          curres=cbind(curres, dim(union_list[[num_pairs]])[1])
                                          curres=cbind(curres, mean(as.numeric(data_m$median)))
                                          
                                          
                                          curres=cbind(curres, mean(as.numeric(data_m$Qscore)))
                                          # curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$Qscore))))
                                          
                                          
                                           curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$median))))
                                         
                                          if(curscore>bestscore){
                                          	
                                          	bestscore<-curscore
                                          	best_i<-num_pairs
                                          	best_pair<-p1_p2
                                          }
                                          
                                          
                                          curres=cbind(curres,curscore)
                                          
					  curres<-as.data.frame(curres)
                                          
                                          finalres=rbind(finalres,curres)
					  
                                             finalname=paste("XCMS_feature_list_at_", p1_p2,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
					  
                                          #Output merge results
                                          write.table(union_list[[num_pairs]],file=paste(subdir3,"/",finalname,sep=""), sep="\t",row.names=FALSE)
					  
		
				
					  num_pairs=num_pairs+1
                                     }
				  }
                          }
                          
                          }else{
                          	
                          
				for(i in 1:length(alignmentresults))
                          {
                                  file_name=sapply(strsplit(alignmentresults[i],".txt"),head)
                                
                                 
                                  j=i
                                  {
                                      if(i==j)
				      				{
				          file_name=sapply(strsplit(alignmentresults[j],".txt"),head)
					
                                         
                                        #  if(i!=j)
                                          {
                                                  p1_p2=paste("p",i,"_U_","p",j,sep="")
                                          }
                                         # else
                                          #{
                                           #       p1_p2="p1"
                                          #}
					  
					    feat.eval.A<-feateval_list[[i]]
					 feat.eval.B<-feateval_list[[j]]
					 feat.eval.A<-feat.eval.A[which(as.numeric(feat.eval.A$median)<=feat.filt.thresh),]
					 feat.eval.B<-feat.eval.B[which(as.numeric(feat.eval.B$median)<=feat.filt.thresh),]
					 
					
					#print(p1_p2)
					 
					 print(paste("Number of good quality features from setting ",i,":", dim(data_rpd_all[[i]])[1],sep=": "))
					 #print(paste("Number of good quality features from setting ",j,":",dim(data_rpd_all[[j]])[1],sep=": "))
	
					  data_m=merge.Results(data_rpd_all[[i]],data_rpd_all[[j]],feat.eval.A,feat.eval.B,max.mz.diff,max.rt.diff, merge.eval.pvalue,alignment.tool="XCMS",
					numnodes=numnodes, mult.test.cor,mergecorthresh,missingvalue)
                                         
					  numcols<-dim(data_m)[2]
					
					

					 data_m_int<-data_m[,-c(1:8,(numcols-8):numcols)]
			
					numsamps<-dim(data_m_int)[2]/num_replicates
					 
					  
					  
					
					
					maxint<-apply(data_m_int,1,function(x){max(x,na.rm=TRUE)})

				
					
					#numpeaks<-lapply(1:dim(data_m_int)[1],function(j){length(which(is.na(data_m_int[j,])==FALSE))})
					numpeaks<-apply(data_m_int,1,countpeaks)		
                
                    data_m[,8]<-numpeaks
                    
                    
                if(is.na(minexp.pct)==FALSE){
                    minexp<-round(minexp.pct*dim(data_m_int)[2])
                    
                    if(length(which(data_m[,8]>=minexp))>0){
                        data_m<-data_m[which(data_m[,8]>=minexp),]
                    }else{
                        stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
                    }
                    
                    
                }
                
				
					 
					  union_list[[num_pairs]]<-data_m[,c(1:8)]
					
				        #union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],numpeaks)	
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$numgoodsamples)
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$median)
					#Qscore<-(as.numeric(data_m$numgoodsamples)/as.numeric(data_m$median+0.1))
					#Qscore<-100*(Qscore/numsamps)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$Qscore)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],maxint)
					
			
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m_int)
					
					featinfo<-colnames(data_m[,c(1:7)])
					cnames<-colnames(data_m_int)
					merge.res.colnames<-c(featinfo,"NumPres.All.Samples","NumPres.Biological.Samples",paste("median",fileroot,sep=""),"Qscore","Max.Intensity",cnames)
					colnames(union_list[[num_pairs]])<-as.character(merge.res.colnames)
					
					
					
					  	  curres={}
					   curres=cbind(curres, p1_p2)
                                          curres=cbind(curres, dim(union_list[[num_pairs]])[1])
                                          curres=cbind(curres, mean(as.numeric(data_m$median)))
                                          curres=cbind(curres, mean(as.numeric(data_m$Qscore)))
                                          
                                          # curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$Qscore))))
                                          
                                           curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$median))))
                                         
                                          if(curscore>bestscore){
                                          	
                                          	bestscore<-curscore
                                          	best_i<-num_pairs
                                          	best_pair<-p1_p2
                                          }
                                          
                                          
                                          curres=cbind(curres,curscore)
					  curres<-as.data.frame(curres)
                                          
                                          finalres=rbind(finalres,curres)

					  
                                             finalname=paste("XCMS_feature_list_at_", p1_p2,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
					  
                                          #Output merge results
                                          write.table(union_list[[num_pairs]],file=paste(subdir3,"/",finalname,sep=""), sep="\t",row.names=FALSE)
					  
					  metlin.res={}
                                          kegg.res={}
					 #length(union_list[[num_pairs]]$mz
					
					  num_pairs=num_pairs+1
                                     }
				  }
                          }

                          	
                          	}
                      
                        
			}else{
			
						data_rpd_all_parameters<-{}
							 feat.eval.all<-{}
							 p1_p2<-"pall"
								
								for(i in 1:length(alignmentresults))
								{
									 
										
													
									data_rpd_all_parameters<-rbind(data_rpd_all_parameters,data_rpd_all[[i]])
									
									feat.eval.A<-feateval_list[[i]]
									 feat.eval.A<-feat.eval.A[which(as.numeric(feat.eval.A$median)<=feat.filt.thresh),]
								
									 feat.eval.all<-rbind(feat.eval.all,feat.eval.A)
									
								}
								
								
								cnames1<-colnames(feateval_list[[i]])
								
								
								
								#feat.eval.all<-unique(feat.eval.all)
								 #data_rpd_all_parameters<-unique(data_rpd_all_parameters)
								 
								 feat.eval.all<-as.data.frame(feat.eval.all)
								 data_rpd_all_parameters<-as.data.frame(data_rpd_all_parameters)
								 
								 
								data_m=merge.Results(data_rpd_all_parameters,data_rpd_all_parameters, feat.eval.all,feat.eval.all,max.mz.diff,max.rt.diff,merge.eval.pvalue,alignment.tool="XCMS",
								 numnodes=numnodes,mult.test.cor,mergecorthresh,missingvalue)
								 
								 
								 numcols<-dim(data_m)[2]



					 data_m_int<-data_m[,-c(1:8,(numcols-8):numcols)]
					numsamps<-dim(data_m_int)[2]/num_replicates
					
                                           numcols<-dim(data_m)[2]

if(is.na(minexp.pct)==FALSE)
{
    minexp<-round(minexp.pct*dim(data_m_int)[2])
    
    if(length(which(data_m[,8]>=minexp))>0){
        data_m<-data_m[which(data_m[,8]>=minexp),]
    }else{
        stop(paste("No features have non-missing value in ",minexp, " samples",sep=""))
    }
    
    
}

					 data_m_int<-data_m[,-c(1:8,(numcols-8):numcols)]
					numsamps<-dim(data_m_int)[2]/num_replicates
					 
					maxint<-apply(data_m_int,1,function(x){max(x,na.rm=TRUE)})
					#numpeaks<-lapply(1:dim(data_m_int)[1],function(j){length(which(is.na(data_m_int[j,])==FALSE))})
					numpeaks<-apply(data_m_int,1,countpeaks)		
					 
					 
					
					 
					  union_list[[num_pairs]]<-data_m[,c(1:7)]
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],numpeaks)
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$numgoodsamples)
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$median)
					#Qscore<-(as.numeric(data_m$numgoodsamples)/as.numeric(data_m$median+0.1))
					#Qscore<-100*(Qscore/numsamps)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m$Qscore)
					
					
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],maxint)
					
			
					union_list[[num_pairs]]<-cbind(union_list[[num_pairs]],data_m_int)
					
					featinfo<-colnames(data_m[,c(1:7)])
					cnames<-colnames(data_m_int)
					merge.res.colnames<-c(featinfo,"NumPres.All.Samples","NumPres.Biological.Samples",paste("median",fileroot,sep=""),"Qscore","Max.Intensity",cnames)
					colnames(union_list[[num_pairs]])<-as.character(merge.res.colnames)
					
					
					
					  curres={}
					   curres=cbind(curres, p1_p2)
                                         
                      curres=cbind(curres, dim(union_list[[num_pairs]])[1])
                      curres=cbind(curres, mean(as.numeric(data_m$median)))
                      curres=cbind(curres, mean(as.numeric(data_m$Qscore)))
                      
                       #  curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$Qscore))))
                         curscore<-(dim(union_list[[num_pairs]])[1]-(scoreweight*mean(as.numeric(data_m$median))))
                                                           
                                          if(curscore>bestscore){
                                          	
                                          	bestscore<-curscore
                                          	best_i<-num_pairs
                                          	best_pair<-p1_p2
                                          }
                                          
                                          
                                          curres=cbind(curres,curscore)
                      
					  curres<-as.data.frame(curres)
                                          
                                          finalres=rbind(finalres,curres)
					  
					
                                             finalname=paste("XCMS_feature_list_at_", p1_p2,"cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
					  
                                          #Output merge results
                                          write.table(union_list[[num_pairs]],file=paste(subdir3,"/",finalname,sep=""), sep="\t",row.names=FALSE)
					  
			}
			      finalres<-as.data.frame(finalres)
                  
                  if(merge.pairwise==TRUE){
			  colnames(finalres)<-c("Parameter Combination", "Number of Features", "median PID/CV between sample replicates", "mean Qscore (Quality score)", "Parameter score")
                          write.table(finalres,file=paste(xMSanalyzer.outloc,"/xcms_with_xMSanalyzer_merge_summary.txt",sep=""), sep="\t", row.names=FALSE)
                  }
                  
	print("Most optimal feature setting:")
    print(best_pair)
    cat("\n")
    #########################################################################
    
    
    #rawQCeval/
     cat("\n")
      print(paste("********Stage 3b: Generating final (pre-batcheffect correction) untargeted and targeted feature tables using ",best_pair," results******",sep=""))
      cat("\n")
    
   # pdf("Stage4a_QE_plots.pdf")
    #pdf(file=paste("subdir4a/Stage4a_QE_plots.pdf",sep=""))
    
    pdf(file=paste(subdir4a,"/QE_plots_All_samples_RAW.pdf",sep=""))
    
    #pdf(file=paste(xMSanalyzer.outloc,"/QE_plots_RAW.pdf",sep=""))
    #most optimal set after merger

    
 	d1<-union_list[[best_i]]
    
    d1<-unique(d1)
    
   if(is.na(refMZ)==FALSE){
			stddata<-read.table(refMZ,sep="\t",header=TRUE)
			print(refMZ)
			print(head(stddata))
			
		}else{
			if(charge_type=="pos"){
			data(example_target_list_pos)
			stddata<-example_target_list_pos
			}else{
				
				if(charge_type=="neg"){
						data(example_target_list_neg)
						stddata<-example_target_list_neg
					}else{
						stop("Invalid option. \'charge_type\' should be \'pos\' or \'neg\'.")
						}
				}
		}


d1_int<-round(d1[,-c(1:12)],0)


rsqres_list<-evaluate.Samples(d1_int, num_replicates, alignment.tool=NA, cormethod,missingvalue,ignore.missing)

rsqres<-as.data.frame(rsqres_list$cor.matrix)

rsqres<-as.data.frame(rsqres)
snames<-colnames(d1_int)
snames_1<-snames[seq(1,length(snames),num_replicates)]
rownames(rsqres)<-snames_1
pcor_outfile=paste(subdir4a,"/Pairwise_Pearson_correlation_technical_replicates.txt",sep="")

write.table(rsqres, pcor_outfile,sep="\t",row.names=TRUE)


if(num_replicates>2)
{
    max_rep_mean_cor<-apply(rsqres,1,max)
    
    max_rep_mean_cor<-as.data.frame(max_rep_mean_cor)
}else{
    max_rep_mean_cor<-rsqres
    
    
}


#if(is.na(sample_info_file)==FALSE)
if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
{
    #sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)
    
    sampleid_mapping[,1]<-gsub(sampleid_mapping[,1],pattern=filepattern,replacement="",ignore.case=TRUE)
    cnames<-colnames(data_m)
    
    sampleid_mapping[,2]<-tolower(sampleid_mapping[,2])
    
    qc_label<-tolower(qc_label)
    
    cnames<-gsub(cnames,pattern=filepattern,replacement="",ignore.case=TRUE)
    colnames(data_m)<-cnames
    
    class_labels<-as.factor(sampleid_mapping[,2])
    
    if(is.na(qc_label)==FALSE){
    
        qc_label<-tolower(qc_label)
    qc_label_check<-gregexpr(pattern=qc_label,text=sampleid_mapping[,2])
    
    qc_index<-which(qc_label_check>0)
    
    if(length(qc_index)>0){
  
    print("QC file index")
    print(qc_index)
    
	d1_int_qc<-d1_int[,qc_index]

    
    rsqres_listQC<-evaluate.Samples(d1_int_qc, numreplicates=length(qc_index), alignment.tool=NA, cormethod,missingvalue,ignore.missing)
   
    rsqres_listQC$cor.matrix<-round(rsqres_listQC$cor.matrix,2)
    
    r1<-matrix(0,nrow=length(qc_index),ncol=length(qc_index))
    
    #r1[upper.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
    diag(r1)<-1
    #r1[lower.tri(r1)]<-rsqres_listQC$cor.matrix #[-length(rsqres_listQC$cor.matrix)]
    
    start_count=0
    for(r1_row in 1:(length(qc_index)-1))
    {
        
          
            r1[r1_row,c((r1_row+1):length(qc_index))]<-rsqres_listQC$cor.matrix[,(start_count+1):(start_count+length(qc_index)-r1_row)]
           
           
            start_count=(start_count+length(qc_index)-r1_row) #(r1_row*length(qc_index))-r1_row-1 #(length(qc_index)-1)


    }
    
    rsqresQC<-as.data.frame(rsqres_listQC$cor.matrix)

pcor_outfile=paste(subdir4a,"/Pairwise_Pearson_correlation_QC1.txt",sep="")

write.table(r1, pcor_outfile,sep="\t",row.names=TRUE)
pcor_outfile=paste(subdir4a,"/Pairwise_Pearson_correlation_QC2.txt",sep="")

write.table(rsqresQC, pcor_outfile,sep="\t",row.names=TRUE)


    
    r1[lower.tri(r1)]<-NA #r1[upper.tri(r1,diag=TRUE)] #[-length(rsqres_listQC$cor.matrix)]
    
    
   
   TIC <- vector("list",length(qc_index))
   
   for(qc_i in 1:length(qc_index)){
       
       
       #plot( d1$time, d1_int_qc[,qc_i],xlab="Retention Time",ylab="TIC",main=mainlab1)
       TIC[[qc_i]] <-cbind(d1$time, d1_int_qc[,qc_i])
   }
    }else{
        
        print("No matches found for the QC samples.")
        qc_label=NA
    }
   
    }
   
  
}



d1<-cbind(d1[,c(1:12)],d1_int)
rm(d1_int)


  
   			Sys.sleep(1)
            
            num_features<-nrow(d1)
            max_ylim<-nrow(d1)+100
            
            num_features<-nrow(d1)
            
            
            par(mfrow=c(1,2))
            if(num_replicates>2){
                
                rep_cv<-d1$median_CV
                
                rep_cv_ylab<-"CV"
                
                h1<-hist(d1$median_CV,breaks=seq(0,max(d1$median_CV,na.rm=TRUE)+10,10),main=paste("Histogram median CV \n (using all ",num_features," features)",sep=""),col="brown",xlab="median CV%", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
                lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
            }else{
                
                rep_cv<-d1$median_PID
                
                rep_cv_ylab<-"PID"
                
                h1<-hist(d1$median_PID,breaks=seq(0,max(d1$median_PID,na.rm=TRUE)+10,10),main=paste("Histogram median CV \n (using all ",num_features," features)",sep=""),col="brown",xlab="median PID%", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
                lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
            }
            lab_text<-paste(h1$breaks,"-",h1$breaks+10,sep="")
            #pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main="Pie chart of median CVs (%) \n using all features",cex.main=0.7)
            
            if(num_replicates>2){
                pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main=paste("Pie chart of median CVs (%) \n using all ",num_features," features\n; average=",round(mean(d1$median_CV),2),sep=""),cex.main=0.7)
            }else{
                pie(h1$counts,col=rainbow(length(h1$counts)),labels=lab_text,main=paste("Pie chart of median PIDs (%) \n using all ",num_features," features\n; average=",round(mean(d1$median_PID),2),sep=""),cex.main=0.7)
                
            }
            
            par(mfrow=c(1,1))
            
            
		
			hist(d1$NumPres.All.Samples,main="Histogram NumPres.All.Samples",col="brown",xlab="Number of samples (including replicates) \n with non-zero intensity values", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
			
			hist(d1$NumPres.Biological.Samples,main="Histogram NumPres.Biological.Samples",col="brown",xlab="Number of biological samples \n with non-zero intensity values", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
			#h1<-hist(d1$median_CV,main="Histogram median CV \n (using all data)",col="brown",xlab="median CV%", ylab="Number of features",cex.main=0.7)
			

			par(mfrow=c(1,1))
			
			#HEREDEBUG
            if(num_replicates>1){
                
                min_rep_pcor<-round(min(max_rep_mean_cor[,1],na.rm=TRUE),2)
                max_rep_pcor<-round(max(max_rep_mean_cor[,1],na.rm=TRUE),2)
                
                if(is.na(samp.filt.thresh)==FALSE){
                    
                    
                    parent_bad_list<-which(max_rep_mean_cor[,1]<samp.filt.thresh)
                    
                    if(length(parent_bad_list)>0){
                        
                        if(is.na(sample_info_file)==FALSE && sample_info_file!="NA"){
                            
                            filt_samp_names<-sampleid_mapping[parent_bad_list,1]
                            filt_samp_names<-paste(filt_samp_names,collapse=";")
                            
                            filt_samp_names<-length(parent_bad_list)
                        }
                    }else{
                        filt_samp_names<-"None"
                        
                    }
                }else{
                    filt_samp_names<-"None"
                    
                }

                
                
                hist(max_rep_mean_cor[,1],breaks=seq(0,1,0.1),main=paste("Histogram for mean Pearson correlation (min: ",min_rep_pcor,"; max: ",max_rep_pcor,") \n within technical replicates after filtering \n # of files filtered at threshold ",samp.filt.thresh,": ",filt_samp_names,sep=""),col="brown",xlab="mean replicate Correlation", ylab="Number of samples",cex.main=0.7)
                
                
                
            }
            

            

par(mfrow=c(2,2))
#5ppm
data_a<-find.Unique.mzs.sameset(dataA=d1[,c("mz","time")],dataB=d1[,c("mz","time")],mz.thresh=5,time.thresh=NA,alignment.tool=NA)
num_unique_features<-nrow(data_a$uniqueA)
total_features<-nrow(d1)

xlab2<-paste("overlapping m/z features \n based on +/- ",5,"ppm overlap criteria",sep="")

temp_df<-cbind(total_features,num_unique_features)

pie_v1<-round(temp_df[1,2]/temp_df[1,1],2)
pie_v2<-1-pie_v1

temp_dfA<-cbind(temp_df[1,1],temp_df[1,2], (temp_df[1,1]-temp_df[1,2]))#cbind(100*pie_v1,100*pie_v2)

colnames(temp_dfA)<-c("Total","Unique","Overlapping")

#barplot(temp_dfA,col="brown",ylab="Percentage (%) of total features", main=paste("Overlapping vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,100))

barplot(temp_dfA,col="brown",ylab="Number of features", main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,temp_df[1,1]+100))
#non-overlapping m/z features")

# pie(c(pie_v2,pie_v1),labels=c("Overlapping","Unique"),col=c("orange","brown"),main=paste("Total vs ",xlab2,sep=""),cex.main=0.8)

pie_label1<-paste("Overlapping \n(",100*pie_v2,"%)",sep="")
pie_label2<-paste("Unique \n(",100*pie_v1,"%)",sep="")

pie(c(pie_v2,pie_v1),labels=c(pie_label1,pie_label2),col=c("orange","brown"),main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8)

#10ppm
data_a<-find.Unique.mzs.sameset(dataA=d1[,c("mz","time")],dataB=d1[,c("mz","time")],mz.thresh=10,time.thresh=NA,alignment.tool=NA)
num_unique_features<-nrow(data_a$uniqueA)
total_features<-nrow(d1)

#xlab2<-paste("Unique (non-overlapping) m/z features \n based on +/- ",10,"ppm overlap criteria",sep="")

xlab2<-paste("overlapping m/z features \n based on +/- ",10,"ppm overlap criteria",sep="")


temp_df<-cbind(total_features,num_unique_features)
#colnames(temp_df)<-c("Total","Unique (non-overlapping)")

pie_v1<-round(temp_df[1,2]/temp_df[1,1],2)
pie_v2<-1-pie_v1

temp_dfB<-cbind(temp_df[1,1],temp_df[1,2], (temp_df[1,1]-temp_df[1,2])) #cbind(100*pie_v1,100*pie_v2)

colnames(temp_dfB)<-c("Total","Unique","Overlapping")

#barplot(temp_dfB,col="brown",ylab="Percentage (%) of total features", main=paste("Total vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,100)) #non-overlapping m/z features")

barplot(temp_dfB,col="brown",ylab="Number of features", main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8,ylim=c(0,temp_df[1,1]+100))


pie_label1<-paste("Overlapping \n(",100*pie_v2,"%)",sep="")
pie_label2<-paste("Unique \n(",100*pie_v1,"%)",sep="")
pie(c(pie_v2,pie_v1),labels=c(pie_label1,pie_label2),col=c("orange","brown"),main=paste("Unique (non-overlapping) vs ",xlab2,sep=""),cex.main=0.8)

par(mfrow=c(1,1))

			hist(d1$Qscore,main="Histogram Qscore",col="brown",xlab="Quality score", ylab="Number of features",cex.main=0.7,ylim=c(0,max_ylim))
			
            
            
			#pie(h1$counts,col=rainbow(length(h1$counts)),labels=h1$breaks,main="Pie chart of median CVs (%) \n using all features",cex.main=0.7)
			
			#pie(h1$counts,col=rainbow(length(h1$counts)),labels=h1$breaks,main="Pie chart of median CVs (%) \n using all features",cex.main=0.7)
			
				
			d2<-d1[order(d1$time),]
	
	
		
		
        plot(d2$time,d2$mz,main="m/z vs Time",col="brown",xlab="Time (s)",ylab="m/z",cex.main=0.7)
        
        
        par(mfrow=c(1,1))
        
        plot(d2$mz,d2$Max.Intensity,main="Intensity vs m/z",col="brown",xlab="m/z",ylab="Max intensity \nacross all samples/profiles",cex.main=0.7)
        plot(d2$time,d2$Max.Intensity,main="Intensity vs time",col="brown",xlab="Time (s)",ylab="Max intensity \nacross all samples/profiles",cex.main=0.7)
        
        plot(rep_cv,d2$Max.Intensity,main=paste("Intensity vs ",rep_cv_ylab,sep=""),col="brown",xlab=rep_cv_ylab,ylab="Max intensity \nacross all samples/profiles",cex.main=0.7)
        
        
        par(mfrow=c(1,1))
        
        plot(d2$mz,d2$Qscore,main="Qscore vs m/z",col="brown",xlab="m/z",ylab="Qscore",cex.main=0.7)
        plot(d2$time,d2$Qscore,main="Qscore vs Time",col="brown",xlab="Time (s)",ylab="Qscore",cex.main=0.7)
        plot(rep_cv,d2$Qscore,main=paste("Qscore vs ",rep_cv_ylab,sep=""),col="brown",xlab=rep_cv_ylab,ylab="Qscore",cex.main=0.7)
        
        
        
        par(mfrow=c(1,1))
			
	max_numzeros<-dim(d1)[2]*1

	if(is.na(void.vol.timethresh)==FALSE){
		dfirst15<-d1[which(d1$time<void.vol.timethresh),]
	
	  
        if(nrow(dfirst15)>1){
            
            
            ind1<-which(dfirst15$Max.Intensity==max(dfirst15$Max.Intensity))[1]
            
            time_thresh<-dfirst15$time[ind1]
            
            time_thresh<-time_thresh-(0.30*time_thresh)
            
            time_thresh<-round(time_thresh,1)
            
            plot(dfirst15$time,dfirst15$Max.Intensity,xlab="Time (s)", col="brown",ylab="Max intensity \nacross all samples", main=paste("Estimated void volume time: ",time_thresh," s",sep=""))
            abline(v=time_thresh,col=4,lty=3)
            
            
            
            d1<-d1[which(d1$time>=time_thresh),]
            
            print("Estimated void volume time cutoff")
            print(time_thresh)
            
            
        }else{
            
            print("No features eluting before void volume time cutoff")
        }
		
		
	
		#sfname<-paste(xMSanalyzer.outloc,"/stage5/feature_table_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,"filtered.txt",sep="")
		finalname<-paste("featuretable_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,"_voidtimefilt.txt",sep="")
		
        #write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
	   

	}else{
		
		
		finalname<-paste("featuretable_",best_pair,"_cor",samp.filt.thresh,fileroot,feat.filt.thresh,".txt",sep="")
		
        #write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
	   

		
		}



	cat("\n")
    print(paste("********Stage 4a: Performing QC checks using ",best_pair," results*******",sep=""))
    cat("\n")
    
    Sys.sleep(1)
    
		
		print("Dim data after void time filtering")
		print(dim(d1))
	 
     
     targeted_feat_raw<-eval.target.mz(dataA=d1[,-c(2:3,5:12)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="raw",xMSanalyzer.outloc=subdir4a)
     
   raw_mz_time<-d1[,c("mz","time")]
	    
        write.table(raw_mz_time,file=paste(subdir3b,"/Precalibration_mz_time.txt",sep=""),sep="\t",row.names=FALSE)
        
	print("m/z calibration")
     
     median_error<-median(targeted_feat_raw$delta_ppm_error,na.rm=TRUE)
     
	print("median m/z error in ppm")
	print(median_error)
	
	  median_time_error<-median(targeted_feat_raw$delta_ppm_error,na.rm=TRUE)

     #d1<-get_calibrated_mz_data(dataA=d1,delta_ppm_error=median_error,delta_time_error=median_time_error)
       #here1
       #     d1<-get_calibrated_mz_data(dataA=d1,delta_ppm_error=median_error,delta_time_error=median_time_error,refMZ=stddata,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff)
           d1<-get_calibrated_mz_data(dataA=d1,refMZ=stddata,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,feature.eval.result=feat.eval.result,calibration.method=calibration.method[1])   
     
     
     finalname<-paste("RAW_mzcalibrated_untargeted_featuretable.txt",sep="")
     
     write.table(d1,file=paste(subdir3b,"/",finalname,sep=""),sep="\t",row.names=FALSE)
     
     #write.table(d1,file=paste(xMSanalyzer.outloc,finalname,sep=""),sep="\t",row.names=FALSE)
     


     targeted_calibfeat_raw<-eval.target.calibratedmz(dataA=d1[,-c(2:3,5:12)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="raw",xMSanalyzer.outloc=subdir4a)
        

	data_m<-d1[,-c(1:12)]
   
    if(replacezeroswithNA==TRUE){
		data_m<-replace(as.matrix(data_m),which(data_m==0),NA)
		
		d1<-cbind(d1[,c(1:12)],data_m)
	}
	
	
	counna<-apply(data_m,1,function(x){length(which(is.na(x)==TRUE))})

	maxzeros<-1*dim(data_m)[2]
	
	if(length(which(counna<maxzeros))){
		data_m<-data_m[which(counna<maxzeros),]
	}

    maxint<-apply(data_m,1,function(x){max(x,na.rm=TRUE)})
   	maxint_ord<-order(maxint,decreasing=TRUE)
   	#[maxint_ord[1:5000]
   	
    X<-t(data_m) #[maxint_ord[1:2000],])
    X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
    
    tic.eval(d1[,-c(1:12)],outloc=subdir4a)
    feat.eval.result=evaluate.Features(d1[,-c(9:12)], numreplicates=num_replicates,min.samp.percent=min.samp.percent,alignment.tool="XCMS",impute.bool=impute.bool,numnodes=numnodes)
 							
    cnames=colnames(feat.eval.result)
		
			

    Sys.sleep(1)
    s1<-"Stage 1 results: QC evaluation of invidual parameters from XCMS"
s2<-"Stage 2 results: filtered results from each paramter setting based on sample and feature quality (CV within replicates) checks"
s3<-"Stage 3a results: merged results using stage 2 filtered files"
s3b<-"Stage 3b results: RAW m/z calibrated untargeted and targeted feature tables (averaged and non-averaged)"
s4a<-"Stage 4a results: QC evaluation of targeted and untargeted data before batch-effect correction"
sm<-rbind(s1,s2,s3,s3b,s4a)

#targeted_feat_raw<-eval.target.mz(dataA=d1[,-c(3:12)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4a,folderheader="raw")


if(summarize.replicates==TRUE && data.norm.pipeline=="AC"){
    
    if(num_replicates>1){
  
    sampleid_mapping<-sampleid_mapping[seq(1,dim(d1[,-c(1:12)])[2],num_replicates),]
   
	 minexp=minexp/num_replicates 
    data_matrix<-data_summarize(Xmat=d1[,- c(2:3,5:12)],Ymat=NA,feature_table_file=NA,parentoutput_dir=subdir3b,class_labels_file=NA,num_replicates=num_replicates,summarize.replicates=summarize.replicates,
    summary.method=summary.method,missing.val=missingvalue, rep.num.max.missing.thresh=rep.num.max.missing.thresh,summary.na.replacement=summary.na.replacement)
    
    print(dim(data_matrix))
    
    d1<-cbind(d1[,c(1:12)],data_matrix[,-c(1:2)])
    
    X<-t(data_matrix[,-c(1:2)])
    data_m<-data_matrix[,-c(1:2)]
    }
    
    num_replicates=1
    
}
	

   # if(is.na(sample_info_file)==FALSE)
     if(is.na(sample_info_file)==FALSE && sample_info_file!="NA")
    {
    	
    	
		#sampleid_mapping<-read.table(sample_info_file,sep="\t",header=TRUE)

	    	sampleid_mapping[,1]<-gsub(sampleid_mapping[,1],pattern=filepattern,replacement="",ignore.case=TRUE)
		cnames<-colnames(data_m)
		
		cnames<-gsub(cnames,pattern=filepattern,replacement="",ignore.case=TRUE)
		colnames(data_m)<-cnames
		
    	batch_inf<-sampleid_mapping[,3]
    	
    	batch_labels<-as.factor(sampleid_mapping[,3])
    
    	l1<-levels(batch_labels)
     
    	   batch_levels<-levels(batch_labels)
		   print("Batch-effect evaluation")
		    try(pca.eval(X=X,samplelabels=batch_labels,filename="raw",ncomp=20,center=TRUE,scale=TRUE,legendlocation="bottomleft",legendcex=0.5,outloc=subdir4a),silent=TRUE)
			 #pca.eval(X,samplelabels=batch_labels,filename="Raw",ncomp=20,center=TRUE,scale=TRUE,legendlocation="topright",legendcex=0.5,outloc=outputloc)	
			Sys.sleep(1)
			
			
	
   		   if(dim(sampleid_mapping)[2]<4){
   		    mod<-rep(1,dim(data_m)[2])
   		    }else{
   		    		mod<-sampleid_mapping[,-c(1:3)]
   		    	}
 	
    		dev.off()

    		Sys.sleep(1)
    #centwave
	   if(length(batch_levels)>1){
		     #######################################################################
		     cat("\n")
		    print("Stage 4b: Performing batch-effect correction and post-correction QC evaluation")
		    cat("\n")
		   # pdf("Stage4b_QE_plots.pdf")
		    #pdf(file=paste("subdir4b/Stage4b_QE_plots.pdf",sep=""))
		    
           # pdf(file=paste(subdir4b,"/QE_plots_All_samples_ComBat.pdf",sep=""))
            
            #pdf(file=paste(xMSanalyzer.outloc,"/QE_plots_ComBat.pdf",sep=""))
            
	    
	               if(summarize.replicates==FALSE){
			pdf(file=paste(subdir4b,"/QE_plots_All_samples_ComBat.pdf",sep=""))
			}else{
               
	        if(summarize.replicates==TRUE && data.norm.pipeline=="AC"){
			if(summary.method=="mean"){
			pdf(file=paste(subdir4b,"/QE_plots_All_samples_1averaged_2ComBat.pdf",sep=""))
			}else{
				pdf(file=paste(subdir4b,"/QE_plots_All_samples_1mediansummarized_2ComBat.pdf",sep=""))
			}
		}else{
		
			pdf(file=paste(subdir4b,"/QE_plots_All_samples_ComBat.pdf",sep=""))
		
		}
		
		
               
           }
		    ##################################################################
		
	
		   
		   	  if(is.na(missingvalue)==FALSE){
			
			  if(replacezeroswithNA==TRUE){
					#data_m<-replace(as.matrix(data_m),which(data_m==0),NA)
					data_m_na=replace(as.matrix(data_m),which(data_m==missingvalue),NA)
					d1<-cbind(d1[,c(1:12)],data_m_na)
					
					sum_check<-apply(data_m_na,1,function(x){length(which(is.na(x)==FALSE))})
					
					
					if(length(which(sum_check<minexp))>0){
					data_m_na<-data_m_na[-which(sum_check<minexp),]
					d1<-d1[-which(sum_check<minexp),]
					
					}
					if(summarize.replicates==TRUE){
					 print(paste("Number of features after minexp ",minexp," filtering post-averaging",sep=""))
                                         print(dim(data_m)[1])
					 }
					
					}else{
						data_m_na<-data_m
					}

			}else{
					data_m_na<-data_m
					sum_check<-apply(data_m_na,1,function(x){length(which(is.na(x)==FALSE))})
					
					if(length(which(sum_check<minexp))>0){
						data_m_na<-data_m_na[-which(sum_check<minexp),]
						d1<-d1[-which(sum_check<minexp),]
					}
					
					 print(paste("Number of features after minexp ",minexp," filtering post-averaging",sep=""))
                                         print(dim(data_m)[1])
					
		   }
		   
		    adjdata<-try(sva::ComBat(dat=data_m_na,batch=batch_inf,mod=mod,par.prior=TRUE),silent=TRUE)
		    
		    
		    if(is(adjdata,"try-error")){
		 
		
				data_m1<-cbind(d1[,c(1:4)],data_m_na)
				
	
		    	adjdata<-MetabComBat(dat=data_m1,saminfo=sampleid_mapping,par.prior=T,filter=F,write=F,prior.plots=F)
                adjdata<-adjdata[,-c(1:4)]

#print("Done with MetabComBat")
#				print(dim(adjdata))
#				print(dim(d1))
		    }
    
    				    
		    maxint<-apply(adjdata,1,function(x){max(x,na.rm=TRUE)})
		    
		    adjdata2<-replace(as.matrix(adjdata),which(is.na(adjdata)==TRUE),missingvalue)
		    adjdata2<-replace(as.matrix(adjdata2),which(adjdata2<0),0)
		    

 			adjdata3<-cbind(d1[,c(1:12)],adjdata2)
		    adjdata2<-cbind(d1[,c(1:8)],adjdata2)
		    
		    

			
           if(summarize.replicates==TRUE){
    
     if(summary.method=="mean"){
        expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_averaged_featuretable.txt",sep="")
	}else{
		if(summary.method=="median"){
			expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt",sep="")
		}
	
	}
 }else{
     expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_featuretable.txt",sep="")
 }
            #expression_xls<-paste(subdir4b,"/ComBat_corrected_feature_tableA.txt",sep="")

		   # write.table(adjdata3,file=expression_xls,sep="\t",row.names=FALSE)
		    
		    feat.eval.result=evaluate.Features(adjdata2, numreplicates=num_replicates,min.samp.percent=min.samp.percent,alignment.tool="XCMS",impute.bool=impute.bool)
			cnames=colnames(feat.eval.result)
										
			feat.eval.result<-apply(feat.eval.result,2,as.numeric)
			feat.eval.result<-as.data.frame(feat.eval.result)
			#feat.eval.result.mat=cbind(adjdata2[,c(1:4)],feat.eval.result)  
			
			numpeaks<-apply(adjdata2[,-c(1:8)],1,countpeaks)						
		    
		    maxint<-apply(adjdata2[,-c(1:8)],1,function(x){max(x,na.rm=TRUE)})
		    
		    adjdata2<-cbind(adjdata2[,c(1:7)],numpeaks,feat.eval.result$numgoodsamples,feat.eval.result$median,feat.eval.result$Qscore,maxint,adjdata2[,-c(1:8)])
		    
		    colnames(adjdata2)<-colnames(d1)
		  
		      adjdata2$time<-round(adjdata2$time,1)
                    adjdata2$mz<-round(adjdata2$mz,4)
                    adjdata2$mzmin<-round(adjdata2$mzmin,5)
                    adjdata2$mzmax<-round(adjdata2$mzmax,5)
                    adjdata2$rtmin<-round(adjdata2$rtmin,1)
                    adjdata2$rtmax<-round(adjdata2$rtmax,1)

            #expression_xls<-paste(subdir4b,"/ComBat_corrected_feature_table.txt",sep="")
		    
            #expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_featuretable.txt",sep="")
            
            if(summarize.replicates==TRUE){
    
     if(summary.method=="mean"){
        expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_averaged_featuretable.txt",sep="")
	}else{
		if(summary.method=="median"){
			expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt",sep="")
		}
	
	}
 }else{
     expression_xls<-paste(subdir4b,"/ComBat_mzcalibrated_untargeted_featuretable.txt",sep="")
 }
            
            write.table(adjdata2,file=expression_xls,sep="\t",row.names=FALSE)
		    
		   # X<-t(adjdata2[maxint_ord[1:2000],-c(1:9)])
		   
		  
			 tic.eval(adjdata2[,-c(1:12)],outloc=subdir4b)
    

		    X<-t(adjdata2[,-c(1:12)])
		   
		     
		    X<-as.matrix(X)
		    
		    X<-replace(as.matrix(X),which(is.na(X)==TRUE),0)
		    
		    try(pca.eval(X,samplelabels=batch_labels,filename="ComBat",ncomp=5,center=TRUE,scale=TRUE,legendlocation="topright",legendcex=0.5,outloc=subdir4b),silent=TRUE)
			    
		    
		    #eval.reference.metabs(dataA=adjdata2,stdData=stddata,mzthresh=10,outloc=getwd(),folderheader="ComBat")
			 #eval.target.mz(dataA=adjdata2,refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=10,outloc=subdir4b,folderheader="ComBat")
			
			
			targeted_feat_combat<-eval.target.calibratedmz(dataA=adjdata2[,-c(2:3,5:12)],refMZ=stddata,feature.eval.result=feat.eval.result,mzthresh=refMZ.mz.diff,timethresh=refMZ.time.diff,outloc=subdir4b,folderheader="ComBat",xMSanalyzer.outloc=subdir4b)
			
				s4b<-"Stage 4b results: Batch-effect evaluation post ComBat and QC evaluation of targeted data post batch-effect correction"
		sm<-rbind(sm,s4b)
		
	dev.off()
	
		}else{
			
			adjdata2=d1
			targeted_feat_combat<-targeted_feat_raw
			
			adjdata2<-replace(as.matrix(adjdata2),which(is.na(adjdata2)==TRUE),missingvalue)
			adjdata2<-replace(as.matrix(adjdata2),which(adjdata2<0),0)
		    
			
			
			}
    
	}
			 
			 
			  	 metlin.res={}
                      kegg.res={}
                      
					 #length(union_list[[num_pairs]]$mz
					 if(is.na(adduct.list)==FALSE){
					
						cat("\n")
					  print("*********Stage 5: Mapping m/z values to known metabolites*********")
					  cat("\n")
					 
                     #annot.res<-feat.batch.annotation.KEGG(adjdata2,mz.tolerance.dbmatch,adduct.list,subdir5, numnodes=numnodes,syssleep=syssleep)
                        
                        
                        for(db in db_name){
                            annot.res<-simpleAnnotation(dataA=adjdata2,max.mz.diff=mz.tolerance.dbmatch,num_nodes=numnodes,queryadductlist=adduct.list,
                            gradienttype="Acetonitrile",mode=charge_type,outloc=subdir5,db_name=db)
                            
                            #fname<-paste("DBmatches_",db,".txt",sep="")
                             
                              fname<-paste(subdir5,"/DBmatches_",db,".txt",sep="")
                             #fname<-paste(xMSanalyzer.outloc,"/DBmatches_",db,".txt",sep="")
                            write.table(annot.res,file=fname,sep="\t",row.names=FALSE)
                            
                        }
                        
					    annotres_list<-annot.res
					    
					     s5<-"Stage 5 results: Annotation of features"
					    sm<-rbind(sm,s5)

                     }
				
			 

			
			
									
			
			
			
			  print("*************Processing complete**********")

			  #print("*********Characterizing metabolites*********")
		}             
                else
                {
                        stop(paste("No files exist in",XCMS.outloc, "Please check the input value for cdfloc", sep=""))
                }
                    
                
                
        }
	 
	 suppressWarnings(sink(file=NULL))
     #sm<-"The results include pairwise Pearson correlation within technical replicates; 2) file with raw m/z and time for each feature; 3) raw (non batch-effect corrected) feature table after m/z calibration; 4) ComBat processed feature table after m/z calibration;  5) Simple (m/z based) database searches"

sm<-as.data.frame(sm)
colnames(sm)<-"Output_description"
setwd(xMSanalyzer.outloc)
write.table(sm,file="Readme.txt",sep="\t",row.names=FALSE)
	
	 print(paste("Processing is complete. Program results can be found at: ",xMSanalyzer.outloc,sep=""))
		
	 return(list("XCMS.merged.res"=union_list, "XCMS.ind.res"=data_rpd_all,"XCMS.ind.res.filtered"=data_rpd, "final.feat.table.annot"=annotres_list, "feat.eval.ind"=feateval_list, "sample.eval.ind"=rsqres_list,"final.feat.table.raw"=d1,"final.feat.table.combat"=adjdata2,
	 "final.targeted.feat.table.raw"=targeted_feat_raw,"final.targeted.feat.table.combat"=targeted_feat_combat))
}
                                                                                     
                                     
############################################################
xMSwrapper<-function()
{
    print("Usage: xMSwrapper.apLCMS or xMSwrapper.XCMS.centWave or xMSwrapper.XCMS.matchedFilter")
}


getpeaktable<-function(j,files,n.nodes = 4, min.exp = 2, 

    min.pres = 0.5, min.run = 12, mz.tol = 1e-05, baseline.correct.noise.percentile = 0.25, 

    shape.model = "bi-Gaussian", baseline.correct = NA, peak.estim.method = "moment", 

    min.bw = NA, max.bw = NA, sd.cut = c(1, 60), sigma.ratio.lim = c(0.33, 

        3), subs = NULL,folder)

    {

	setwd(folder)
	    

	    suf.prof <- paste(min.pres, min.run, mz.tol, baseline.correct, 

        sep = "_")

    suf <- paste(suf.prof, shape.model, sd.cut[1], sd.cut[2], 

        sep = "_")

    if (shape.model == "bi-Gaussian") 

        suf <- paste(suf, sigma.ratio.lim[1], sigma.ratio.lim[2], 

            sep = "_")


	    

    		    {

                  this.name <- paste(strsplit(tolower(files[j]), 

                    "\\.")[[1]][1], suf, min.bw, max.bw, ".feature", 

                    sep = "_")

                  this.feature <- NA

                  that.name <- paste(strsplit(tolower(files[j]), 

                    "\\.")[[1]][1], suf.prof, ".profile", sep = "_")

                  processable <- "goodgood"

                  processable <- try(this.prof <- proc.cdf(files[j], 

                    min.pres = min.pres, min.run = min.run, tol = mz.tol, 

                    baseline.correct = baseline.correct, baseline.correct.noise.percentile = baseline.correct.noise.percentile, 

                    do.plot = FALSE))

                  if (substr(processable, 1, 5) == "Error") {

                    file.copy(from = files[j], to = "error_files")

                    file.remove(files[j])

                  }

                  else {

			#print("here")
                    save(this.prof, file = that.name)

                  }

                  if (substr(processable, 1, 5) != "Error") {

                    processable.2 <- "goodgood"

                    processable.2 <- try(this.feature <- prof.to.features(this.prof, 

                      min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, 

                      shape.model = shape.model, estim.method = peak.estim.method, 

                      do.plot = FALSE))

                    if (substr(processable.2, 1, 5) == "Error") {

                      file.copy(from = files[j], to = "error_files")

                      file.remove(files[j])

                      this.feature <- NA

                    }

                    else {

                      save(this.feature, file = this.name)

                    }

                  }

                }

		return(this.feature)

    }

 



    

   
   


    


    do.recover<-function(j,files,min.pres,min.run,aligned,features,f2,mz.tol=1e-5,recover.mz.range = NA, 

    recover.chr.range = NA, use.observed.range = TRUE, min.bw = NA, max.bw = NA,

                  bandwidth = 0.5, recover.min.count = 3,suf=NA)

    {

			  

	    this.name <- paste(strsplit(tolower(files[j]),  "\\.")[[1]][1], suf, ".recover", sep = "_")

                  this.recovered <- try(recover.weaker(filename = files[j], 

                    loc = j, aligned.ftrs = aligned$aligned.ftrs, 

                    pk.times = aligned$pk.times, align.mz.tol = aligned$mz.tol, 

                    align.chr.tol = aligned$chr.tol, this.f1 = features[[j]], 

                    this.f2 = f2[[j]], mz.range = recover.mz.range, 

                    chr.range = recover.chr.range, use.observed.range = use.observed.range, 

                    orig.tol = mz.tol, min.bw = min.bw, max.bw = max.bw, 

                    bandwidth = 0.5, recover.min.count = recover.min.count),silent=TRUE)

			if(is(this.recovered,"try-error")){
	
				print("Error in recovery. Skipping to next file")
				this.recovered<-{}
				
				save(this.recovered,file=this.fname)
		
			}
                  save(this.recovered, file = this.name)

	    return(this.recovered)

    }


 cdf.to.ftr.linux<-function (folder, file.pattern = ".cdf", n.nodes = 4,min.exp = 2, 

    min.pres = 0.5, min.run = 12, mz.tol = 1e-05, baseline.correct.noise.percentile = 0.25, 

    shape.model = "bi-Gaussian", baseline.correct = NA, peak.estim.method = "moment", 

    min.bw = NA, max.bw = NA, sd.cut = c(1, 60), sigma.ratio.lim = c(0.33, 

        3), subs = NULL, align.mz.tol = NA, align.chr.tol = NA, 

    max.align.mz.diff = 0.01, pre.process = FALSE, recover.mz.range = NA, 

    recover.chr.range = NA, use.observed.range = TRUE, recover.min.count = 3,reference_sample=NA) 

{

   
    setwd(folder)


    files <- dir(pattern = file.pattern, ignore.case = TRUE)

    files <- files[order(files)]

    if (!is.null(subs)) {

        if (!is.na(subs[1])) 

            files <- files[subs]

    }

    

    s1<-1:length(files)

    dir.create("error_files")

    message("***************************** prifiles --> feature lists *****************************")

    suf.prof <- paste(min.pres, min.run, mz.tol, baseline.correct, 

        sep = "_")

    suf <- paste(suf.prof, shape.model, sd.cut[1], sd.cut[2], 

        sep = "_")

    if (shape.model == "bi-Gaussian") 

        suf <- paste(suf, sigma.ratio.lim[1], sigma.ratio.lim[2], 

            sep = "_")

    
    to.do <- paste(matrix(unlist(strsplit(tolower(files), "\\.")), 

        nrow = 2)[1, ], suf, min.bw, max.bw, ".feature", sep = "_")

    to.do <- which(!(to.do %in% dir()))
    
    rawprof.names <- paste(unlist(strsplit(tolower(files), "\\."))[seq(1, 
            2 * length(files), by = 2)], "_", min.run, "_", min.pres, 
            "_", mz.tol, ".rawprof", sep = "")

   to.do2 <- which(!(rawprof.names %in% dir()))
   
   to.do<-c(to.do,to.do2)
   to.do<-unique(to.do)
   
    message(c("number of files to process: ", length(to.do)))

    if (length(to.do) > 0) {

	s1_temp<-s1[to.do]

	features<-mclapply(s1_temp,getpeaktable,files,n.nodes = n.nodes, min.exp = min.exp, 

    min.pres = min.pres, min.run = min.run, mz.tol = mz.tol, baseline.correct.noise.percentile = baseline.correct.noise.percentile, 

    shape.model =  shape.model, baseline.correct = baseline.correct, peak.estim.method = peak.estim.method, 

    min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, sigma.ratio.lim = sigma.ratio.lim, subs = subs,folder=folder)
	
           

    }

    all.files <- dir()

    sel <- which(files %in% all.files)

    files <- files[sel]

    features <- new("list")

    for (i in 1:length(files)) {

        this.name <- paste(strsplit(tolower(files[i]), "\\.")[[1]][1], 

            suf, min.bw, max.bw, ".feature", sep = "_")

        cat(this.name, " ")

        load(this.name)

        features[[i]] <- this.feature

    }

    gc()

    if (!pre.process) {

        message("****************************** time correction ***************************************")

        suf <- paste(suf, align.mz.tol, align.chr.tol, subs[1], 

            subs[length(subs)], sep = "_")

        this.name <- paste("time_correct_done_", suf, ".bin", 

            sep = "")

        all.files <- dir()

        is.done <- all.files[which(all.files == this.name)]

        if (length(is.done) == 0) {

            message(c("***** correcting time, CPU time (seconds) ", 

                as.vector(system.time(f2 <- adjust.time(features, 

                  mz.tol = align.mz.tol, chr.tol = align.chr.tol, 

                  find.tol.max.d = 10 * mz.tol, max.align.mz.diff = max.align.mz.diff)))[1]))

            save(f2, file = this.name)

        }

        else {

            load(this.name)

        }

        gc()

        message("****************************  aligning features **************************************")

        suf <- paste(suf, min.exp, sep = "_")

        this.name <- paste("aligned_done_", suf, ".bin", sep = "")

        all.files <- dir()

        is.done <- all.files[which(all.files == this.name)]

        if (length(is.done) == 0) {

            message(c("***** aligning features, CPU time (seconds): ", 

                as.vector(system.time(aligned <- feature.align(f2, 

                  min.exp = min.exp, mz.tol = align.mz.tol, chr.tol = align.chr.tol, 

                  find.tol.max.d = 10 * mz.tol, max.align.mz.diff = max.align.mz.diff)))[1]))

            save(aligned, file = this.name)

        }

        else {

            load(this.name)

        }

        gc()

        message("**************************** recovering weaker signals *******************************")

        suf <- paste(suf, recover.mz.range, recover.chr.range, 

            use.observed.range, sep = "_")

	#vt_150912_128_NA_.recover
       
	#suf="NA"
	 worklist <- paste(matrix(unlist(strsplit(tolower(files), "\\.")), nrow = 2)[1, ], suf, ".recover",  sep = "_")

        to.do <- which(!(worklist %in% dir()))

        grps <- round(seq(0, length(to.do), length = n.nodes +  1))

        grps <- unique(grps)

        message(c("number of files to process: ", length(to.do)))

	if(length(to.do)>0){
	s1_temp<-s1[to.do]

        features.recov <-mclapply(s1_temp,do.recover,files,min.pres,min.run,aligned,features,f2,mz.tol=mz.tol,recover.mz.range=recover.mz.range,

                  recover.chr.range = recover.chr.range, use.observed.range = use.observed.range,

                  min.bw = min.bw, max.bw = max.bw,

                  bandwidth = 0.5, recover.min.count = recover.min.count,suf=suf,mc.cores=n.nodes,mc.preschedule=FALSE)

	}

        gc()

        new.aligned <- aligned

          for (i in 1:length(files)) {

            this.name <- paste(strsplit(tolower(files[i]), "\\.")[[1]][1], 

                suf, ".recover", sep = "_")
		load(this.name)
	  
	    
	    #print(length(this.recovered))
	if(length(this.recovered)>1){
            new.aligned$aligned.ftrs[, i + 4] <- this.recovered$this.ftrs

            new.aligned$pk.times[, i + 4] <- this.recovered$this.times

            new.aligned$features[[i]] <- this.recovered$this.f1

            new.aligned$f2[[i]] <- this.recovered$this.f2
	}else{
		print("File did not pass the recovery step.")
		print(i)
		print(files[i])
	}
            gc()

        }

        rec <- new("list")

        colnames(aligned$aligned.ftrs) <- colnames(aligned$pk.times) <- colnames(new.aligned$aligned.ftrs) <- colnames(new.aligned$pk.times) <- c("mz", 

            "time", "mz.min", "mz.max", files)

        rec$features <- new.aligned$features

        rec$features2 <- new.aligned$f2

        rec$aligned.ftrs <- aligned$aligned.ftrs

        rec$pk.times <- aligned$pk.times

        rec$final.ftrs <- new.aligned$aligned.ftrs

        rec$final.times <- new.aligned$pk.times

        rec$align.mz.tol <- new.aligned$mz.tol

        rec$align.chr.tol <- new.aligned$chr.tol

        rec$mz.tol <- mz.tol

       aligned<-rec
	rm(rec)
	rm(new.aligned)


        return(aligned)

    }

}

cdf.to.ftr.windows<-function (folder, file.pattern = ".cdf", n.nodes = 4, min.exp = 2, 
    min.pres = 0.5, min.run = 12, mz.tol = 1e-05, baseline.correct.noise.percentile = 0.25, 
    shape.model = "bi-Gaussian", baseline.correct = NA, peak.estim.method = "moment", 
    min.bw = NA, max.bw = NA, sd.cut = c(1, 60), sigma.ratio.lim = c(0.33, 
        3), subs = NULL, align.mz.tol = NA, align.chr.tol = NA, 
    max.align.mz.diff = 0.01, pre.process = FALSE, recover.mz.range = NA, 
    recover.chr.range = NA, use.observed.range = TRUE, recover.min.count = 3,reference_sample=NA) 
{
    library(mzR)
    library(doParallel)
    setwd(folder)
    cl <- makePSOCKcluster(n.nodes, error = recover)
    registerDoParallel(cl)
    clusterEvalQ(cl, library(apLCMS))
    files <- dir(pattern = file.pattern, ignore.case = TRUE)
    files <- files[order(files)]
    if (!is.null(subs)) {
        if (!is.na(subs[1])) 
            files <- files[subs]
    }
    dir.create("error_files")
    message("***************************** prifiles --> feature lists *****************************")
    suf.prof <- paste(min.pres, min.run, mz.tol, baseline.correct, 
        sep = "_")
    suf <- paste(suf.prof, shape.model, sd.cut[1], sd.cut[2], 
        sep = "_")
    if (shape.model == "bi-Gaussian") 
        suf <- paste(suf, sigma.ratio.lim[1], sigma.ratio.lim[2], 
            sep = "_")
    to.do <- paste(matrix(unlist(strsplit(tolower(files), "\\.")), 
        nrow = 2)[1, ], suf, min.bw, max.bw, ".feature", sep = "_")
	
    to.do <- which(!(to.do %in% dir()))
      rawprof.names <- paste(unlist(strsplit(tolower(files), "\\."))[seq(1, 
            2 * length(files), by = 2)], "_", min.run, "_", min.pres, 
            "_", mz.tol, ".rawprof", sep = "")

   to.do2 <- which(!(rawprof.names %in% dir()))
   
   to.do<-c(to.do,to.do2)
   to.do<-unique(to.do)
    
    message(c("number of files to process: ", length(to.do)))
    if (length(to.do) > 0) {
        grps <- round(seq(0, length(to.do), length = n.nodes + 
            1))
        grps <- unique(grps)
        features <- foreach(i = 2:length(grps), .export = ls(envir = globalenv())) %dopar% 
            {
                this.subset <- to.do[(grps[i - 1] + 1):grps[i]]
                for (j in this.subset) {
                  this.name <- paste(strsplit(tolower(files[j]), 
                    "\\.")[[1]][1], suf, min.bw, max.bw, ".feature", 
                    sep = "_")
                  this.feature <- NA
                  that.name <- paste(strsplit(tolower(files[j]), 
                    "\\.")[[1]][1], suf.prof, ".profile", sep = "_")
                  processable <- "goodgood"
                  processable <- try(this.prof <- proc.cdf(files[j], 
                    min.pres = min.pres, min.run = min.run, tol = mz.tol, 
                    baseline.correct = baseline.correct, baseline.correct.noise.percentile = baseline.correct.noise.percentile, 
                    do.plot = FALSE))
                  if (substr(processable, 1, 5) == "Error") {
                    file.copy(from = files[j], to = "error_files")
                    file.remove(files[j])
                  }
                  else {
                    save(this.prof, file = that.name)
                  }
                  if (substr(processable, 1, 5) != "Error") {
                    processable.2 <- "goodgood"
                    processable.2 <- try(this.feature <- prof.to.features(this.prof, 
                      min.bw = min.bw, max.bw = max.bw, sd.cut = sd.cut, 
                      shape.model = shape.model, estim.method = peak.estim.method, 
                      do.plot = FALSE))
                    if (substr(processable.2, 1, 5) == "Error") {
                      file.copy(from = files[j], to = "error_files")
                      file.remove(files[j])
                      this.feature <- NA
                    }
                    else {
                      save(this.feature, file = this.name)
                    }
                  }
                }
                1
            }
    }
    all.files <- dir()
    sel <- which(files %in% all.files)
    files <- files[sel]
    features <- new("list")
    for (i in 1:length(files)) {
        this.name <- paste(strsplit(tolower(files[i]), "\\.")[[1]][1], 
            suf, min.bw, max.bw, ".feature", sep = "_")
        cat(this.name, " ")
        load(this.name)
        features[[i]] <- this.feature
    }
    gc()
    if (!pre.process) {
        message("****************************** time correction ***************************************")
        suf <- paste(suf, align.mz.tol, align.chr.tol, subs[1], 
            subs[length(subs)], sep = "_")
        this.name <- paste("time_correct_done_", suf, ".bin", 
            sep = "")
        all.files <- dir()
        is.done <- all.files[which(all.files == this.name)]
        if (length(is.done) == 0) {
            message(c("***** correcting time, CPU time (seconds) ", 
                as.vector(system.time(f2 <- adjust.time(features, 
                  mz.tol = align.mz.tol, chr.tol = align.chr.tol, 
                  find.tol.max.d = 10 * mz.tol, max.align.mz.diff = max.align.mz.diff)))[1]))
            save(f2, file = this.name)
        }
        else {
            load(this.name)
        }
        gc()
        message("****************************  aligning features **************************************")
        suf <- paste(suf, min.exp, sep = "_")
        this.name <- paste("aligned_done_", suf, ".bin", sep = "")
        all.files <- dir()
        is.done <- all.files[which(all.files == this.name)]
        if (length(is.done) == 0) {
            message(c("***** aligning features, CPU time (seconds): ", 
                as.vector(system.time(aligned <- feature.align(f2, 
                  min.exp = min.exp, mz.tol = align.mz.tol, chr.tol = align.chr.tol, 
                  find.tol.max.d = 10 * mz.tol, max.align.mz.diff = max.align.mz.diff)))[1]))
            save(aligned, file = this.name)
        }
        else {
            load(this.name)
        }
        gc()
        message("**************************** recovering weaker signals *******************************")
        suf <- paste(suf, recover.mz.range, recover.chr.range, 
            use.observed.range, sep = "_")
        worklist <- paste(matrix(unlist(strsplit(tolower(files), 
            "\\.")), nrow = 2)[1, ], suf, ".recover", 
            sep = "_")
        to.do <- which(!(worklist %in% dir()))
        grps <- round(seq(0, length(to.do), length = n.nodes + 
            1))
        grps <- unique(grps)
        message(c("number of files to process: ", length(to.do)))
	if(length(to.do)>0){
        features.recov <- foreach(i = 2:length(grps), .export = ls(envir = globalenv())) %dopar% 
            {
                this.subset <- to.do[(grps[i - 1] + 1):grps[i]]
                for (j in this.subset) {
                  this.name <- paste(strsplit(tolower(files[j]), 
                    "\\.")[[1]][1], suf, ".recover", sep = "_")
                  this.recovered <- try(recover.weaker(filename = files[j], 
                    loc = j, aligned.ftrs = aligned$aligned.ftrs, 
                    pk.times = aligned$pk.times, align.mz.tol = aligned$mz.tol, 
                    align.chr.tol = aligned$chr.tol, this.f1 = features[[i]], 
                    this.f2 = f2[[i]], mz.range = recover.mz.range, 
                    chr.range = recover.chr.range, use.observed.range = use.observed.range, 
                    orig.tol = mz.tol, min.bw = min.bw, max.bw = max.bw, 
                    bandwidth = 0.5, recover.min.count = recover.min.count),silent=TRUE)
		    
		    if(is(this.recovered,"try-error")){
	
				print("Error in recovery. Skipping to next file")
				this.recovered<-{}
			}
			
                  save(this.recovered, file = this.name)
                }
            }
	    }
        gc()
        new.aligned <- aligned
        for (i in 1:length(files)) {
            this.name <- paste(strsplit(tolower(files[i]), "\\.")[[1]][1], 
                suf, ".recover", sep = "_")
            load(this.name)
  
	    if(length(this.recovered)>1){
            new.aligned$aligned.ftrs[, i + 4] <- this.recovered$this.ftrs

            new.aligned$pk.times[, i + 4] <- this.recovered$this.times

            new.aligned$features[[i]] <- this.recovered$this.f1

            new.aligned$f2[[i]] <- this.recovered$this.f2
	}else{
		print("File did not pass the recovery step.")
		print(i)
		print(files[i])
	}
            gc()
	    
	    
        }
        rec <- new("list")
        colnames(aligned$aligned.ftrs) <- colnames(aligned$pk.times) <- colnames(new.aligned$aligned.ftrs) <- colnames(new.aligned$pk.times) <- c("mz", 
            "time", "mz.min", "mz.max", files)
        rec$features <- new.aligned$features
        rec$features2 <- new.aligned$f2
        rec$aligned.ftrs <- aligned$aligned.ftrs
        rec$pk.times <- aligned$pk.times
        rec$final.ftrs <- new.aligned$aligned.ftrs
        rec$final.times <- new.aligned$pk.times
        rec$align.mz.tol <- new.aligned$mz.tol
        rec$align.chr.tol <- new.aligned$chr.tol
        rec$mz.tol <- mz.tol
        stopImplicitCluster()
        #return(rec)
	
	aligned<-rec
	rm(rec)
	rm(new.aligned)


        return(aligned)
    }
}

apLCMS.EIC.plot<-function(aligned, rows = NA, colors = NA, transform = "none",
    subset = NA, mz.list=NA,time.list=NA,minrt=NA, maxrt=NA, min.run=12, 
    min.pres=0.5, max.spline.time.points = 1000,
    chem.names=NA,rawprofileloc,cvvec=NA,plotEIC=TRUE,numnodes=NA)
{
    if (!(transform %in% c("none", "log", "sqrt", "cuberoot"))) {
        message("Invalid transformation. It has to be from: none, log, sqrt, and cuberoot")
        break
    }
	setwd(rawprofileloc)

	if(is.na(numnodes)==TRUE){
		numnodes<-detectCores()-1
	}
	 
    peak_score_all<-{}
    #print(rows)

	print("length rows")
	print(length(rows))
    if (!is.na(rows[1])) {

        library(splines)
        num.exp <- nrow(summary(aligned$features))
	
	  files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
	 if(is.na(mz.list)==TRUE){
	 
		mz.list<-aligned$final.ftrs[rows,1]
	 }
	  if(is.na(time.list)==TRUE){
	 
		time.list<-aligned$final.ftrs[rows,2]
	 }
	
	if(length(files)<1){
		aligned$final.ftrs<-aligned$aligned.ftrs
		 files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
		aligned$aligned.ftrs<-NULL
	}
	if(is.na(minrt)==TRUE){
		minrt=min(aligned$final.ftrs[,2],na.rm=TRUE)
	}
	if(is.na(maxrt)==TRUE){
                maxrt=max(aligned$final.ftrs[,2],na.rm=TRUE)
        }

        if (is.na(subset[1]))
            subset <- 1:num.exp
        if (is.na(colors))
            colors <- rainbow(length(subset)) #
	   
     
        rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1,
            2 * num.exp, by = 2)], "_", min.run, "_", min.pres,
            "_", aligned$mz.tol, ".rawprof", sep = "")
	rawprof.names<-tolower(rawprof.names)
        adj.times <- new("list")
        #for (i in subset) 
	adj.times<-lapply(1:length(subset),function(j)
	{
		i<-subset[j]
            load(rawprof.names[i])
            times <- unique(raw.prof$labels)
            times <- times[order(times)]
            orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][,
                3]), 2]
            adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][,
                3]), 2]
            orig.time <- round(orig.time, 8)
            adjusted.time <- round(adjusted.time, 8)
            ttt <- table(adjusted.time)
            ttt2 <- table(orig.time)
            to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt ==
                1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 ==
                1]))
            if (length(to.use) > max.spline.time.points)
                to.use <- sample(to.use, max.spline.time.points,
                  replace = FALSE)
            orig.time <- orig.time[to.use]
            adjusted.time <- adjusted.time[to.use]
            sp <- interpSpline(adjusted.time ~ orig.time)
            #adj.times[[i]] <- 
	    return(predict(sp, times)$y)
        })
	
	#print("length rows")
	#print(length(rows))
	
#for (n in 1:length(rows)) {
peak_score_all<-lapply(1:length(rows),function(n){	

	#print("New feature")
	if(FALSE){
		cl<-parallel::makeCluster(numnodes)
	
		
			clusterEvalQ(cl, "get_peakscore")
			clusterExport(cl, "get_peakscore")
		clusterEvalQ(cl, "compare_intensities_ttest")
	}			
#peak_score_all<-parLapply(cl,1:length(rows),function(n){ 
            if (n%%4 == 1) {
                par(mfrow = c(2, 2))
            }
            this.row <- rows[n]
            this.intensi <- aligned$final.ftrs[this.row, 5:(num.exp +
                4)]
            this.times <- aligned$final.times[this.row, 5:(num.exp +
                4)]
		#print("this.times")
		#print(this.times)
		med_rt<-median(this.times[this.times>0],na.rm=TRUE)
		range_rt<-range(this.times[this.times>0],na.rm=TRUE)
		
            this.mz <- aligned$final.ftrs[this.row, 1]
           max_int_index<-order(this.intensi,decreasing=TRUE)
	  # print(aligned$final.ftrs[this.row, 1:2])
	if(length(max_int_index)>5)
	  {
	  
	        set.seed(555)
		rand_set<-sample(x=max_int_index[-c(1:3)],size=3)
		subset<-c(max_int_index[1:3],rand_set) #,max_int_index[length(max_int_index)])
		subset<-unique(subset)
	  }
	
	   #print(this.row)
	    to.plot <- new("list")
            y.max <- 0
            x.max <- 0
            x.min <- Inf
	    good_i<-{}
	    sel_l<-new("list")
                      for (iii in 1:length(subset)) {
                i <- subset[iii]
		to.plot[[i]]<-{}
		sel_l[[i]]<-{}
                #if (this.intensi[i] != 0) 
		{
                  load(rawprof.names[i])
                  times <- unique(raw.prof$labels)
                  times <- times[order(times)]
                  mz.diff <- abs(aligned$features2[[i]][, 1] -
                    this.mz)
                  time.diff <- abs(aligned$features2[[i]][, 2] -
                    this.times[i])
                  sel <- which(time.diff < 1e-17)
                  if (length(sel) > 1)
                    sel <- sel[which(mz.diff[sel] == min(mz.diff[sel]))[1]]
                  if (length(sel) == 0) {
                    mz.lim <- aligned$final.ftrs[this.row, c(3,
                      4)]
                    sel <- which(aligned$features2[[i]][, 1] >=
                      mz.lim[1] & aligned$features2[[i]][, 1] <=
                      mz.lim[2] & !is.na(aligned$features2[[i]][,
                      6]))
                    sub.features <- aligned$features2[[i]][sel,
                      ]
		      print(length(sub.features))
		      if(length(sub.features)>7){
			sub.time.diff <- abs(sub.features[, 2] - this.times[i])
                    }else{
		    sub.time.diff <- abs(sub.features[2] - this.times[i])
		    }
		    sel <- sel[which(sub.time.diff == min(sub.time.diff))][1]
                  }
                  target.mz <- aligned$features2[[i]][sel, 1]
                  target.time <- aligned$features[[i]][sel, 2]
                  time.adjust <- aligned$features2[[i]][sel,
                    2] - aligned$features[[i]][sel, 2]
                  mz.diff <- abs(raw.prof$masses - target.mz)
                  sel.slice <- raw.prof$grps[which(mz.diff ==
                    min(mz.diff))[1]]
                  sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                    sel.slice])
			#print("sel.time.range is")
		   # print(sel.time.range)
		    if(is.na(sel.time.range[1])==TRUE){
			next;
		    }
                  while (target.time < sel.time.range[1] | target.time >
                    sel.time.range[2]) {
                    mz.diff[raw.prof$grps == sel.slice] <- 100
                    sel.slice <- raw.prof$grps[which(mz.diff ==
                      min(mz.diff))[1]]
                    sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                      sel.slice])
                  }
                  sel.time <- raw.prof$labels[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- raw.prof$intensi[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- sel.intensi[order(sel.time)]
                  sel.time <- sel.time[order(sel.time)]
                  all.intensi <- times * 0
                  all.intensi[times %in% sel.time] <- sel.intensi
                  if(i<length(adj.times)){
		  all.times <- adj.times[[i]]
			if(length(which(is.na(all.times)==TRUE))>0){
				sn=NA
			}
		  }else{
			all.times <-{}
			sn=NA
                                #return(sn)
		  }
		 
                  if (transform == "log")
                    all.intensi <- log10(all.intensi + 1)
                  if (transform == "sqrt")
                    all.intensi <- sqrt(all.intensi)
                  if (transform == "cuberoot")
                    all.intensi <- all.intensi^(1/3)
		   #if (max(all.intensi,na.rm=TRUE) > y.max)
		    #	y.max<-max(all.intensi)
                  if (max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE) > y.max)
			y.max <- max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE)
		#if (max(all.intensi[which(all.times>(target.time-20) & all.times<(target.time+20))],na.rm=TRUE) > y.max){
			#y.max <- max(all.intensi,na.rm=TRUE) #max(all.intensi[which(all.times<(target.time+30) & all.times>(target.time-30))],na.rm=TRUE)
			
			#y.max <-max(all.intensi[which(all.times>(target.time-20) & all.times<(target.time+20))],na.rm=TRUE)
		# }
                  if (max(all.times[all.intensi > 0],na.rm=TRUE) > x.max)
                    x.max <- max(all.times[all.intensi > 0],na.rm=TRUE)
                  if (min(all.times[all.intensi > 0],na.rm=TRUE) < x.min)
                    x.min <- min(all.times[all.intensi > 0],na.rm=TRUE)
                  to.plot[[i]] <- cbind(all.times, all.intensi)
		  good_i<-c(good_i,i)
                  sel_l[[i]]<-sel
		}
            }
	    
	    max_int_vec<-{}
	    peak_score_vec<-{}
	    bool_plot=0
	   # print("ymax is")
	   # print(y.max)
	 
	    y_max_final<-y.max+10000
	   # print( nrow(summary(to.plot)))
	    #print(to.plot)
	  #  print(length(to.plot))
	#    print(minrt)
	#    print(maxrt)
	    minrt<-min(c(0,minrt),na.rm=TRUE)
	    
	     peak_score<-(-1)
	    peak_score_vec<-c(-1)
	    
	    if(length(to.plot)>0){
	    
	      peak_score_vec<-{}
            for (iii in 1:min(length(subset), nrow(summary(to.plot)))) 
		{
		peak_score<-0
                i <- subset[iii]
		curtime<-NA
		#if(i<length(sel_l))
		#if (this.intensi[i] != 0)
		if(is.na(this.times[i])==FALSE)
		{
		
	#	print("length of to.plot")
		#print(length(to.plot))
		#print(head(to.plot))
		#print(good_i)
		#print(i)
				#
		#if(length((to.plot[[i]]))>1)
		#if(i%in%good_i)
		test_file<-try(to.plot[[i]],silent=TRUE)
		
		   if(is(test_file,"try-error") || (is(test_file,"NULL")))
		{
			next;
		}else{
		sel<-sel_l[[i]]
		if(dim(to.plot[[i]])[2]>1){
		time_v1<-to.plot[[i]][,1]
			int_v1<-to.plot[[i]][,2]
			if(length(which(sel>dim(aligned$features[[i]])[1]))<1){
			target.time <- aligned$features[[i]][sel, 2]
			time.adjust <- aligned$features2[[i]][sel,2]
			curtime<-time.list[which(mz.list==this.mz)[1]]

			cur_cv<-round(cvvec[which(mz.list==this.mz)],1)
			max_int<-max(to.plot[[i]][,2],na.rm=TRUE)
			d1<-density(time_v1,na.rm=TRUE)
			max_int_index<-which(int_v1==max_int)
			peak_time<-time_v1[max_int_index]
			baseline_range<-20*peak_time
			#curtime<-time.adjust
			#print("getting peak score")
			#print(summary(int_v1))
			#print(summary(time_v1))
			#print(time.adjust)
			#print(target.time)
			#print(curtime)
			peak_score<-get_peakscore(int_v1=int_v1,time_v1=time_v1,cur.time=curtime,bw=d1$bw,rt_range=range_rt)
			
			
			}else{
				peak_score<-NA
				curtime<-time.list[which(mz.list==this.mz)[1]]

			cur_cv<-round(cvvec[which(mz.list==this.mz)],1)
			max_int<-max(to.plot[[i]][,2],na.rm=TRUE)
			d1<-density(time_v1,na.rm=TRUE)
			max_int_index<-which(int_v1==max_int)
			peak_time<-time_v1[max_int_index]
			baseline_range<-20*peak_time
			}
			if(is.na(chem.names[1])==FALSE){
                        curname<-chem.names[which(mz.list==this.mz)]
			}else{
				curname<-""
			}
			if(length(peak_score_vec)<1){
				bool_plot=0
			}
			
			max_int_vec<-c(max_int_vec,max_int)	
			peak_score_vec<-c(peak_score_vec,peak_score)	
			minrt=max(0,(curtime-2*d1$bw))
			maxrt=min(600,(curtime+2*d1$bw))
			
			colors[iii][which(is.na(colors[iii])==TRUE)]<-"gray"
			#if(FALSE)
			y.maxt<-y.max #+10000
			if(plotEIC==TRUE)
			{
		
	    y.max<-y.maxt+10000
                if (bool_plot == 0) {
		bool_plot=1
		  plot(to.plot[[i]], xlim = c(minrt, maxrt),
                    ylim = c(0,  y_max_final), type = "l", col = colors[iii],
                    xlab = "retention time", ylab = "intensity",
                    main = paste("m/z",
                      round(this.mz, 5)," time:",round(curtime,1),"\n", " RT range: ",paste(round(range_rt,0),collapse=" to "),"\n", curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
			}else {
			if(bool_plot==1){
				lines(to.plot[[i]], col = colors[iii])
			}else{
				bool_plot=1
				  plot(to.plot[[i]], xlim = c(minrt, maxrt),
                    ylim = c(0,   y_max_final), type = "l", col = colors[iii],
                    xlab = "retention time", ylab = "intensity",
                    main = paste("EIC {",min.run,",",min.pres,"} m/z",
                      #round(this.mz, 5)," time:", round(curtime,1),"\n",curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
		        round(this.mz, 5)," time:",round(curtime,1),"\n", " RT range: ",paste(round(range_rt,0),collapse=" to "),"\n", curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
				}
                	}
			
			}
			
			
		}
                }
		}
	    }
	    }
	  
	    peak_score_vec<-max(peak_score_vec,na.rm=TRUE) #dim(aligned$features[[i]])[1] #
	
	 max_int_all<-max(max_int_vec,na.rm=TRUE)
	    peak_score_all<-c(peak_score_all,peak_score_vec)
	    if(bool_plot==1){
	    	#mtext(paste("peak score:",round(peak_score_vec,2),";","max int (log10): ",log10(round(max_int_all)+1),sep=""),cex=0.7)
		mtext(paste("peak score:",round(peak_score_vec,2),sep=""),cex=0.7)
		}
	
	if(peak_score_vec==(-Inf)){
	
		peak_score_vec=0.5
	}
	
	if(peak_score_vec==(Inf)){
	
		peak_score_vec=10
	}
      
	})   
 }
return(peak_score_all)
}


apLCMS.get.peakscore<-function(aligned, rows = NA, colors = NA, transform = "none",
    subset = NA, mz.list=NA,time.list=NA,minrt=NA, maxrt=NA, min.run=12, min.pres=0.5,
    max.spline.time.points = 1000,chem.names=NA,rawprofileloc,cvvec=NA,plotEIC=TRUE,numnodes=NA,cvthresh=100)
{
    if (!(transform %in% c("none", "log", "sqrt", "cuberoot"))) {
        message("Invalid transformation. It has to be from: none, log, sqrt, and cuberoot")
        break
    }
	setwd(rawprofileloc)
    
    if(is.na(cvvec)==TRUE){
    
	cvvec<-rep(1,dim(aligned$final.ftrs)[1])
    }
    
    sys_name<-Sys.info()['sysname']

    if(is.na(numnodes)==TRUE){
		numnodes<-detectCores()*0.5
	}
    if(sys_name=="windows" || sys_name=="Windows"){
        
        cl<-parallel::makeCluster(numnodes)
        
        clusterEvalQ(cl, library(splines))
    }

	
	 
    peak_score_all<-{}
    #print(rows)

	#print("length rows")
	#print(length(rows))
    if (!is.na(rows[1])) {

        library(splines)
        num.exp <- nrow(summary(aligned$features))
	
	  files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
	 if(is.na(mz.list)==TRUE){
	 
		mz.list<-aligned$final.ftrs[rows,1]
	 }
	  if(is.na(time.list)==TRUE){
	 
		time.list<-aligned$final.ftrs[rows,2]
	 }
	
	if(length(files)<1){
		aligned$final.ftrs<-aligned$aligned.ftrs
		 files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
		aligned$aligned.ftrs<-NULL
	}
	if(is.na(minrt)==TRUE){
		minrt=min(aligned$final.ftrs[,2],na.rm=TRUE)
	}
	if(is.na(maxrt)==TRUE){
                maxrt=max(aligned$final.ftrs[,2],na.rm=TRUE)
        }

        if (is.na(subset[1]))
            subset <- 1:num.exp
        if (is.na(colors))
            colors <- rainbow(length(subset)) #
	   
     
        rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1,
            2 * num.exp, by = 2)], "_", min.run, "_", min.pres,
            "_", aligned$mz.tol, ".rawprof", sep = "")
	rawprof.names<-tolower(rawprof.names)
        adj.times <- new("list")
        #for (i in subset)
        if(sys_name=="windows" || sys_name=="Windows"){
	adj.times<-parLapply(cl,1:length(subset),function(j)
	{
		i<-subset[j]
            load(rawprof.names[i])
            times <- unique(raw.prof$labels)
            times <- times[order(times)]
            orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][,
                3]), 2]
            adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][,
                3]), 2]
            orig.time <- round(orig.time, 8)
            adjusted.time <- round(adjusted.time, 8)
            ttt <- table(adjusted.time)
            ttt2 <- table(orig.time)
            to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt ==
                1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 ==
                1]))
            if (length(to.use) > max.spline.time.points)
                to.use <- sample(to.use, max.spline.time.points,
                  replace = FALSE)
            orig.time <- orig.time[to.use]
            adjusted.time <- adjusted.time[to.use]
            sp <- interpSpline(adjusted.time ~ orig.time)
            #adj.times[[i]] <- 
	    return(predict(sp, times)$y)
        })
    stopCluster(cl)
        }else{
            
            adj.times<-mclapply(1:length(subset),function(j)
            {
                i<-subset[j]
                load(rawprof.names[i])
                times <- unique(raw.prof$labels)
                times <- times[order(times)]
                orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][,
                3]), 2]
                adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][,
                3]), 2]
                orig.time <- round(orig.time, 8)
                adjusted.time <- round(adjusted.time, 8)
                ttt <- table(adjusted.time)
                ttt2 <- table(orig.time)
                to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt ==
                1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 ==
                1]))
                if (length(to.use) > max.spline.time.points)
                to.use <- sample(to.use, max.spline.time.points,
                replace = FALSE)
                orig.time <- orig.time[to.use]
                adjusted.time <- adjusted.time[to.use]
                sp <- interpSpline(adjusted.time ~ orig.time)
                #adj.times[[i]] <- 
                return(predict(sp, times)$y)
            },mc.cores=numnodes,mc.preschedule=FALSE)
            
        }
	#print("length rows")
	#print(length(rows))
	
#
#peak_score_all<-lapply(1:length(rows),function(n){	

	#print("New feature")
	#if(FALSE){
		cl<-parallel::makeCluster(numnodes)
	
		
			clusterEvalQ(cl, "get_peakscore")
			clusterExport(cl, "get_peakscore")
		clusterEvalQ(cl, "compare_intensities_ttest")
	#}			
peak_score_all<-parLapply(cl,1:length(rows),function(n){
    
#for (n in 1:length(rows)) {
    if(cvvec[n]>cvthresh){
        
            return(-1)
    }
            if (n%%4 == 1) {
                par(mfrow = c(2, 2))
            }
            this.row <- rows[n]
            this.intensi <- aligned$final.ftrs[this.row, 5:(num.exp +
                4)]
            this.times <- aligned$final.times[this.row, 5:(num.exp +
                4)]
		#print("this.times")
		#print(this.times)
		med_rt<-median(this.times[this.times>0],na.rm=TRUE)
		range_rt<-range(this.times[this.times>0],na.rm=TRUE)
		
            this.mz <- aligned$final.ftrs[this.row, 1]
           max_int_index<-order(this.intensi,decreasing=TRUE)
	  
	   
	  
	   
	if(length(max_int_index)>5)
	  {
		set.seed(555)
		rand_set<-sample(x=max_int_index[-c(1:3)],size=3)
		subset<-c(max_int_index[1:3],rand_set)
		subset<-unique(subset)
	  }
	   
	 # print(this.row)
	    to.plot <- new("list")
            y.max <- 0
            x.max <- 0
            x.min <- Inf
	    good_i<-{}
	    sel_l<-new("list")
                      for (iii in 1:length(subset)) {
                i <- subset[iii]
		to.plot[[i]]<-{}
		sel_l[[i]]<-{}
                #if (this.intensi[i] != 0) 
		{
                  load(rawprof.names[i])
                  times <- unique(raw.prof$labels)
                  times <- times[order(times)]
                  mz.diff <- abs(aligned$features2[[i]][, 1] -
                    this.mz)
                  time.diff <- abs(aligned$features2[[i]][, 2] -
                    this.times[i])
                  sel <- which(time.diff < 1e-17)
                  if (length(sel) > 1)
                    sel <- sel[which(mz.diff[sel] == min(mz.diff[sel]))[1]]
                  if (length(sel) == 0) {
                    mz.lim <- aligned$final.ftrs[this.row, c(3,
                      4)]
                    sel <- which(aligned$features2[[i]][, 1] >=
                      mz.lim[1] & aligned$features2[[i]][, 1] <=
                      mz.lim[2] & !is.na(aligned$features2[[i]][,
                      6]))
                    sub.features <- aligned$features2[[i]][sel,
                      ]
		      print(length(sub.features))
		      if(length(sub.features)>7){
			sub.time.diff <- abs(sub.features[, 2] - this.times[i])
                    }else{
		    sub.time.diff <- abs(sub.features[2] - this.times[i])
		    }
		    sel <- sel[which(sub.time.diff == min(sub.time.diff))][1]
                  }
                  target.mz <- aligned$features2[[i]][sel, 1]
                  target.time <- aligned$features[[i]][sel, 2]
                  time.adjust <- aligned$features2[[i]][sel,
                    2] - aligned$features[[i]][sel, 2]
                  mz.diff <- abs(raw.prof$masses - target.mz)
                  sel.slice <- raw.prof$grps[which(mz.diff ==
                    min(mz.diff))[1]]
                  sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                    sel.slice])
			print("sel.time.range is")
		    print(sel.time.range)
		    if(is.na(sel.time.range[1])==TRUE){
			next;
		    }
                  while (target.time < sel.time.range[1] | target.time >
                    sel.time.range[2]) {
                    mz.diff[raw.prof$grps == sel.slice] <- 100
                    sel.slice <- raw.prof$grps[which(mz.diff ==
                      min(mz.diff))[1]]
                    sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                      sel.slice])
                  }
                  sel.time <- raw.prof$labels[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- raw.prof$intensi[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- sel.intensi[order(sel.time)]
                  sel.time <- sel.time[order(sel.time)]
                  all.intensi <- times * 0
                  all.intensi[times %in% sel.time] <- sel.intensi
                  if(i<length(adj.times)){
		  all.times <- adj.times[[i]]
			if(length(which(is.na(all.times)==TRUE))>0){
				sn=NA
			}
		  }else{
			all.times <-{}
			sn=NA
                                #return(sn)
		  }
		 
                  if (transform == "log")
                    all.intensi <- log10(all.intensi + 1)
                  if (transform == "sqrt")
                    all.intensi <- sqrt(all.intensi)
                  if (transform == "cuberoot")
                    all.intensi <- all.intensi^(1/3)
		   #if (max(all.intensi,na.rm=TRUE) > y.max)
		    #	y.max<-max(all.intensi)
                  if (max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE) > y.max)
			y.max <- max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE)
		#if (max(all.intensi[which(all.times>(target.time-20) & all.times<(target.time+20))],na.rm=TRUE) > y.max){
			#y.max <- max(all.intensi,na.rm=TRUE) #max(all.intensi[which(all.times<(target.time+30) & all.times>(target.time-30))],na.rm=TRUE)
			
			#y.max <-max(all.intensi[which(all.times>(target.time-20) & all.times<(target.time+20))],na.rm=TRUE)
		# }
                  if (max(all.times[all.intensi > 0],na.rm=TRUE) > x.max)
                    x.max <- max(all.times[all.intensi > 0],na.rm=TRUE)
                  if (min(all.times[all.intensi > 0],na.rm=TRUE) < x.min)
                    x.min <- min(all.times[all.intensi > 0],na.rm=TRUE)
                  to.plot[[i]] <- cbind(all.times, all.intensi)
		  good_i<-c(good_i,i)
                  sel_l[[i]]<-sel
		}
            }
	    
	    max_int_vec<-{}
	    peak_score_vec<-{}
	    bool_plot=0
	  #  print("ymax is")
	   # print(y.max)
	 
	    y_max_final<-y.max+10000
	   # print( nrow(summary(to.plot)))
	    #print(to.plot)
	  #  print(length(to.plot))
	#    print(minrt)
	#    print(maxrt)
	    minrt<-min(c(0,minrt),na.rm=TRUE)
	    
	    peak_score<-(-1)
	    peak_score_vec<-c(-1)
	    if(length(to.plot)>0){
	    
	     peak_score_vec<-{}
            for (iii in 1:min(length(subset), nrow(summary(to.plot)))) 
		{
		peak_score<-0
                i <- subset[iii]
		curtime<-NA
		#if(i<length(sel_l))
		#if (this.intensi[i] != 0)
		if(is.na(this.times[i])==FALSE)
		{
		
	#	print("length of to.plot")
		#print(length(to.plot))
		#print(head(to.plot))
		#print(good_i)
		#print(i)
				#
		#if(length((to.plot[[i]]))>1)
		#if(i%in%good_i)
		test_file<-try(to.plot[[i]],silent=TRUE)
		
		   if(is(test_file,"try-error") || (is(test_file,"NULL")))
		{
			next;
		}else{
		sel<-sel_l[[i]]
		if(dim(to.plot[[i]])[2]>1){
		time_v1<-to.plot[[i]][,1]
			int_v1<-to.plot[[i]][,2]
			if(length(which(sel>dim(aligned$features[[i]])[1]))<1){
			target.time <- aligned$features[[i]][sel, 2]
			time.adjust <- aligned$features2[[i]][sel,2]
			curtime<-time.list[which(mz.list==this.mz)[1]]

			cur_cv<-round(cvvec[which(mz.list==this.mz)],1)
			max_int<-max(to.plot[[i]][,2],na.rm=TRUE)
			d1<-density(time_v1,na.rm=TRUE)
			max_int_index<-which(int_v1==max_int)
			peak_time<-time_v1[max_int_index]
			baseline_range<-20*peak_time
			#curtime<-time.adjust
			print("getting peak score")
			#print(summary(int_v1))
			#print(summary(time_v1))
			#print(time.adjust)
			#print(target.time)
			#print(curtime)
			peak_score<-get_peakscore(int_v1=int_v1,time_v1=time_v1,cur.time=curtime,bw=d1$bw,rt_range=range_rt)
			
			
			}else{
				peak_score<-NA
					curtime<-time.list[which(mz.list==this.mz)[1]]

			cur_cv<-round(cvvec[which(mz.list==this.mz)],1)
			max_int<-max(to.plot[[i]][,2],na.rm=TRUE)
			d1<-density(time_v1,na.rm=TRUE)
			max_int_index<-which(int_v1==max_int)
			peak_time<-time_v1[max_int_index]
			baseline_range<-20*peak_time
			}
			if(is.na(chem.names[1])==FALSE){
                        curname<-chem.names[which(mz.list==this.mz)]
			}else{
				curname<-""
			}
			if(length(peak_score_vec)<1){
				bool_plot=0
			}
			
			max_int_vec<-c(max_int_vec,max_int)	
			peak_score_vec<-c(peak_score_vec,peak_score)	
			minrt=max(0,(curtime-2*d1$bw))
			maxrt=min(600,(curtime+2*d1$bw))
			
			colors[iii][which(is.na(colors[iii])==TRUE)]<-"gray"
			#if(FALSE)
			y.maxt<-y.max #+10000
			if(plotEIC==TRUE)
			{
		
	    y.max<-y.maxt+10000
                if (bool_plot == 0) {
		bool_plot=1
		  plot(to.plot[[i]], xlim = c(minrt, maxrt),
                    ylim = c(0,  y_max_final), type = "l", col = colors[iii],
                    xlab = "retention time", ylab = "intensity",
                    main = paste("m/z",
                      round(this.mz, 5)," time:",round(curtime,1),"\n", " RT range: ",paste(round(range_rt,0),collapse=" to "),"\n", curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
			}else {
			if(bool_plot==1){
				lines(to.plot[[i]], col = colors[iii])
			}else{
				bool_plot=1
				  plot(to.plot[[i]], xlim = c(minrt, maxrt),
                    ylim = c(0,   y_max_final), type = "l", col = colors[iii],
                    xlab = "retention time", ylab = "intensity",
                    main = paste("EIC m/z",
                      #round(this.mz, 5)," time:", round(curtime,1),"\n",curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
		        round(this.mz, 5)," time:",round(curtime,1),"\n", " RT range: ",paste(round(range_rt,0),collapse=" to "),"\n", curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
				}

                	}
			
			}
			
			
		}
                }
		}
	    }
	    }
	  
	    peak_score_vec<-max(peak_score_vec,na.rm=TRUE) #dim(aligned$features[[i]])[1] #
	
	 max_int_all<-max(max_int_vec,na.rm=TRUE)
	    peak_score_all<-c(peak_score_all,peak_score_vec)
	    if(bool_plot==1){
	    	#mtext(paste("peak score:",round(peak_score_vec,2),";","max int (log10): ",log10(round(max_int_all)+1),sep=""),cex=0.7)
		mtext(paste("peak score:",round(peak_score_vec,2),sep=""),cex=0.7)
		}
	
	if(peak_score_vec==(-Inf)){
	
		peak_score_vec=0.5
	}
	
	if(peak_score_vec==(Inf)){
	
		peak_score_vec=10
	}

	#peak_score_all<-c(peak_score_all,peak_score_vec)
	#}
	return(round(peak_score_vec,2))
	})   
 }
 
return(peak_score_all)
}
apLCMS.get.peakscorevold<-function(aligned, rows = NA, colors = NA, transform = "none",
    subset = NA, mz.list=NA,time.list=NA,minrt=NA, maxrt=NA, min.run=12, min.pres=0.5, max.spline.time.points = 1000,chem.names=NA,rawprofileloc,cvvec=NA,plotEIC=TRUE,numnodes=NA)
{
    if (!(transform %in% c("none", "log", "sqrt", "cuberoot"))) {
        message("Invalid transformation. It has to be from: none, log, sqrt, and cuberoot")
        break
    }
	setwd(rawprofileloc)

	if(is.na(numnodes)==TRUE){
		numnodes<-detectCores()-1
	}
	 
    peak_score_all<-{}

	
	
	feat_names<-paste(aligned$final.ftrs[,1],aligned$final.ftrs[,2],sep="_")
    if (!is.na(rows[1])) 
{

        library(splines)
        num.exp <- nrow(summary(aligned$features))
	
	  files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
	 if(is.na(mz.list)==TRUE){
	 
		mz.list<-aligned$final.ftrs[rows,1]
	 }
	  if(is.na(time.list)==TRUE){
	 
		time.list<-aligned$final.ftrs[rows,2]
	 }
	
	feat_names<-paste(aligned$final.ftrs[rows,1],aligned$final.ftrs[rows,2],sep="_")
	
	if(length(files)<1){
		aligned$final.ftrs<-aligned$aligned.ftrs
		 files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
		aligned$aligned.ftrs<-NULL
	}
	if(is.na(minrt)==TRUE){
		minrt=min(aligned$final.ftrs[,2],na.rm=TRUE)
	}
	if(is.na(maxrt)==TRUE){
                maxrt=max(aligned$final.ftrs[,2],na.rm=TRUE)
        }

        if (is.na(subset[1]))
            subset <- 1:num.exp
        if (is.na(colors))
            colors <- rainbow(length(subset)) #
	
	
        rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1, 
            2 * num.exp, by = 2)], "_", min.run, "_", min.pres, 
            "_", aligned$mz.tol, ".rawprof", sep = "")
	    
        #adj.times <- new("list")

     
      #  rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1,
       #     2 * num.exp, by = 2)], "_", min.run, "_", min.pres,
       #     "_", aligned$mz.tol, ".rawprof", sep = "")
	
	rawprof.names<-tolower(rawprof.names)
        adj.times <- new("list")
        #for (i in subset) 
	adj.times<-lapply(1:length(subset),function(j)
	{
		i<-subset[j]
            load(rawprof.names[i])
            times <- unique(raw.prof$labels)
            times <- times[order(times)]
            orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][,
                3]), 2]
            adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][,
                3]), 2]
            orig.time <- round(orig.time, 8)
            adjusted.time <- round(adjusted.time, 8)
            ttt <- table(adjusted.time)
            ttt2 <- table(orig.time)
            to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt ==
                1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 ==
                1]))
            if (length(to.use) > max.spline.time.points)
                to.use <- sample(to.use, max.spline.time.points,
                  replace = FALSE)
            orig.time <- orig.time[to.use]
            adjusted.time <- adjusted.time[to.use]
            sp <- interpSpline(adjusted.time ~ orig.time)
            #adj.times[[i]] <- 
	    return(predict(sp, times)$y)
        })
	

	
#for (n in 1:length(rows)) {
#peak_score_all<-lapply(1:length(rows),function(n){	
#if(FALSE){

		cl<-parallel::makeCluster(numnodes)
	
		
			clusterEvalQ(cl, "get_peakscore")
			clusterExport(cl, "get_peakscore")
		clusterEvalQ(cl, "compare_intensities_ttest")
	#}	

print("number of rows")
print(length(rows))
peak_score_all<-parLapply(cl,1:length(rows),function(n){ 

		print(n)
            if (n%%4 == 1) {
                par(mfrow = c(2, 2))
            }
            this.row <- rows[n]
            this.intensi <- aligned$final.ftrs[this.row, 5:(num.exp +
                4)]
            this.times <- aligned$final.times[this.row, 5:(num.exp +
                4)]
		
		med_rt<-median(this.times[this.times>0],na.rm=TRUE)
		range_rt<-range(this.times[this.times>0],na.rm=TRUE)
		
            this.mz <- aligned$final.ftrs[this.row, 1]
           max_int_index<-order(this.intensi,decreasing=TRUE)
	   
	  if(length(max_int_index)>6)
	  {
		rand_set<-sample(x=max_int_index[-c(1:3)],size=3)
		subset<-c(max_int_index[1:3],rand_set)
		
	  }
	   
	   #print(this.row)
	    to.plot <- new("list")
            y.max <- 0
            x.max <- 0
            x.min <- Inf
	    good_i<-{}
	    sel_l<-new("list")
                    
	    
	    max_int_vec<-{}
	    peak_score_vec<-{}
	    bool_plot=0

	    y_max_final<-y.max+10000

	    minrt<-min(c(0,minrt),na.rm=TRUE)
			    if(length(to.plot)>0){
         #   for (iii in 1:min(length(subset), nrow(summary(to.plot))))
	   peak_score_vec<-lapply(1:length(rows),function(n)
		{
		peak_score<-0
                i <- subset[iii]
		curtime<-NA
		
		 i <- subset[iii]
		to.plot[[i]]<-{}
		sel_l[[i]]<-{}
		
                #if (this.intensi[i] != 0) 
		if(is.na(this.times[i])==FALSE)
		{
                  load(rawprof.names[i])
                  times <- unique(raw.prof$labels)
                  times <- times[order(times)]
                  mz.diff <- abs(aligned$features2[[i]][, 1] -
                    this.mz)
                  time.diff <- abs(aligned$features2[[i]][, 2] -
                    this.times[i])
                  sel <- which(time.diff < 1e-17)
                  if (length(sel) > 1)
                    sel <- sel[which(mz.diff[sel] == min(mz.diff[sel]))[1]]
                  if (length(sel) == 0) {
                    mz.lim <- aligned$final.ftrs[this.row, c(3,
                      4)]
                    sel <- which(aligned$features2[[i]][, 1] >=
                      mz.lim[1] & aligned$features2[[i]][, 1] <=
                      mz.lim[2] & !is.na(aligned$features2[[i]][,
                      6]))
                    sub.features <- aligned$features2[[i]][sel,
                      ]
		      print(length(sub.features))
		      if(length(sub.features)>7){
			sub.time.diff <- abs(sub.features[, 2] - this.times[i])
                    }else{
		    sub.time.diff <- abs(sub.features[2] - this.times[i])
		    }
		    sel <- sel[which(sub.time.diff == min(sub.time.diff))][1]
                  }
                  target.mz <- aligned$features2[[i]][sel, 1]
                  target.time <- aligned$features[[i]][sel, 2]
                  time.adjust <- aligned$features2[[i]][sel,
                    2] - aligned$features[[i]][sel, 2]
                  mz.diff <- abs(raw.prof$masses - target.mz)
                  sel.slice <- raw.prof$grps[which(mz.diff ==
                    min(mz.diff))[1]]
                  sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                    sel.slice])
			print("sel.time.range is")
		    print(sel.time.range)
		    if(is.na(sel.time.range[1])==TRUE){
			return(0)
		    }
                  while (target.time < sel.time.range[1] | target.time >
                    sel.time.range[2]) {
                    mz.diff[raw.prof$grps == sel.slice] <- 100
                    sel.slice <- raw.prof$grps[which(mz.diff ==
                      min(mz.diff))[1]]
                    sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                      sel.slice])
                  }
                  sel.time <- raw.prof$labels[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- raw.prof$intensi[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- sel.intensi[order(sel.time)]
                  sel.time <- sel.time[order(sel.time)]
                  all.intensi <- times * 0
                  all.intensi[times %in% sel.time] <- sel.intensi
                  if(i<length(adj.times)){
		  all.times <- adj.times[[i]]
			if(length(which(is.na(all.times)==TRUE))>0){
				sn=NA
			}
		  }else{
			all.times <-{}
			sn=NA
                                #return(sn)
		  }
		 
                  if (transform == "log")
                    all.intensi <- log10(all.intensi + 1)
                  if (transform == "sqrt")
                    all.intensi <- sqrt(all.intensi)
                  if (transform == "cuberoot")
                    all.intensi <- all.intensi^(1/3)
		   #if (max(all.intensi,na.rm=TRUE) > y.max)
		    #	y.max<-max(all.intensi)
                  if (max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE) > y.max)
			y.max <- max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE)
	
                  if (max(all.times[all.intensi > 0],na.rm=TRUE) > x.max)
                    x.max <- max(all.times[all.intensi > 0],na.rm=TRUE)
                  if (min(all.times[all.intensi > 0],na.rm=TRUE) < x.min)
                    x.min <- min(all.times[all.intensi > 0],na.rm=TRUE)
                  to.plot[[i]] <- cbind(all.times, all.intensi)
		  good_i<-c(good_i,i)
                  sel_l[[i]]<-sel
		}
		
		#if(i<length(sel_l))
		#if (this.intensi[i] != 0)
		{
		
		#print("length of to.plot")
		#print(length(to.plot))
		
				#
		#if(length((to.plot[[i]]))>1)
		#if(i%in%good_i)
		test_file<-try(to.plot[[i]],silent=TRUE)
		print(test_file)
		   if(is(test_file,"try-error") || (is(test_file,"NULL")))
		{
			next;
		}else{
			sel<-sel_l[[i]]
				if(dim(to.plot[[i]])[2]>1){
				time_v1<-to.plot[[i]][,1]
				int_v1<-to.plot[[i]][,2]
				if(length(which(sel>dim(aligned$features[[i]])[1]))<1){
					target.time <- aligned$features[[i]][sel, 2]
					time.adjust <- aligned$features2[[i]][sel,2]
					curtime<-time.list[which(mz.list==this.mz)[1]]

					cur_cv<-round(cvvec[which(mz.list==this.mz)],1)
					max_int<-max(to.plot[[i]][,2],na.rm=TRUE)
					d1<-density(time_v1,na.rm=TRUE)
					max_int_index<-which(int_v1==max_int)
					peak_time<-time_v1[max_int_index]
					baseline_range<-20*peak_time
					#curtime<-time.adjust
					print("getting peak score")
					print(summary(int_v1))
					print(summary(time_v1))
					print("mz")
					print(this.mz)
					
					
					#print(time.adjust)
					#print(target.time)
					#print(curtime)
					peak_score<-get_peakscore(int_v1=int_v1,time_v1=time_v1,cur.time=curtime,bw=d1$bw,rt_range=range_rt)
					
					
					}else{
						peak_score<-NA
					}
					if(is.na(chem.names[1])==FALSE){
					curname<-chem.names[which(mz.list==this.mz)]
					}else{
						curname<-""
					}
					if(length(peak_score_vec)<1){
						bool_plot=0
					}
					peak_score_vec<-c(peak_score_vec,peak_score)	
					max_int_vec<-c(max_int_vec,max_int)
					minrt=max(0,(curtime-2*d1$bw))
					maxrt=min(600,(curtime+2*d1$bw))
					
					colors[iii][which(is.na(colors[iii])==TRUE)]<-"gray"
					#if(FALSE)
					y.maxt<-y.max #+10000

					
			
			}
                }
		}
		return(peak_score)
	    })
	    }
	  
	  
	    peak_score_vec<-max(peak_score_vec,na.rm=TRUE) #dim(aligned$features[[i]])[1] #
	   max_int_all<-max(max_int_vec,na.rm=TRUE)
	    peak_score_all<-c(peak_score_all,peak_score_vec)
	    if(bool_plot==1){
	    	#mtext(paste("peak score:",round(peak_score_vec,2)," ","max int: ",max_int_all,sep=""),cex=0.7)
		mtext(paste("peak score:",round(peak_score_vec,2),";","max int (log10): ",log10(round(max_int_all)+1),sep=""),cex=0.7)
		}
	if(peak_score_vec==(-Inf)){
	
		peak_score_vec=1
	}
	
       if(is.na(peak_score_vec)==TRUE){
       
	peak_score_vec=0.1
       }
	return(peak_score_vec)
	})   
 }
    peak_score_all<-unlist(peak_score_all)
    
	names(peak_score_all)<-as.character(feat_names)
    return(peak_score_all)
}

apLCMS.get.peakscorev3<-function(aligned, rows = NA, colors = NA, transform = "none",
    subset = NA, mz.list=NA,time.list=NA,minrt=NA, maxrt=NA, min.run=12, min.pres=0.5, max.spline.time.points = 1000,chem.names=NA,rawprofileloc,cvvec=NA,plotEIC=TRUE,numnodes=NA)
{
    if (!(transform %in% c("none", "log", "sqrt", "cuberoot"))) {
        message("Invalid transformation. It has to be from: none, log, sqrt, and cuberoot")
        break
    }
	setwd(rawprofileloc)

	if(is.na(numnodes)==TRUE){
		numnodes<-detectCores()-1
	}
	 
    peak_score_all<-{}


	
    if (!is.na(rows[1])) {

        library(splines)
        num.exp <- nrow(summary(aligned$features))
	
	  files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
	 if(is.na(mz.list)==TRUE){
	 
		mz.list<-aligned$final.ftrs[rows,1]
	 }
	  if(is.na(time.list)==TRUE){
	 
		time.list<-aligned$final.ftrs[rows,2]
	 }
	
	if(length(files)<1){
		aligned$final.ftrs<-aligned$aligned.ftrs
		 files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
		aligned$aligned.ftrs<-NULL
	}
	if(is.na(minrt)==TRUE){
		minrt=min(aligned$final.ftrs[,2],na.rm=TRUE)
	}
	if(is.na(maxrt)==TRUE){
                maxrt=max(aligned$final.ftrs[,2],na.rm=TRUE)
        }

        if (is.na(subset[1]))
            subset <- 1:num.exp
        if (is.na(colors))
            colors <- rainbow(length(subset)) #
	
	
        rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1, 
            2 * num.exp, by = 2)], "_", min.run, "_", min.pres, 
            "_", aligned$mz.tol, ".rawprof", sep = "")
	    
        #adj.times <- new("list")

     
      #  rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1,
       #     2 * num.exp, by = 2)], "_", min.run, "_", min.pres,
       #     "_", aligned$mz.tol, ".rawprof", sep = "")
	
	rawprof.names<-tolower(rawprof.names)
        adj.times <- new("list")
        #for (i in subset) 
	adj.times<-lapply(1:length(subset),function(j)
	{
		i<-subset[j]
            load(rawprof.names[i])
            times <- unique(raw.prof$labels)
            times <- times[order(times)]
            orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][,
                3]), 2]
            adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][,
                3]), 2]
            orig.time <- round(orig.time, 8)
            adjusted.time <- round(adjusted.time, 8)
            ttt <- table(adjusted.time)
            ttt2 <- table(orig.time)
            to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt ==
                1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 ==
                1]))
            if (length(to.use) > max.spline.time.points)
                to.use <- sample(to.use, max.spline.time.points,
                  replace = FALSE)
            orig.time <- orig.time[to.use]
            adjusted.time <- adjusted.time[to.use]
            sp <- interpSpline(adjusted.time ~ orig.time)
            #adj.times[[i]] <- 
	    return(predict(sp, times)$y)
        })
	

	
#for (n in 1:length(rows)) {
#peak_score_all<-lapply(1:length(rows),function(n){	

	#if(FALSE){
		cl<-parallel::makeCluster(numnodes)
	
		
			clusterEvalQ(cl, "get_peakscore")
			clusterExport(cl, "get_peakscore")
		clusterEvalQ(cl, "compare_intensities_ttest")
	#}			
peak_score_all<-parLapply(cl,1:length(rows),function(n){ 

		print(n)
            if (n%%4 == 1) {
                par(mfrow = c(2, 2))
            }
            this.row <- rows[n]
            this.intensi <- aligned$final.ftrs[this.row, 5:(num.exp +
                4)]
            this.times <- aligned$final.times[this.row, 5:(num.exp +
                4)]
		
		med_rt<-median(this.times[this.times>0],na.rm=TRUE)
		range_rt<-range(this.times[this.times>0],na.rm=TRUE)
		
            this.mz <- aligned$final.ftrs[this.row, 1]
           max_int_index<-order(this.intensi,decreasing=TRUE)
	   
	  if(length(max_int_index)>20)
	  {
		rand_set<-sample(x=max_int_index[-c(1:10)],size=10)
		subset<-c(max_int_index[1:10],rand_set)
		
	  }
	   
	   #print(this.row)
	    to.plot <- new("list")
            y.max <- 0
            x.max <- 0
            x.min <- Inf
	    good_i<-{}
	    sel_l<-new("list")
                      for (iii in 1:length(subset)) {
                i <- subset[iii]
		to.plot[[i]]<-{}
		sel_l[[i]]<-{}
		
                #if (this.intensi[i] != 0) 
		if(is.na(this.times[i])==FALSE)
		{
                  load(rawprof.names[i])
                  times <- unique(raw.prof$labels)
                  times <- times[order(times)]
                  mz.diff <- abs(aligned$features2[[i]][, 1] -
                    this.mz)
                  time.diff <- abs(aligned$features2[[i]][, 2] -
                    this.times[i])
                  sel <- which(time.diff < 1e-17)
                  if (length(sel) > 1)
                    sel <- sel[which(mz.diff[sel] == min(mz.diff[sel]))[1]]
                  if (length(sel) == 0) {
                    mz.lim <- aligned$final.ftrs[this.row, c(3,
                      4)]
                    sel <- which(aligned$features2[[i]][, 1] >=
                      mz.lim[1] & aligned$features2[[i]][, 1] <=
                      mz.lim[2] & !is.na(aligned$features2[[i]][,
                      6]))
                    sub.features <- aligned$features2[[i]][sel,
                      ]
		      print(length(sub.features))
		      if(length(sub.features)>7){
			sub.time.diff <- abs(sub.features[, 2] - this.times[i])
                    }else{
		    sub.time.diff <- abs(sub.features[2] - this.times[i])
		    }
		    sel <- sel[which(sub.time.diff == min(sub.time.diff))][1]
                  }
                  target.mz <- aligned$features2[[i]][sel, 1]
                  target.time <- aligned$features[[i]][sel, 2]
                  time.adjust <- aligned$features2[[i]][sel,
                    2] - aligned$features[[i]][sel, 2]
                  mz.diff <- abs(raw.prof$masses - target.mz)
                  sel.slice <- raw.prof$grps[which(mz.diff ==
                    min(mz.diff))[1]]
                  sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                    sel.slice])
			print("sel.time.range is")
		    print(sel.time.range)
		    if(is.na(sel.time.range[1])==TRUE){
			next;
		    }
                  while (target.time < sel.time.range[1] | target.time >
                    sel.time.range[2]) {
                    mz.diff[raw.prof$grps == sel.slice] <- 100
                    sel.slice <- raw.prof$grps[which(mz.diff ==
                      min(mz.diff))[1]]
                    sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                      sel.slice])
                  }
                  sel.time <- raw.prof$labels[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- raw.prof$intensi[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- sel.intensi[order(sel.time)]
                  sel.time <- sel.time[order(sel.time)]
                  all.intensi <- times * 0
                  all.intensi[times %in% sel.time] <- sel.intensi
                  if(i<length(adj.times)){
		  all.times <- adj.times[[i]]
			if(length(which(is.na(all.times)==TRUE))>0){
				sn=NA
			}
		  }else{
			all.times <-{}
			sn=NA
                                #return(sn)
		  }
		 
                  if (transform == "log")
                    all.intensi <- log10(all.intensi + 1)
                  if (transform == "sqrt")
                    all.intensi <- sqrt(all.intensi)
                  if (transform == "cuberoot")
                    all.intensi <- all.intensi^(1/3)
		   #if (max(all.intensi,na.rm=TRUE) > y.max)
		    #	y.max<-max(all.intensi)
                  if (max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE) > y.max)
			y.max <- max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE)
		#if (max(all.intensi[which(all.times>(target.time-20) & all.times<(target.time+20))],na.rm=TRUE) > y.max){
			#y.max <- max(all.intensi,na.rm=TRUE) #max(all.intensi[which(all.times<(target.time+30) & all.times>(target.time-30))],na.rm=TRUE)
			
			#y.max <-max(all.intensi[which(all.times>(target.time-20) & all.times<(target.time+20))],na.rm=TRUE)
		# }
                  if (max(all.times[all.intensi > 0],na.rm=TRUE) > x.max)
                    x.max <- max(all.times[all.intensi > 0],na.rm=TRUE)
                  if (min(all.times[all.intensi > 0],na.rm=TRUE) < x.min)
                    x.min <- min(all.times[all.intensi > 0],na.rm=TRUE)
                  to.plot[[i]] <- cbind(all.times, all.intensi)
		  good_i<-c(good_i,i)
                  sel_l[[i]]<-sel
		}
            }
	    
	    max_int_vec<-{}
	    peak_score_vec<-{}
	    bool_plot=0

	    y_max_final<-y.max+10000

	    minrt<-min(c(0,minrt),na.rm=TRUE)
			    if(length(to.plot)>0){
            for (iii in 1:min(length(subset), nrow(summary(to.plot)))) 
		{
		peak_score<-0
                i <- subset[iii]
		curtime<-NA
		#if(i<length(sel_l))
		#if (this.intensi[i] != 0)
		{
		
		print("length of to.plot")
		print(length(to.plot))
		
				#
		#if(length((to.plot[[i]]))>1)
		#if(i%in%good_i)
		test_file<-try(to.plot[[i]],silent=TRUE)
		print(test_file)
		   if(is(test_file,"try-error") || (is(test_file,"NULL")))
		{
			next;
		}else{
		sel<-sel_l[[i]]
		if(dim(to.plot[[i]])[2]>1){
		time_v1<-to.plot[[i]][,1]
			int_v1<-to.plot[[i]][,2]
			if(length(which(sel>dim(aligned$features[[i]])[1]))<1){
			target.time <- aligned$features[[i]][sel, 2]
			time.adjust <- aligned$features2[[i]][sel,2]
			curtime<-time.list[which(mz.list==this.mz)[1]]

			cur_cv<-round(cvvec[which(mz.list==this.mz)],1)
			max_int<-max(to.plot[[i]][,2],na.rm=TRUE)
			d1<-density(time_v1,na.rm=TRUE)
			max_int_index<-which(int_v1==max_int)
			peak_time<-time_v1[max_int_index]
			baseline_range<-20*peak_time
			#curtime<-time.adjust
			print("getting peak score")
			print(summary(int_v1))
			print(summary(time_v1))
			#print(time.adjust)
			#print(target.time)
			#print(curtime)
			peak_score<-get_peakscore(int_v1=int_v1,time_v1=time_v1,cur.time=curtime,bw=d1$bw,rt_range=range_rt)
			
			
			}else{
				peak_score<-NA
			}
			if(is.na(chem.names[1])==FALSE){
                        curname<-chem.names[which(mz.list==this.mz)]
			}else{
				curname<-""
			}
			if(length(peak_score_vec)<1){
				bool_plot=0
			}
			peak_score_vec<-c(peak_score_vec,peak_score)	
			max_int_vec<-c(max_int_vec,max_int)
			minrt=max(0,(curtime-2*d1$bw))
			maxrt=min(600,(curtime+2*d1$bw))
			
			colors[iii][which(is.na(colors[iii])==TRUE)]<-"gray"
			#if(FALSE)
			y.maxt<-y.max #+10000
			if(plotEIC==TRUE)
			{
		
	    y.max<-y.maxt+10000
                if (bool_plot == 0) {
		bool_plot=1
		  plot(to.plot[[i]], xlim = c(minrt, maxrt),
                    ylim = c(0,  y_max_final), type = "l", col = colors[iii],
                    xlab = "retention time", ylab = "intensity",
                    main = paste("m/z",
                      round(this.mz, 5)," time:",round(curtime,1),"\n", " RT range: ",paste(round(range_rt,0),collapse=" to "),"\n", curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
			}else {
			if(bool_plot==1){
				lines(to.plot[[i]], col = colors[iii])
			}else{
				bool_plot=1
				  plot(to.plot[[i]], xlim = c(minrt, maxrt),
                    ylim = c(0,   y_max_final), type = "l", col = colors[iii],
                    xlab = "retention time", ylab = "intensity",
                    main = paste("EIC m/z",
                      #round(this.mz, 5)," time:", round(curtime,1),"\n",curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
		        round(this.mz, 5)," time:",round(curtime,1),"\n", " RT range: ",paste(round(range_rt,0),collapse=" to "),"\n", curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
				}
                	}
			
			}
			
			
		}
                }
		}
	    }
	    }
	  
	  
	    peak_score_vec<-max(peak_score_vec,na.rm=TRUE) #dim(aligned$features[[i]])[1] #
	   max_int_all<-max(max_int_vec,na.rm=TRUE)
	    peak_score_all<-c(peak_score_all,peak_score_vec)
	    if(bool_plot==1){
	    	#mtext(paste("peak score:",round(peak_score_vec,2)," ","max int: ",max_int_all,sep=""),cex=0.7)
		mtext(paste("peak score:",round(peak_score_vec,2),";","max int (log10): ",log10(round(max_int_all)+1),sep=""),cex=0.7)
		}
	if(peak_score_vec==(-Inf)){
	
		peak_score_vec=1
	}
	
       if(is.na(peak_score_vec)==TRUE){
       
	peak_score_vec=0.1
       }
	return(peak_score_vec)
	})   
 }
    peak_score_all<-unlist(peak_score_all)
    return(peak_score_all)
}


apLCMS.get.peakscore2<-function(aligned, rows = NA, colors = NA, transform = "none",
    subset = NA, mz.list=NA,time.list=NA,minrt=NA, maxrt=NA, min.run=12, min.pres=0.5, max.spline.time.points = 1000,chem.names=NA,rawprofileloc,cvvec=NA,plotEIC=TRUE,numnodes=NA)
{
if (!(transform %in% c("none", "log", "sqrt", "cuberoot"))) {
        message("Invalid transformation. It has to be from: none, log, sqrt, and cuberoot")
        break
    }
	setwd(rawprofileloc)

	if(is.na(numnodes)==TRUE){
		numnodes<-detectCores()-1
	}
	 
    peak_score_all<-{}
    #print(rows)

	print("length rows")
	print(length(rows))
    if (!is.na(rows[1])) {

        library(splines)
        num.exp <- nrow(summary(aligned$features))
	
	  files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
	 if(is.na(mz.list)==TRUE){
	 
		mz.list<-aligned$final.ftrs[rows,1]
	 }
	  if(is.na(time.list)==TRUE){
	 
		time.list<-aligned$final.ftrs[rows,2]
	 }
	
	if(length(files)<1){
		aligned$final.ftrs<-aligned$aligned.ftrs
		 files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
		aligned$aligned.ftrs<-NULL
	}
	if(is.na(minrt)==TRUE){
		minrt=min(aligned$final.ftrs[,2],na.rm=TRUE)
	}
	if(is.na(maxrt)==TRUE){
                maxrt=max(aligned$final.ftrs[,2],na.rm=TRUE)
        }

        if (is.na(subset[1]))
            subset <- 1:num.exp
        if (is.na(colors))
            colors <- rainbow(length(subset)) #
	   
     
        rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1,
            2 * num.exp, by = 2)], "_", min.run, "_", min.pres,
            "_", aligned$mz.tol, ".rawprof", sep = "")
	rawprof.names<-tolower(rawprof.names)
        adj.times <- new("list")
        #for (i in subset) 
	adj.times<-lapply(1:length(subset),function(j)
	{
		i<-subset[j]
            load(rawprof.names[i])
            times <- unique(raw.prof$labels)
            times <- times[order(times)]
            orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][,
                3]), 2]
            adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][,
                3]), 2]
            orig.time <- round(orig.time, 8)
            adjusted.time <- round(adjusted.time, 8)
            ttt <- table(adjusted.time)
            ttt2 <- table(orig.time)
            to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt ==
                1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 ==
                1]))
            if (length(to.use) > max.spline.time.points)
                to.use <- sample(to.use, max.spline.time.points,
                  replace = FALSE)
            orig.time <- orig.time[to.use]
            adjusted.time <- adjusted.time[to.use]
            sp <- interpSpline(adjusted.time ~ orig.time)
            #adj.times[[i]] <- 
	    return(predict(sp, times)$y)
        })
	
	print("length rows")
	print(length(rows))
	
#for (n in 1:length(rows)) {
#peak_score_all<-lapply(1:length(rows),function(n){	

	#if(FALSE){
		cl<-parallel::makeCluster(numnodes)
	
		
			clusterEvalQ(cl, "get_peakscore")
			clusterExport(cl, "get_peakscore")
		clusterEvalQ(cl, "compare_intensities_ttest")
	#}			
peak_score_all<-parLapply(cl,1:length(rows),function(n){ 
            if (n%%4 == 1) {
                par(mfrow = c(2, 2))
            }
            this.row <- rows[n]
            this.intensi <- aligned$final.ftrs[this.row, 5:(num.exp +
                4)]
            this.times <- aligned$final.times[this.row, 5:(num.exp +
                4)]
		#print("this.times")
		#print(this.times)
		med_rt<-median(this.times[this.times>0],na.rm=TRUE)
		range_rt<-range(this.times[this.times>0],na.rm=TRUE)
		print(med_rt)
            this.mz <- aligned$final.ftrs[this.row, 1]
           max_int_index<-order(this.intensi,decreasing=TRUE)
	   
	  if(length(max_int_index)>20)
	  {
		rand_set<-sample(x=max_int_index[-c(1:10)],size=10)
		subset<-c(max_int_index[1:10],rand_set)
		
	  }
	   
	   #print(this.row)
	    to.plot <- new("list")
            y.max <- 0
            x.max <- 0
            x.min <- Inf
	    good_i<-{}
	    sel_l<-new("list")
                      for (iii in 1:length(subset)) {
                i <- subset[iii]
		to.plot[[i]]<-{}
		sel_l[[i]]<-{}
                if (this.intensi[i] != 0) {
                  load(rawprof.names[i])
                  times <- unique(raw.prof$labels)
                  times <- times[order(times)]
                  mz.diff <- abs(aligned$features2[[i]][, 1] -
                    this.mz)
                  time.diff <- abs(aligned$features2[[i]][, 2] -
                    this.times[i])
                  sel <- which(time.diff < 1e-17)
                  if (length(sel) > 1)
                    sel <- sel[which(mz.diff[sel] == min(mz.diff[sel]))[1]]
                  if (length(sel) == 0) {
                    mz.lim <- aligned$final.ftrs[this.row, c(3,
                      4)]
                    sel <- which(aligned$features2[[i]][, 1] >=
                      mz.lim[1] & aligned$features2[[i]][, 1] <=
                      mz.lim[2] & !is.na(aligned$features2[[i]][,
                      6]))
                    sub.features <- aligned$features2[[i]][sel,
                      ]
		      print(length(sub.features))
		      if(length(sub.features)>7){
			sub.time.diff <- abs(sub.features[, 2] - this.times[i])
                    }else{
		    sub.time.diff <- abs(sub.features[2] - this.times[i])
		    }
		    sel <- sel[which(sub.time.diff == min(sub.time.diff))][1]
                  }
                  target.mz <- aligned$features2[[i]][sel, 1]
                  target.time <- aligned$features[[i]][sel, 2]
                  time.adjust <- aligned$features2[[i]][sel,
                    2] - aligned$features[[i]][sel, 2]
                  mz.diff <- abs(raw.prof$masses - target.mz)
                  sel.slice <- raw.prof$grps[which(mz.diff ==
                    min(mz.diff))[1]]
                  sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                    sel.slice])
			print("sel.time.range is")
		    print(sel.time.range)
		    if(is.na(sel.time.range[1])==TRUE){
			next;
		    }
                  while (target.time < sel.time.range[1] | target.time >
                    sel.time.range[2]) {
                    mz.diff[raw.prof$grps == sel.slice] <- 100
                    sel.slice <- raw.prof$grps[which(mz.diff ==
                      min(mz.diff))[1]]
                    sel.time.range <- range(raw.prof$labels[raw.prof$grps ==
                      sel.slice])
                  }
                  sel.time <- raw.prof$labels[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- raw.prof$intensi[raw.prof$grps ==
                    sel.slice]
                  sel.intensi <- sel.intensi[order(sel.time)]
                  sel.time <- sel.time[order(sel.time)]
                  all.intensi <- times * 0
                  all.intensi[times %in% sel.time] <- sel.intensi

                  if(i<length(adj.times)){
		  all.times <- adj.times[[i]]
			if(length(which(is.na(all.times)==TRUE))>0){
				sn=NA
			}
		  }else{
			all.times <-{}
			sn=NA
                                #return(sn)
		  }
		 
                  if (transform == "log")
                    all.intensi <- log10(all.intensi + 1)
                  if (transform == "sqrt")
                    all.intensi <- sqrt(all.intensi)
                  if (transform == "cuberoot")
                    all.intensi <- all.intensi^(1/3)
                  if (max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE) > y.max)
			y.max <- max(all.intensi[which(all.times<(target.time+60) & all.times>(target.time-60))],na.rm=TRUE)

                  if (max(all.times[all.intensi > 0],na.rm=TRUE) > x.max)
                    x.max <- max(all.times[all.intensi > 0],na.rm=TRUE)
                  if (min(all.times[all.intensi > 0],na.rm=TRUE) < x.min)
                    x.min <- min(all.times[all.intensi > 0],na.rm=TRUE)
                  to.plot[[i]] <- cbind(all.times, all.intensi)
		  good_i<-c(good_i,i)
                  sel_l[[i]]<-sel
		}
            }
	    
	    
	    peak_score_vec<-{}
	    bool_plot=0
	    print("ymax is")
	    print(y.max)

	    y_max_final<-y.max+10000

	    minrt<-min(c(0,minrt),na.rm=TRUE)
	    if(length(to.plot)>0){
            for (iii in 1:min(length(subset), nrow(summary(to.plot)))) 
		{
		peak_score<-0
                i <- subset[iii]
		#if(i<length(sel_l))
		if (this.intensi[i] != 0)
		{
		
		print("length of to.plot")
		print(length(to.plot[[i]]))
		
		#if(i%in%good_i)
		if(length((to.plot[[i]]))>1)
		{
		
		
		sel<-sel_l[[i]]
		if(dim(to.plot[[i]])[2]>1){
		time_v1<-to.plot[[i]][,1]
			int_v1<-to.plot[[i]][,2]
			if(length(which(sel>dim(aligned$features[[i]])[1]))<1){
			target.time <- aligned$features[[i]][sel, 2]
			time.adjust <- aligned$features2[[i]][sel,2]
			curtime<-time.list[which(mz.list==this.mz)[1]]

			cur_cv<-round(cvvec[which(mz.list==this.mz)],1)
			max_int<-max(to.plot[[i]][,2],na.rm=TRUE)
			d1<-density(time_v1,na.rm=TRUE)
			max_int_index<-which(int_v1==max_int)
			peak_time<-time_v1[max_int_index]
			baseline_range<-20*peak_time
			#curtime<-time.adjust
			print("getting peak score")
			print(summary(int_v1))
			print(summary(time_v1))
			#print(time.adjust)
			#print(target.time)
			#print(curtime)
			peak_score<-get_peakscore(int_v1=int_v1,time_v1=time_v1,cur.time=curtime,bw=d1$bw,rt_range=range_rt)
			
			
			}else{
				peak_score<-NA
			}
			if(is.na(chem.names[1])==FALSE){
                        curname<-chem.names[which(mz.list==this.mz)]
			}else{
				curname<-""
			}
			if(length(peak_score_vec)<1){
				bool_plot=0
			}
			peak_score_vec<-c(peak_score_vec,peak_score)	
			minrt=max(0,(curtime-2*d1$bw))
			maxrt=min(600,(curtime+2*d1$bw))
			
			colors[iii][which(is.na(colors[iii])==TRUE)]<-"gray"
			#if(FALSE)
			y.maxt<-y.max #+10000
			if(plotEIC==TRUE)
			{
			 print("ymax is")
	    print(y.max)
	    y.max<-y.maxt+10000
                if (bool_plot == 0) {
		bool_plot=1
		  plot(to.plot[[i]], xlim = c(minrt, maxrt),
                    ylim = c(0,  y_max_final), type = "l", col = colors[iii],
                    xlab = "retention time", ylab = "intensity",
                    main = paste("m/z",
                      round(this.mz, 5)," time:",round(curtime,1),"\n", " RT range: ",paste(round(range_rt,0),collapse=" to "),"\n", curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
			}else {
			if(bool_plot==1){
				lines(to.plot[[i]], col = colors[iii])
			}else{
				bool_plot=1
				  plot(to.plot[[i]], xlim = c(minrt, maxrt),
                    ylim = c(0,   y_max_final), type = "l", col = colors[iii],
                    xlab = "retention time", ylab = "intensity",
                    main = paste("EIC m/z",
                      #round(this.mz, 5)," time:", round(curtime,1),"\n",curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
		        round(this.mz, 5)," time:",round(curtime,1),"\n", " RT range: ",paste(round(range_rt,0),collapse=" to "),"\n", curname, "\n","medCV:",cur_cv,"%"),cex.main=0.7)
				}
                	}
			
			}
			
			
		}
                }
		}
	    }
	    }
	   # print(to.plot[[i]])
	    print(peak_score_vec)
	  # print(bool_plot)
	    peak_score_vec<-max(peak_score_vec,na.rm=TRUE) #dim(aligned$features[[i]])[1] #
	    #print("next")
	    peak_score_all<-c(peak_score_all,peak_score_vec)
	    if(bool_plot==1){
	    	mtext(paste("peak score:",round(peak_score_vec,2),sep=""),cex=0.7)
		}
	
	
        message("the colors used are:")
        print(paste(subset, ": ", colors[1:length(subset)], sep = ""))
	 print(peak_score_vec)
	return(peak_score_vec)
	})   
 }
    peak_score_all<-unlist(peak_score_all)
    return(peak_score_all)
}



get_peakscore<-function(int_v1,time_v1,cur.time,bw,rt_range){
	
			max_int_overall<-max(int_v1,na.rm=TRUE)
                        max_int_index<-which(int_v1==max_int_overall)
                        peak_time<-cur.time #time_v1[which(time_v1>(cur.time-1) & time_v1<(cur.time+1))] #[max_int_index]
			
	
			sn<-(0)
			rt_window<-abs(rt_range[2]-rt_range[1])
			
			
			if(length(time_v1)<1 || length(cur.time)<1){
				
				
				return(sn)
			}
			if(max(time_v1,na.rm=TRUE)<(cur.time+rt_window)){
				maxtimelimit<-max(time_v1,na.rm=TRUE)
			}else{
				maxtimelimit<-(cur.time+rt_window)
			}
			 if(min(time_v1,na.rm=TRUE)>(cur.time-rt_window)){
                                mintimelimit<-min(time_v1,na.rm=TRUE)
                        }else{

				mintimelimit<-(cur.time-rt_window)
			}
                        
			max_int<-max(int_v1[which(time_v1>(cur.time-rt_window) & time_v1<(cur.time+rt_window))],na.rm=TRUE)
			print("max int")
			#print(int_v1)
			print(max_int)
			
			if(is.na(max_int)==TRUE){
			
				return(0)
			}else{
			half_int<-0.5*max_int
			
			low_end_half_int<-time_v1[which(int_v1<=half_int & time_v1<(cur.time-rt_window))]
			low_end_half_int<-low_end_half_int[length(low_end_half_int)]
			
			
			
			high_end_half_int<-time_v1[which(int_v1<=half_int & time_v1>(cur.time+rt_window))]
			high_end_half_int<-high_end_half_int[1]
			
			baseline_range<-max(rt_window*4,bw) #min(rt_window,(high_end_half_int-low_end_half_int)*0.5) #min(10,bw) #*peak_time
                        peak_time_ind<-which(time_v1>(peak_time-baseline_range) & time_v1<(peak_time+baseline_range))
			
			peak_time_ind_nonzero<-which(time_v1>(peak_time-baseline_range) & time_v1<(peak_time+baseline_range) & int_v1>0)
                        base_time_ind<-seq(1,length(time_v1))
			
                        base_time_ind<-which(time_v1<(peak_time-baseline_range) | time_v1>(peak_time+baseline_range))

			#get 80th percentile threshold
			time_lowpeaks<-which(int_v1<=quantile(int_v1[which(int_v1>0)],0.8,na.rm=TRUE))
			
			base_time_ind<-intersect(base_time_ind,time_lowpeaks)
			

                        if(length(base_time_ind)<1){
                                low_end<-0
                                high_end<-100
                                base_o<-1000
                                base_o<-log10(base_o)
                                max_int1<-max_int
                                max_int<-log10(max_int)
				int_v1<-log10(int_v1+1)

                            
				               sn<-(max_int-3)/(sd(int_v1,na.rm=TRUE)+0.01)
                        
					max_int2<-max(int_v1[peak_time_ind],na.rm=TRUE)

					int_factor<-1*(max_int2-3)

					if(int_factor<0){
			
						int_factor<-0.001
					}

			
			max_peak_int<-max(int_v1[peak_time_ind_nonzero],na.rm=TRUE)
			min_peak_int<-min(1,int_v1[peak_time_ind_nonzero],na.rm=TRUE)
			int_check<-1*(max_peak_int-3)+0.01
			
		
				if(int_factor>2){

					time_diff<-min(100,(high_end_half_int-low_end_half_int),na.rm=TRUE)
					
					#time_diff<-min(100,(high_end-low_end),na.rm=TRUE)
					time_diff<-0.25*sqrt(time_diff)
				}else{
					time_diff<-min(100,(high_end_half_int-low_end_half_int),na.rm=TRUE)
					
					#time_diff<-min(100,(high_end-low_end),na.rm=TRUE)
					time_diff<-0.25*sqrt(time_diff)
				}

		
				
			#sn<-1*(sn*(int_factor)*(int_check))*(1/(time_diff))*(1/(0.25*sqrt(rt_range[2]-rt_range[1]+1)))
			sn<-1*(sn*(int_factor)*(int_check))*(1/(time_diff))*(1/(sqrt(rt_range[2]-rt_range[1]+1)))	
			
			print("eval")
			print(sn)
			print(max_int)
			print(base_o)
			print(sd(int_v1,na.rm=TRUE))
			print(int_factor)
			print(int_check)
			print(high_end_half_int)
			print(low_end_half_int)
			print(high_end)
			print(low_end)
			print(rt_range)
			print(baseline_range)
			
                        }else{

			low_end<-time_v1[base_time_ind[1]]
                        high_end<-time_v1[base_time_ind[length(base_time_ind)]]
               
			base_int_vec<-int_v1[base_time_ind]
			
			base_int_vec_zero_count<-length(which(base_int_vec==0))
			
			if(base_int_vec_zero_count>0){
			
				base_o<-min(base_int_vec[-c(which(base_int_vec==0))],na.rm=TRUE)
				
			}else{
				base_o<-max(base_int_vec,na.rm=TRUE,3)
					
			}
			base_o<-mean(int_v1[base_time_ind],na.rm=TRUE)
			
			#base_o<-min(c(int_v1[base_time_ind[1]],int_v1[base_time_ind[length(base_time_ind)]]),na.rm=TRUE)

			  #base_o<-median(int_v1[base_time_ind],na.rm=TRUE)+0.01
			
			max_int1<-max_int
                        max_int<-log10(max_int+1)
			
			
			if(is.na(max_int)==FALSE){


					#base_o<-max(3,base_o)
			
                        base_o<-log10(base_o+1.01)
                        #int_v1[base_time_ind]<-log10(int_v1[base_time_ind])
                        int_v1<-log10(int_v1+1)
			base_int_vec<-log10(base_int_vec+1.01)
			

                        sn<-(max_int-base_o)/(sd(int_v1,na.rm=TRUE)+0.01)
                        
                        max_int2<-max(int_v1[peak_time_ind],na.rm=TRUE)

			int_factor<-1*(max_int2-3)

			if(int_factor<0){
			
						int_factor<-0.001
					}

			
			max_peak_int<-max(int_v1[peak_time_ind_nonzero],na.rm=TRUE)
			min_peak_int<-min(1,int_v1[peak_time_ind_nonzero],na.rm=TRUE)
			int_check<-1*(max_peak_int-3)+0.01
			
		
				if(int_factor>2){

					time_diff<-min(100,(high_end_half_int-low_end_half_int),na.rm=TRUE)
					
					#time_diff<-min(100,(high_end-low_end),na.rm=TRUE)
					time_diff<-0.25*sqrt(time_diff)
				}else{
					time_diff<-min(100,(high_end_half_int-low_end_half_int),na.rm=TRUE)
					
					#time_diff<-min(100,(high_end-low_end),na.rm=TRUE)
					time_diff<-0.25*sqrt(time_diff)
				}

		
				
			#sn<-1*(sn*(int_factor)*(int_check))*(1/(time_diff))*(1/(0.25*sqrt(rt_range[2]-rt_range[1]+1)))
				
			sn<-1*(sn*(int_factor)*(int_check))*(1/(time_diff))*(1/(sqrt(rt_range[2]-rt_range[1]+1)))
			
			
			print("eval")
			print(sn)
			print(max_int)
			print(base_o)
			print(sd(int_v1,na.rm=TRUE))
			print(sd(base_int_vec,na.rm=TRUE))
			print(int_factor)
			print(int_check)
			print(high_end_half_int)
			print(low_end_half_int)
			print(high_end)
			print(low_end)
			print(rt_range)
			print(baseline_range)
			
			#print(base_time_ind)
			
                        #sn<-((max_int-base_o))*sn

				#sn<-(max_int-base_o)/sd(int_v1,na.rm=TRUE)

                       
				                        

			}else{
				  #sn<-NA
				
				  sn<-0
				
			}
			
		}
		
		return(sn)
		}
	}
	
	
	fwhm <- function(int_v1,time_v1,i) {

	i<-which(int_v1==max(int_v1,na.rm=TRUE))
	
    n <- length(int_v1)
    left <- ifelse(i <= 1, 1, i)
    right <- ifelse(i >= n, n, i)
    
    hm <- max(int_v1)/2
    
    while (left > 1 && int_v1[left] > hm) {
      left <- left-1
    }
    while (right < n && int_v1[right] > hm) {
      right <- right+1
    }
    
    ## interpolate x values
    xleft <- approx(x=time_v1[left:(left+1)],
                    y=int_v1[left:(left+1)], xout=hm)$y
    xright <- approx(x=time_v1[(right-1):right],
                    y=int_v1[(right-1):right], xout=hm)$y
                    
    return(abs(xleft-xright))
  }
  
 
get_peakscorev2<-function(int_v1,time_v1,cur.time,bw){
	
			max_int_overall<-max(int_v1,na.rm=TRUE)
                        max_int_index<-which(int_v1==max_int_overall)
			peak_int<-int_v1[which(time_v1>(cur.time-3) & time_v1<(cur.time+3))]
			
			time_diff_vec<-abs(time_v1-cur.time)
			
			print("new mz")
			print(cur.time)
			
			#print(time_diff_vec)
			
			peak_time<-time_v1[which(time_diff_vec==min(time_diff_vec,na.rm=TRUE))]
			
			peak_int<-int_v1[which(time_diff_vec==min(time_diff_vec,na.rm=TRUE))]
			
			peak_ind<-which(time_diff_vec==min(time_diff_vec,na.rm=TRUE))
			
			
			#print(which(time_diff_vec==min(time_diff_vec,na.rm=TRUE)))
			
			#print(int_v1[(peak_ind-5):(peak_ind+5)])
			#print(cur.time)
			#print(peak_time)
			#print(peak_int)
			
                       # peak_time<-time_v1[which(time_v1>(cur.time-3) & time_v1<(cur.time+3))] #cur.time # #[max_int_index]
			sn<-NA
			rt_window<-5
#			print(cur.time)
#			print(max(time_v1,na.rm=TRUE))


			half_peak_int<-log10((0.5*peak_int)+1)
			
			
		
			if(length(time_v1)<1){
				return(0)

			}
			if(length(cur.time)<1){
				return(0)
							
			}	
			if(max(time_v1,na.rm=TRUE)<(cur.time+rt_window)){
				maxtimelimit<-max(time_v1,na.rm=TRUE)
			}else{
				maxtimelimit<-(cur.time+rt_window)
			}
			 if(min(time_v1,na.rm=TRUE)>(cur.time-rt_window)){
                                mintimelimit<-min(time_v1,na.rm=TRUE)
                        }else{

				mintimelimit<-(cur.time-rt_window)
			}
                        #
			
			if(length(which(time_v1>mintimelimit & time_v1<maxtimelimit))>0){
				max_int<-max(int_v1[which(time_v1>mintimelimit & time_v1<maxtimelimit)],na.rm=TRUE)
			}else{
				max_int<-100 #max_int_overall
			}

			#max_int<-int_v1[which(time_v1==cur.time)]
			#max_int<-max(int_v1[which(time_v1>(cur.time-5) & time_v1<(cur.time+5))],na.rm=TRUE)
			
			
			
			if(max_int==0){

					max_int<-max(int_v1[which(time_v1>(cur.time-10) & time_v1<(cur.time+10))],na.rm=TRUE)	
			}
			baseline_range<-max(10,bw) #*peak_time
                        peak_time_ind<-which(time_v1>(peak_time-baseline_range) & time_v1<(peak_time+baseline_range))
                        base_time_ind<-seq(1,length(time_v1))
                        base_time_ind<-which(time_v1<(peak_time-baseline_range) | time_v1>(peak_time+baseline_range))
                        
				time_nonzero<-which(is.na(int_v1)==FALSE)
                        base_time_ind<-intersect(base_time_ind,time_nonzero)

                        if(length(base_time_ind)<1){
                                low_end<-0
                                high_end<-100
                                base_o<-1000
                                base_o<-log10(base_o)
                                max_int1<-max_int
                                max_int<-log10(max_int)
					  int_v1<-log10(int_v1+1)

                                sn<-max_int #(max_int-base_o)/(0.5*max_int) #sd(int_v1,na.rm=TRUE)
                        }else{

			low_end<-time_v1[base_time_ind[1]]
                        high_end<-time_v1[base_time_ind[length(base_time_ind)]]
               

                        #base_o<-max(c(to.plot[[i]][which(time_v1==low_end),2],to.plot[[i]][which(time_v1==high_end),2]),na.rm=TRUE)

                        #base_o<-median(int_v1[base_time_ind],na.rm=TRUE)
                        
			base_o<-median(c(int_v1[base_time_ind[1]],int_v1[base_time_ind[length(base_time_ind)]]),na.rm=TRUE)+0.001

			#base_o<-max(3,base_o)
			max_int1<-max_int
                        max_int<-log10(max_int+1)
                        base_o<-log10(base_o+1)
                        #int_v1[base_time_ind]<-log10(int_v1[base_time_ind])
                        int_v1<-log10(int_v1+1)

                        sn<-(peak_int-base_o)/sd(int_v1,na.rm=TRUE)
                        
                        max_int2<-max(int_v1[base_time_ind],na.rm=TRUE)

			int_factor<-(peak_int-3)

			if(int_factor<0){
			
				int_factor<-0
			}

			int_check<-0.5*max(int_v1[peak_time_ind],na.rm=TRUE)-min(int_v1[base_time_ind],na.rm=TRUE)
                        
			if(is.na(max_int)==FALSE){

					
				if(int_factor>2){

					time_diff<-max(10,1*sqrt(bw))
				}else{
                        time_diff<-max(10,(high_end-low_end))
                        time_diff<-1*sqrt(time_diff)
				}
                        #time_diff<-1
		
			#time_diff<-max(10,(high_end-low_end))	
			if(FALSE)
			{
			
			print("here")
			print(sn)
			print(peak_int)
			print(base_o)
			print(sd(int_v1,na.rm=TRUE))
			print(peak_time_ind)
			#print(base_time_ind)
			print(int_factor)
			print(time_diff)
			print(max_int)
			print(high_end)
			print(bw)
			print(low_end)
			} 
	              # 	sn<-0.001+(sn*(int_factor/((time_diff*0.1)))*int_check)

			sn<-sn*(int_factor)*int_check 
			                        
                        #sn<-((max_int-base_o))*sn

			#sn<-(max_int-base_o)/sd(int_v1,na.rm=TRUE)

                       
				                   print("sn")
							print(sn)
							print(peak_int)
							print(base_o)
							print(sd(int_v1,na.rm=TRUE))

			}else{
				  sn<-NA
				

				
			}
			
		}
		
		return(sn)
	}


custom.EIC.plot<-function(aligned, rows = NA, colors = NA, transform = "none", 
    subset = NA, mz.list=NA,time.list=NA,chem.names=NA,minrt=NA, maxrt=NA, min.run, min.pres, max.spline.time.points = 1000,rawprofileloc) 
{
    if (!(transform %in% c("none", "log", "sqrt", "cuberoot"))) {
        message("Invalid transformation. It has to be from: none, log, sqrt, and cuberoot")
        break
    }
        setwd(rawprofileloc)

    if (!is.na(rows[1])) {
        library(splines)
        num.exp <- nrow(summary(aligned$features))
        if (is.na(subset[1])) 
            subset <- 1:num.exp
        colors<-topo.colors(length(subset))
        #print(colors)
        files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
        #files<-gsub(files,pattern=".cdf",replacement="")
        #print(files)
        rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1, 
            2 * num.exp, by = 2)], "_", min.run, "_", min.pres, 
            "_", aligned$mz.tol, ".rawprof", sep = "")
        adj.times <- new("list")
        #print(rawprof.names)
         rawprof.names <-tolower(rawprof.names)
        
        adj.times<-lapply(subset,function(i){
            load(rawprof.names[i])
            times <- unique(raw.prof$labels)
            times <- times[order(times)]
            orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][, 
                3]), 2]
            adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][, 
                3]), 2]
            orig.time <- round(orig.time, 8)
            adjusted.time <- round(adjusted.time, 8)
            ttt <- table(adjusted.time)
            ttt2 <- table(orig.time)
            to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt == 
                1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 == 
                1]))
            if (length(to.use) > max.spline.time.points) 
                to.use <- sample(to.use, max.spline.time.points, 
                  replace = FALSE)
            orig.time <- orig.time[to.use]
            adjusted.time <- adjusted.time[to.use]
            sp <- interpSpline(adjusted.time ~ orig.time)
            adj.times1 <- predict(sp, times)$y
           return(adj.times1)
         })
        
        #cl<-parallel::makeCluster(numcluster)
        #for (n in 1:length(rows)) 

	#print("length adj times")
	#print(length(adj.times))	
        mz.list<-round(mz.list,5)
        l1<-lapply(1:length(rows),function(n)
	{
	
		
            if (n%%4 == 1) {
                par(mfrow = c(2, 2))
            }
            this.row <- rows[n]
            this.intensi <- aligned$final.ftrs[this.row, 5:(num.exp + 
                4)]
            this.times <- aligned$final.times[this.row, 5:(num.exp + 
                4)]
            this.mz <- round(aligned$final.ftrs[this.row, 1],5)
	    this.rowtime<-aligned$final.ftrs[this.row, 2] 
            to.plot <- new("list")
            y.max <- 0
            x.max <- 0
            x.min <- Inf
                sn_vec<-{}
            for (iii in 1:length(subset)) 
            {
                i <- subset[iii]
                if (this.intensi[i] != 0) 
		{
                  load(rawprof.names[i])
                 #raw.prof<-as.data.frame(raw.prof)
                  times <- unique(raw.prof$labels)
                  times <- times[order(times)]
                  mz.diff <- abs(aligned$features2[[i]][, 1] - 
                    this.mz)
                  time.diff <- abs(aligned$features2[[i]][, 2] - 
                    this.times[i])
                  sel <- which(time.diff < 1e-17)
                  if (length(sel) > 1) 
                    sel <- sel[which(mz.diff[sel] == min(mz.diff[sel]))[1]]
                  if (length(sel) == 0) {
                    mz.lim <- aligned$final.ftrs[this.row, c(3, 
                      4)]
                    sel <- which(aligned$features2[[i]][, 1] >= 
                      mz.lim[1] & aligned$features2[[i]][, 1] <= 
                      mz.lim[2] & !is.na(aligned$features2[[i]][, 
                      6]))
                    sub.features <- aligned$features2[[i]][sel, 
                      ]
                    sub.time.diff <- abs(sub.features[, 2] - 
                      this.times[i])
                    sel <- sel[which(sub.time.diff == min(sub.time.diff))][1]
                  }
                  target.mz <- aligned$features2[[i]][sel, 1]
                  target.time <- aligned$features[[i]][sel, 2]
                  time.adjust <- aligned$features2[[i]][sel, 
                    2] - aligned$features[[i]][sel, 2]
                  mz.diff <- abs(raw.prof$masses - target.mz)
                  sel.slice <- raw.prof$grps[which(mz.diff == 
                    min(mz.diff))[1]]
                  sel.time.range <- range(raw.prof$labels[raw.prof$grps == 
                    sel.slice])
                  while (target.time < sel.time.range[1] | target.time > 
                    sel.time.range[2]) {
                    mz.diff[raw.prof$grps == sel.slice] <- 100
                    sel.slice <- raw.prof$grps[which(mz.diff == 
                      min(mz.diff))[1]]
                    sel.time.range <- range(raw.prof$labels[raw.prof$grps == 
                      sel.slice])
                  }
                  sel.time <- raw.prof$labels[raw.prof$grps == 
                    sel.slice]
                  sel.intensi <- raw.prof$intensi[raw.prof$grps == 
                    sel.slice]
                  sel.intensi <- sel.intensi[order(sel.time)]
                  sel.time <- sel.time[order(sel.time)]
                  all.intensi <- times * 0
                  all.intensi[times %in% sel.time] <- sel.intensi
                  #all.times <- adj.times[[i]]
		  sn<-0

		  if(i<length(adj.times)){
                  all.times <- adj.times[[i]]
                        if(length(which(is.na(all.times)==TRUE))>0){
                                sn=NA
                                #return(sn)
                        	sn_vec<-c(sn_vec,sn)
			}
                  }else{

			all.times<-{}
                        sn=NA
                                #return(sn)
                		sn_vec<-c(sn_vec,sn)  
		  }

                  if(is.na(sn)==FALSE && is.na(x.max)==FALSE){
		  if (transform == "log") 
                    all.intensi <- log10(all.intensi + 1)
                  if (transform == "sqrt") 
                    all.intensi <- sqrt(all.intensi)
                  if (transform == "cuberoot") 
                    all.intensi <- all.intensi^(1/3)
                  if (max(all.intensi) > y.max) 
                    y.max <- max(all.intensi)
                  if (max(all.times[all.intensi > 0],na.rm=TRUE) > x.max) 
                    x.max <- max(all.times[all.intensi > 0])
                  if (min(all.times[all.intensi > 0],na.rm=TRUE) < x.min) 
                    x.min <- min(all.times[all.intensi > 0])
                   int_mat1<-cbind(all.times, all.intensi)
                   colnames(int_mat1)<-c("time","intensity")
                   int_mat1<-as.data.frame(int_mat1)    

		  to.plot[[i]] <- int_mat1 #cbind(all.times, all.intensi)
                
                         if(length(to.plot[[i]])>0)
			 {
                         time_v1<-as.numeric(to.plot[[i]][,1])

				
                        int_v1<-as.numeric(to.plot[[i]][,2])
								if(length(which(int_v1==0))>0){
                        int_v1<-replace(int_v1,which(int_v1==0),NA)
				}
				
				d1<-density(time_v1,na.rm=TRUE)
				curtime<- aligned$aligned.ftrs[n,2]
				sn<-get_peakscore(int_v1=int_v1,time_v1=time_v1,cur.time=curtime,bw=d1$bw)
				sn_vec<-c(sn_vec,sn)
                    
			}

				}
			}
        
		}
            
        #print(sn_vec)       
              if(is.na(minrt)==FALSE)
                  {
                        x.min=minrt
                  }else{

                                if(is.na(maxrt)==FALSE){

                                        x.max=maxrt
                                }
                  }
        sn<-max(sn_vec,na.rm=TRUE)

	if(sn=="Inf"){
		sn=10
	}else{
		if(sn=="-Inf"){
			sn=1
		}		
	}	
        #sn<-summary(sn_vec,na.rm=TRUE)
        #sn<-sn[5]
  
	bool_plot=0
	
	  this.row <- rows[n]
        min_v1<-min(length(subset), nrow(summary(to.plot))) 
        for (iii in 1:min_v1) #lapply(1:min_v1,function(iii)
        {
                i <- subset[iii]
                
              
                #sn<-"NA"
               	print(length(to.plot))
		if(i<length(to.plot))
		{

			 
                if (iii == 1 && i<=length(to.plot) && length(to.plot[[i]])>0 && length(to.plot)>0) {
                        test1<-try(to.plot[[i]])
                        if (is(test1, "try-error")){
                        break;
                        }else{
                        
                        curtime<-aligned$aligned.ftrs[this.row,2] #time.list[which(mz.list==this.mz)]
                        curname<-chem.names[which(mz.list==this.mz)]
			
			if(is.na(curname)==TRUE){
				curname<-{}
			}
			if(is.na(x.min)==TRUE){
				x.min=0
			}
			if(is.na(x.max)==TRUE){

				x.max=max(to.plot[[i]][,1],na.rm=TRUE)
			}
			
			bool_plot=1
			plot(to.plot[[i]], xlim = c(x.min, x.max), 
                    ylim = c(0, y.max), type = "l", col = colors[iii], 
                    xlab = "retention time", ylab = "intensity", 
                    main = paste("row", rows[n], ",m/z", 
                      round(this.mz, 5), "time", round(curtime,2),"\n ",curname,"\n ","peak score:",round(sn,2)))
                        
                        }

                }
                else {
			if(i<=length(to.plot) && length(to.plot)>0  && length(to.plot[[i]])>0){
                        test1<-try(to.plot[[i]])
                        if (is(test1, "try-error")){
                        x=1 #break;
                        }else{
				
			curtime<-aligned$aligned.ftrs[this.row,2] # time.list[which(mz.list==this.mz)]
                        curname<-chem.names[which(mz.list==this.mz)]
			
			if(is.na(curname)==TRUE){
				curname<-{}
			}
				if(bool_plot==1){
                                lines(to.plot[[i]], col = colors[iii])
				}else{
				bool_plot=1
			plot(to.plot[[i]], xlim = c(x.min, x.max), 
                    ylim = c(0, y.max), type = "l", col = colors[iii], 
                    xlab = "retention time", ylab = "intensity", 
                    main = paste("row", rows[n], ",m/z", 
                      round(this.mz, 5), "time", round(curtime,2),"\n ",curname,"\n ","peak score:",round(sn,2)))
				
				
				}
			
			}
			}
			
                }

		}
            }
                return(sn)
                #return(to.plot)        #return(sn)
        })
        message("the colors used are:")
        print(paste(subset, ": ", colors[1:length(subset)], sep = ""))
        return(l1)
        }
        
        
}


custom.EIC.plot_vold<-function (aligned, rows = NA, colors = NA, transform = "none", 
    subset = NA, mz.list=NA,time.list=NA,chem.names=NA,minrt=NA, maxrt=NA, min.run, min.pres, max.spline.time.points = 1000,cdfloc) 
{
    if (!(transform %in% c("none", "log", "sqrt", "cuberoot"))) {
        message("Invalid transformation. It has to be from: none, log, sqrt, and cuberoot")
        break
    }
	setwd(cdfloc)

	print("here")
	print(getwd())

    if (!is.na(rows[1])) {
        library(splines)
        num.exp <- nrow(summary(aligned$features))
        if (is.na(subset[1])) 
            subset <- 1:num.exp
	colors<-topo.colors(length(subset))
	print(colors)
        files <- colnames(aligned$final.ftrs)[5:(num.exp + 4)]
	#files<-gsub(files,pattern=".cdf",replacement="")
	print(files)
        rawprof.names <- paste(unlist(strsplit(files, "\\."))[seq(1, 
            2 * num.exp, by = 2)], "_", min.run, "_", min.pres, 
            "_", aligned$mz.tol, ".rawprof", sep = "")
        adj.times <- new("list")
	print(rawprof.names)
	 rawprof.names <-tolower(rawprof.names)
	
	adj.times<-lapply(subset,function(i){
            load(rawprof.names[i])
            times <- unique(raw.prof$labels)
            times <- times[order(times)]
            orig.time <- aligned$features[[i]][!is.na(aligned$features[[i]][, 
                3]), 2]
            adjusted.time <- aligned$features2[[i]][!is.na(aligned$features2[[i]][, 
                3]), 2]
            orig.time <- round(orig.time, 8)
            adjusted.time <- round(adjusted.time, 8)
            ttt <- table(adjusted.time)
            ttt2 <- table(orig.time)
            to.use <- which(adjusted.time %in% as.numeric(names(ttt)[ttt == 
                1]) & orig.time %in% as.numeric(names(ttt2)[ttt2 == 
                1]))
            if (length(to.use) > max.spline.time.points) 
                to.use <- sample(to.use, max.spline.time.points, 
                  replace = FALSE)
            orig.time <- orig.time[to.use]
            adjusted.time <- adjusted.time[to.use]
            sp <- interpSpline(adjusted.time ~ orig.time)
            adj.times1 <- predict(sp, times)$y
       	   return(adj.times1)
	 })
	
	#cl<-parallel::makeCluster(numcluster)
	
        #for (n in 1:length(rows)) 
	l1<-lapply(1:length(rows),function(n){
            if (n%%4 == 1) {
                par(mfrow = c(2, 2))
            }
            this.row <- rows[n]
            this.intensi <- aligned$final.ftrs[this.row, 5:(num.exp + 
                4)]
            this.times <- aligned$final.times[this.row, 5:(num.exp + 
                4)]
            this.mz <- aligned$final.ftrs[this.row, 1]
            to.plot <- new("list")
            y.max <- 0
            x.max <- 0
            x.min <- Inf
		sn_vec<-{}
            for (iii in 1:length(subset)) 
	    {
                i <- subset[iii]
                if (this.intensi[i] != 0) {
                  load(rawprof.names[i])
                 #raw.prof<-as.data.frame(raw.prof)
		  times <- unique(raw.prof$labels)
                  times <- times[order(times)]
                  mz.diff <- abs(aligned$features2[[i]][, 1] - 
                    this.mz)
                  time.diff <- abs(aligned$features2[[i]][, 2] - 
                    this.times[i])
                  sel <- which(time.diff < 1e-17)
                  if (length(sel) > 1) 
                    sel <- sel[which(mz.diff[sel] == min(mz.diff[sel]))[1]]
                  if (length(sel) == 0) {
                    mz.lim <- aligned$final.ftrs[this.row, c(3, 
                      4)]
                    sel <- which(aligned$features2[[i]][, 1] >= 
                      mz.lim[1] & aligned$features2[[i]][, 1] <= 
                      mz.lim[2] & !is.na(aligned$features2[[i]][, 
                      6]))
                    sub.features <- aligned$features2[[i]][sel, 
                      ]
                    sub.time.diff <- abs(sub.features[, 2] - 
                      this.times[i])
                    sel <- sel[which(sub.time.diff == min(sub.time.diff))][1]
                  }
                  target.mz <- aligned$features2[[i]][sel, 1]
                  target.time <- aligned$features[[i]][sel, 2]
                  time.adjust <- aligned$features2[[i]][sel, 
                    2] - aligned$features[[i]][sel, 2]
                  mz.diff <- abs(raw.prof$masses - target.mz)
                  sel.slice <- raw.prof$grps[which(mz.diff == 
                    min(mz.diff))[1]]
                  sel.time.range <- range(raw.prof$labels[raw.prof$grps == 
                    sel.slice])
                  while (target.time < sel.time.range[1] | target.time > 
                    sel.time.range[2]) {
                    mz.diff[raw.prof$grps == sel.slice] <- 100
                    sel.slice <- raw.prof$grps[which(mz.diff == 
                      min(mz.diff))[1]]
                    sel.time.range <- range(raw.prof$labels[raw.prof$grps == 
                      sel.slice])
                  }
                  sel.time <- raw.prof$labels[raw.prof$grps == 
                    sel.slice]
                  sel.intensi <- raw.prof$intensi[raw.prof$grps == 
                    sel.slice]
                  sel.intensi <- sel.intensi[order(sel.time)]
                  sel.time <- sel.time[order(sel.time)]
                  all.intensi <- times * 0
                  all.intensi[times %in% sel.time] <- sel.intensi
                  #all.times <- adj.times[[i]]
                  if(i<length(adj.times)){
                  all.times <- adj.times[[i]]
                        if(length(which(is.na(all.times)==TRUE))>0){
          
	                        sn=NA
                                #return(sn)
                        	 sn_vec<-c(sn_vec,sn)
			}
                  }else{

                        sn=NA
                                #return(sn)
			 sn_vec<-c(sn_vec,sn)
                  }
		  if (transform == "log") 
                    all.intensi <- log10(all.intensi + 1)
                  if (transform == "sqrt") 
                    all.intensi <- sqrt(all.intensi)
                  if (transform == "cuberoot") 
                    all.intensi <- all.intensi^(1/3)
                  if (max(all.intensi) > y.max) 
                    y.max <- max(all.intensi)
                  if (max(all.times[all.intensi > 0]) > x.max) 
                    x.max <- max(all.times[all.intensi > 0])
                  if (min(all.times[all.intensi > 0]) < x.min) 
                    x.min <- min(all.times[all.intensi > 0])
               	   int_mat1<-cbind(all.times, all.intensi)
	           colnames(int_mat1)<-c("time","intensity")
		   int_mat1<-as.data.frame(int_mat1)	
		   to.plot[[i]] <- int_mat1 #cbind(all.times, all.intensi)
                
			 if(length(to.plot[[i]])>0){
                         time_v1<-as.numeric(to.plot[[i]][,1])
                        int_v1<-as.numeric(to.plot[[i]][,2])
			int_v1<-replace(int_v1,which(int_v1==0),NA)

                        #print(time_v1)
                        d1<-density(time_v1)
                        max_int<-max(int_v1,na.rm=TRUE)
                        max_int_index<-which(int_v1==max_int)
                        peak_time<-time_v1[max_int_index]
                        baseline_range<-max(60,d1$bw) #*peak_time
                        peak_time_ind<-which(time_v1>(peak_time-baseline_range) & time_v1<(peak_time+baseline_range))
			base_time_ind<-seq(1,length(time_v1))
			base_time_ind<-which(time_v1<(peak_time-baseline_range) | time_v1>(peak_time+baseline_range))
			time_nonzero<-which(int_v1>0)
			base_time_ind<-intersect(base_time_ind,time_nonzero)
			if(length(base_time_ind)<1){
				low_end<-0
				high_end<-100
				base_o<-1000
				base_o<-log10(base_o)
				max_int1<-max_int
                        	max_int<-log10(max_int)
				sn<-max_int #(max_int-base_o)/(0.5*max_int) #sd(int_v1,na.rm=TRUE)
			}else{

				low_end<-time_v1[base_time_ind[1]]
                        high_end<-time_v1[base_time_ind[length(base_time_ind)]]
			

                      	#base_o<-max(c(to.plot[[i]][which(time_v1==low_end),2],to.plot[[i]][which(time_v1==high_end),2]),na.rm=TRUE)

			base_o<-max(int_v1[base_time_ind],na.rm=TRUE)
			max_int1<-max_int
			max_int<-log10(max_int)
			base_o<-log10(base_o)
			#int_v1[base_time_ind]<-log10(int_v1[base_time_ind])
                        int_v1<-log10(int_v1)

			sn<-(max_int-base_o)/sd(int_v1,na.rm=TRUE)
			}
			int_factor<-2^(max_int-3)
			#int_factor<-1
			time_diff<-max(10,(high_end-low_end))
			time_diff<-1*sqrt(time_diff)
			
               		sn<-sn*(int_factor/((time_diff*0.1))) 
			
			sn<-((max_int-base_o))*sn
			sn_vec<-c(sn_vec,sn)
			}


	
			}
	    }
	      if(is.na(minrt)==FALSE)
                  {
                        x.min=minrt
                  }else{

                                if(is.na(maxrt)==FALSE){

                                        x.max=maxrt
                                }
                  }
	sn<-max(sn_vec,na.rm=TRUE)
     #sn<-summary(sn_vec,na.rm=TRUE)
	#sn<-sn[5]
  
	min_v1<-min(length(subset), nrow(summary(to.plot))) 
	#for (iii in 1:min_v1) 
	lapply(1:min_v1,function(iii)
	{
                i <- subset[iii]
		
		print("i is ")
		print(i)
	        #sn<-"NA"
		
                if (iii == 1) {
			test1<-try(to.plot[[i]])
			if (is(test1, "try-error")){
			break;
			}else{
			
			curtime<-time.list[which(mz.list==this.mz)]
			curname<-chem.names[which(mz.list==this.mz)]
                  plot(to.plot[[i]], xlim = c(x.min, x.max), 
                    ylim = c(0, y.max), type = "l", col = colors[iii], 
                    xlab = "retention time", ylab = "intensity", 
                    main = paste("row", rows[n], ",m/z", 
                      round(this.mz, 5), "time", round(curtime,2),"\n ",curname,"\n ","peak score:",round(sn,2)))
			
			}

		}
                else {
			test1<-try(to.plot[[i]])
			if (is(test1, "try-error")){
			break;
			}else{
				lines(to.plot[[i]], col = colors[iii])
			}
		}
            })
		return(sn)
        	#return(to.plot)	#return(sn)
	})
        message("the colors used are:")
        print(paste(subset, ": ", colors[1:length(subset)], sep = ""))
    	return(l1)
	}
	
	
}


adjust.time.custom<-function(features, mz.tol = NA, chr.tol = NA, colors = NA, find.tol.max.d = 1e-04, 
    max.align.mz.diff = 0.01,template=NA,reference_sample=NA) 
{


    num.exp <- nrow(summary(features))
    if (num.exp > 1) {
        par(mfrow = c(2, 2))
        plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", 
            main = "", axes = FALSE)
        text(x = 0, y = 0, "Retention time \n adjustment", cex = 2)
        a <- summary(features)
        sizes <- as.numeric(a[, 1])/ncol(features[[1]])
        sizes <- cumsum(sizes)
        sel <- length(sizes)
        mz <- chr <- lab <- rep(0, sizes[sel])
        sizes <- c(0, sizes)
        for (i in 1:sel) {
            mz[(sizes[i] + 1):sizes[i + 1]] <- features[[i]][, 
                1]
            chr[(sizes[i] + 1):sizes[i + 1]] <- features[[i]][, 
                2]
            lab[(sizes[i] + 1):sizes[i + 1]] <- i
        }
        if (is.na(mz.tol)) {
            mz.tol <- find.tol(mz, uppermost = find.tol.max.d)
        }
        else {
            plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", 
                main = "m/z tolerance level given", axes = FALSE)
            text(x = 0, y = 0, mz.tol, cex = 1.2)
        }
        if (!is.na(chr.tol)) {
            plot(c(-1, 1), c(-1, 1), type = "n", xlab = "", ylab = "", 
                main = "retention time \n tolerance level given", 
                axes = FALSE)
            text(x = 0, y = 0, chr.tol, cex = 1.2)
        }
        all.ft <- find.tol.time(mz, chr, lab, num.exp = num.exp, 
            mz.tol = mz.tol, chr.tol = chr.tol, max.mz.diff = max.align.mz.diff)
        chr.tol <- all.ft$chr.tol
        message("**** performing time correction ****")
        message(paste("m/z tolerance level: ", mz.tol))
        message(paste("time tolerance level:", chr.tol))
        for (i in 1:num.exp) {
            this <- features[[i]]
            sel <- which(all.ft$lab == i)
            that <- cbind(all.ft$mz[sel], all.ft$chr[sel], all.ft$grps[sel])
            this <- this[order(this[, 1], this[, 2]), ]
            that <- that[order(that[, 1], that[, 2]), ]
            this <- cbind(this, rep(i, nrow(this)), that[, 3])
            features[[i]] <- this
        }
        num.ftrs <- as.vector(table(all.ft$lab))
	
	if(is.na(template)==TRUE){
		template <- which(num.ftrs == max(num.ftrs))[1]
	}
        message(paste("the template is sample", template))
        if (is.na(colors[1])) 
            colors <- c("red", "blue", "dark blue", "orange", 
                "green", "yellow", "cyan", "pink", "violet", 
                "bisque", "azure", "brown", "chocolate", rep("grey", 
                  num.exp))
	
		if(is.na(reference_sample)==FALSE){
	
			cnames<-c("mz","pos","sd1","sd2","area")
			
			temp_mat<-cbind(reference_sample[,1:2],rep(3,nrow(reference_sample)),rep(5,nrow(reference_sample)),reference_sample[,3])
			
			colnames(temp_mat)<-as.character(cnames)
			
			features[[num.exp+1]]<-temp_mat
			
			template<-num.exp+1
		}
        candi <- features[[template]][, 1:2]
	cor.method="apLCMS"
        features.2 <- foreach(j = 1:num.exp, .export = ls(envir = globalenv())) %dopar% 
            {
                this.feature <- features[[j]]
                if (j != template) {
                  this.comb <- rbind(cbind(candi, rep(template, 
                    nrow(candi))), cbind(this.feature[, 1:2], 
                    rep(j, nrow(this.feature))))
                  this.comb <- this.comb[order(this.comb[, 1]), 
                    ]
                  l <- nrow(this.comb)
                  sel <- which(this.comb[2:l, 1] - this.comb[1:(l - 
                    1), 1] < mz.tol * this.comb[1:(l - 1), 1] * 
                    2 & abs(this.comb[2:l, 2] - this.comb[1:(l - 
                    1), 2]) < chr.tol & this.comb[2:l, 3] != 
                    this.comb[1:(l - 1), 3])
                  if (length(sel) < 20) {
                    cat("too few, aborted")
                  }
                  else {
		  
		      if(cor.method=="dtw"){
		      
			reference<-candi[order(candi[, 2]), 2]
			query<-this.feature[order(this.feature[, 2]), 2]
			
			this.corrected<-dtw(query,reference,keep=TRUE)
			
			this.corrected<-warp(this.corrected)
			
		      }else{
                    all.ftr.table <- cbind(this.comb[sel, 2], 
                      this.comb[sel + 1, 2])
                    to.flip <- which(this.comb[sel, 3] == j)
                    temp <- all.ftr.table[to.flip, 2]
                    all.ftr.table[to.flip, 2] <- all.ftr.table[to.flip, 
                      1]
                    all.ftr.table[to.flip, 1] <- temp
                    cat(c("sample", j, "using", nrow(all.ftr.table), 
                      ","))
                    if (j%%3 == 0) 
                      cat("\n")
                    all.ftr.table <- all.ftr.table[order(all.ftr.table[, 
                      2]), ]
		  
                    this.dev <- all.ftr.table[, 2]
                    aver.time <- all.ftr.table[, 1] - this.dev
                    this.feature <- this.feature[order(this.feature[, 
                      2], this.feature[, 1]), ]
                    this.corrected <- this.old <- this.feature[, 
                      2]
                    to.correct <- this.old[this.old >= min(this.dev) & 
                      this.old <= max(this.dev)]
                    this.smooth <- ksmooth(this.dev, aver.time, 
                      kernel = "normal", bandwidth = (max(this.dev) - 
                        min(this.dev))/5, x.points = to.correct)
                    this.corrected[this.old >= min(this.dev) & 
                      this.old <= max(this.dev)] <- this.smooth$y + 
                      to.correct
                    this.corrected[this.old < min(this.dev)] <- this.corrected[this.old < 
                      min(this.dev)] + mean(this.smooth$y[this.smooth$x == 
                      min(this.smooth$x)])
                    this.corrected[this.old > max(this.dev)] <- this.corrected[this.old > 
                      max(this.dev)] + mean(this.smooth$y[this.smooth$x == 
                      max(this.smooth$x)])
		      }
                    this.feature[, 2] <- this.corrected
                    this.feature <- this.feature[order(this.feature[, 
                      1], this.feature[, 2]), ]
                  }
                }
                this.feature
            }
    }
    else {
        message("Only one sample.  No need to correct for time.")
    }
    plot(range(features[[1]][, 2]), c(-chr.tol, chr.tol), type = "n", 
        ylab = "Retention time deviation", xlab = "Original Retention time")
    for (i in 1:num.exp) {
        features[[i]] <- features[[i]][order(features[[i]][, 
            1], features[[i]][, 2]), ]
        points(features[[i]][, 2], features.2[[i]][, 2] - features[[i]][, 
            2], col = colors[i], cex = 0.2)
    }
    return(features.2)
}


