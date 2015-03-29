normalizeRatios=function(allIDs=UV2all,targetIDs=UV2digly){
	logRatiosAll<-log(allIDs@data$xpress,base=2)
	allIDs@data<-cbind(allIDs@data,logRatiosAll)
	logRatiosTarget<-log(targetIDs@data$xpress,base=2)
	targetIDs@data<-cbind(targetIDs@data,logRatiosTarget)
	#### plot the distributions
	par(mfrow=c(2,2),cex.axis=c(1.5))
	#### plot all IDs
	hist(logRatiosAll,breaks=seq(-11,10,by=0.1),main="All IDs Dist")
	upperlim<-mean(logRatiosAll)+2*sd(logRatiosAll)
	lowerlim<-mean(logRatiosAll)-2*sd(logRatiosAll)
	oneSDup<-mean(logRatiosAll)+1*sd(logRatiosAll)
	oneSDdown<-mean(logRatiosAll)-1*sd(logRatiosAll)
	abline(v=mean(logRatiosAll),col="blue")	
	abline(v=upperlim,col="red")
	abline(v=lowerlim,col="red")
	abline(v=oneSDup,col="red")
	abline(v=oneSDdown,col="red")
	#### plot diGly IDs
	hist(logRatiosTarget,breaks=seq(-11,10,by=0.1),main="DiGly Dist")
	abline(v=mean(logRatiosAll),col="blue")	
	abline(v=upperlim,col="red")
	abline(v=lowerlim,col="red")
	abline(v=oneSDup,col="red")
	abline(v=oneSDdown,col="red")
	#### correct the distributions
	correction<-mean(logRatiosAll)
	allIDs@data$logRatiosAll=allIDs@data$logRatiosAll-correction
	#### addd plot of corrected all ID distribution
	hist(allIDs@data$logRatiosAll,breaks=seq(-11,10,by=0.1),main="DiGly Dist")
	abline(v=mean(allIDs@data$logRatiosAll),col="blue")
	upperlim<-mean(allIDs@data$logRatiosAll)+2*sd(allIDs@data$logRatiosAll)
	lowerlim<-mean(allIDs@data$logRatiosAll)-2*sd(allIDs@data$logRatiosAll)
	oneSDup<-mean(allIDs@data$logRatiosAll)+1*sd(allIDs@data$logRatiosAll)
	oneSDdown<-mean(allIDs@data$logRatiosAll)-1*sd(allIDs@data$logRatiosAll)
	abline(v=upperlim,col="red")
	abline(v=lowerlim,col="red")
	abline(v=oneSDup,col="red")
	abline(v=oneSDdown,col="red")
	### correct target IDs
	targetIDs@data$logRatiosTarget=targetIDs@data$logRatiosTarget-correction
	hist(targetIDs@data$logRatiosTarget,breaks=seq(-11,10,by=0.1),main="DiGly Dist")
	abline(v=upperlim,col="red")
	abline(v=lowerlim,col="red")	
	abline(v=mean(allIDs@data$logRatiosAll),col="blue")
	abline(v=oneSDup,col="red")
	abline(v=oneSDdown,col="red")
	#mean(targetIDs@data$logRatiosTarget)
	print("mean")
	print(correction)
	print("st.dev")
	print(sd(logRatiosAll))
	allIDs@data$xpress
	return(list(allIDs,targetIDs))
	}
