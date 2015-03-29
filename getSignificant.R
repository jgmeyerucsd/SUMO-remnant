### update this to output three tables - 
### all PSMs, unique sites only significant changers, and unique sites only nonsignificant changers

getSignificant=function(allIDs=mgpl.all,targetIDs=mgpl.ave,name="MG132plus.filtered",
		writetsv=F,
		stdev=0.9404181){
	logRatiosAll<-log(allIDs@data$xpress,base=2)
	allIDs@data<-cbind(allIDs@data,logRatiosAll)
	meanAll<-mean(logRatiosAll)
	medianAll<-median(logRatiosAll)
	logRatiosAllnorm<-logRatiosAll-medianAll
	norm.median.all<-median(logRatiosAllnorm)
	### get distribution of weighted averages for unique sites
	#targetIDs@modindex
	numuniq<-length(unique(unlist(targetIDs@modindex)))-1
	unique.position.singleline<-rep(0,times=numuniq)
	unique.position.singleweight<-rep(0,times=numuniq)
	#targetIDs@data$weighted.ratios
	for(i in 1:numuniq){
		unique.position.singleline[i]<-which(targetIDs@modindex==i)[1]
		unique.position.singleweight[i]<-log(targetIDs@data[which(targetIDs@modindex==i)[1],"weighted.ratios"],base=2)
		}
	median(unique.position.singleweight)
	range(unique.position.singleweight)

	norm.uniq.pos.singleweight<-unique.position.singleweight-medianAll
	#2*sd(norm.uniq.pos.singleweight)
	if(length(stdev)==1){
		upperlim<-norm.median.all+2*stdev
		lowerlim<-norm.median.all-2*stdev
		onesdup<-norm.median.all+1*stdev
		onesddown<-norm.median.all-1*stdev
		}
	if(length(stdev)==0){
		upperlim<-norm.median.all+2*sd(logRatiosAllnorm)
		lowerlim<-norm.median.all-2*sd(logRatiosAllnorm)
		onesdup<-norm.median.all+1*sd(logRatiosAllnorm)
		onesddown<-norm.median.all-1*sd(logRatiosAllnorm)
		}
	#?pairlist
	#log(10)
	
	#unique.position.singleline[1]
	#norm.uniq.pos.singleweight[1]
	names(norm.uniq.pos.singleweight)<-unique.position.singleline
	sort.norm.uniq.pos.singleweight<-sort(norm.uniq.pos.singleweight*-1)
	#median(sort.norm.uniq.pos.singleweight)
	
	plot(x=1:numuniq,y=sort.norm.uniq.pos.singleweight)
	abline(h=upperlim,col="red")
	abline(h=lowerlim,col="red")
	###
	print("one standard deviation is")
	print(onesdup)
	abline(h=onesdup,col="red")
	abline(h=onesddown,col="red")
	#### now flag those IDs that are outside 2 std. dev.
	idlen<-length(targetIDs@data$xpress)

	oversigma<-rep(0,times=idlen)
	count=0
	#### loop through and test which are outside 2 std. dev
	for(i in 1:numuniq){
		if(norm.uniq.pos.singleweight[i]>=onesdup | norm.uniq.pos.singleweight[i]<=onesddown){
			#print(norm.uniq.pos.singleweight[i])
			count=count+1
			#print(count)
			oversigma[unique.position.singleline[i]]<-1
			}
		}
	over2sigma<-rep(0,times=idlen)
	count=0
	#### loop through and test which are outside 2 std. dev
	for(i in 1:numuniq){
		if(norm.uniq.pos.singleweight[i]>=upperlim | norm.uniq.pos.singleweight[i]<=lowerlim){
			#print(norm.uniq.pos.singleweight[i])
			count=count+1
			#print(count)
			oversigma[unique.position.singleline[i]]<-2
			}
		}
	# for those missing their weighted average value, put in xpress value
	for(i in 1:length(targetIDs@data[,"weighted.ratios"])){
		if(targetIDs@data[i,"weighted.ratios"]==0){
			targetIDs@data[i,"weighted.ratios"]<-targetIDs@data[i,"xpress"]
			}
		}

	targetIDs@data<-cbind(targetIDs@data,log.ratios=log(targetIDs@data[,"weighted.ratios"],base=2)-medianAll)
	targetIDs@data<-cbind(targetIDs@data,oversigma)
	#targetIDs@data<-cbind(targetIDs@data,over1sigma)

	### part to output only unique significant changers
	#unique.position.singleline
	#targetIDs@data[,"log.ratios"]
	uniquelines<-targetIDs@data[unique.position.singleline,]
	changelines<-uniquelines[uniquelines[,"oversigma"]>=1,]
	nochangelines<-uniquelines[uniquelines[,"oversigma"]==0,]


	### part to output only unique non changers


	### now have new object with binary whether outside 2*sigma
	print(i)
	#### maybe write new table?
	print("unique sites")
	print(numuniq)
	print("over 2 stdev")
	print(length(oversigma[oversigma==2]))
	print("over 1 stdev")
	print(length(oversigma[oversigma==1]))
	if(writetsv){
		write.table(file=paste(name,".all.tsv",sep=""),targetIDs@data,quote=F,sep="\t",col.names=T,row.names=F)
		write.table(file=paste(name,".changing.tsv",sep=""),changelines,quote=F,sep="\t",col.names=T,row.names=F)
		write.table(file=paste(name,".nochange.tsv",sep=""),nochangelines,quote=F,sep="\t",col.names=T,row.names=F)
		}
	return(list(allIDs,targetIDs,nochangelines,changelines))
	}
