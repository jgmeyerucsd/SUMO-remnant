flagRatios=function(allIDs=UV2all,targetIDs=UV2digly,name="MGneg.norm.tsv"){
	logRatiosAll<-allIDs@data$logRatiosAll
	logRatiosTarget<-targetIDs@data$logRatiosTarget

	twoSDabove<-mean(logRatiosAll)+2*sd(logRatiosAll)
	twoSDbelow<-mean(logRatiosAll)-2*sd(logRatiosAll)
	oneSDabove<-mean(logRatiosAll)+1*sd(logRatiosAll)
	oneSDbelow<-mean(logRatiosAll)-1*sd(logRatiosAll)

	idlen<-length(targetIDs@data$xpress)
	oversigma<-rep(0,times=idlen)
	#over2sigma<-rep(0,times=idlen)

	for(i in 1:idlen){
		if(targetIDs@data$logRatiosTarget[i]>=oneSDabove | targetIDs@data$logRatiosTarget[i]<=oneSDbelow){
			#print(i)
			oversigma[i]=1
			}
		if(targetIDs@data$logRatiosTarget[i]>=twoSDabove | targetIDs@data$logRatiosTarget[i]<=twoSDbelow){
			#print(i)
			oversigma[i]=2
			}	
		}
	#as.data.frame(over2sigma)
	#diglyIDobject@data[1:100,]
	#newobject@data[1:100,]
	targetIDs@data<-cbind(x=targetIDs@data,y=as.data.frame(oversigma))
	### now have new object with binary whether outside 2*sigma
	print(i)
	#### maybe write new table?
	print("over 2 stdev")
	print(length(oversigma[oversigma==2]))	
	print("over 1 stdev")
	print(length(oversigma[oversigma==1]))	
	write.table(file=name,targetIDs@data,quote=F,sep="\t",col.names=T,row.names=F)
	return(list(allIDs,targetIDs))
	}
