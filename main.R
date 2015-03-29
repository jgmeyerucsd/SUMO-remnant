setwd("C:/Users/Jesse/my documents/SUMO/2014/silacquant")
setwd("C:/Users/JGmeyer/my documents/R/silacquant")

files<-list.files()
files


setClass("pepsum",representation(summary="matrix",scans="ANY",specdir="character",
	specfiles="list",data="ANY", totalheat="matrix",residues="character",
	sequence="ANY",score="character",filetype="character",scanvec="ANY",
	modposition.protein="ANY", chargevec="ANY",proteinmodlist="ANY",index="ANY", 
	fraction="character",pepvec="ANY",provec="ANY",filenames="ANY",p1.index="ANY",
	modposition.peptide="ANY",countlist="list",pepraw="ANY",ionCoverage="list",
	ionCovSum="list",modindex="list",modsummary="list"))

#### basic script to read ptmprophet table into object@data
###########
read.PTMproph=function(file=files[1]){
	object<-new("pepsum")
	object@filetype="PTMprophet"
	object@data<-read.delim(file[1],fill=TRUE,header=TRUE,row.names=NULL)
	object@fraction<-file
	return(object)
	}

MGpl.all<-read.PTMproph(file=files[3])
MGpl.target<-read.PTMproph(file=files[14])


normalizeRatios(allIDs=MGpl.all,targetIDs=wall)->MGpl.norm
flagRatios(allIDs=MGpl.norm[[1]],targetIDs=MGpl.norm[[2]],name="MGpl.norm.tsv")->MGpl.norm.flag

MGneg.all<-read.PTMproph(file=files[1])
MGneg.target<-read.PTMproph(file=files[2])
normalizeRatios(allIDs=MGneg.all,targetIDs=mgneg)->mgneg.norm

MGneg.target.norm[[2]]

flagRatios(allIDs=mgneg.norm[[1]],targetIDs=mgneg.norm[[2]])->mgneg.norm.flag


UV2all<-read.PTMproph()
UV2digly<-read.PTMproph(file=files[2])
#### plot all ratios


??log
?hist
object<-UV2digly


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

median(logRatiosAll)

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


