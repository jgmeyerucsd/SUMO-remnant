
#### change working directory
setwd("C:/Users/Jesse/my documents/SUMO/2014/silacquant")
setwd("C:/Users/JGmeyer/my documents/R/silacquant")

files<-list.files()
files



### create class for peptide identification data
###

setClass("pepsum",representation(summary="matrix",scans="ANY",specdir="character",
	specfiles="list",data="ANY", totalheat="matrix",residues="character",
	sequence="ANY",score="character",filetype="character",scanvec="ANY",
	modposition.protein="ANY", chargevec="ANY",proteinmodlist="ANY",index="ANY", 
	fraction="character",pepvec="ANY",provec="ANY",filenames="ANY",p1.index="ANY",
	modposition.peptide="ANY",countlist="list",pepraw="ANY",ionCoverage="list",
	ionCovSum="list",modindex="list",modsummary="list"))


mgpl.pos<-proteinPositions(object=mgpl,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=FALSE,
		name="mgplus.localized.tsv")

#### first read the peptide prophet results for all IDs and SUMO-remnant-only IDs into pepsum objects
mgpl.s<-read.PepProph(input=files[85])
mgpl.s.all<-read.PepProph(input=files[79])

#### remove those IDs with localization scores below an arbitrary value
mgpl.s<-removeLowLocalization(object=mgpl.s,minscore=0.75)

#### determine to position of each site ID in their protein
mgpl.s.pos<-proteinPositions(object=mgpl.s,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=FALSE,
		name="mgplus.localized.tsv")

#### summarize and combine the unique protein site IDs
mgpl.s.ave<-summarizeProtPositions(object=mgpl.s.pos)

#### compute the weighted average of sites identified by with multiple sequences,
#### then determine which sites are outside at least one standard deviation
#### will produce a plot of quantification values as given in figure 3a
mgpl.s.sig<-getSignificant(allIDs=mgpl.s.all,targetIDs=mgpl.s.ave)

#### correlate the sites between two different replicates, plot their quant values, and determine a linear fit 
newobject<-correlatepositions(object1=MG132.single, object2=MG132.hprp)





uv0<-read.PepProph(input=files[39])
uv0.all<-read.PepProph(input=files[38])
uv0<-removeLowLocalization(object=uv0,minscore=0.75)
uv0<-proteinPositions(object=uv0,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=FALSE,
		name="mgplus.localized.tsv")

uv0<-summarizeProtPositions(object=uv0)
uv0.sig<-getSignificant(allIDs=uv0.all,targetIDs=uv0,writetsv=T,name="uv0.norm.sig.tsv")

uv2<-read.PepProph(input=files[44])
uv2.all<-read.PepProph(input=files[40])
uv2<-removeLowLocalization(object=uv2,minscore=0.75)
uv2<-proteinPositions(object=uv2,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=F,
		name="UV2.localized.tsv")

uv2<-summarizeProtPositions(object=uv2)
uv2.sig<-getSignificant(allIDs=uv2.all,targetIDs=uv2,writetsv=T,name="uv2.norm.sig.tsv")

uv8<-read.PepProph(input=files[46])
uv8.all<-read.PepProph(input=files[45])
uv8<-removeLowLocalization(object=uv8,minscore=0.75)
uv8<-proteinPositions(object=uv8,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=F,
		name="UV8.localized.tsv")

uv8<-summarizeProtPositions(object=uv8)
uv8.sig<-getSignificant(allIDs=uv8.all,targetIDs=uv8,writetsv=T,name="uv8.norm.sig.tsv")



length(object@data[1,])
mgpl@data[1,]


#### part to generate histograms as shown in figure 3a
?hist
par(mfcol=c(2,1))
hist(log(object@data[,26],2),breaks=200,main="weighted average")
upperlim2sd<-mean(log(object@data[,26],2))+2*sd(log(object@data[,26],2))
lowerlim2sd<-mean(log(object@data[,26],2))-2*sd(log(object@data[,26],2))
upperlim1sd<-mean(log(object@data[,26],2))+1*sd(log(object@data[,26],2))
lowerlim1sd<-mean(log(object@data[,26],2))-1*sd(log(object@data[,26],2))
	abline(v=mean(log(object@data[,26],2)),col="blue")	
	abline(v=upperlim2sd,col="red")
	abline(v=lowerlim2sd,col="red")
	abline(v=upperlim1sd,col="red")
	abline(v=lowerlim1sd,col="red")


hist(log(object@data$xpress,2),breaks=200,main="raw xpress scores")
upperlim2sd<-mean(log(object@data$xpress,2))+2*sd(log(object@data$xpress,2))
lowerlim2sd<-mean(log(object@data$xpress,2))-2*sd(log(object@data$xpress,2))
upperlim1sd<-mean(log(object@data$xpress,2))+1*sd(log(object@data$xpress,2))
lowerlim1sd<-mean(log(object@data$xpress,2))-1*sd(log(object@data$xpress,2))
	abline(v=mean(log(object@data$xpress,2)),col="blue")	
	abline(v=upperlim2sd,col="red")
	abline(v=lowerlim2sd,col="red")
	abline(v=upperlim1sd,col="red")
