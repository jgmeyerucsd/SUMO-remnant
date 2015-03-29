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
mgpl.pos<-proteinPositions(object=mgpl,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=FALSE,
		name="mgplus.localized.tsv")


mgpl.s<-read.PepProph(input=files[85])
mgpl.s.all<-read.PepProph(input=files[79])
mgpl.s<-removeLowLocalization(object=mgpl.s,minscore=0.75)
mgpl.s.pos<-proteinPositions(object=mgpl.s,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=FALSE,
		name="mgplus.localized.tsv")
mgpl.s.ave<-summarizeProtPositions(object=mgpl.s.pos)
mgpl.s.sig<-getSignificant(allIDs=mgpl.s.all,targetIDs=mgpl.s.ave)


mgpl<-read.PepProph(input=files[100])
mgpl.all<-read.PepProph(input=files[98])
mgpl<-removeLowLocalization(object=mgpl,minscore=0.75)

mgpl<-proteinPositions(object=mgpl,
		fasta="F:/MSGFplus.20140716/database/110712_human.cc.fasta",
		writetsv=FALSE,
		name="mgplus.localized.tsv")
mgpl<-summarizeProtPositions(mgpl)

mgpl.pos<-proteinPositions(object=mgpl,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=FALSE,
		name="mgplus.localized.tsv")

mgpl.ave<-summarizeProtPositions(object=mgpl.pos)
mgpl.sig<-getSignificant(stdev=0.9404181,allIDs=mgpl.all,targetIDs=mgpl.ave)

mgneg<-read.PepProph(input=files[11])
mgneg.all<-read.PepProph(input=files[8])
mgneg@modposition.peptide
mgneg<-removeLowLocalization(object=mgneg,minscore=0.75)
mgneg.pos<-proteinPositions(object=mgneg,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=FALSE,
		name="mgneg.localized.tsv")
mgneg.pos@modposition.protein
mgneg.sum<-summarizeProtPositions(object=mgneg.pos)
mgneg.sum@modindex
mgneg.sig<-getSignificant(allIDs=mgneg.all,targetIDs=mgneg.sum,name="MG132neg.norm.significant.tsv")
mgneg@modindex
mgneg@modposition.peptide


targetIDs@modindex
mgpl@modindex
#### need to fix this to ignore positions where weighted ave==0


files<-list.files()
files

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
