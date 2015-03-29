setClass("pepsum",representation(summary="matrix",scans="ANY",specdir="character",
	specfiles="list",data="ANY", totalheat="matrix",residues="character",
	sequence="ANY",score="character",filetype="character",scanvec="ANY",
	modposition.protein="ANY", chargevec="ANY",proteinmodlist="ANY",index="ANY", 
	fraction="character",pepvec="ANY",provec="ANY",filenames="ANY",p1.index="ANY",
	modposition.peptide="ANY",countlist="list",pepraw="ANY",ionCoverage="list",
	ionCovSum="list",modindex="list"))


getwd()
setwd("C:/Users/Jgmeyer/Documents/R/silacquant")
files<-list.files()
files


mgpl<-read.PepProph(input=files[14])
mgpl.all<-read.PepProph(input=files[12])
mgpl<-removeLowLocalization(object=mgpl,minscore=0.75)
mgpl<-proteinPositions(object=mgpl,
		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
		writetsv=FALSE,
		name="mgplus.localized.tsv")

mgneg@pepvec

### 1 
### 2 remove peptides with ptm location <0.9
### 3 assign postition in protein
### 4 report motif
### 5 compare positions in proteins
### 





object@data[1,]


mgneg<-removeLowLocalization(object=mgneg)
test2@data[,"peptide"]
test2@data[1,]

fasta="C:/MSGFplus/database/110712_human.cc.fasta"
proteinPositions(object=mgpl,name="mgpl.positions.tsv")->mgpl


summarizeProtPositions()




#	Tall<-object



reportMotif=function(object){
	datamat<-object@data
	
	
	
	for(i in 1:peptempvec.l){
		#print(peptempvec[i])
		#print(datamat[i,"protein"])		
		currentpro<-substr(start=4,stop=9,datamat[i,"protein"])
		currentprot<-fastaobj[[which(fastaacc==substr(start=4,stop=9,datamat[i,"protein"]))]][1]
		currentpep<-peptempvec[i]
		#print(currentpro)
		#print(currentpep)
		matched<-matchPattern(currentpep,currentprot)
		proteinposition[i]<-matched@ranges[[1]][1]-1+as.numeric(as.vector(datamat[i,1]))
		print(i)
		}
	
	
	}


getModPositions=function(object){
	datamat<-object@data
	proteins<-substr(start=4,stop=9,as.character(datamat[,"protein"]))
	uniqproteins<-substr(start=4,stop=9,as.character(unique(datamat[,"protein"])))
	positionlist<-list()
	x<-uniqproteins[1]
	for(x in uniqproteins){
		print(x)
		grep(x, proteins)
		#positionlist[[x]]<-grep(x, proteins)
		positionlist[[x]]<-as.numeric(as.character(unique(datamat[grep(x, proteins),1])))
		}
	object@proteinmodlist<-positionlist
	return(object)
	}

walltest<-getModPositions(object=wall)
wcsttest<-getModPositions(object=wcst)
length(wcsttest@proteinmodlist)
length(unlist(wcsttest@proteinmodlist))


whprptest<-getModPositions(object=whprp)
length(whprptest@proteinmodlist)
length(unlist(whprptest@proteinmodlist))

walltest@proteinmodlist
talltest<-getModPositions(object=tall)
length(talltest@proteinmodlist)

tcsttest<-getModPositions(object=tcst)
length(tcsttest@proteinmodlist)
length(unlist(tcsttest@proteinmodlist))

thprptest<-getModPositions(object=thprp)
length(thprptest@proteinmodlist)
length(unlist(thprptest@proteinmodlist))

### need to fix peptides with n-terminal mods before actually running this
### 


#object2<-object

### proteins object1 must be shorter
modPositionOverlap=function(object1=wall,object2=tall){
	#require(made4)
	### get mod lists
	object1<-getModPositions(object1)
	object2<-getModPositions(object2)
	modlist1<-object1@proteinmodlist
	modlist2<-object2@proteinmodlist
	
	### compare mod lists
	samelist<-list()
	uniqlist<-list()
	proteins1.len <- length(modlist1)
	proteins2.len <- length(modlist2)
	protein.names1<-names(modlist1)
	protein.names2<-names(modlist2)
	
	
	for(x in protein.names1){
		print(x)
		
		if(length(modlist2[[x]])==0){
			uniqlist[x]<-modlist1[x]
			}
	
		### if that protein is present in the second set
		if(length(modlist2[[x]])>=1){
			## add the intersection to samelist
			if(length(intersect(modlist1[[x]],modlist2[[x]]))>0){
				samelist[[x]]<-intersect(modlist1[[x]],modlist2[[x]])
				}
			#print(samelist[[x]])
			if(length(samelist[[x]])<length(modlist1[[x]])){
				uniqlist[[x]]<-setdiff(modlist1[[x]],modlist2[[x]])
				}
			}
		
		}
		
		
	modlist <- 
	
	
	#### output the overlap and unique object1
	}

write.table(uniqlist,file="testuniq.txt")
lapply(uniqlist, cat, "\n", file="test.txt", append=TRUE)
lapply(names(uniqlist), write, "names.txt", append=TRUE, ncolumns=1000)
length(uniqlist)
length(samelist)

files<-list.files()
files

phosphositeplus<-object
<-object1


readknownsumo=function(input=files[15]){
	object<-new("pepsum")
	object@filetype="knownsites"
	
	object@data<-read.delim(input,header=T,skip=0)
	

	positionvec<-as.numeric(object@data[,"Position"])
	#object@data[,"MOD_RSD"]<-as.numeric(positionvec)
	datamat<-object@data
	accvec<-substr(as.character(object@data[,2]),start=1,stop=6)
	uniqproteins<-unique(accvec)
	x<-uniqproteins[1]
	positionlist<-list()
	for(x in uniqproteins){
		#print(x)
		#grep(x, proteins)
		#positionlist[[x]]<-grep(x, proteins)
		positionlist[[x]]<-datamat[grep(x, accvec),"Position"]
		}
	object@proteinmodlist<-positionlist

	}


readknownsumo=function(input=files[15]){
	object<-new("pepsum")
	object@filetype="knownsites"
	
	object@data<-read.delim(input,header=T,skip=3)
	

	positionvec<-substr(object@data[,"MOD_RSD"],start=2,stop=nchar(as.character(object@data[,"MOD_RSD"])))
	object@data[,"MOD_RSD"]<-as.numeric(positionvec)
	datamat<-object@data
	accvec<-substr(as.character(object@data[,2]),start=1,stop=6)
	uniqproteins<-unique(accvec)
	x<-uniqproteins[1]
	positionlist<-list()
	for(x in uniqproteins){
		#print(x)
		grep(x, proteins)
		#positionlist[[x]]<-grep(x, proteins)
		positionlist[[x]]<-datamat[grep(x, accvec),"MOD_RSD"]
		}
	object@proteinmodlist<-positionlist

	}



