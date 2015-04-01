
#### function to determine the mod postion within a protein using the fasta DB
####  usage:
### 	#### determine to position of each site ID in their protein

### mgpl.s.pos<-proteinPositions(object=mgpl.s,
###		fasta="C:/MSGFplus/database/110712_human.cc.fasta",
###		writetsv=FALSE,
###		name="mgplus.localized.tsv")
		
		

proteinPositions=function(fasta="C:/Users/JgMeyer/Documents/R/pepsum/pepsum/110712_human.cc.fasta",
	object=mgpl.single.pos,
	writetsv=FALSE,
	name="MG132.all.filtered.tsv"){

	### read in fasta file
	require(seqinr)
	require(Biostrings)
	fastaobj<-read.fasta(fasta,seqtype="AA",as.string=TRUE)
	fastaacc<-substr(names(fastaobj),start=4,stop=9)
	datamat<-object@data
	#proteins<-levels(datamat[,"protein"])
	proteins<-as.character(unique(datamat[,"protein"]))
	uniqproteins<-substr(start=4,stop=9,proteins)
	allprotacc<-substr(as.character(datamat[,"protein"]),start=4,stop=9)
	### replace object@pepvec with cleaned peptides
	peptempvec<-as.character(object@pepvec)
	peptempvec.l<-length(peptempvec)
	### make empty vector for protein positions
	proteinposition<-rep(0,times=peptempvec.l)
	####  loop through each line in the datamat
	proteinposition<-list()
	#### fix to correctly assign n-terminal position
	for(i in 1:peptempvec.l){
		print(i)
		#print(peptempvec[i])
		#print(datamat[i,"protein"])		
		#currentpro<-substr(start=4,stop=9,datamat[i,"protein"])
		currentprot<-fastaobj[[which(fastaacc==substr(start=4,stop=9,datamat[i,"protein"]))]][1]
		currentpep<-peptempvec[i]
		#print(currentpro)
		#print(currentpep)
		matched<-matchPattern(currentpep,currentprot)
		proteinposition[[i]]<-matched@ranges[[1]][1]+(unlist(object@modposition.peptide[i])-1)
		}
	object@data<-cbind(protein.position=as.character(proteinposition),object@data)
	#object@data[1,]
	object@modposition.protein<-proteinposition
	#name="testMGneg.tsv"
	if(writetsv==TRUE){
		write.table(object@data,file=name,quote=FALSE,sep="\t",row.names=F)
		}
	return(object)
	}
